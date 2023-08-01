
using HDF5
using Rotations
using WignerD

"Sample provider which yields initial electron samples through WFAT formula."
struct WFATSampler <: ElectronSampleProvider
    laser           ::Laser;
    target          ::Molecule;     # WFAT only supports [Molecule]
    t_samples       ::AbstractVector;
    ss_kd_samples   ::AbstractVector;
    ss_kz_samples   ::AbstractVector;
    cutoff_limit    ::Real;
    tun_exit        ::Symbol;       # :Para for tunneling, :IpF for over-barrier, automatically specified.
    ion_orbit_idx   ::Integer;

    function WFATSampler(;  laser               ::Laser,
                            target              ::Molecule,
                            sample_t_intv       ::Tuple{<:Real,<:Real},
                            sample_t_num        ::Int,
                            ss_kd_max           ::Real,
                            ss_kd_num           ::Int,
                            ss_kz_max           ::Real,
                            ss_kz_num           ::Int,
                            sample_cutoff_limit ::Real,
                            mol_orbit_idx       ::Int,
                            kwargs...   # kwargs are surplus params.
                            )
        # check sampling parameters.
        @assert (sample_t_num>0) "[WFATSampler] Invalid time sample number $sample_t_num."
        @assert (sample_cutoff_limit≥0) "[WFATSampler] Invalid cut-off limit $sample_cutoff_limit."
        @assert (ss_kd_num>0 && ss_kz_num>0) "[WFATSampler] Invalid kd/kz sample number $ss_kd_num/$ss_kz_num."
        # load WFAT IntData.
        if ! (mol_orbit_idx in MolWFATAvailableIndices(target))
            MolCalcWFATData!(target, mol_orbit_idx)
        end
        # check Keldysh parameter & over-barrier condition.
        Ip = IonPotential(target, mol_orbit_idx)
        F0 = LaserF0(laser)
        γ0 = AngFreq(laser) * sqrt(2Ip) / F0
        if γ0 ≥ 0.5
            @warn "[WFATSampler] Keldysh parameter γ=$γ0, adiabatic (tunneling) condition [γ<<1] not sufficiently satisfied."
        elseif γ0 ≥ 1.0
            @warn "[WFATSampler] Keldysh parameter γ=$γ0, adiabatic (tunneling) condition [γ<<1] unsatisfied."
        end
        F_crit = Ip^2/4/(1-sqrt(Ip/2))
        tun_exit = :Para
        if F0 ≥ F_crit
            @warn "[WFATSampler] Peak electric field strength F0=$F0, reaching the over-barrier critical value, weak-field condition unsatisfied. Tunneling exit method switched from [Para] to [IpF]."
            tun_exit = :IpF
        elseif F0 ≥ F_crit*2/3
            @warn "[WFATSampler] Peak electric field strength F0=$F0, reaching 2/3 of over-barrier critical value, weak-field condition not sufficiently satisfied."
        end
        # finish initialization.
        return new( laser, target,
                    range(sample_t_intv[1],sample_t_intv[2];length=sample_t_num),
                    range(-abs(ss_kd_max),abs(ss_kd_max);length=ss_kd_num), range(-abs(ss_kz_max),abs(ss_kz_max);length=ss_kz_num),
                    sample_cutoff_limit, tun_exit,
                    mol_orbit_idx
                    )
    end
end

"Gets the total number of batches."
function batch_num(sp::WFATSampler)
    return length(sp.t_samples)
end

"Generates a batch of electrons of `batchId` from `sp` using WFAT method."
function gen_electron_batch(sp::WFATSampler, batchId::Int)
    t = sp.t_samples[batchId]
    Fx::Function = LaserFx(sp.laser)
    Fy::Function = LaserFy(sp.laser)
    Fxt = Fx(t)
    Fyt = Fy(t)
    Ft = hypot(Fxt,Fyt)
    φ_field = atan( Fyt, Fxt)   # direction of field vector F.
    φ_exit  = atan(-Fyt,-Fxt)   # direction of tunneling exit, which is opposite to F.
    Ip = IonPotential(sp.target, sp.ion_orbit_idx)
    κ  = sqrt(2Ip)
    Z  = AsympNuclCharge(sp.target)
    cutoff_limit = sp.cutoff_limit
    if Ft == 0 || ADKRateExp(sp.target)(Ip,Ft,0.0,0.0) < cutoff_limit
        return nothing
    end

    # determining tunneling exit position (using ADK's parabolic tunneling exit method if tunExit=:Para)
    r_exit =
        if sp.tun_exit == :Para
            Ip_eff -> (Ip_eff + sqrt(Ip_eff^2 - 4*(1-sqrt(Ip_eff/2))*Ft)) / 2Ft
        else
            Ip_eff -> Ip_eff / Ft
        end

    # determining Euler angles (β,γ) (α contributes nothing to Γ, thus is neglected)
    mα,mβ,mγ = MolRotation(sp.target)
    RotMol = RotZYZ(mγ, mβ, mα)  # the molecule's rotation
    RotLaser = RotMatrix3([[ 0  -sin(φ_field)  cos(φ_field)];
                           [ 0   cos(φ_field)  sin(φ_field)];
                           [-1              0             0]])          # the laser's rotation (directions of x&y can be arbitrary)
    α,β,γ = Rotations.params(RotZYZ(inv(RotLaser)*RotMol))      # the ZYZ Euler angles of the rotations from laser frame (F points to Z axis) to molecular frame.

    # determining tunneling rate Γ.
    # The total rate Γ consists of partial rates of different channels ν=(nξ,m): Γ = ∑ Γ_ν
    # The partial rate consists of structure factor part |G_ν(β,γ)|² and field factor W_ν(F): Γ_ν = |G_ν(β,γ)|²*W_ν(F)
    nξMax, mMax = MolWFATMaxChannels(sp.target, sp.ion_orbit_idx)
    G_data = zeros(ComplexF64, nξMax+1, 2mMax+1)    # to obtain G_nξ,m, call G_data[nξ+1, m+mMax+1]
    for nξ in 0:nξMax, m in -mMax:mMax
        G_data[nξ+1, m+mMax+1] = MolWFATStructureFactor_G(sp.target,sp.ion_orbit_idx,nξ,m,β,γ)
    end
    rate::Function =
        function (F,kd,kz)
            Γsum = 0.
            for nξ in 0:nξMax, m in -mMax:mMax
                G2 = abs2(G_data[nξ+1, m+mMax+1])  # G², structural part
                WF = (κ/2) * (4κ^2/F)^(2Z/κ-2nξ-abs(m)-1) * exp(-2(κ^2+kd^2+kz^2)^1.5/3F)   # W_F, field part
                Γsum += G2*WF
            end
            return Γsum
        end

    # generating samples
    dim = 8
    kd_samples = sp.ss_kd_samples
    kz_samples = sp.ss_kz_samples
    kdNum, kzNum = length(kd_samples), length(kz_samples)

    sample_count_thread = zeros(Int,nthreads())
    init_thread = zeros(Float64, dim, nthreads(), kdNum*kzNum) # initial condition (support for multi-threading)

    @threads for ikd in 1:kdNum
        for ikz in 1:kzNum
            kd0, kz0 = kd_samples[ikd], kz_samples[ikz]
            kx0 = kd0*-sin(φ_exit)
            ky0 = kd0* cos(φ_exit)
            r0 = r_exit(Ip+(kd0^2+kz0^2)/2)
            x0 = r0*cos(φ_exit)
            y0 = r0*sin(φ_exit)
            z0 = 0.0
            rate_ = rate(Ft,kd0,kz0)
            if rate_ < cutoff_limit
                continue    # discard the sample
            end
            sample_count_thread[threadid()] += 1
            init_thread[1:8,threadid(),sample_count_thread[threadid()]] = [x0,y0,z0,kx0,ky0,kz0,t,rate_]
        end
    end
    if sum(sample_count_thread) == 0
        # @warn "[WFATSampler] All sampled electrons are discarded in batch #$(batchId), corresponding to t=$t."
        return nothing
    end
    # collect electron samples from different threads.
    init = zeros(Float64, dim, sum(sample_count_thread))
    N = 0
    for i in 1:nthreads()
        init[:,N+1:N+sample_count_thread[i]] = init_thread[:,i,1:sample_count_thread[i]]
        N += sample_count_thread[i]
    end
    return init
end

