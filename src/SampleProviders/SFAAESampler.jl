
using Base.Threads
using ForwardDiff

"Sample provider which yields initial electron samples through SFA-AE formula."
struct SFAAESampler <: ElectronSampleProvider
    laser           ::Laser;
    target          ::SAEAtomBase;      # SFA-AE only supports [SAEAtomBase].
    t_samples       ::AbstractVector;
    ss_kd_samples   ::AbstractVector;
    ss_kz_samples   ::AbstractVector;
    cutoff_limit    ::Real;
    phase_method    ::Symbol;           # currently supports :CTMC, :QTMC, :SCTS.
    rate_prefix     ::Symbol;           # currently supports :ExpRate, :ExpPre, :ExpJac, :Full.
    function SFAAESampler(; laser               ::Laser,
                            target              ::SAEAtomBase,
                            sample_t_intv       ::Tuple{<:Real,<:Real},
                            sample_t_num        ::Integer,
                            sample_cutoff_limit ::Real,
                            traj_phase_method   ::Symbol,
                            rate_prefix         ::Symbol,
                            ss_kd_max           ::Real,
                            ss_kd_num           ::Integer,
                            ss_kz_max           ::Real,
                            ss_kz_num           ::Integer,
                            kwargs...   # kwargs are surplus params.
                            )
        F0 = LaserF0(laser)
        Ip = IonPotential(target)
        γ0 = AngFreq(laser) * sqrt(2Ip) / F0
        # check phase method support.
        if ! (traj_phase_method in [:CTMC, :QTMC, :SCTS])
            error("[SFAAESampler] Undefined phase method [$traj_phase_method].")
            return
        end
        # check rate prefix support.
        if ! (rate_prefix in [:ExpRate, :ExpPre, :ExpJac, :Full])
            error("[SFAAESampler] Undefined tunneling rate prefix [$rate_prefix].")
            return
        end
        # check Keldysh paramater.
        if γ0 ≥ 1.0
            @warn "[SFAAESampler] Keldysh parameter γ=$γ0, adiabatic (tunneling) condition [γ<<1] unsatisfied."
        end
        # check sampling parameters.
        @assert (sample_t_num>0) "[SFAAESampler] Invalid time sample number $sample_t_num."
        @assert (sample_cutoff_limit≥0) "[SFAAESampler] Invalid cut-off limit $sample_cutoff_limit."
        @assert (ss_kd_num>0 && ss_kz_num>0) "[SFAAESampler] Invalid kd/kz sample number $ss_kd_num/$ss_kz_num."
        # finish initialization.
        return new(laser, target,
                range(sample_t_intv[1],sample_t_intv[2];length=sample_t_num),
                range(-abs(ss_kd_max),abs(ss_kd_max);length=ss_kd_num), range(-abs(ss_kz_max),abs(ss_kz_max);length=ss_kz_num),
                sample_cutoff_limit, traj_phase_method, rate_prefix
                )
    end
end

"Gets the total number of batches."
function batch_num(sp::SFAAESampler)
    return length(sp.t_samples)
end

"Generates a batch of electrons of `batchId` from `sp` using SFA-AE method."
function gen_electron_batch(sp::SFAAESampler, batchId::Int)
    t = sp.t_samples[batchId]
    Fx::Function = LaserFx(sp.laser)
    Fy::Function = LaserFy(sp.laser)
    dFxt = ForwardDiff.derivative(Fx,t)
    dFyt = ForwardDiff.derivative(Fy,t)
    Fxt = Fx(t)
    Fyt = Fy(t)
    Ft  = hypot(Fxt,Fyt)
    Ax::Function = LaserAx(sp.laser)
    Ay::Function = LaserAy(sp.laser)
    φ   = atan(-Fyt,-Fxt)
    Z   = AsympNuclCharge(sp.target)
    Ip  = IonPotential(sp.target)
    cutoff_limit = sp.cutoff_limit
    if Ft == 0 || ADKRateExp(sp.target)(Ft,0.0,0.0) < cutoff_limit
        return nothing
    end
    rate::Function =
        if  sp.rate_prefix == :ExpRate
            ADKRateExp(sp.target)
        else
            ADKRate_Exp = ADKRateExp(sp.target)
            α = 1.0+Z/sqrt(2Ip)
            if  sp.rate_prefix == :ExpPre
                @inline function ADKRate_ExpPre(F,kd,kz)
                    return ADKRate_Exp(F,kd,kz) / ((kd^2+kz^2+2Ip)*Ft^2)^(0.5α)
                end
            elseif  sp.rate_prefix == :ExpJac
                @inline function ADKRate_ExpJac(F,kd,kz)
                    transform((tr,kd_)) = SVector(kd_*Fy(tr)/sqrt(Fx(tr)^2+Fy(tr)^2) - Ax(tr), -kd_*Fx(tr)/sqrt(Fx(tr)^2+Fy(tr)^2) - Ay(tr))
                    jac = abs(det(ForwardDiff.jacobian(transform,SVector(t,kd))))
                    return ADKRate_Exp(F,kd,kz) * jac
                end
            elseif  sp.rate_prefix == :Full
                @inline function ADKRate_Full(F,kd,kz)
                    transform((tr,kd_)) = SVector(kd_*Fy(tr)/sqrt(Fx(tr)^2+Fy(tr)^2) - Ax(tr), -kd_*Fx(tr)/sqrt(Fx(tr)^2+Fy(tr)^2) - Ay(tr))
                    jac = abs(det(ForwardDiff.jacobian(transform,SVector(t,kd))))
                    return ADKRate_Exp(F,kd,kz) * jac / ((kd^2+kz^2+2Ip)*Ft^2)^(0.5α)
                end
            end
        end
    @inline F2eff(kx,ky) = Ft^2 - (kx*dFxt+ky*dFyt)  # F2eff=F²-p⟂⋅F'
    @inline r_exit(kx,ky) = Ft/2*(kx^2+ky^2+2Ip)/F2eff(kx,ky)
    phase_method = sp.phase_method
    dim = (phase_method == :CTMC) ? 8 : 9 # x,y,z,kx,ky,kz,t0,rate[,phase]
    kd_samples = sp.ss_kd_samples
    kz_samples = sp.ss_kz_samples
    kdNum, kzNum = length(kd_samples), length(kz_samples)

    sample_count_thread = zeros(Int,nthreads())
    init_thread = zeros(Float64, dim, nthreads(), kdNum*kzNum) # initial condition (support for multi-threading)

    @threads for ikd in 1:kdNum
        for ikz in 1:kzNum
            kd0, kz0 = kd_samples[ikd], kz_samples[ikz]
            kx0 = kd0*-sin(φ)
            ky0 = kd0* cos(φ)
            F2eff_ = F2eff(kx0,ky0)
            if F2eff_ ≤ 0
                continue
            end
            r0 = r_exit(kx0,ky0)
            x0 = r0*cos(φ)
            y0 = r0*sin(φ)
            z0 = 0.0
            rate_ = rate(sqrt(F2eff_),kd0,kz0)
            if rate_ < cutoff_limit
                continue
            end
            sample_count_thread[threadid()] += 1
            init_thread[1:8,threadid(),sample_count_thread[threadid()]] = [x0,y0,z0,kx0,ky0,kz0,t,rate_]
            if phase_method != :CTMC
                init_thread[9,threadid(),sample_count_thread[threadid()]] = 0.0
            end
        end
    end
    if sum(sample_count_thread) == 0
        # @warn "[SFAAESampler] All sampled electrons are discarded in batch #$(batchId), corresponding to t=$t."
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