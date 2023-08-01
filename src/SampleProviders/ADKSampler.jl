
using LinearAlgebra

"Sample provider which yields initial electron samples through ADK rate formula."
struct ADKSampler <: ElectronSampleProvider
    laser           ::Laser;
    target          ::SAEAtomBase;      # ADK only supports [SAEAtomBase].
    monte_carlo     ::Bool;
    t_samples       ::AbstractVector;
    ss_kd_samples   ::AbstractVector;
    ss_kz_samples   ::AbstractVector;
    mc_kp_num       ::Integer;
    mc_kp_max       ::Real;
    cutoff_limit    ::Real;
    phase_method    ::Symbol;           # currently supports :CTMC, :QTMC, :SCTS.
    rate_prefix     ::Symbol;           # currently supports :ExpRate, :ExpPre, :ExpJac, :Full.
    ADK_tun_exit    ::Symbol;           # currently supports :IpF, :FDM, :Para.
    function ADKSampler(;   laser               ::Laser,
                            target              ::SAEAtomBase,
                            sample_t_intv       ::Tuple{<:Real,<:Real},
                            sample_t_num        ::Int,
                            sample_cutoff_limit ::Real,
                            sample_monte_carlo  ::Bool,
                            traj_phase_method   ::Symbol,
                            rate_prefix         ::Symbol,
                            adk_tun_exit        ::Symbol,
                                #* for step-sampling (!rate_monteCarlo)
                            ss_kd_max           ::Real,
                            ss_kd_num           ::Int,
                            ss_kz_max           ::Real,
                            ss_kz_num           ::Int,
                                #* for Monte-Carlo-sampling (rate_monteCarlo)
                            mc_kp_num           ::Int,
                            mc_kp_max           ::Real,
                            kwargs...   # kwargs are surplus params.
                            )
        F0 = LaserF0(laser)
        Ip = IonPotential(target)
        γ0 = AngFreq(laser) * sqrt(2Ip) / F0
        # check phase method support.
        if ! (traj_phase_method in [:CTMC, :QTMC, :SCTS])
            error("[ADKSampler] Undefined phase method [$traj_phase_method].")
            return
        end
        # check tunneling exit support.
        if ! (adk_tun_exit in [:IpF, :FDM, :Para])
            error("[ADKSampler] Undefined tunneling exit method [$adk_tun_exit].")
            return
        end
        # switch tunneling exit method if inappropriate.
        if adk_tun_exit == :FDM && Ip^2 < 4F0
            adk_tun_exit == :IpF
            @warn "[ADKSampler] Tunneling exit method has changed from [FDM] to [IpF] due to failure in tunneling exit determination of [FDM] model."
        elseif adk_tun_exit == :Para && Ip^2 < 4*(1-sqrt(Ip/2))*F0
            adk_tun_exit == :IpF
            @warn "[ADKSampler] Tunneling exit method has changed from [Para] to [IpF] due to failure in tunneling exit determination of [Para] model."
        end
        # check rate prefix support.
        if ! (rate_prefix in [:ExpRate, :ExpPre, :ExpJac, :Full])
            error("[ADKSampler] Undefined tunneling rate prefix [$rate_prefix].")
            return
        end
        # check Keldysh parameter.
        if γ0 ≥ 0.5
            @warn "[ADKSampler] Keldysh parameter γ=$γ0, adiabatic (tunneling) condition [γ<<1] not sufficiently satisfied."
        elseif γ0 ≥ 1.0
            @warn "[ADKSampler] Keldysh parameter γ=$γ0, adiabatic (tunneling) condition [γ<<1] unsatisfied."
        end
        # check sampling parameters.
        @assert (sample_t_num>0) "[ADKSampler] Invalid time sample number $sample_t_num."
        @assert (sample_cutoff_limit≥0) "[ADKSampler] Invalid cut-off limit $sample_cutoff_limit."
        if ! sample_monte_carlo    # check SS sampling parameters.
            @assert (ss_kd_num>0 && ss_kz_num>0) "[ADKSampler] Invalid kd/kz sample number $ss_kd_num/$ss_kz_num."
        else                    # check MC sampling parameters.
            @assert (sample_t_intv[1] < sample_t_intv[2]) "[ADKSampler] Invalid sampling time interval $sample_t_intv."
            @assert (mc_kp_num>0) "[ADKSampler] Invalid sampling kp_num $mc_kp_num."
            @assert (mc_kp_max>0) "[ADKSampler] Invalid sampling kp_max $mc_kp_max."
        end
        # finish initialization.
        return if ! sample_monte_carlo
            new(laser,target,
                sample_monte_carlo,
                range(sample_t_intv[1],sample_t_intv[2];length=sample_t_num),
                range(-abs(ss_kd_max),abs(ss_kd_max);length=ss_kd_num), range(-abs(ss_kz_max),abs(ss_kz_max);length=ss_kz_num),
                0,0,    # for MC params. pass meaningless values
                sample_cutoff_limit,traj_phase_method,rate_prefix,adk_tun_exit)
        else
            t_samples = rand(sample_t_num) .* (sample_t_intv[2]-sample_t_intv[1]) .+ sample_t_intv[1]
            new(laser,target,
                sample_monte_carlo,
                t_samples,
                0:0,0:0,    # for SS params. pass meaningless values
                mc_kp_num, mc_kp_max,
                sample_cutoff_limit,traj_phase_method,rate_prefix,adk_tun_exit)
        end
    end
end
"Gets the total number of batches."
function batch_num(sp::ADKSampler)
    return length(sp.t_samples)
end
"Generates a batch of electrons of `batchId` from `sp` using ADK method."
function gen_electron_batch(sp::ADKSampler, batchId::Int)
    t = sp.t_samples[batchId]
    Fx::Function = LaserFx(sp.laser)
    Fy::Function = LaserFy(sp.laser)
    Fxt = Fx(t)
    Fyt = Fy(t)
    Ft = hypot(Fxt,Fyt)
    Ax::Function = LaserAx(sp.laser)
    Ay::Function = LaserAy(sp.laser)
    φ  = atan(-Fyt,-Fxt)
    Z  = AsympNuclCharge(sp.target)
    Ip = IonPotential(sp.target)
    cutoff_limit = sp.cutoff_limit
    if Ft == 0 || ADKRateExp(sp.target)(Ft,0.0,0.0) < cutoff_limit
        return nothing
    end
    # determining tunneling exit position. here Ip_eff=Ip+k0^2/2
    r_exit =
        if      sp.ADK_tun_exit == :IpF
            Ip_eff -> Ip_eff / Ft
        elseif  sp.ADK_tun_exit == :FDM
            Ip_eff -> (Ip_eff + sqrt(Ip_eff^2 - 4Ft)) / 2Ft
        elseif  sp.ADK_tun_exit == :Para
            Ip_eff -> (Ip_eff + sqrt(Ip_eff^2 - 4*(1-sqrt(Ip_eff/2))*Ft)) / 2Ft
        end
    # determining ADK rate. ADKRate(F,φ,kd,kz)
    ADKRate::Function =
        if      sp.rate_prefix == :ExpRate
            ADKRateExp(sp.target)
        else
            ADKRate_Exp = ADKRateExp(sp.target)
            α = 1.0+Z/sqrt(2Ip)
            if  sp.rate_prefix == :ExpPre
                function ADKRate_ExpPre(F,kd,kz)
                    return ADKRate_Exp(F,kd,kz) / ((kd^2+kz^2+2Ip)*Ft^2)^(0.5α)
                end
            elseif  sp.rate_prefix == :ExpJac
                function ADKRate_ExpJac(F,kd,kz)
                    transform((tr,kd)) = SVector(kd*Fy(tr)/sqrt(Fx(tr)^2+Fy(tr)^2) - Ax(tr), -kd*Fx(tr)/sqrt(Fx(tr)^2+Fy(tr)^2) - Ay(tr))
                    jac = abs(det(ForwardDiff.jacobian(transform,SVector(t,kd))))
                    return ADKRate_Exp(F,kd,kz) * jac
                end
            elseif  sp.rate_prefix == :Full
                function ADKRate_Full(F,kd,kz)
                    transform((tr,kd)) = SVector(kd*Fy(tr)/sqrt(Fx(tr)^2+Fy(tr)^2) - Ax(tr), -kd*Fx(tr)/sqrt(Fx(tr)^2+Fy(tr)^2) - Ay(tr))
                    jac = abs(det(ForwardDiff.jacobian(transform,SVector(t,kd))))
                    return ADKRate_Exp(F,kd,kz) * jac / ((kd^2+kz^2+2Ip)*Ft^2)^(0.5α)
                end
            end
        end
    phase_method = sp.phase_method
    dim = (phase_method == :CTMC) ? 8 : 9 # x,y,z,kx,ky,kz,t0,rate[,phase]

    sample_count_thread = zeros(Int,nthreads())
    init_thread = if ! sp.monte_carlo
        zeros(Float64, dim, nthreads(), length(sp.ss_kd_samples)*length(sp.ss_kz_samples)) # initial condition (support for multi-threading)
    else
        zeros(Float64, dim, nthreads(), sp.mc_kp_num)
    end

    if ! sp.monte_carlo
        kd_samples = sp.ss_kd_samples
        kz_samples = sp.ss_kz_samples
        kdNum, kzNum = length(kd_samples), length(kz_samples)
        @threads for ikd in 1:kdNum
            for ikz in 1:kzNum
                kd0, kz0 = kd_samples[ikd], kz_samples[ikz]
                kx0 = kd0*-sin(φ)
                ky0 = kd0* cos(φ)
                r0 = r_exit(Ip+(kd0^2+kz0^2)/2)
                x0 = r0*cos(φ)
                y0 = r0*sin(φ)
                z0 = 0.0
                rate = ADKRate(Ft,kd0,kz0)
                if rate < cutoff_limit
                    continue    # discard the sample
                end
                sample_count_thread[threadid()] += 1
                init_thread[1:8,threadid(),sample_count_thread[threadid()]] = [x0,y0,z0,kx0,ky0,kz0,t,rate]
                if phase_method != :CTMC
                    init_thread[9,threadid(),sample_count_thread[threadid()]] = 0.0
                end
            end
        end
    else
        kpMax = sp.mc_kp_max
        @threads for i in 1:sp.mc_kp_num
            # generates random (kd0,kz0) inside circle kd0^2+kz0^2=ktMax^2.
            kd0, kz0 = (rand()-0.5)*2kpMax, (rand()-0.5)*2kpMax
            while kd0^2+kz0^2 > kpMax^2
                kd0, kz0 = (rand()-0.5)*2kpMax, (rand()-0.5)*2kpMax
            end
            kx0 = kd0*-sin(φ)
            ky0 = kd0* cos(φ)
            r0 = r_exit(Ip+(kd0^2+kz0^2)/2)
            x0 = r0*cos(φ)
            y0 = r0*sin(φ)
            z0 = 0.0
            rate = ADKRate(Ft,kd0,kz0)
            if rate < cutoff_limit
                continue    # discard the sample
            end
            sample_count_thread[threadid()] += 1
            init_thread[1:8,threadid(),sample_count_thread[threadid()]] = [x0,y0,z0,kx0,ky0,kz0,t,rate]
            if phase_method != :CTMC
                init_thread[9,threadid(),sample_count_thread[threadid()]] = 0.0
            end
        end
    end
    if sum(sample_count_thread) == 0
        # @warn "[ADKSampler] All sampled electrons are discarded in batch #$(batchId), corresponding to t=$t."
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