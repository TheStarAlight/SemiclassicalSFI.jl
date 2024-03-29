
"Represents a monochromatic elliptically polarized laser field with Cos2-shape envelope propagating in z direction."
struct Cos2Laser <: MonochromaticLaser
    "Peak intensity of the laser field (in W/cm^2)."
    peak_int;
    "Wavelength of the laser field (in NANOMETER)."
    wave_len;
    "Cycle number of the laser field."
    cyc_num;
    "Ellipticity of the laser field."
    ellip;
    "Azimuth angle of the laser's polarization's principle axis relative to x axis (in radians)."
    azi;
    "Carrier-Envelope-Phase (CEP) of the laser field."
    cep;
    "Time shift of the laser relative to the peak (in a.u.)."
    t_shift;
    """
    Constructs a new monochromatic elliptically polarized laser field with Cos2-shape envelope.
    # Parameters
    - `peakInt`     : Peak intensity of the laser field (in W/cm²).
    - `WaveLen`     : Wavelength of the laser field (in nm).
    - `cycNum`      : Number of cycles of the laser field.
    - `ellip`       : Ellipticity of the laser field [-1≤e≤1, 0 indicates linear polarization and ±1 indicates circular polarization].
    - `azi`         : Azimuth angle of the laser's polarization's principle axis relative to x axis (in radians) (optional, default 0).
    - `cep`         : Carrier-Envelope-Phase of the laser field (optional, default 0).
    - `t_shift`     : Time shift of the laser (in a.u.) relative to the peak (optional, default 0).
    """
    function Cos2Laser(peakInt, waveLen, cycNum, ellip, azi=0., cep=0., t_shift=0.)
        @assert peakInt>0   "[Cos2Laser] Peak intensity must be positive."
        @assert waveLen>0   "[Cos2Laser] Wavelength must be positive."
        @assert cycNum>0    "[Cos2Laser] Cycle number must be positive."
        @assert -1≤ellip≤1  "[Cos2Laser] Ellipticity must be in [-1,1]."
        new(peakInt,waveLen,cycNum,ellip,azi,cep,t_shift)
    end
    """
    Constructs a new monochromatic elliptically polarized laser field with Cos2-shape envelope.
    # Parameters
    - `peakInt`     : Peak intensity of the laser field (in W/cm²).
    - `WaveLen`     : Wave length of the laser field (in nm). Must specify either `waveLen` or `angFreq`.
    - `angFreq`     : Angular frequency of the laser field (in a.u.). Must specify either `waveLen` or `angFreq`.
    - `cycNum`      : Number of cycles of the laser field. Must specify either `cycNum` or `duration`.
    - `duration`    : Duration of the laser field (in a.u.). Must specify either `cycNum` or `duration`.
    - `ellip`       : Ellipticity of the laser field [-1≤e≤1, 0 indicates linear polarization and ±1 indicates circular polarization].
    - `azi`         : Azimuth angle of the laser's polarization's principle axis relative to x axis (in radians) (optional, default 0).
    - `cep`         : Carrier-Envelope-Phase of the laser field (optional, default 0).
    - `t_shift`     : Time shift of the laser (in a.u.) relative to the peak (optional, default 0).
    """
    function Cos2Laser(;peakInt,
                        waveLen=-1, angFreq=-1,     # must specify either waveLen or angFreq.
                        cycNum=-1,  duration=-1,    # must specify either cycNum or duration.
                        ellip, azi=0., cep=0., t_shift=0.)
        @assert waveLen>0 || angFreq>0  "[Cos2Laser] Must specify either waveLen or angFreq."
        @assert cycNum>0 || duration>0  "[Cos2Laser] Must specify either cycNum or duration."
        if waveLen>0 && angFreq>0
            @warn "[Cos2Laser] Both waveLen & angFreq are specified, will use waveLen."
        end
        if cycNum>0 && duration>0
            @warn "[Cos2Laser] Both cycNum & duration are specified, will use cycNum."
        end
        if waveLen==-1
            waveLen = 45.563352525 / angFreq
        end
        if cycNum==-1
            cycNum = duration / (2π/angFreq)
        end
        Cos2Laser(peakInt,waveLen,cycNum,ellip,azi,cep,t_shift)
    end
end
"Gets the peak intensity of the laser field (in W/cm²)."
PeakInt(l::Cos2Laser) = l.peak_int
"Gets the wave length of the laser field (in nm)."
WaveLen(l::Cos2Laser) = l.wave_len
"Gets the cycle number of the laser field."
CycNum(l::Cos2Laser) = l.cyc_num
"Gets the ellipticity of the laser field."
Ellipticity(l::Cos2Laser) = l.ellip
"Gets the azimuth angle of the laser's polarization's principle axis relative to x axis (in radians)."
Azimuth(l::Cos2Laser) = l.azi
"Gets the angular frequency (ω) of the laser field (in a.u.)."
AngFreq(l::Cos2Laser) = 45.563352525 / l.wave_len
"Gets the period of the laser field (in a.u.)."
Period(l::Cos2Laser) = 2π / AngFreq(l)
"Gets the time shift relative to the peak (in a.u.)."
TimeShift(l::Cos2Laser) = l.t_shift
"Gets the peak electric field intensity of the laser field (in a.u.)."
LaserF0(l::Cos2Laser) = sqrt(l.peak_int/(1.0+l.ellip^2)/3.50944521e16)
"Gets the peak vector potential intensity of the laser field (in a.u.)."
LaserA0(l::Cos2Laser) = LaserF0(l) / AngFreq(l)

"Gets the time-dependent x component of the vector potential under dipole approximation."
function LaserAx(l::Cos2Laser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local N = l.cycNum; local φ = l.cep; local Δt = l.t_shift; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            A0 * cos(ω*t/(2N))^2 * (abs(ω*real(t))<N*π) * cos(ω*t+φ)
        end
    else
        function(t)
            t -= Δt
            A0 * cos(ω*t/(2N))^2 * (abs(ω*real(t))<N*π) * (cos(ω*t+φ)*cos(ϕ)+sin(ω*t+φ)*ε*sin(ϕ))
        end
    end
end
"Gets the time-dependent y component of the vector potential under dipole approximation."
function LaserAy(l::Cos2Laser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local N = l.cycNum; local φ = l.cep; local Δt = l.t_shift; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            A0 * cos(ω*t/(2N))^2 * (abs(ω*real(t))<N*π) * sin(ω*t+φ) * ε
        end
    else
        function(t)
            t -= Δt
            A0 * cos(ω*t/(2N))^2 * (abs(ω*real(t))<N*π) * (cos(ω*t+φ)*-sin(ϕ)+sin(ω*t+φ)*ε*cos(ϕ))
        end
    end
end
"Gets the time-dependent x component of the electric field strength under dipole approximation."
function LaserFx(l::Cos2Laser)
    local F0 = LaserF0(l); local ω = AngFreq(l); local N = l.cycNum; local φ = l.cep; local Δt = l.t_shift; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            F0 * cos(ω*t/(2N)) * (abs(ω*real(t))<N*π) * ( cos(ω*t/(2N))*sin(ω*t+φ) + 1/N*sin(ω*t/(2N))*cos(ω*t+φ))
        end
    else
        function(t)
            t -= Δt
            F0 * cos(ω*t/(2N)) * (abs(ω*real(t))<N*π) * ( (cos(ω*t/(2N))*sin(ω*t+φ) + 1/N*sin(ω*t/(2N))*cos(ω*t+φ))*cos(ϕ) - (cos(ω*t/(2N))*cos(ω*t+φ) - 1/N*sin(ω*t/(2N))*sin(ω*t+φ))*ε*sin(ϕ) )
        end
    end
end
"Gets the time-dependent y component of the electric field strength under dipole approximation."
function LaserFy(l::Cos2Laser)
    local F0 = LaserF0(l); local ω = AngFreq(l); local N = l.cycNum; local φ = l.cep; local Δt = l.t_shift; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            F0 * cos(ω*t/(2N)) * (abs(ω*real(t))<N*π) * ( -cos(ω*t/(2N))*cos(ω*t+φ) + 1/N*sin(ω*t/(2N))*sin(ω*t+φ)) * ε
        end
    else
        function(t)
            t -= Δt
            F0 * cos(ω*t/(2N)) * (abs(ω*real(t))<N*π) * ( (cos(ω*t/(2N))*sin(ω*t+φ) + 1/N*sin(ω*t/(2N))*cos(ω*t+φ))*-sin(ϕ) - (cos(ω*t/(2N))*cos(ω*t+φ) - 1/N*sin(ω*t/(2N))*sin(ω*t+φ))*ε*cos(ϕ) )
        end
    end
end

"Prints the information about the laser."
Base.show(io::IO, l::Cos2Laser) = print(io,"[MonochromaticLaser] Envelope cos², Wavelength=$(l.waveLen) nm, $(l.cycNum) cycle(s), e=$(l.ellip)"
                                           * (l.ellip==0 ? " [Linearly polarized]" : "") * (abs(l.ellip)==1 ? " [Circularly polarized]" : "")
                                           * ", PrincipleAxisAzimuth=$(l.azi/π*180)°" * (l.t_shift==0 ? "" : ", Peaks at t₀=$(l.t_shift) a.u.") * (l.cep==0 ? "" : ", CEP=$(l.cep)"))