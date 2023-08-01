
"Represents a monochromatic elliptically polarized laser field with Cos4-shape envelope propagating in z direction."
struct Cos4Laser <: MonochromaticLaser
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
    Constructs a new monochromatic elliptically polarized laser field with Cos4-shape envelope.
    # Parameters
    - `peak_int`    : Peak intensity of the laser field (in W/cm²).
    - `wave_len`    : Wavelength of the laser field (in nm).
    - `cyc_num`     : Number of cycles of the laser field.
    - `ellip`       : Ellipticity of the laser field [-1≤e≤1, 0 indicates linear polarization and ±1 indicates circular polarization].
    - `azi`         : Azimuth angle of the laser's polarization's principle axis relative to x axis (in radians) (optional, default 0).
    - `cep`         : Carrier-Envelope-Phase of the laser field (optional, default 0).
    - `t_shift`     : Time shift of the laser (in a.u.) relative to the peak (optional, default 0).
    """
    function Cos4Laser(peak_int, wave_len, cyc_num, ellip, azi=0., cep=0., t_shift=0.)
        @assert peak_int>0  "[Cos4Laser] Peak intensity must be positive."
        @assert wave_len>0  "[Cos4Laser] Wavelength must be positive."
        @assert cyc_num>0   "[Cos4Laser] Cycle number must be positive."
        @assert -1≤ellip≤1  "[Cos4Laser] Ellipticity must be in [-1,1]."
        new(peak_int, wave_len, cyc_num, ellip, azi, cep, t_shift)
    end
    """
    Constructs a new monochromatic elliptically polarized laser field with Cos4-shape envelope.
    # Parameters
    - `peak_int`    : Peak intensity of the laser field (in W/cm²).
    - `wave_len`    : Wave length of the laser field (in nm). Must specify either `wave_len` or `ang_freq`.
    - `ang_freq`    : Angular frequency of the laser field (in a.u.). Must specify either `wave_len` or `ang_freq`.
    - `cyc_num`     : Number of cycles of the laser field. Must specify either `cyc_num` or `duration`.
    - `duration`    : Duration of the laser field (in a.u.). Must specify either `cyc_num` or `duration`.
    - `ellip`       : Ellipticity of the laser field [-1≤ε≤1, 0 indicates linear polarization and ±1 indicates circular polarization].
    - `azi`         : Azimuth angle of the laser's polarization's principle axis relative to x axis (in radians) (optional, default 0).
    - `cep`         : Carrier-Envelope-Phase of the laser field (optional, default 0).
    - `t_shift`     : Time shift of the laser (in a.u.) relative to the peak (optional, default 0).
    """
    function Cos4Laser(;peak_int,
                        wave_len=-1, ang_freq=-1,     # must specify either wave_len or ang_freq.
                        cyc_num=-1,  duration=-1,    # must specify either cyc_num or duration.
                        ellip, azi=0., cep=0., t_shift=0.)
        @assert wave_len>0 || ang_freq>0    "[Cos4Laser] Must specify either wave_len or ang_freq."
        @assert cyc_num>0 || duration>0     "[Cos4Laser] Must specify either cyc_num or duration."
        if wave_len>0 && ang_freq>0
            @warn "[Cos4Laser] Both wave_len & ang_freq are specified, will use wave_len."
        end
        if cyc_num>0 && duration>0
            @warn "[Cos4Laser] Both cyc_num & duration are specified, will use cyc_num."
        end
        if wave_len==-1
            wave_len = 45.563352525 / ang_freq
        else
            ang_freq = 45.563352525 / wave_len
        end
        if cyc_num==-1
            cyc_num = duration / (2π/ang_freq)
        end
        Cos4Laser(peak_int, wave_len, cyc_num, ellip, azi, cep, t_shift)
    end
end
"Gets the peak intensity of the laser field (in W/cm²)."
PeakInt(l::Cos4Laser) = l.peak_int
"Gets the wave length of the laser field (in nm)."
WaveLen(l::Cos4Laser) = l.wave_len
"Gets the cycle number of the laser field."
CycNum(l::Cos4Laser) = l.cyc_num
"Gets the ellipticity of the laser field."
Ellipticity(l::Cos4Laser) = l.ellip
"Gets the azimuth angle of the laser's polarization's principle axis relative to x axis (in radians)."
Azimuth(l::Cos4Laser) = l.azi
"Gets the angular frequency (ω) of the laser field (in a.u.)."
AngFreq(l::Cos4Laser) = 45.563352525 / l.wave_len
"Gets the period of the laser field (in a.u.)."
Period(l::Cos4Laser) = 2π / AngFreq(l)
"Gets the Carrier-Envelope Phase (CEP) of the laser field."
CEP(l::Cos4Laser) = l.cep
"Gets the time shift relative to the peak (in a.u.)."
TimeShift(l::Cos4Laser) = l.t_shift
"Gets the peak electric field intensity of the laser field (in a.u.)."
LaserF0(l::Cos4Laser) = sqrt(l.peak_int/(1.0+l.ellip^2)/3.50944521e16)
"Gets the peak vector potential intensity of the laser field (in a.u.)."
LaserA0(l::Cos4Laser) = LaserF0(l) / AngFreq(l)

"Gets the unit envelope function (the peak value is 1) of the laser field."
function UnitEnvelope(l::Cos4Laser)
    local ω = AngFreq(l); local N = CycNum(l); local Δt = l.t_shift;
    function (t)
        t -= Δt
        cos(ω*t/(2N))^4 * (abs(ω*real(t))<N*π)
    end
end

"Gets the time-dependent x component of the vector potential under dipole approximation."
function LaserAx(l::Cos4Laser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local N = CycNum(l); local φ = l.cep; local Δt = l.t_shift; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            A0 * cos(ω*t/(2N))^4 * (abs(ω*real(t))<N*π) * cos(ω*t+φ)
        end
    else
        function(t)
            t -= Δt
            A0 * cos(ω*t/(2N))^4 * (abs(ω*real(t))<N*π) * (cos(ω*t+φ)*cos(ϕ)+sin(ω*t+φ)*ε*sin(ϕ))
        end
    end
end
"Gets the time-dependent y component of the vector potential under dipole approximation."
function LaserAy(l::Cos4Laser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local N = CycNum(l); local φ = l.cep; local Δt = l.t_shift; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            A0 * cos(ω*t/(2N))^4 * (abs(ω*real(t))<N*π) * sin(ω*t+φ) * ε
        end
    else
        function(t)
            t -= Δt
            A0 * cos(ω*t/(2N))^4 * (abs(ω*real(t))<N*π) * (cos(ω*t+φ)*-sin(ϕ)+sin(ω*t+φ)*ε*cos(ϕ))
        end
    end
end
"Gets the time-dependent x component of the electric field strength under dipole approximation."
function LaserFx(l::Cos4Laser)
    local F0 = LaserF0(l); local ω = AngFreq(l); local N = CycNum(l); local φ = l.cep; local Δt = l.t_shift; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            F0 * cos(ω*t/(2N))^3 * (abs(ω*real(t))<N*π) * ( cos(ω*t/(2N))*sin(ω*t+φ) + 2/N*sin(ω*t/(2N))*cos(ω*t+φ))
        end
    else
        function(t)
            t -= Δt
            F0 * cos(ω*t/(2N))^3 * (abs(ω*real(t))<N*π) * ( (cos(ω*t/(2N))*sin(ω*t+φ) + 2/N*sin(ω*t/(2N))*cos(ω*t+φ))*cos(ϕ) + (-cos(ω*t/(2N))*cos(ω*t+φ) + 2/N*sin(ω*t/(2N))*sin(ω*t+φ))*ε*sin(ϕ) )
        end
    end
end
"Gets the time-dependent y component of the electric field strength under dipole approximation."
function LaserFy(l::Cos4Laser)
    local F0 = LaserF0(l); local ω = AngFreq(l); local N = CycNum(l); local φ = l.cep; local Δt = l.t_shift; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            F0 * cos(ω*t/(2N))^3 * (abs(ω*real(t))<N*π) * (-cos(ω*t/(2N))*cos(ω*t+φ) + 2/N*sin(ω*t/(2N))*sin(ω*t+φ)) * ε
        end
    else
        function(t)
            t -= Δt
            F0 * cos(ω*t/(2N))^3 * (abs(ω*real(t))<N*π) * ( (cos(ω*t/(2N))*sin(ω*t+φ) + 2/N*sin(ω*t/(2N))*cos(ω*t+φ))*-sin(ϕ) + (-cos(ω*t/(2N))*cos(ω*t+φ) + 2/N*sin(ω*t/(2N))*sin(ω*t+φ))*ε*cos(ϕ) )
        end
    end
end

"Prints the information about the laser."
Base.show(io::IO, l::Cos4Laser) = println(io,"[MonochromaticLaser] Envelope cos⁴, Wavelength=$(l.wave_len) nm, $(l.cyc_num) cycle(s), ε=$(l.ellip)"
                                           * (l.ellip==0 ? " [Linearly polarized]" : "") * (abs(l.ellip)==1 ? " [Circularly polarized]" : "")
                                           * ", PrincipleAxisAzimuth=$(l.azi/π*180)°" * (l.t_shift==0 ? "" : ", Peaks at t₀=$(l.t_shift) a.u.") * (l.cep==0 ? "" : ", CEP=$(l.cep)"))

using Parameters, OrderedCollections
"Returns a `Dict{Symbol,Any}` containing properties of the object."
function Serialize(l::Cos4Laser)
    dict = OrderedDict{Symbol,Any}()
    type        = typeof(l)
    peak_int    = l.peak_int
    wave_len    = l.wave_len
    cyc_num     = l.cyc_num
    ellip       = l.ellip
    azi         = l.azi
    cep         = l.cep
    t_shift     = l.t_shift
    @pack! dict = (type, peak_int, wave_len, cyc_num, ellip, azi, cep, t_shift)
    return dict
end
