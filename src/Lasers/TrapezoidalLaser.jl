
"Represents a monochromatic elliptically polarized laser field with Trapezoidal-shape envelope propagating in z direction."
struct TrapezoidalLaser <: MonochromaticLaser
    "Peak intensity of the laser field (in W/cm^2)."
    peak_int;
    "Wavelength of the laser field (in NANOMETER)."
    wave_len;
    "Cycle number of the laser field in the turn-on."
    cyc_num_turn_on;
    "Cycle number of the laser field in the turn-off."
    cyc_num_turn_off;
    "Cycle number of the laser field in the constant-intensity."
    cyc_num_const;
    "Ellipticity of the laser field."
    ellip;
    "Azimuth angle of the laser's polarization's principle axis relative to x axis (in radians)."
    azi;
    "Carrier-Envelope-Phase (CEP) of the laser field."
    cep;
    "Time shift of the laser relative to the beginning of TURN-ON (in a.u.)."
    t_shift;
    """
    Constructs a new monochromatic elliptically polarized laser field with Trapezoidal-shape envelope.
    # Parameters
    - `peak_int`        : Peak intensity of the laser field (in W/cm²).
    - `wave_len`        : Wavelength of the laser field (in nm).
    - `cyc_num_turn_on` : Number of cycles of the laser field in the turn-on.
    - `cyc_num_turn_off`: Number of cycles of the laser field in the turn-off.
    - `cyc_num_const`   : Number of cycles of the laser field in the constant-intensity.
    - `ellip`           : Ellipticity of the laser field [-1≤ε≤1, 0 indicates linear polarization and ±1 indicates circular polarization].
    - `azi`             : Azimuth angle of the laser's polarization's principle axis relative to x axis (in radians) (optional, default 0).
    - `t_shift`         : Time shift of the laser (in a.u.) relative to the beginning of TURN-ON (optional, default 0).
    """
    function TrapezoidalLaser(peak_int, wave_len, cyc_num_turn_on, cyc_num_turn_off, cyc_num_const, ellip, azi=0., cep=0., t_shift=0.)
        @assert peak_int>0  "[TrapezoidalLaser] Peak intensity must be positive."
        @assert wave_len>0  "[TrapezoidalLaser] Wavelength must be positive."
        @assert cyc_num_turn_on>0 && cyc_num_turn_off>0 && cyc_num_const≥0  "[TrapezoidalLaser] Cycle number must be positive."
        @assert -1≤ellip≤1  "[TrapezoidalLaser] Ellipticity must be in [-1,1]."
        new(peak_int, wave_len, cyc_num_turn_on, cyc_num_turn_off, cyc_num_const, ellip, azi, cep, t_shift)
    end
    """
    Constructs a new monochromatic elliptically polarized laser field with Trapezoidal-shape envelope.
    # Parameters
    - `peak_int`        : Peak intensity of the laser field (in W/cm²).
    - `wave_len`        : Wavelength of the laser field (in nm). Must specify either `wave_len` or `ang_freq`.
    - `ang_freq`        : Angular frequency of the laser field (in a.u.). Must specify either `wave_len` or `ang_freq`.
    - `cyc_num_turn_on` : Number of cycles of the laser field in the turn-on.
    - `cyc_num_turn_off`: Number of cycles of the laser field in the turn-off.
    - `cyc_num_const`   : Number of cycles of the laser field in the constant-intensity.
    - `ellip`           : Ellipticity of the laser field [-1≤ε≤1, 0 indicates linear polarization and ±1 indicates circular polarization].
    - `azi`             : Azimuth angle of the laser's polarization's principle axis relative to x axis (in radians) (optional, default 0).
    - `t_shift`         : Time shift of the laser (in a.u.) relative to the beginning of TURN-ON (optional, default 0).
    """
    function TrapezoidalLaser(; peak_int,
                                wave_len=-1, ang_freq=-1,     # must specify either waveLen or angFreq.
                                cyc_num_turn_on, cyc_num_turn_off, cyc_num_const,
                                ellip, azi=0., cep=0., t_shift=0.)
        @assert wave_len>0 || ang_freq>0    "[TrapezoidalLaser] Must specify either wave_len or ang_freq."
        @assert cyc_num_turn_on>0 && cyc_num_turn_off>0 && cyc_num_const≥0  "[TrapezoidalLaser] Cycle number must be positive."
        if wave_len>0 && ang_freq>0
            @warn "[TrapezoidalLaser] Both wave_len & ang_freq are specified, will use wave_len."
        end
        if wave_len==-1
            wave_len = 45.563352525 / ang_freq
        end
        TrapezoidalLaser(peak_int, wave_len, cyc_num_turn_on, cyc_num_turn_off, cyc_num_const, ellip, azi, cep, t_shift)
    end
end
"Gets the peak intensity of the laser field (in W/cm²)."
PeakInt(l::TrapezoidalLaser) = l.peak_int
"Gets the wave length of the laser field (in nm)."
WaveLen(l::TrapezoidalLaser) = l.wave_len
"Gets the total cycle number of the laser field."
CycNumTotal(l::TrapezoidalLaser) = l.cyc_num_turn_on + l.cyc_num_turn_off + l.cyc_num_const
"Gets the cycle number of the laser field in the turn-on."
CycNumTurnOn(l::TrapezoidalLaser) = l.cyc_num_turn_on
"Gets the cycle number of the laser field in the turn-off."
CycNumTurnOff(l::TrapezoidalLaser) = l.cyc_num_turn_off
"Gets the cycle number of the laser field in the constant-intensity."
CycNumConst(l::TrapezoidalLaser) = l.cyc_num_const
"Gets the ellipticity of the laser field."
Ellipticity(l::TrapezoidalLaser) = l.ellip
"Gets the azimuth angle of the laser's polarization's principle axis relative to x axis (in radians)."
Azimuth(l::TrapezoidalLaser) = l.azi
"Gets the angular frequency (ω) of the laser field (in a.u.)."
AngFreq(l::TrapezoidalLaser) = 45.563352525 / l.wave_len
"Gets the period of the laser field (in a.u.)."
Period(l::TrapezoidalLaser) = 2π / AngFreq(l)
"Gets the Carrier-Envelope Phase (CEP) of the laser field."
CEP(l::TrapezoidalLaser) = l.cep
"Gets the time shift relative to the beginning of TURN-ON (in a.u.)."
TimeShift(l::TrapezoidalLaser) = l.t_shift
"Gets the peak electric field intensity of the laser field (in a.u.)."
LaserF0(l::TrapezoidalLaser) = sqrt(l.peak_int/(1.0+l.ellip^2)/3.50944521e16)
"Gets the peak vector potential intensity of the laser field (in a.u.)."
LaserA0(l::TrapezoidalLaser) = LaserF0(l) / AngFreq(l)

"Gets the unit envelope function (the peak value is 1) of the laser field."
function UnitEnvelope(l::TrapezoidalLaser)
    local T = Period(l); local Δt = TimeShift(l);
    local N_on = CycNumTurnOn(l); local N_const = CycNumConst(l); local N_off = CycNumTurnOff(l);
    local t_on = N_on*T; local t_const = N_const*T; local t_off = N_off*T;
    function (t)
        t -= Δt
        ( (0<t<t_on)*(t/t_on) + (t_on≤t≤(t_on+t_const))*1.0 + ((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off) )
    end
end

"Gets the time-dependent x component of the vector potential under dipole approximation."
function LaserAx(l::TrapezoidalLaser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local T = Period(l); local Δt = TimeShift(l);
    local N_on = CycNumTurnOn(l); local N_const = CycNumConst(l); local N_off = CycNumTurnOff(l);
    local t_on = N_on*T; local t_const = N_const*T; local t_off = N_off*T;
    local φ = l.cep; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)  # this function might be accessed by GPU, thus the conditional expressions should be avoided.
            t -= Δt
            A0 * ( (0<t<t_on)*(t/t_on) + (t_on≤t≤(t_on+t_const))*1.0 + ((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off) ) * cos(ω*t+φ)
        end
    else
        function(t)
            t -= Δt
            A0 * ( (0<t<t_on)*(t/t_on) + (t_on≤t≤(t_on+t_const))*1.0 + ((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off) )*(cos(ω*t+φ)*cos(ϕ) + sin(ω*t+φ)*ε*sin(ϕ))
        end
    end
end
"Gets the time-dependent y component of the vector potential under dipole approximation."
function LaserAy(l::TrapezoidalLaser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local T = Period(l); local Δt = TimeShift(l);
    local N_on = CycNumTurnOn(l); local N_const = CycNumConst(l); local N_off = CycNumTurnOff(l);
    local t_on = N_on*T; local t_const = N_const*T; local t_off = N_off*T;
    local φ = l.cep; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            A0 * ( (0<t<t_on)*(t/t_on) + (t_on≤t≤(t_on+t_const))*1.0 + ((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off) ) * sin(ω*t+φ) * ε
        end
    else
        function(t)
            t -= Δt
            A0 * ( (0<t<t_on)*(t/t_on) + (t_on≤t≤(t_on+t_const))*1.0 + ((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off) )*(cos(ω*t+φ)*-sin(ϕ) + sin(ω*t+φ)*ε*cos(ϕ))
        end
    end
end
"Gets the time-dependent x component of the electric field strength under dipole approximation."
function LaserFx(l::TrapezoidalLaser)
    local F0 = LaserF0(l); local ω = AngFreq(l); local T = Period(l); local Δt = TimeShift(l);
    local N_on = CycNumTurnOn(l); local N_const = CycNumConst(l); local N_off = CycNumTurnOff(l);
    local t_on = N_on*T; local t_const = N_const*T; local t_off = N_off*T;
    local φ = l.cep; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            F0 * ( ((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*sin(ω*t+φ)*ω - ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*cos(ω*t+φ))
        end
    else
        function(t)
            t -= Δt
            F0 * ( ((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*sin(ω*t+φ)*ω + ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*cos(ω*t+φ)*cos(ϕ) + (((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*cos(ω*t+φ)*ω + ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*sin(ω*t+φ))*-ε*sin(ϕ) )
        end
    end
end
"Gets the time-dependent y component of the electric field strength under dipole approximation."
function LaserFy(l::TrapezoidalLaser)
    local F0 = LaserF0(l); local ω = AngFreq(l); local T = Period(l); local Δt = TimeShift(l);
    local N_on = CycNumTurnOn(l); local N_const = CycNumConst(l); local N_off = CycNumTurnOff(l);
    local t_on = N_on*T; local t_const = N_const*T; local t_off = N_off*T;
    local φ = l.cep; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            -F0 * ε * ( ((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*cos(ω*t+φ)*ω + ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*sin(ω*t+φ))
        end
    else
        function(t)
            t -= Δt
            F0 * ( ((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*sin(ω*t+φ)*ω + ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*cos(ω*t+φ)*-sin(ϕ) + (((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*cos(ω*t+φ)*ω + ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*sin(ω*t+φ))*-ε*cos(ϕ) )
        end
    end
end

"Prints the information about the laser."
Base.show(io::IO, l::TrapezoidalLaser) = println(io,"[MonochromaticLaser] Envelope Trapezoidal, Wavelength=$(l.wave_len) nm, TurnOn/Constant/TurnOff: $(l.cyc_num_turn_on)/$(l.cyc_num_const)/$(l.cyc_num_turn_off) cycle(s), ε=$(l.ellip)"
                                                * (l.ellip==0 ? " [Linearly polarized]" : "") * (abs(l.ellip)==1 ? " [Circularly polarized]" : "")
                                                * ", PrincipleAxisAzimuth=$(l.azi/π*180)°" * (l.t_shift==0 ? "" : ", Rises at t₀=$(l.t_shift) a.u.") * (l.cep==0 ? "" : ", CEP=$(l.cep)"))

using Parameters, OrderedCollections
"Returns a `Dict{Symbol,Any}` containing properties of the object."
function Serialize(l::TrapezoidalLaser)
    dict = OrderedDict{Symbol,Any}()
    type                = typeof(l)
    peak_int            = l.peak_int
    wave_len            = l.wave_len
    cyc_num_turn_on     = l.cyc_num_turn_on
    cyc_num_turn_off    = l.cyc_num_turn_off
    cyc_num_turn_const  = l.cyc_num_const
    ellip               = l.ellip
    azi                 = l.azi
    cep                 = l.cep
    t_shift             = l.t_shift
    @pack! dict = (type, peak_int, wave_len, cyc_num_turn_on, cyc_num_turn_off, cyc_num_turn_const, ellip, azi, cep, t_shift)
    return dict
end
