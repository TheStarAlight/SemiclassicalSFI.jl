
"""
The `Lasers` module provides information about laser fields.
"""
module Lasers

# Lasers
export Laser, Cos4Laser, Cos2Laser, TrapezoidalLaser, GaussianLaser
# Properties
export PeakInt, WaveLen, Ellipticity, AngFreq, Period, CEP, Azimuth, TimeShift, UnitEnvelope   # General
export LaserF0, LaserA0, LaserAx, LaserAy, LaserFx, LaserFy     # General
export Serialize
export CycNum                                                   # Cos2, Cos4
export CycNumTurnOn, CycNumConst, CycNumTurnOff, CycNumTotal    # Trapezoidal
export SpreadCycNum, SpreadDuration, FWHM_Duration              # Gaussian

abstract type Laser end
abstract type MonochromaticLaser <: Laser end
abstract type BichromaticLaser <: Laser end

include("Cos4Laser.jl")
include("Cos2Laser.jl")
include("TrapezoidalLaser.jl")
include("GaussianLaser.jl")

end