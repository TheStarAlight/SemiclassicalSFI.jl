using SemiclassicalSFI
using SemiclassicalSFI.Lasers
using Test

@info "# Testing Lasers ..."

@testset verbose=true "Lasers" begin

    @info "Testing Cos2Laser ..."
    @testset verbose=true "Cos2Laser" begin
        l1 = Cos2Laser(peak_int=4e14, wave_len=800., cyc_num=2., ellip=1., azi=π/2, cep=π, t_shift=10.)
        l2 = Cos2Laser(peak_int=4e14, ang_freq=0.05695419065625, cyc_num=2., ellip=1., azi=π/2, cep=π, t_shift=10.)
        l3 = Cos2Laser(peak_int=4e14, wave_len=800., duration=220.63996467273427, ellip=1., azi=π/2, cep=π, t_shift=10.)
        l = Cos2Laser(4e14,800.,2.,1.,π/2,π,10.)
        @test l == l1
        @test l == l2
        @test l == l3
        @test begin
            show(l1)
            true
        end
        @test PeakInt(l1)       == 4e14
        @test WaveLen(l1)       == 800.
        @test CycNum(l1)        == 2.
        @test Ellipticity(l1)   == 1.
        @test Azimuth(l1)       == π/2
        @test AngFreq(l1)       == 45.563352525 / WaveLen(l1)
        @test Period(l1)        == 2π / AngFreq(l1)
        @test CEP(l1)           == π
        @test TimeShift(l1)     == 10.0
        @test LaserF0(l1)       == sqrt(PeakInt(l1)/(1.0+Ellipticity(l1)^2)/3.50944521e16)
        @test LaserA0(l1)       == LaserF0(l1) / AngFreq(l1)
    end

    @info "Testing Cos4Laser ..."
    @testset verbose=true "Cos4Laser" begin
        l1 = Cos4Laser(peak_int=4e14, wave_len=800., cyc_num=2., ellip=1., azi=π/2, cep=π, t_shift=10.)
        l2 = Cos4Laser(peak_int=4e14, ang_freq=0.05695419065625, cyc_num=2., ellip=1., azi=π/2, cep=π, t_shift=10.)
        l3 = Cos4Laser(peak_int=4e14, wave_len=800., duration=220.63996467273427, ellip=1., azi=π/2, cep=π, t_shift=10.)
        l = Cos4Laser(4e14,800.,2.,1.,π/2,π,10.)
        @test l == l1
        @test l == l2
        @test l == l3
        @test begin
            show(l1)
            true
        end
        @test PeakInt(l1)       == 4e14
        @test WaveLen(l1)       == 800.
        @test CycNum(l1)        == 2.
        @test Ellipticity(l1)   == 1.
        @test Azimuth(l1)       == π/2
        @test AngFreq(l1)       == 45.563352525 / WaveLen(l1)
        @test Period(l1)        == 2π / AngFreq(l1)
        @test CEP(l1)           == π
        @test TimeShift(l1)     == 10.0
        @test LaserF0(l1)       == sqrt(PeakInt(l1)/(1.0+Ellipticity(l1)^2)/3.50944521e16)
        @test LaserA0(l1)       == LaserF0(l1) / AngFreq(l1)
    end

    @info "Testing GaussianLaser ..."
    @testset verbose=true "GaussianLaser" begin
        l1 = GaussianLaser(peak_int=4e14, wave_len=800., spread_cyc_num=2., ellip=1., azi=π/2, cep=π, t_shift=10.)
        l2 = GaussianLaser(peak_int=4e14, ang_freq=0.05695419065625, spread_cyc_num=2., ellip=1., azi=π/2, cep=π, t_shift=10.)
        l3 = GaussianLaser(peak_int=4e14, wave_len=800., spread_duration=220.63996467273427, ellip=1., azi=π/2, cep=π, t_shift=10.)
        l4 = GaussianLaser(peak_int=4e14, wave_len=800., FWHM_duration=519.5674115462751, ellip=1., azi=π/2, cep=π, t_shift=10.)
        l = GaussianLaser(4e14, 800., 2., 1., π/2, π, 10.)
        @test l == l1
        @test l == l2
        @test l == l3
        @test l == l4
        @test begin
            show(l1)
            true
        end
        @test PeakInt(l1)       == 4e14
        @test WaveLen(l1)       == 800.
        @test SpreadCycNum(l1)  == 2.
        @test SpreadDuration(l1)== SpreadCycNum(l1) * Period(l1)
        @test FWHM_Duration(l1) == SpreadDuration(l1) * (2*sqrt(2*log(2)))
        @test Ellipticity(l1)   == 1.
        @test Azimuth(l1)       == π/2
        @test AngFreq(l1)       == 45.563352525 / WaveLen(l1)
        @test Period(l1)        == 2π / AngFreq(l1)
        @test CEP(l1)           == π
        @test TimeShift(l1)     == 10.0
        @test LaserF0(l1)       == sqrt(PeakInt(l1)/(1.0+Ellipticity(l1)^2)/3.50944521e16)
        @test LaserA0(l1)       == LaserF0(l1) / AngFreq(l1)
    end

    @info "Testing TrapezoidalLaser ..."
    @testset verbose=true "TrapezoidalLaser" begin
        l1 = TrapezoidalLaser(peak_int=4e14, wave_len=800., cyc_num_turn_on=2., cyc_num_turn_off=2., cyc_num_const=6., ellip=1., azi=π/2, cep=π, t_shift=10.)
        l2  = TrapezoidalLaser(peak_int=4e14, ang_freq=0.05695419065625, cyc_num_turn_on=2., cyc_num_turn_off=2., cyc_num_const=6., ellip=1., azi=π/2, cep=π, t_shift=10.)
        l = TrapezoidalLaser(4e14, 800., 2., 2., 6., 1., π/2, π, 10.)
        @test l == l1
        @test l == l2
        @test begin
            show(l1)
            true
        end
        @test PeakInt(l1)       == 4e14
        @test WaveLen(l1)       == 800.
        @test CycNumTurnOn(l1)  == 2.
        @test CycNumTurnOff(l1) == 2.
        @test CycNumConst(l1)   == 6.
        @test CycNumTotal(l1)   == 10.
        @test Ellipticity(l1)   == 1.
        @test Azimuth(l1)       == π/2
        @test AngFreq(l1)       == 45.563352525 / WaveLen(l1)
        @test Period(l1)        == 2π / AngFreq(l1)
        @test CEP(l1)           == π
        @test TimeShift(l1)     == 10.0
        @test LaserF0(l1)       == sqrt(PeakInt(l1)/(1.0+Ellipticity(l1)^2)/3.50944521e16)
        @test LaserA0(l1)       == LaserF0(l1) / AngFreq(l1)
    end

end