using SemiclassicalSFI
using SemiclassicalSFI.Targets
using SemiclassicalSFI.Lasers
using Test
using Base.Threads

@info "# Testing main method performSFI ..."

@testset verbose=true "performSFI" begin

    if Threads.nthreads() == 1
        @warn "The process is running with single available thread, to enable multi-threading, start julia with parameter `-t N` to use N threads."
    end

    rm("./performSFI_test_output/", recursive=true, force=true)
    mkdir("./performSFI_test_output/")

    @info "Testing ADK ..."
    @testset "ADK" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,2.0),
                final_p_num         = (200,200,1),
                ss_kd_max           = 2.0,
                ss_kd_num           = 100,
                ss_kz_max           = 2.0,
                ss_kz_num           = 100,
                save_path           = "./performSFI_test_output/test_ADK_4e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing ADK (GPU) ..."
    @testset "ADK (GPU)" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,2.0),
                final_p_num         = (200,200,1),
                ss_kd_max           = 2.0,
                ss_kd_num           = 100,
                ss_kz_max           = 2.0,
                ss_kz_num           = 100,
                save_path           = "./performSFI_test_output/test_ADK_GPU_4e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full,
                traj_GPU            = true
            )
            true
        end
    end

    @info "Testing ADK (Monte-Carlo) ..."
    @testset "ADK (Monte-Carlo)" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,2.0),
                final_p_num         = (200,200,1),
                mc_kp_num           = 10000,
                mc_kp_max           = 2.0,
                save_path           = "./performSFI_test_output/test_ADK_MC_4e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full,
                sample_monte_carlo  = true
            )
            true
        end
    end

    @info "Testing ADK (QTMC) ..."
    @testset "ADK (QTMC)" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,2.0),
                final_p_num         = (200,200,1),
                ss_kd_max           = 2.0,
                ss_kd_num           = 100,
                ss_kz_max           = 2.0,
                ss_kz_num           = 100,
                save_path           = "./performSFI_test_output/test_ADK_QTMC_4e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :QTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing ADK (SCTS) ..."
    @testset "ADK (SCTS)" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,2.0),
                final_p_num         = (200,200,1),
                ss_kd_max           = 2.0,
                ss_kd_num           = 100,
                ss_kz_max           = 2.0,
                ss_kz_num           = 100,
                save_path           = "./performSFI_test_output/test_ADK_SCTS_4e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :SCTS,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing SFA-AE ..."
    @testset "SFA-AE" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :SFAAE,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,2.0),
                final_p_num         = (200,200,1),
                ss_kd_max           = 2.0,
                ss_kd_num           = 100,
                ss_kz_max           = 2.0,
                ss_kz_num           = 100,
                save_path           = "./performSFI_test_output/test_SFAAE_4e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing SFA ..."
    @testset "SFA" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :SFA,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,2.0),
                final_p_num         = (200,200,1),
                ss_kd_max           = 2.0,
                ss_kd_num           = 100,
                ss_kz_max           = 2.0,
                ss_kz_num           = 100,
                save_path           = "./performSFI_test_output/test_SFA_4e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MO-ADK ..."
    @testset "MO-ADK" begin
        t = Molecule("Molecule_Hydrogen.h5")
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :MOADK,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,2.0),
                final_p_num         = (200,200,1),
                ss_kd_max           = 2.0,
                ss_kd_num           = 100,
                ss_kz_max           = 2.0,
                ss_kz_num           = 100,
                save_path           = "./performSFI_test_output/test_MOADK_Hydrogen_4e14_800nm_2cyc_CP.h5",
                traj_rtol           = 1e-6
            )
            true
        end
    end

    @info "Testing WFAT ..."
    @testset "WFAT" begin
        t = Molecule("Molecule_Hydrogen.h5")
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :WFAT,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,2.0),
                final_p_num         = (200,200,1),
                ss_kd_max           = 2.0,
                ss_kd_num           = 100,
                ss_kz_max           = 2.0,
                ss_kz_num           = 100,
                save_path           = "./performSFI_test_output/test_WFAT_Hydrogen_4e14_800nm_2cyc_CP.h5",
                traj_rtol           = 1e-6
            )
            true
        end
    end

end