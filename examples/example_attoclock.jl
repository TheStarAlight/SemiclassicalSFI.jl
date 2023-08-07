# === example_attoclock.jl ===
#
# With the aim of studying the influence of non-adiabatic effects on the attoclock signal, we employ ADK, SFA-AE and SFA to provide the initial conditions of the electrons, and perform trajectory simulations.
# The parameters are taken from [J. Phys. B 54, 144001 (2020)].
#
# To run the example, a julia intepreter (version ≥ 1.7) with the *SemiclassicalSFI* package (and its dependencies) installed is required.
#
# Usage: place example_attoclock.jl in the working directory, and execute:
#   julia -t <n_threads> example_attoclock.jl
# where the `n_threads` should be replaced with the desired number of threads allocated for julia.


@info "SemiclassicalSFI Example - Attoclock"

@info "Loading Packages..."
using SemiclassicalSFI
using SemiclassicalSFI.Targets
using SemiclassicalSFI.Lasers

t = HeAtom()
l = Cos4Laser(peak_int=4e14, wave_len=800, cyc_num=2, ellip=1.0)

init_cond_list = [:ADK, :SFAAE, :SFA]

for init_cond in init_cond_list
    @info "Running $(init_cond)..."
    filename = "SCSFI_He_4e14_800nm_2cyc_CP_$(init_cond)_CTMC.h5"
    performSFI(
        target = t,
        laser = l,
        init_cond_method = init_cond,
        sample_t_intv = (-80,80),
        sample_t_num = 20000,
        ss_kd_max = 1.5,
        ss_kd_num = 1000,
        ss_kz_max = 1.5,
        ss_kz_num = 50,
        traj_phase_method = :CTMC,
        traj_t_final = 120,
        final_p_max = (2,2,2),
        final_p_num = (500,500,1),
        save_path = filename
    )
end
