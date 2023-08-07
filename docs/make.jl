using Documenter, SemiclassicalSFI

makedocs(
        sitename = "SemiclassicalSFI.jl",
        pages = [
            "Home"      => "index.md",
            "Theory"    => [
                "Initial Conditions"                        => "theory1_initial_conditions.md",
                "Trajectory Simulation and Phase Methods"   => "theory2_trajectory_simulation_phase_methods.md"
                ],
            "Manual"    => [
                "Targets"       => "manual1_targets.md",
                "Lasers"        => "manual2_lasers.md",
                "Main Method"   => "manual3_main_method.md"
                ],
            "Examples"  => [
                "Attoclock and Initial Condition Methods"   => "example1_attoclock.md",
                "Phase Methods"                             => "example2_phase_methods.md"
            ]
        ],
        format = Documenter.HTML(
            ansicolor=true,
            edit_link=nothing,
            footer="· *SemiclassicalSFI.jl* Documentation · by *Mingyu Zhu* and other contributors"
            )
        )