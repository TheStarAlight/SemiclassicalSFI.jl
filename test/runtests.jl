using SemiclassicalSFI
using Test

@testset verbose=true "SemiclassicalSFI" begin
    include("Targets_test.jl")
    include("Lasers_test.jl")
    include("PySCFMolecularCalculator_test.jl")
    include("TrajectorySimulation_test.jl")
    include("performSFI_test.jl")
end
