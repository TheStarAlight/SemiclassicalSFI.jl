
"Contains interfaces relating to molecular calculation."
module MolecularCalculators

using ..Targets

export MolecularCalculatorBase
export PySCFMolecularCalculator

abstract type MolecularCalculatorBase end

include("PySCFMolecularCalculator.jl")

end
