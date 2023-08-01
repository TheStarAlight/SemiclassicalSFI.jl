
using SemiclassicalSFI
using SemiclassicalSFI.Targets
using Test

@info "# Testing PySCFMolecularCalculator ..."

if Sys.iswindows()
    @warn "The OS is Windows, the PySCFMolecularCalculator is inoperable."
    @testset verbose=true "PySCFMolecularCalculator" begin
        @test_skip 0==0
    end
else

    @testset verbose=true "PySCFMolecularCalculator" begin
        mol = Molecule(["H","H"], [0 0 -0.375; 0 0 0.375], 0, "Hydrogen")
        #* initialization and energy calculation
        @testset verbose=true "Energy" begin
            MolCalcEnergyData!(mol, MCType=MolecularCalculators.PySCFMolecularCalculator, basis="sto-3g")
            @test MolHOMOEnergy(mol) ≈ -0.57443656
            @test MolEnergyLevels(mol)[MolHOMOIndex(mol)+1] ≈ 0.6609100513
        end
        #* WFAT structure factor calculation
        @testset verbose=true "WFAT" begin
            MolCalcWFATData!(mol, MCType=MolecularCalculators.PySCFMolecularCalculator, grid_rNum=100, grid_θNum=30, grid_ϕNum=30, sf_nξMax=1, sf_mMax=1, sf_lMax=6)
            @test MolWFATStructureFactor_G(mol,0,0,0,0.0,0.0) ≈ 2.274804156 && MolWFATStructureFactor_G(mol,0,0,0,π/2,0.0) ≈ 2.053653020 && MolWFATStructureFactor_G(mol,0,0,1,π/4,0.0) ≈ -0.211109323
        end
        #* MOADK coefficients calculation
        @testset verbose=true "MOADK" begin
            MolCalcMOADKCoeff!(mol, MCType=MolecularCalculators.PySCFMolecularCalculator, grid_rNum=50, grid_rReg=(3,8), grid_θNum=30, grid_ϕNum=30, l_max=6, m_max=3)
            @test MolMOADKCoeffs(mol)[1,1] ≈ 2.07865748426 && MolMOADKCoeffs(mol)[3,1] ≈ 0.252536756346
            @test abs2(MolMOADKStructureFactor_B(mol,0,0,0,0.0,0.0)) ≈ 3.4500046097 && abs2(MolMOADKStructureFactor_B(mol,0,0,0,π/2,0.0)) ≈ 1.642343941
        end
    end

end
