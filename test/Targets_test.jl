using SemiclassicalSFI
using SemiclassicalSFI.Targets
using Test

@info "# Testing Targets ..."

@testset verbose=true "Targets" begin

    @info "Testing HydrogenLikeAtom ..."
    @testset verbose=true "HydrogenLikeAtom" begin
        t1 = HydrogenLikeAtom(Ip=0.5,Z=1.0,soft_core=0.2,name="H")
        t2 = HydrogenLikeAtom(0.5,1.0,0.2,"H")
        t3 = HAtom()
        @test t1 == t2 == t3
        @test begin
            show(t1)
            true
        end
        @test IonPotential(t1)       == 0.5
        @test AsympNuclCharge(t1)    == 1
        @test SoftCore(t1)           == 0.2
        @test TargetName(t1)         == "H"
        @test TargetPotential(t1)(1.0,1.0,1.0) ≈ -1/sqrt(3.2)
        @test reduce(*, TargetForce(t1)(1.0,1.0,1.0) .≈ (-1.0,-1.0,-1.0) .* (3.2)^(-1.5))
        #@test TrajectoryFunction(t) #TODO: add this test after revising the code.
        @test ADKRateExp(t1)(0.1,1.0,1.0) ≈ exp(-2*3.0^1.5/0.3)
    end

    @info "Testing SAEAtom ..."
    @testset verbose=true "SAEAtom" begin
        t1 = SAEAtom(Ip=0.9035698802, Z=1.0, a1= 1.230723 , b1=0.6620055, a2=-1.325040 , b2=1.236224 , a3=-0.2307230 , b3=0.4804286, name="He")
        t2 = SAEAtom(0.9035698802, 1.0, 1.230723, 0.6620055, -1.325040, 1.236224, -0.2307230, 0.4804286, "He")
        t3 = HeAtom()
        @test t1 == t2 == t3
        @test begin
            show(t1)
            true
        end
        @test IonPotential(t1)       ≈ 0.9035698802
        @test AsympNuclCharge(t1)    ≈ 1.0
        @test TargetName(t1)         == "He"
        @test TargetPotential(t1)(1.0,1.0,1.0) ≈ - (1.0 + 1.230723*exp(-0.6620055*1.73205081) + -1.325040*1.73205081*exp(-1.236224*1.73205081) + -0.2307230*exp(-0.4804286*1.73205081)) / 1.73205081
        @test reduce(*, TargetForce(t1)(1.0,1.0,1.0) .≈ (-1.0,-1.0,-1.0) .* (1.73205081^(-3) * (1.0 + 1.230723*(1+0.6620055*1.73205081)*exp(-0.6620055*1.73205081) + -0.2307230*(1+0.4804286*1.73205081)*exp(-0.4804286*1.73205081)) + -1.325040*1.236224/1.73205081 * exp(-1.236224*1.73205081)))
        #@test TrajectoryFunction(t) #TODO: add this test after revising the code.
        @test ADKRateExp(t1)(0.1,1.0,1.0) ≈ exp(-2*(2.0+2*0.9035698802)^1.5/0.3)
    end

    @info "Testing Molecule ..."
    @testset verbose=true "Molecule" begin
        # fresh init
        m1 = Molecule(atoms=["C","C","C","C","C","C","H","H","H","H","H","H"],
                    atom_coords=[0 0 1.4 ; 1.212 0 0.7 ; 1.212 0 -0.7 ; 0 0 -1.4 ; -1.212 0 -0.7 ; -1.212 0 0.7 ; 0 0 2.48; 2.148 0 1.24; 2.148 0 -1.24; 0 0 -2.48; -2.148 0 -1.24; -2.148 0 1.24],
                    charge=0, name="Benzene")
        m2 = Molecule(["C","C","C","C","C","C","H","H","H","H","H","H"],
                    [0 0 1.4 ; 1.212 0 0.7 ; 1.212 0 -0.7 ; 0 0 -1.4 ; -1.212 0 -0.7 ; -1.212 0 0.7 ; 0 0 2.48; 2.148 0 1.24; 2.148 0 -1.24; 0 0 -2.48; -2.148 0 -1.24; -2.148 0 1.24],
                    0,"Benzene")
        @test MolAtoms(m1) == MolAtoms(m2)
        @test MolAtomCoords(m1) == MolAtomCoords(m2)
        @test MolCharge(m1) == MolCharge(m2)
        @test TargetName(m1) == TargetName(m2)
        # init from existing file
        m3 = Molecule("Molecule_Hydrogen.h5")
        @test begin
            show(m3)
            true
        end
        @test MolAtoms(m3)      == ["H","H"]
        @test MolAtomCoords(m3) == [0.0 0.0 -0.375; 0.0 0.0 0.375]
        @test MolCharge(m3)     == 0
        @test MolEnergyDataAvailable(m3) == true
        @test MolEnergyLevels(m3)[1] ≈ -0.5909757276151718
        @test MolHOMOIndex(m3)  == 1
        @test MolHOMOEnergy(m3) == MolEnergyLevels(m3)[MolHOMOIndex(m3)]
        @test MolWFATAvailableIndices(m3) == Set([0])
        begin
            μ,IntData = MolWFATData(m3,0)
            @test reduce(*, abs.(μ) .<= [1e-10,1e-10,1e-10])
            @test IntData[1,1,3,1] ≈ -4.240521134419741e-5 - 2.6425097195784973e-19im
        end
        @test reduce(*, abs2.(MolWFATStructureFactor_G(m3,0,0,0,[0.0,π/2],[0.0,0.0])) .≈ [3.662226956020088, 2.761095365444768])
        @test MolMOADKAvailableIndices(m3) == Set([0])
        @test MolMOADKCoeffs(m3, 0)[1,1] ≈ 2.47486649
        @test reduce(*, abs2.(MolMOADKStructureFactor_B(m3,0,0,0,[0.0,π/2],[0.0,0.0])) .≈ [4.90607965289471, 2.36850621912187])
        begin
            @test begin
                SetMolRotation(m3,π,π/2,π/3)
                reduce(*, MolRotation(m3) .≈ (π,π/2,π/3))
            end
            @test begin
                SetMolRotation(m3,(π/2,π/3,π/4))
                reduce(*, MolRotation(m3) .≈ (π/2,π/3,π/4))
            end
        end
        @test MolExportAtomInfo(m3) == "H [0.0, 0.0, -0.375]; H [0.0, 0.0, 0.375]"
        # MolCalcEnergyData!, _MolSaveEnergyData, MolCalcWFATData!, _MolSaveWFATData, MolSaveDataAs are not included in the test
        @test IonPotential(m3)      == -MolHOMOEnergy(m3)
        @test IonPotential(m3,0)    == -MolHOMOEnergy(m3)
        @test AsympNuclCharge(m3)   == 1
        @test TargetName(m3)        == "Hydrogen"
        @test TargetPotential(m3)(1.0,1.0,1.0) == -0.5
        @test reduce(*, TargetForce(m3)(1.0,1.0,1.0) .== (-1/8,-1/8,-1/8))
        #@test TrajectoryFunction(t) #TODO: add this test after revising the code.
    end
end