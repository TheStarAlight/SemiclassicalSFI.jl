using ..Targets
using SpecialFunctions
using SphericalHarmonics
using Folds
using Einsum
using LsqFit
using PyCall
using libcint_jll
using HDF5
using Dates
using ProgressMeter

"""
```
struct PySCFMolecularCalculator <: MolecularCalculatorBase
```
An interface of molecular calculation using PySCF.

An instance of `PySCFMolecularCalculator` can be initialized via the following constructor method:

```
PySCFMolecularCalculator(; mol, basis::String="pc-1", kwargs...)
```

# Parameters
- `mol::Molecule`   : The `Molecule` to be calculated.
- `basis::String`   : Basis set used for calculation (default `pc-1`).
"""
mutable struct PySCFMolecularCalculator <: MolecularCalculatorBase
    "The molecule to be calculated."
    # mol::Molecule; # linter reported error when specifying type Molecule.
    mol;
    "The basis function used for calculation."
    basis::String;

    "PySCF library."
    _pyscf;
    "The PySCF molecule object."
    _pymol;
    "The PySCF computation task object."
    _pytask;

    "HOMO energy."
    HOMO_energy;
    "Total dipole momentum vector in the MF."
    dip_momentum;

    """
    Initializes an instance of `PySCFMolecularCalculator` with given parameter.

    # Parameters
    - `mol::Molecule`   : The molecule to be calculated.
    - `basis::String`   : Basis set used for calculation (default `pc-1`).
    """
    function PySCFMolecularCalculator(; mol, basis::String="pc-1", kwargs...)
        mc::PySCFMolecularCalculator = new(mol,basis)
        @info "[PySCFMolecularCalculator] Running molecular calculation..."
        time = [0.0]
        try
            mc._pyscf  = pyimport("pyscf")
            mc._pymol  = mc._pyscf.gto.M(atom=MolExportAtomInfo(mol), charge=MolCharge(mol), basis=basis, verbose=0)
            mc._pytask = mc._pyscf.scf.RHF(mc._pymol)
            mc._pytask.chkfile = nothing
            time[1] = @elapsed mc._pytask.run()
        catch
            @error "[PySCFMolecularCalculator] Encountered error when calling pyscf."
            rethrow()
        end
        if ! mc._pytask.converged
            error("[PySCFMolecularCalculator] SCF calculation unable to converge.")
        end
        mc.HOMO_energy = mc._pytask.mo_energy[HOMOIndex(mc)]
        mc.dip_momentum = mc._pytask.dip_moment(unit="AU", verbose=0)
        @info "Finished initialization [taking $(time[1]) second(s)]."
        return mc
    end
end

"Gets the index of the HOMO."
function HOMOIndex(mc::PySCFMolecularCalculator)
    task     = mc._pytask
    mo_occ   = Int.(task.mo_occ)
    num_AO   = size(task.mo_coeff,1)
    idx = 1
    for i in 1:num_AO
        if mo_occ[i]==0
            idx = i; break
        end
    end
    return idx-1    # idx is LUMO, while idx-1 is HOMO
end

"Gets the energy level of the molecule's HOMO (in a.u.)."
function HOMOEnergy(mc::PySCFMolecularCalculator)
    return mc.HOMO_energy
end

"""
Gets the energy level of a specific molecular orbital (MO) of the molecule (in a.u.).
- `orbitIdx::Int`   : Index of the orbital.
"""
function EnergyLevel(mc::PySCFMolecularCalculator, orbitIdx::Int)
    return mc._pytask.mo_energy[orbitIdx]
end

"Gets the energy levels of all the molecular orbitals (MO) of the molecule (in a.u.)."
function EnergyLevels(mc::PySCFMolecularCalculator)
    return mc._pytask.mo_energy
end

"Gets the permanent dipole momentum vector of the molecule in the molecular frame (MF) (in a.u.)."
function DipoleMomentum(mc::PySCFMolecularCalculator)
    return mc.dip_momentum
end

"""
```
calc_WFAT_data(;
    mc::PySCFMolecularCalculator,
    orbitIdx_relHOMO::Int = 0,
    grid_rNum::Int  = 200,
    grid_rMax::Real = 10.,
    grid_θNum::Int  = 60,
    grid_ϕNum::Int  = 60,
    sf_nξMax ::Int  = 5,
    sf_mMax  ::Int  = 5,
    sf_lMax  ::Int  = 10,
    kwargs...)
```

Calculates the DATA used in structure factor calculation in WFAT of the given molecule.

# Returns
`(μ, int_data)` : Orbital dipole momentum and the array which stores the integrals.

# Parameters
- `mc`              : The molecular calculator.
- `orbitIdx_relHOMO`: Index of selected orbit relative to the HOMO (e.g., 0 indicates HOMO and -1 indicates HOMO-1) (default 0).
- `grid_rNum`       : The number of radial grid (default 200).
- `grid_rMax`       : The maximum radius of the radial grid (default 10.0).
- `grid_θNum`       : The number of angular grid in the ``θ`` direction (default 60).
- `grid_ϕNum`       : The number of angular grid in the ``ϕ`` direction (default 60).
- `sf_nξMax`        : The maximum number of ``n_ξ`` used in calculation (default 5).
- `sf_mMax`         : The maximum number of ``|m|`` used in calculation (default 5).
- `sf_lMax`         : The maximum angular quantum number ``l`` used in calculation (default 10).
"""
function calc_WFAT_data(;
                        mc::PySCFMolecularCalculator,
                        orbitIdx_relHOMO::Int = 0,
                        grid_rNum::Int  = 200,
                        grid_rMax::Real = 10.,
                        grid_θNum::Int  = 60,
                        grid_ϕNum::Int  = 60,
                        sf_nξMax ::Int  = 5,
                        sf_mMax  ::Int  = 5,
                        sf_lMax  ::Int  = 10,
                        kwargs...)
    # == PROCEDURE ==
    # 0. Obtain the coefficients (finished in the initialization).
    # 1. Calculate the effective core potential.
    # 2. Calculate the integrals and save them as output.

    @info "[PySCFMolecularCalculator] Running calculation of WFAT structure factor data... (ionizing orbital $orbitIdx_relHOMO relative to HOMO)"

    #* Preprocess molecular information

    pyscf = mc._pyscf
    pyscf_df = pyimport("pyscf.df")
    pymol = mc._pymol   # storing the molecule's info and the basis's info.
    task  = mc._pytask  # storing the calculation result.

    mo_coeff    = task.mo_coeff     # Linear combination coefficients of AO to make up MO.
    num_atom    = pymol.natm        # Total number of atoms.
    num_AO      = size(mo_coeff,1)  # Number of atomic orbits (aka AO) or Gaussian basis. (pymol.nao doesn't return an interger!)
    HOMO_orbitIdx = HOMOIndex(mc)
    orbitIdx = HOMO_orbitIdx + orbitIdx_relHOMO
    @assert 0<orbitIdx<num_AO "[PySCFMolecularCalculator] Orbit index $(orbitIdx_relHOMO) out of range."

    Z   = MolCharge(mc.mol)+1
    Ip  = -task.mo_energy[orbitIdx]
    if Ip ≤ 0
        error("[PySCFMolecularCalculator] The energy of the selected molecular orbital is positive.")
    end
    κ   = sqrt(2Ip)

    #* Define the spherical grids.
    grid_rMin = 0.001
    grid_dr = (grid_rMax-grid_rMin)/(grid_rNum-1)
    r_grid = range(start=grid_rMin, stop=grid_rMax, length=grid_rNum)
    θ_grid = range(start=0., stop=π, length=grid_θNum)
    ϕ_grid = range(start=0., stop=2π, length=grid_ϕNum)
    N = grid_rNum*grid_θNum*grid_ϕNum

    "Returns the spherical grid indices of the given index of the point."
    function ptIdx2sphCoordIdx(i::Int)
        iϕ = (i-1) % grid_ϕNum + 1
        iθ = (ceil(Int, i/grid_ϕNum)-1) % grid_θNum + 1
        ir = ceil(Int, i/grid_ϕNum/grid_θNum)
        return (ir,iθ,iϕ)
    end

    # converting to cartesian grid points.
    # x = r sinθ cosϕ, y = r sinθ sinϕ, z = r cosθ.
    pt_x = zeros(N)
    pt_y = zeros(N)
    pt_z = zeros(N)
    dV   = zeros(N) # dV = r²sinθ drdθdϕ, used in the integration.
    dr   = r_grid[2]-r_grid[1]
    dθ   = θ_grid[2]-θ_grid[1]
    dϕ   = ϕ_grid[2]-ϕ_grid[1]
    Threads.@threads for i in 1:N
        ir,iθ,iϕ = ptIdx2sphCoordIdx(i)
        r,θ,ϕ = r_grid[ir],θ_grid[iθ],ϕ_grid[iϕ]
        pt_x[i] = r*sin(θ)*cos(ϕ)
        pt_y[i] = r*sin(θ)*sin(ϕ)
        pt_z[i] = r*cos(θ)
        dV[i]   = r^2*sin(θ)*dr*dθ*dϕ
    end
    pt_xyz = hcat(pt_x,pt_y,pt_z)

    #* 1. Calculate the effective core potential
    # Vc is composed of:
    # 1. Asymptotic Coulomb potential Z/r
    # 2. Nuclear potential Vnuc
    # 3. Inter-electron interaction: direct part Vd
    # 4. Inter-electron interaction: exchange part Vex
    # The direct expression of Vex couldn't be obtained, but Vex*ψ0 is obtainable.
    # In the following scheme, Z/r, Vnuc, Vd would be first calculated and directly added to Vc_ψ0,
    # and then would be multiplied by ψ0.
    # Vc_ψ0 = (Z/r + Vnuc + Vd) * ψ0 + Vex_ψ0.

    Vc_ψ0 = zeros(N)

    #*  1.1 Calculate the wavefunction ψ0 & dip_moment μ
    # wavefunction ψ0
    χi = pymol.eval_gto("GTOval",pt_xyz)        # Size: N×Num_AO. Wavefunction of all AOs by calling eval_gto.
    orbit_coeff = @view mo_coeff[:,orbitIdx]    # Select the coefficients related to the interested MO.
    ψ0 = χi * orbit_coeff                       # Calculate wavefunction of the interested MO (matmul operation).

    # selected orbit's dipole momentum μ (shouldn't be mixed with the total dip. moment. D)
    dip_int = -1 .* pymol.intor_symmetric("int1e_r")  # I[i,α,β] = ∫(ψα*)(ri)(ψβ)dr, where r1=x, r2=y, r3=z.
    μ = zeros(3); i,α,β=0,0,0
    @einsum μ[i] := orbit_coeff[α] * orbit_coeff[β] * dip_int[i,α,β]

    #*  1.2 Calculate the asymptotic Coulomb potential Z/r
    Threads.@threads for ir in 1:grid_rNum
        Vc_ψ0[(ir-1)*grid_ϕNum*grid_θNum+1:ir*grid_ϕNum*grid_θNum] .+= Z / r_grid[ir]
    end

    #*  1.3 Calculate the nuclear potential Vnuc
    atomCharges = pymol.atom_charges()
    atomCoords  = pymol.atom_coords()   # Size: Num_Atoms×3
    for iatm in 1:num_atom
        Folds.map(
            function (i)
                Vc_ψ0[i] -= atomCharges[iatm]/sqrt((pt_xyz[i,1]-atomCoords[iatm,1])^2+(pt_xyz[i,2]-atomCoords[iatm,2])^2+(pt_xyz[i,3]-atomCoords[iatm,3])^2+1e-5)
            end
        , 1:N)
    end

    #*  1.4 Calculate the inter-electron interaction: Vd & Vex
    batch_size = 200  # the integral takes huge memory and thus needs to be performed in batches.
    batch_num = ceil(Int, N/batch_size)

    # density matrix
    den_mat = zeros(num_AO,num_AO); i,α,β=0,0,0 # den_mat_αβ = ∑_i c[i,α]*c[i,β]
    occupied_mo_coeff = mo_coeff[:,1:HOMO_orbitIdx]
    @einsum den_mat[α,β] := 2 * occupied_mo_coeff[α,i] * occupied_mo_coeff[β,i]

    prog11 = ProgressUnknown(dt=0.2, desc="Calculating the effective potential... ($N pts)", color = :cyan, spinner = true)
    prog12 = Progress(N; dt=0.2, color = :cyan, barlen = 25, barglyphs = BarGlyphs('[', '●', ['◔', '◑', '◕'], '○', ']'), showspeed = true, offset=1)

    # libcint param. preparation

    m_atm  = pymol._atm
    m_nbas = pymol.nbas
    m_bas = pymol._bas
    m_env = pymol._env
    ao_loc = pymol.ao_loc   # example: [0,2,5,7,10] indicates 4 AOs, whose GTOs are of indices 0~1,2~4,5~6,7~9 (indices starting from 0).

    """
    (translated from pyscf) Constructs the data of 'fakemol' (which is actually a delta function) used in `libcint` integration.
    - Returns: `(atm,bas,env)` of the constructed fakemol.
    """
    function construct_fakemol_data(coords)
        PTR_COORD  = 1
        ATM_SLOTS  = 6
        ATOM_OF    = 0
        NPRIM_OF   = 2
        NCTR_OF    = 3
        PTR_EXP    = 5
        PTR_COEFF  = 6
        BAS_SLOTS  = 8
        PTR_ENV_START = 20
        f_nbas = size(coords,1)
        f_atm = zeros(Int32, f_nbas, ATM_SLOTS)
        f_bas = zeros(Int32, f_nbas, BAS_SLOTS)
        f_env = vcat(zeros(PTR_ENV_START), reshape(transpose(coords),:), [1e16, 6.366197723675814e23])
        pos = PTR_ENV_START
        #* Note: in julia the indices begin at ONE!
        f_atm[:,PTR_COORD+1] = pos:3:(pos+3*f_nbas-3)
        pos += 3*f_nbas
        f_bas[:,ATOM_OF+1  ]  = 0:(f_nbas-1)
        f_bas[:,NPRIM_OF+1 ] .= 1
        f_bas[:,NCTR_OF+1  ] .= 1
        f_bas[:,PTR_EXP+1  ] .= pos
        f_bas[:,PTR_COEFF+1] .= pos+1
        return f_atm,f_bas,f_env
    end
    """
    Concatenates the two molecules' integral parameters.
    - Returns: `(atm,bas,env)` of the concatenated param.
    """
    function conc_mol_param(atm1,bas1,env1,atm2,bas2,env2)
        PTR_COORD  = 1
        PTR_ZETA   = 3
        ATOM_OF    = 0
        PTR_EXP    = 5
        PTR_COEFF  = 6
        offset = size(env1,1)
        natm_offset = size(atm1,1)
        atm2_ = copy(atm2)
        bas2_ = copy(bas2)
        atm2_[:,PTR_COORD+1] .+= offset
        atm2_[:,PTR_ZETA+1 ] .+= offset
        bas2_[:,ATOM_OF+1  ] .+= natm_offset
        bas2_[:,PTR_EXP+1  ] .+= offset
        bas2_[:,PTR_COEFF+1] .+= offset
        return Int32.(vcat(atm1,atm2_)), Int32.(vcat(bas1,bas2_)), vcat(env1,env2)
    end
    "Wraps the `int3c2e_sph` function in the `libcint`. Def: int3c2e_sph(double *out, int *dims, int *shls, int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt, double *cache)."
    function int3c2e_sph!(out,dims,shls,atm,natm,bas,nbas,env,opt=C_NULL,cache=C_NULL)
        @ccall libcint.int3c2e_sph( out  :: Ptr{Cdouble},
                                    dims :: Ptr{Cint},
                                    shls :: Ptr{Cint},
                                    atm  :: Ptr{Cint},
                                    natm :: Cint,
                                    bas  :: Ptr{Cint},
                                    nbas :: Cint,
                                    env  :: Ptr{Cdouble},
                                    opt  :: Ptr{Cvoid},
                                    cache:: Ptr{Cvoid}
                                    )::Cvoid
    end
    #TODO: add support for cintopt.
    Threads.@threads for i in 1:batch_num
        pt_idx = if i < batch_num
            CartesianIndices((((i-1)*batch_size+1): i*batch_size,))
        else
            CartesianIndices((((i-1)*batch_size+1): N,))
        end
        batch_size_this = size(pt_idx,1)
        I = zeros(num_AO,num_AO,batch_size_this)   # Size: num_AO × num_AO × num_pts, I_αβ = ∫dr' (χα(r'-Rα)*χβ(r'-Rβ))/|r-r'|
        begin
            # here we denote α,β as indices of AOs (in julia, starting from 1);
            # i,j as shells (used in libcint API, starting from 0);
            # ipt as indices of the fakemol's delta function.
            f_atm,f_bas,f_env = construct_fakemol_data(pt_xyz[pt_idx,:])
            atm,bas,env = conc_mol_param(m_atm,m_bas,m_env,f_atm,f_bas,f_env)
            atm = collect(atm)
            bas = collect(bas)
            env = collect(env)
            natm = size(atm,1)
            nbas = size(bas,1)
            atm = permutedims(atm,(2,1))    # julia stores arrays in column order, calling C requires arrays in row order.
            bas = permutedims(bas,(2,1))
            for ipt in 1:batch_size_this
                for i in 0:m_nbas-1, j in 0:i     # I[i,j]=I[j,i], thus only half of the matrix needs to be calculated.
                    k = ipt-1
                    shls = Int32[i,j,m_nbas+k]      # denoting the shells to integrate (e.g. H1's 1s and H2's 2px)
                    di = ao_loc[i+2] - ao_loc[i+1]  # each shell i contains di AOs, these AOs' integrations would be evaluated in a single call.
                    dj = ao_loc[j+2] - ao_loc[j+1]
                    dims = Int32[di,dj,1]           # denoting the count of AOs (in each dimension) in this evaluation.
                    out = zeros(di,dj)
                    int3c2e_sph!(out,dims,shls,atm,natm,bas,nbas,env)
                    out = permutedims(out,(2,1))
                    I[ao_loc[j+1]+1:ao_loc[j+2], ao_loc[i+1]+1:ao_loc[i+2], ipt] = out
                end
            end
            if num_AO != 1
                # filling the remaining part according to the symmetry.
                for α in 2:num_AO, β in 1:α-1
                    I[α,β,:] = I[β,α,:]
                end
            end
        end
        #* Calculate Vd
        Vd = zeros(batch_size_this); α,β,i_pt = 0,0,0     # α,β are indices of the basis, i_pt is the index of the points.
        @einsimd Vd[i_pt] := I[α,β,i_pt] * den_mat[α,β]     # Einstein's summation notation is used.
        Vc_ψ0[pt_idx] .+= Vd            # now it is (Z/r + Vnuc + Vd).
        Vc_ψ0[pt_idx] .*= ψ0[pt_idx]    # now it is (Z/r + Vnuc + Vd) * ψ0.
        #* Calculate Vex_ψ0
        χα = χi[pt_idx,:]   # basis function values on the selected points.
        Vex_ψ0 = zeros(batch_size_this); α,β,γ,i_pt = 0,0,0,0
        @einsimd Vex_ψ0[i_pt] := -1/2 * den_mat[α,β] * χα[i_pt,α] * orbit_coeff[γ] * I[β,γ,i_pt]
        Vc_ψ0[pt_idx] .+= Vex_ψ0    # finished building Vc_ψ0.

        next!(prog11,spinner=raw"-\|/"); next!(prog12,step=size(pt_idx,1));
    end
    finish!(prog11); finish!(prog12); println()

    #* 2. Calculate the integral
    #*  2.1 Define some special functions
    "Kummer's confluent hypergeometric function M(a,b,z) = ₁F₁(a;b;z)."
    function M(a,b,z)
        S₀, S₁, j = 1, 1+a*z/b, 1
        while abs(S₀-S₁) > 1e-10 || j ≤ 1
            rⱼ = (a+j)/((b+j)*(j+1))
            S₀, S₁ = S₁, S₁+(S₁-S₀)*rⱼ*z
            j += 1
        end
        return S₁
    end
    "Factorial."
    function fact(N::Int)
        F = 1; n = 2
        while n≤N
            F*=n; n+=1
        end
        return F
    end
    "Normalization coefficient ω_l^ν for radial function R_l^ν of Ω_{lm'}^ν."
    function ω(nξ,m,l,Z,κ)
        F1 = ((-1)^(l+(abs(m)-m)/2+1)) * (2^(l+3/2)) * (κ^(Z/κ-(abs(m)+1)/2-nξ))
        F2 = sqrt(1.0*(2l+1)*fact(l+m)*fact(l-m)*fact(abs(m)+nξ)*fact(nξ)) * fact(l)/fact(2l+1)     # 1.0 to avoid overflow
        F3 = 0  # Factor3 is a sum over k from 0 to min(n,l-|m|).
        for k in 0:min(nξ,l-abs(m))
            F3 += gamma(l+1-Z/κ+nξ-k) / (fact(k)*fact(l-k)*fact(abs(m)+k)*fact(l-abs(m)-k)*fact(nξ-k))
        end
        return F1*F2*F3
    end
    "Radial function R_l^ν of Ω_{lm'}^ν EXCLUDING the normalization constant ω."
    function R_(l,Z,κ,r)
        return (κ*r)^l*exp(-κ*r)*M(l+1-Z/κ,2l+2,2κ*r)
    end

    #*  2.2 Create pre-computation data to accelerate.
    R_precomp_data = zeros(sf_nξMax+1, 2*sf_mMax+1, sf_lMax+1, grid_rNum)
    # obtain the R_{nξ,m,l}(r) by indexing [nξ+1,m+mMax+1,l+1,r_idx].
    Threads.@threads for l in 0:sf_lMax
        R_precomp_data[1,1,l+1,:] = map(r->R_(l,Z,κ,r), r_grid)
        for nξ in 0:sf_nξMax
        for m  in -sf_mMax:sf_mMax
            if nξ==0 && m==-sf_mMax
                continue
            end
            R_precomp_data[nξ+1,m+sf_mMax+1,l+1,:] = R_precomp_data[1,1,l+1,:] .* ω(nξ,m,l,Z,κ)
        end
        end
        R_precomp_data[1,1,l+1,:] .*= ω(0,0,l,Z,κ)
    end
    Y_precomp_data = zeros(ComplexF64, sf_lMax+1, 2*sf_lMax+1, grid_θNum, grid_ϕNum)  # for given l, -l ≤ m' ≤ l.
    # obtain the Y_lm'(θ,ϕ) by indexing [l+1,m'+l+1,θ_idx,ϕ_idx].
    Threads.@threads for l in 0:sf_lMax
    for m_ in -l:l
        for iθ in eachindex(θ_grid)
        for iϕ in eachindex(ϕ_grid)
            Y_precomp_data[l+1,m_+l+1,iθ,iϕ] = SphericalHarmonics.sphericalharmonic(θ_grid[iθ], ϕ_grid[iϕ]; l=l, m=m_) / (-1)^m_
        end; end
    end; end
    "Utilizes the pre-computed data to calculate the Ω_{lm'}^{ν}, where ν=(nξ,m)."
    function Ω_precomp(nξ,m,l,m_, i_pt)
        ir,iθ,iϕ = ptIdx2sphCoordIdx(i_pt)
        abs(m_)>l && return 0.0
        return R_precomp_data[nξ+1,m+sf_mMax+1,l+1,ir] * Y_precomp_data[l+1,m_+l+1,iθ,iϕ]
    end

    #*  2.3 Calculate the integral

    prog21 = ProgressUnknown(dt=0.2, desc="Calculating the integrals... ($((sf_nξMax+1)*(2*sf_mMax+1)*(sf_lMax+1)^2) integrals)", color = :cyan, spinner = true)
    prog22 = Progress((sf_nξMax+1)*(2*sf_mMax+1)*(sf_lMax+1)^2; dt=0.2, color = :cyan, barlen = 25, barglyphs = BarGlyphs('[', '●', ['◔', '◑', '◕'], '○', ']'), showspeed = true, offset=1)

    # `int_data` would store the final data: The integral I(nξ,m,l,m')=∫Ω(nξ,m,l,m')*Vc_ψ0(r)*dV.
    # nξ=0,1,⋯,nξMax;  m=0,±1,⋯,±mMax;  l=0,1,⋯,lMax;  m'=-l,-l+1,⋯,0,1,⋯,l.
    # Obtain I(nξ,m,l,m') by indexing [nξ+1, m+mMax+1, l+1, m'+l+1]
    int_data = zeros(ComplexF64, sf_nξMax+1, 2*sf_mMax+1, sf_lMax+1, 2*sf_lMax+1)
    Vc_ψ0_dV = Vc_ψ0 .* dV
    for nξ in 0:sf_nξMax
    for m in -sf_mMax:sf_mMax
    for l in 0:sf_lMax
    for m_ in -l:l
        int_data[nξ+1, m+sf_mMax+1, l+1, m_+l+1] = Folds.mapreduce(i->conj(Ω_precomp(nξ,m,l,m_,i))*Vc_ψ0_dV[i], +, 1:N)
        next!(prog21,spinner=raw"-\|/"); next!(prog22)
    end; end; end; end
    finish!(prog21); finish!(prog22); println()

    return μ, int_data
end

"""
```
calc_asymp_coeff(;
    mc::PySCFMolecularCalculator,
    orbitIdx_relHOMO::Int = 0,
    grid_rNum::Int  = 200,
    grid_rReg::Tuple{<:Real,<:Real} = (3,8),
    grid_θNum::Int  = 60,
    grid_ϕNum::Int  = 60,
    l_max::Int      = 6,
    kwargs...)
```

Calculates the asymptotic coefficients (used in MO-ADK and MO-SFA) of the given molecule.

# Parameters
- `mc`              : The molecular calculator.
- `orbitIdx_relHOMO`: Index of selected orbit relative to the HOMO (e.g., 0 indicates HOMO and -1 indicates HOMO-1) (default 0).
- `grid_rNum`       : The number of radial grid (default 200).
- `grid_rReg`       : The region of radial distance to fit the wavefunction to obtain the coefficients (default (3,8)).
- `grid_θNum`       : The number of angular grid in the θ direction (default 60).
- `grid_ϕNum`       : The number of angular grid in the ϕ direction (default 60).
- `l_max`           : The maximum number of l calculated (default 6).
"""
function calc_asymp_coeff(;
                            mc::PySCFMolecularCalculator,
                            orbitIdx_relHOMO::Int = 0,
                            grid_rNum::Int  = 200,
                            grid_rReg::Tuple{<:Real,<:Real} = (3,8),
                            grid_θNum::Int  = 60,
                            grid_ϕNum::Int  = 60,
                            l_max::Int      = 6,
                            kwargs...)

    @info "[PySCFMolecularCalculator] Running calculation of asymptotic coefficients... (ionizing orbital $orbitIdx_relHOMO relative to HOMO)"

    #* Preprocess molecular information

    pymol = mc._pymol   # storing the molecule's info and the basis's info.
    task  = mc._pytask  # storing the calculation result.

    mo_coeff    = task.mo_coeff     # Linear combination coefficients of AO to make up MO.
    num_AO      = size(mo_coeff,1)  # Number of atomic orbits (aka AO) or Gaussian basis. (pymol.nao doesn't return an interger!)
    HOMO_orbitIdx = HOMOIndex(mc)
    orbitIdx = HOMO_orbitIdx + orbitIdx_relHOMO
    @assert 0<orbitIdx<num_AO "[PySCFMolecularCalculator] Orbit index $(orbitIdx_relHOMO) out of range."

    Z   = MolCharge(mc.mol)+1
    Ip  = -task.mo_energy[orbitIdx]
    if Ip ≤ 0
        error("[PySCFMolecularCalculator] The energy of the selected molecular orbital is positive.")
    end
    κ   = sqrt(2Ip)

    #* Define the spherical grids.
    r_grid = range(start=grid_rReg[1], stop=grid_rReg[2], length=grid_rNum)
    θ_grid = range(start=0., stop=π, length=grid_θNum)
    ϕ_grid = range(start=0., stop=2π, length=grid_ϕNum)
    N = grid_rNum*grid_θNum*grid_ϕNum

    "Returns the spherical grid indices of the given index of the point."
    function ptIdx2sphCoordIdx(i::Int)
        iϕ = (i-1) % grid_ϕNum + 1
        iθ = (ceil(Int, i/grid_ϕNum)-1) % grid_θNum + 1
        ir = ceil(Int, i/grid_ϕNum/grid_θNum)
        return (ir,iθ,iϕ)
    end

    # converting to cartesian grid points.
    # x = r sinθ cosϕ, y = r sinθ sinϕ, z = r cosθ.
    pt_x = zeros(N)
    pt_y = zeros(N)
    pt_z = zeros(N)
    dΩ   = zeros(N) # dΩ = sinθ dθdϕ, used in the integration.
    dr   = r_grid[2]-r_grid[1]
    dθ   = θ_grid[2]-θ_grid[1]
    dϕ   = ϕ_grid[2]-ϕ_grid[1]
    Threads.@threads for i in 1:N
        ir,iθ,iϕ = ptIdx2sphCoordIdx(i)
        r,θ,ϕ = r_grid[ir],θ_grid[iθ],ϕ_grid[iϕ]
        pt_x[i] = r*sin(θ)*cos(ϕ)
        pt_y[i] = r*sin(θ)*sin(ϕ)
        pt_z[i] = r*cos(θ)
        # dV[i]   = r^2*sin(θ)*dr*dθ*dϕ
        dΩ[i]   = sin(θ)*dθ*dϕ
    end
    pt_xyz = hcat(pt_x,pt_y,pt_z)

    #* defines output
    # C_lm stores the asymptotic coefficients, l=0,⋯,lMax; m=-l,⋯,l.
    # to obtain C_lm, refer to index [l+1,m+l+1].
    C_lm = zeros(l_max+1, 2l_max+1)

    #* Calculate the wavefunction ψ0
    χi = pymol.eval_gto("GTOval",pt_xyz)        # Size: N×Num_AO. Wavefunction of all AOs by calling eval_gto.
    orbit_coeff = @view mo_coeff[:,orbitIdx]    # Select the coefficients related to the interested MO.
    ψ0 = χi * orbit_coeff                       # Calculate wavefunction of the interested MO (matmul operation).

    #* Calculate F_lm(r) and fit C_lm.
    #*  project ψ0 to spherical harmonics Y_lm, obtaining F_lm(r).
    Threads.@threads for l in 0:l_max
        for m in -l:l
            F_lm = zeros(ComplexF64, grid_rNum)
            Y_lm = SphericalHarmonics.sphericalharmonic.(θ_grid, ϕ_grid'; l=l, m=m)
            # obtain F_lm(r)
            for i in 1:N
                ir,iθ,iϕ = ptIdx2sphCoordIdx(i)
                F_lm[ir] += conj(Y_lm[iθ,iϕ]) * ψ0[i] * dΩ[i]
            end
            # fit C_lm
            @. model(x,p) = p[1] * x^(Z/κ-1) * exp(-κ*x)
            p0 = [1.0]
            if sum(abs.(real.(F_lm))) > 1e-10
                fit_re = curve_fit(model, r_grid, real.(F_lm), p0)
                coeff = coef(fit_re)[1]     # the fitted C_lm
                conf_int = confidence_interval(fit_re)[1]   # confidence interval (95%)
                if coeff == 1.0   # returning the original guess means the fit is unsuccessful.
                    @warn "[PySCFMolecularCalculator] The fit of molecular wavefunction (l=$l, m=$m) is unsuccessful, try a more precise basis set or adjust the `grid_rReg` (the upper limit shouldn't be to large (<10 a.u.) !)."
                elseif abs((conf_int[2]-conf_int[1])/coeff) > 0.5   # the error is too large
                    @warn "[PySCFMolecularCalculator] The fit result of molecular wavefunction (l=$l, m=$m) is unsuccessful due to unacceptable error, try a more precise basis set or adjust the `grid_rReg` (the upper limit shouldn't be to large (<10 a.u.) !)."
                else
                    C_lm[l+1,m+l+1] += coeff
                end
            end
        end
    end
    return C_lm
end
