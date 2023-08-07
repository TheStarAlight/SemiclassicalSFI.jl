# Main Method `performSFI`

*This section introduces the main method* [`performSFI`](@ref)*, which performs the trajectory simulation and saves the final electron momentum spectrum.*

```@contents
Pages = ["manual3_main_method.md"]
Depth = 3
```

```@meta
CurrentModule = SemiclassicalSFI
```


## Brief Documentation

```@docs
performSFI
```



## Lasers & Targets

### Lasers

- `laser::Laser`

A `Lasers.Laser` object containing information of the laser field.
Cf. the documentation for [Lasers](@ref lasers_doc).

### Targets

- `target::Target`

A `Targets.Target` object containing information of the atom/molecule target.
Cf. the documentation for [Targets](@ref targets_doc).

### Workflow for Preparation of the `Molecule` Target

To use the [`Molecule`](@ref) as the target and the relating initial condition methods [MO-ADK](@ref MOADK) or [WFAT](@ref WFAT), some coefficients need to be calculated in advance.
Here we present a workflow for preparation of the `Molecule` target before invoking the `performSFI` method, taking the carbon monoxide molecule as an example.

!!! note "Note"
    Currently only the PySCF is implemented as the library's molecular computation interface, which only supports the Linux platform.
    What's more, in current version, only *close-shell* molecules are supported!


#### Initialization

First of all, initialize a `Molecule` object, and provide the necessary information of the molecule. Cf. the documentation of [`Molecule`](@ref).

```julia
using SemiclassicalSFI.Targets
mol = Molecule(atoms=["C","O"], atom_coords=[0 0 -0.180; 0 0 0.950],
               charge=0, name="Carbon Monoxide",
               data_path="./Molecule_CarbonMonoxide.h5")
```

!!! warning "Unit of the atom coordinates"
    The input of the atom coordinates is in **Angstrom (Å)**, which is not atomic unit (Bohr). 1 Bohr = 0.53 Å.

!!! tip "Data saving of the Molecule object"
    If the user specifies `data_path` in the constructor method of `Molecule`, the data would be automatically saved each time the user invokes the [`MolCalcAsympCoeff!`](@ref) and [`MolCalcWFATData!`](@ref).
    However, if doesn't specify (in case the user does not wish to save the data), the data would not be saved, and the user has to manually invoke [`MolSaveDataAs`](@ref) to save the data afterwards.


#### Calculate Asymptotic Coefficients

To calculate the asymptotic coefficients of the `Molecule` object which is used in MO-SFA and [MO-ADK](@ref MOADK), invoke the [`MolCalcAsympCoeff!`](@ref) method:
```julia
MolCalcAsympCoeff!(mol)
```

For typical small molecules, using default parameters usually gives satisfactory results.
However, for special demands, the user may refer to the documentation of [`MolCalcAsympCoeff!`](@ref) and [`Targets.MolecularCalculators.calc_asymp_coeff`](@ref) for more calculation parameters.

#### Calculate WFAT Data

To calculate the data necessary for the [WFAT](@ref WFAT) of the `Molecule` object, invoke the [`MolCalcWFATData!`](@ref) method:
```julia
MolCalcWFATData!(mol, orbitIdx_relHOMO = 0)
```

To obtain the data of other orbitals besides HOMO, the user may alter the `orbitIdx_relHOMO` parameter.
For custom calculation parameters, refer to the documentation of [`MolCalcWFATData!`](@ref) and [`Targets.MolecularCalculators.calc_WFAT_data`](@ref).

#### Setting the Molecule's Orientation

The orientation of the molecule can be specified by invoking the [`SetMolRotation`](@ref) method:
```julia
SetMolRotation(mol, 0.0,π/2,π/3)
```
and can be obtained through the [`MolRotation`](@ref) method:
```julia
MolRotation(mol)
```

!!! tip "Testing molecular calculation setup"
    The test sets of this library include molecular calculations, the user may run the tests of this library to check if the environment is correctly set up for molecular calculation:
    ```julia
    using Pkg
    Pkg.test("SemiclassicalSFI")
    # In pkg mode of REPL:
    # (@v1.8) pkg> test SemiclassicalSFI
    ```


## Initial Condition Methods

- `init_cond_method = <:ADK|:SFA|:SFAAE|:WFAT|:MOADK>`

Method of electrons' initial conditions.
Currently supports `:ADK`, `:SFA`, `:SFAAE` for atoms, and `:WFAT`, `:MOADK` for molecules.

For more information about the theories, cf. [Theory - Initial Conditions](@ref theory_init_cond).

!!! note "Note"
    For initial condition methods [ADK](@ref ADK), [SFA](@ref SFA) or [SFA-AE](@ref SFAAE), you must specify an atom target of types [`HydrogenLikeAtom`](@ref) or [`SAEAtom`](@ref);
    For initial condition methods [WFAT](@ref WFAT) or [MO-ADK](@ref MOADK), you must specify a molecule target of type [`Molecule`](@ref).

    |  | [ADK](@ref ADK) | [SFA](@ref SFA) | [SFA-AE](@ref SFAAE) | [WFAT](@ref WFAT) | [MO-ADK](@ref MOADK) |
    | :------------------------- |:-:|:-:|:-:|:-:|:-:|
    | [`HydrogenLikeAtom`](@ref) | ✔ | ✔ | ✔ |   |   |
    | [`SAEAtom`](@ref)          | ✔ | ✔ | ✔ |   |   |
    | [`Molecule`](@ref)         |   |   |   | ✔ | ✔ |

### Atomic Rate Prefix

- `rate_prefix = <:ExpRate|:ExpPre|:ExpJac|:Full>`

Prefix of the exponential term in the ionization rate (default `:ExpRate`).

For atomic [ADK](@ref ADK), [SFA](@ref SFA) and [SFA-AE](@ref SFAAE), we obtained the ionization probability in the following form:
```math
\mathrm{d}W/\mathrm{d}\bm{k}_\perp \mathrm{d}t_{\mathrm{r}} = \bm{J}(k_d,t_{\mathrm{r}}) \lvert P_{\bm{p}}(t_{\mathrm{s}}) \rvert^2 \exp(2\ \mathrm{Im}\ \Phi_{\mathrm{tun}}),
```
where the ``\lvert P_{\bm{p}}(t_{\mathrm{s}}) \rvert^2`` denotes the prefactor, which has different expressions for different theories:

| [SFA](@ref SFA) | [SFA-AE](@ref SFAAE) | [ADK](@ref ADK) |
| :-: | :-: | :-: |
| ``\{ [\bm{p}+\bm{A}(t_{\mathrm{s}})] \cdot \bm{F}(t_{\mathrm{s}}) \}^{-\alpha}`` | ``\left[ (k_\perp^2+2I_{\mathrm{p}})(F^2-\bm{k}_\perp \cdot \bm{F}') \right]^{-\alpha/2}`` | ``\left[ (k_\perp^2+2I_{\mathrm{p}})F^2\right]^{-\alpha/2}`` |

The ``\bm{J}(k_d,t_{\mathrm{r}})`` denotes the Jacobian which arises from the coordinate transformation.

For initial condition methods [ADK](@ref ADK), [SFA](@ref SFA) and [SFA-AE](@ref SFAAE),
specifying `ExpRate` would not add any prefix besides the exponential term;
`ExpPre` would include the prefactor ``\lvert P_{\bm{p}}(t_{\mathrm{s}}) \rvert^2``;
`ExpJac` would include the Jacobian ``\bm{J}(k_d,t_{\mathrm{r}})``;
`Full` indicates inclusion of both the prefactor and Jacobian.

### Tunneling Exit Methods For Atomic ADK

- `adk_tun_exit = <:IpF|:FDM|:Para>`

Tunneling exit method for atomic ADK methods (when `init_cond_method==:ADK`) (default `:IpF`).
Cf. [Tunneling Exit Methods For Atomic ADK](@ref tun_exit_atomic_adk).



## Sampling Methods and Parameters

In *SemiclassicalSFI.jl* the initial electrons are sampled in the ``(t_{\mathrm{r}},k_d,k_z)`` coordinate, where ``k_d`` denotes the initial momentum's component in the ``xy`` plane (which is perpendicular to the electric field).
There are two ways of sampling in these coordinates, namely *step sampling* and *Monte-Carlo sampling*.

- `sample_monte_carlo = false` : Determines whether Monte-Carlo sampling is used when generating electron samples (default `false`).

- `sample_t_intv = (start,stop)` : Time interval in which the initial electrons are sampled.

- `sample_t_num` : Number of time samples.

- `sample_cutoff_limit = 1e-16` : The cut-off limit of the probability of the sampled electron, electrons with probabilities lower than the limit would be discarded (default `1e-16`).

### Step Sampling

In the step sampling scheme, the `sample_t_num` time samples are uniformly distributed in the interval `sample_t_intv`.
In each time sample, a batch of electrons of different initial conditions are launched and collected, whose initial momenta ``\bm{k}_\perp`` are distributed on a Cartesian grid ``(k_d,k_z)``.
The Cartesian grid is defined by the following parameters:

- `ss_kd_max`, `ss_kd_num`, `ss_kz_max`, `ss_kz_num`

In the ``k_d`` dimension, `ss_kd_num` samples distribute uniformly in the interval (-`ss_kd_max`,`ss_kd_max`);
and in the ``k_z`` dimension, there are `ss_kz_num` equidistant samples in the interval (-`ss_kz_max`,`ss_kz_max`).

The step sampling method is supported for all initial condition methods.

!!! warning "ss_kz_num should be an even number"
    Under the dipole case, electrons with zero initial ``k_z`` would remain zero ``k_z``, which may induce some problems in the step-length choice of the ODE solver and result in "final electron with anomalously large momentum" warning. Thus the sampling on the ``k_z`` direction should avoid zero. Setting `ss_kz_num` to an even number is a good choice.

### Monte-Carlo Sampling

In the Monte-Carlo sampling scheme, the `sample_t_num` time samples are randomly chosen in the `sample_t_intv`;
A batch containing `mc_kp_num` electrons would be sampled in a single time sample, the electrons' initial momenta ``\bm{k}_\perp`` are also randomly sampled inside a circle ``k_d^2+k_z^2 \leq k_{\perp\mathrm{max}}^2``, where the ``k_{\perp\mathrm{max}}`` is defined in the parameter as `mc_kp_max`.

- `mc_kp_num` : Number of kp (initial momentum which is perpendicular to field direction, two dimensional) samples in a single time sample.
- `mc_kp_max` : Maximum value of momentum's transversal component (perpendicular to field direction).

Currently the Monte-Carlo sampling method is only supported for the ADK initial condition method.



## Trajectory Simulation

After preparation of the initial electrons, the electrons evolve classically in the combined potential of the nucleus and laser field, and the trajectory simulation terminates at `traj_t_final`.

- `traj_t_final` : Time when every trajectory simulation ends.

### Phase Methods

- `traj_phase_method = <:CTMC|:QTMC|:SCTS>`

Method of classical trajectories' phase (default `CTMC`).
Currently `:QTMC` and `:SCTS` only support atomic cases.

|                   | [ADK](@ref ADK) | [SFA](@ref SFA) | [SFA-AE](@ref SFAAE) | [WFAT](@ref WFAT) | [MO-ADK](@ref MOADK) |
| :---------------- |:-:|:-:|:-:|:-:|:-:|
| [CTMC](@ref CTMC) | ✔ | ✔ | ✔ | ✔ | ✔ |
| [QTMC](@ref QTMC) | ✔ | ✔ | ✔ |   |   |
| [SCTS](@ref SCTS) | ✔ | ✔ | ✔ |   |   |

For more information about phase methods, cf. [Theory - Trajectory Simulation and Phase Methods](@ref theory_traj_phase).

### Non-dipole Effect

- `traj_nondipole = false`

Determines whether the non-dipole effect is taken account in the simulation (default `false`).

For more information about the non-dipole effects, cf. [Theory - Non-dipole Effects on the Trajectory Motion](@ref traj_nondipole).

Currently, all targets support the inclusion of non-dipole effects.

### Accuracy Control

- `traj_rtol = 1e-6`

Relative error tolerance when solving classical trajectories using adaptive methods (default `1e-6`).

The classical trajectories of the electrons are obtained by solving ordinary differential equations using the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) package, where the solver would adjust the time-step length according to the relative error.
Stricter error tolerance is favorable if the electron moves close to the nucleus, which usually happens in a linearly polarized laser field.

### GPU Acceleration (Experimental)

- `traj_GPU = false`

Determines whether to enable GPU acceleration in trajectory simulation (default `false`).

The ordinary differential equations related to the electrons' classical trajectories can also be solved by GPU via the [DiffEqGPU.jl](https://github.com/SciML/DiffEqGPU.jl) package.
The GPU acceleration now only supports the NVIDIA graphic cards, which also requires the CUDA driver to be correctly installed.

!!! tip "Testing GPU capability"
    The test sets of this library include GPU tests, the user may run the tests of this library to check if the graphic card is ready for GPU acceleration:
    ```julia
    using Pkg
    Pkg.test("SemiclassicalSFI")
    # In pkg mode of REPL:
    # (@v1.8) pkg> test SemiclassicalSFI
    ```

!!! compat "Note: Experimental feature"
    GPU acceleration is an experimental feature of the package and the API may change in the near future.

!!! note "Note: GPU acceleration efficiency"
    The GPU's advantage of speed over CPU may not be so obvious in the trajectory simulation.
    Sometimes the program is even slower when using GPU.
    Thus it's better to benchmark for specific hardware and choose the most efficient option.

## Final Electron Collecting & Saving

After the trajectory simulation ends, the electrons would be analyzed and collected.
Those with positive energies finally become free electrons and would reach the detectors;
while those with negative energies finally fall on Rydberg states.
The collected momentum spectra and Rydberg spectra would be saved in an HDF5 file together with the simulation abstract.
This library provides some parameters to customize the collecting and saving procedure.

### 2D/3D Momentum Spectrum Collecting

- `final_p_max = (pxMax,pyMax,pzMax)` : Boundaries of final momentum spectrum collected in three dimensions.
- `final_p_num = (pxNum,pyNum,pzNum)` : Numbers of final momentum spectrum collected in three dimensions.
- `save_3D_spec = false` : Determines whether the 3D momentum spectrum is saved (if not, will only save 2D by flattening on the xy plane) (default `false`).

Electrons with positive final energies would be collected and placed on the 3D momentum grid determined by `final_p_max` and `final_p_num`.
When setting `save_3D_spec = false`, the three-dimensional final momentum spectrum would be squashed into two-dimensional ones (by summing over the z axis).

The ``p_x``, ``p_y`` and ``p_z`` grids would be saved in entries that are named after `px`, `py` and `pz` respectively in the output file, and the final momentum spectrum would be saved in the `momentum_spec_2D` or `momentum_spec_3D` entry (compressed to reduce the file size).

### Rydberg Final State Collecting

- `final_ryd_collect = false` : Determines whether the rydberg final states are collected (default `false`).
- `final_ryd_n_max` : Determines the maximum principle quantum number n for rydberg final states to be collected.

Rydberg final states would be collected if `final_ryd_collect` is set to `true`.
Only Rydberg states with principle quantum number that below `final_ryd_n_max` would be collected.

The Rydberg spectrum is saved as a three-dimensional array in the entry named after `ryd_spec` in the output file.
To get the probability of Rydberg state at ``(n,l,m)``, index with the indices `(n,l+1,m+final_ryd_n_max)`.

### Output File

- `save_path` : Output HDF5 file path.

The output HDF5 file would be saved in the `save_path`.
If the parameter is left unspecified or an error occurs when trying to write to the specified path, the output path would be set to `./SCSFI-yyyymmdd-hhmmss.h5`.

Apart from the information related to the momentum and Rydberg spectrum, an abstract encoded in YAML, which contains necessary input parameters, is also saved in the `abstract` entry of the output file.

The output file can be opened and accessed using the [HDF5.jl](https://github.com/JuliaIO/HDF5.jl) or [JLD2.jl](https://github.com/JuliaIO/JLD2.jl) packages.
Below shows the structure of typical output files.
```
# saving only 2D momentum spectrum
🗂️ HDF5.File
├─ 🔢 abstract
├─ 🔢 info
├─ 🔢 ion_prob
├─ 🔢 ion_prob_uncollected
├─ 🔢 momentum_spec_2D
├─ 🔢 num_effective_traj
├─ 🔢 px
└─ 🔢 py
```
```
# saving 3D momentum spectrum
🗂️ HDF5.File
├─ 🔢 abstract
├─ 🔢 info
├─ 🔢 ion_prob
├─ 🔢 ion_prob_uncollected
├─ 🔢 momentum_spec_2D
├─ 🔢 momentum_spec_3D
├─ 🔢 num_effective_traj
├─ 🔢 px
├─ 🔢 py
└─ 🔢 pz
```

