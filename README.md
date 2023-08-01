# <p align="center"> 🎆SemiclassicalSFI.jl </p>
<p align="center">
    Implementation of semiclassical methods in strong field ionization of atoms and molecules.
</p>

<p align="center">
    <a href="#background">Background</a> •
    <a href="#installation">Installation</a> •
    <a href="#usage">Usage</a> •
    <a href="#example">Example</a> •
    <a href="#troubleshooting">Troubleshooting</a> •
    <a href="#contributors">Contributors</a> •
    <a href="#license">License</a>
</p>

<p align="center">
    Last updated: May 3, 2023
</p>
<p align="center">
    Version 1.4.0
</p>

---------------------------

## Background

The interaction between laser and matter has attracted widespread interest since the invention of laser technology decades ago.
To study the interaction between an ultrafast and intense laser pulse and atoms/molecules, where the electrons are ionized from the targets through multi-photon or tunneling/over-barrier processes, a time-dependent Schrödinger equation (TDSE) simulation is usually required to be carried out.
However, its high demand in computational resources and limited application scope (atoms and simple molecules) prevent it from its extensive application.

As an alternative, semiclassical/classical electron trajectory simulation is widely used in numerical simulation in studies of strong-field ionization because it is less demanding in computational resources, which, in addition, provides a clear physical picture of strong-field ionization. This library written in julia aims to provide a general, efficient and out-of-box solution to perform trajectory simulations.

---------------------------

## Installation

### Prerequisites

Some prerequisites are listed below:
- [ **Hardware** ] A supported graphic card (optional, if you need to use GPU acceleration)
- [ **Platform** ] *Linux* & Windows (Linux is suggested because `PySCFMolecularCalculator` which is used to calculate molecular structure factors is incompatible with Windows)
- [ **Environment** ] Julia 1.7 & Python 3 (Python 3 is optional, but is required if you need to calculate the molecular structure factors)

### Installing dependency PySCF (optional)

Currently the calculation of molecular structure factors relies on the [PySCF](https://github.com/pyscf/pyscf) python package, which only supports the Linux platform. This library calls the PySCF using the julia [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) package, and there are two ways to set up the Python environment used by PyCall, here we suggest using your local Python environment for convenience:
```
$ julia
julia> ENV["PYTHON"] = "<path_to_python_exec>"
(@v1.7) pkg> add PyCall     # if you haven't installed it yet
(@v1.7) pkg> build PyCall
```
Remember to install PySCF in your python !!

The alternative choice is to install a private Python environment for julia via [Conda.jl](https://github.com/Luthaf/Conda.jl), for more information, please refer to [Conda.jl](https://github.com/JuliaPy/Conda.jl) and [PyCall.jl](https://github.com/JuliaPy/PyCall.jl).

### Installing the package

The package is still under development, to install the package, it's better to use the dev mode in the julia pkg manager (\<path> refers to a location where you want to store the source code, \<branch_name> refers to the branch or tag version of the library you want to install):
```
$ cd <path>
$ git clone https://github.com/TheStarAlight/SemiclassicalSFI.jl.git
$ cd ./SemiclassicalSFI.jl
$ git switch <branch_name/tag_name>
$ julia
(@v1.7) pkg> dev .
```
When you want to update or use another version, switch to the version you want using `git`, and "`dev`" it in julia's pkg manager.

---------------------------

## Usage

The usage of the program is simple and straightforward.
The main function of the program is wrapped in a method `performSFI`, with all simulation parameters provided, the program would automatically generate initial electron samples and perform trajectory simulations, collecting the electrons and mapping them to the momentum spectrum, and finally save the spectrum to the specified location.

The following are parameters of the method `performSFI` (the more detailed documentation is coming soon):

### Required params. for all methods
- `ionRateMethod = <:ADK|:SFA|:SFA_AE|:WFAT>`   : Method of determining ionization rate. Currently supports `:ADK`, `:SFA`, `:SFA_AE` for atoms and `:WFAT` for molecules.
- `laser::Laser`                                : Parameters of the laser field.
- `target::Target`                              : Parameters of the target.
- `sample_tSpan = (start,stop)`                 : Time span in which electrons are sampled.
- `sample_tSampleNum`                           : Number of time samples.
- `simu_tFinal`                                 : Time when every trajectory simulation ends.
- `finalMomentum_pMax = (pxMax,pyMax,pzMax)`    : Boundaries of final momentum spectrum collecting in three dimensions.
- `finalMomentum_pNum = (pxNum,pyNum,pzNum)`    : Numbers of final momentum spectrum collecting in three dimensions.

### Required params. for step-sampling methods
- `ss_kdMax`    : Boundary of kd (momentum's component along transverse direction (in xy plane)) samples.
- `ss_kdNum`    : Number of kd (momentum's component along transverse direction (in xy plane)) samples.
- `ss_kzMax`    : Boundary of kz (momentum's component along propagation direction (z ax.)) samples.
- `ss_kzNum`    : Number of kz (momentum's component along propagation direction (z ax.)) samples.

### Required params. for Monte-Carlo-sampling methods
- `mc_tBatchSize`   : Number of electron samples in a single time sample.
- `mc_ktMax`        : Maximum value of momentum's transversal component (perpendicular to field direction).

### Optional params. for all methods
- `save_fileName`                                           : Output HDF5 file name.
- `save_3D_momentumSpec = false`                            : Determines whether 3D momentum spectrum is saved.
- `simu_phaseMethod = <:CTMC|:QTMC|:SCTS>`                  : Method of classical trajectories' phase.
- `simu_relTol = 1e-6`                                      : Relative error tolerance when solving classical trajectories.
- `simu_nondipole = false`                                  : Determines whether non-dipole effect is taken account in the simulation.
- `simu_GPU = false`                                        : Determines whether GPU acceleration in trajectory simulation is used.
- `rate_monteCarlo = false`                                 : Determines whether Monte-Carlo sampling is used when generating electron samples.
- `rate_ionRatePrefix = <:ExpRate|:ExpPre|:ExpJac|:Full>`   : Prefix of the exponential term in the ionization rate. For MOADK & WFAT, `:ExpRate` & `:ExpPre` are the same.
- `rydberg_collect = false`                                 : Determines whether rydberg final states are collected.
- `rydberg_prinQNMax`                                       : Maximum principle quantum number n to be collected.

### Optional params. for target `Molecule`
- `mol_ionOrbitRelHOMO = 0` : Index of the ionizing orbit relative to the HOMO (e.g., 0 indicates HOMO, and -1 indicates HOMO-1) (default 0).

### Optional params. for MOADK method:
- `moadk_ionOrbit_m = 0`    : Magnetic quantum number m of the ionizing orbital along the z axis. m = 0,1,2 indicate σ, π and δ respectively.

### Optional params. for ADK method
- `adk_ADKTunExit = <:IpF|:FDM|:Para>`      : Tunneling exit method for ADK methods (when `ionRateMethod==:ADK`).


---------------------------

## Example

This section presents a short minimal example.

```julia
using SemiclassicalSFI

filename = "SCSFI_HydLike_Ip_0.5662.h5"  # output file name

# The target is a Hydrogen-like atom with ionization potential of 0.5662 a.u. and single nuclear charge.
t = SemiclassicalSFI.Targets.HydrogenLikeAtom(Ip=0.5662, Z=1)

# The laser is a 800nm (NIR) 6-cycle circularly polarized pulse with a cos^4-shaped envelope (propagating in z direction).
l = SemiclassicalSFI.Lasers.Cos4Laser(peakInt=3e14, waveLen=800, cycNum=6, ellip=1.0)

# Invokes the main method.
SemiclassicalSFI.performSFI(
    # Using ADK as the rate method to give initial conditions of electron samples.
    ionRateMethod = :ADK,
    laser = l,
    target = t,
    # Electrons ejected from the target between time range -300 and 300 a.u. would be sampled, and 20000 time samples would be picked at equal distances in the above time range.
    sample_tSpan = (-300,300),
    sample_tSampleNum = 20000,
    # Trajectory simulation of the ejected electrons would end at t=400 a.u. (a short span after the laser ends)
    simu_tFinal = 400,
    # Parameters of the final momentum spectrum. Electrons with momentum within (±3,±3,±3) a.u. would be collected, and mapped to a 800×800×1 cartesian grid (grid size in z direction is 1 and save_3D is false, which indicates that the simulation only saves 2D momentum spectrum).
    finalMomentum_pMax = (3,3,3),
    finalMomentum_pNum = (800,800,1),
    save_3D_momentumSpec = false,
    # Step-sampling parameters of electrons' initial momentum. `kd` refers to the momentum's component along transverse direction (in xy plane), while `kz` refers to the momentum's component along propagation direction (z ax.).
    ss_kdMax = 2.,
    ss_kdNum = 500,
    ss_kzMax = 2.,
    ss_kzNum = 150,
    # Classical Trajectory Monte-Carlo (CTMC), which indicates no account of phase.
    simu_phaseMethod = :CTMC,
    # GPU acceleration in trajectory simulation is enabled.
    simu_GPU = true,
    # output file name.
    save_fileName = filename
    )
```

---------------------------

## Troubleshooting

### Precompilation failure
Sometimes the precompilation of the package and its dependencies fails, which usually happens on SciML's packages, while no action in the pkg manager works.

Under such circumstances, try to delete the compiled julia code (usually stored in ~/.julia/compiled/\<julia_version>) and precompile again.

If the problem still exists after precompiling from scratch, you may try switching the SciML dependencies' versions in the julia.
As for the author, *OrdinaryDiffEq@6.51* and *DiffEqGPU@1.26* runs well on Windows 10, while for OrdinaryDiffEq@6.20/37/41, the precompilation never succeeded.

### Runtime warning: electrons with anomalous momentum
Sometimes (or usually, for some unlucky users) you may encounter such warning during the simulation if you use GPU acceleration:
```
[Ensemble Simulation] Found electron (#<...> in the batch) with anomalous momentum <...>.
```
This is possibly resulted from the DiffEqGPU package for some unknown reasons.
To solve this problem, try switching to another DiffEqGPU version (v1.26 is suggested).

---------------------------

## Contributors

- [Mingyu Zhu](https://github.com/TheStarAlight) @ ECNU
- Hongcheng Ni @ ECNU

---------------------------

## License

This package is licensed under the Apache 2.0 license, and is copyrighted by Mingyu Zhu and the other contributors.