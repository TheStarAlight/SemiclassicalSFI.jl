# <p align="center"> 🎆SemiclassicalSFI.jl </p>
<p align="center">
    Implementation of semiclassical methods in strong field ionization of atoms and molecules.
</p>

<p align="center">
    <a href="#background">Background</a> •
    <a href="#features">Features</a> •
    <a href="#installation">Installation</a> •
    <a href="#usage">Usage</a> •
    <a href="#example">Example</a> •
    <a href="#troubleshooting">Troubleshooting</a> •
    <a href="#contributors">Contributors</a> •
    <a href="#license">License</a>
</p>

<p align="center">
    • Documentation available at <a href="https://TheStarAlight.github.io/SemiclassicalSFI.jl">TheStarAlight.github.io/SemiclassicalSFI.jl</a> •
</p>

<p align="center">
    Last updated: July 23, 2023
</p>
<p align="center">
    Version 1.4.0
</p>

---------------------------

## Background

The interaction between laser and matter has attracted widespread interest since the invention of laser technology decades ago.
To study the interaction between an ultrafast and intense laser pulse and atoms/molecules, where the electrons are ionized from the targets through multi-photon or tunneling/over-barrier processes, a time-dependent Schrödinger equation (TDSE) simulation is usually required to be carried out.
However, its high demand in computational resources and limited application scope (atoms and simple molecules) prevents it from its extensive application.

To overcome the shortcomings of TDSE, Corkum *et al.* [^Corkum_1989] proposed a scheme, where the electron is first ionized from the target through the tunneling mechanism, and then acts as a classical electron in the laser field.
This scheme was further developed by Hu *et al.* [^Hu_1997], in which the initial conditions of the classical electrons and the Coulomb potential of the parent ion are more appropriatedly taken account.
This scheme is named after the *Classical Trajectory Monte-Carlo (CTMC)* method, which has been widely adopted for research in interaction between high-intensity ultra-fast laser pulses and atoms/molecules.
Compared with TDSE, trajectory simulation schemes including CTMC and its variants, are less demanding in computational resources, which, in addition, provides a clear physical picture of strong-field ionization.

The essence of the trajectory simulation scheme lies in two aspects:
(1) The initial conditions of the classical electron samples at the beginning of the classical trajectories, which consists of initial position $\vec{r}_0$ (i.e., the tunneling exit position), initial momenta $\vec{p}_0$, and the corresponding ionization probability $W$ carried by the electron sample.
(2) The quantum phase property of classical trajectories, while the full classical trajectory (i.e., the CTMC) is widely adopted, there are schemes (e.g., QTMC and SCTS, which would be discussed further in the documentation) which introduce quantum phases in the electron trajectories and develop a semiclassical method for trajectory simulations.

After decades of accumulation of research and development, the trajectory simulation has grown to a complete solution of research on strong-field ionization of atoms and molecules. Developing a library with implementation of existing methods, efficiency of calculation, extensibility for future development and ease of maintenance would provide great convenience for theoretical research on strong-field ionization. With such aim, here we present *SemiclassicalSFI.jl*, a program package written in julia language, which provides a general, efficient and out-of-box solution of performing trajectory simulations.

[^Corkum_1989]: P. B. Corkum *et al.*, Above-Threshold Ionization in the Long-Wavelength Limit. *Phys. Rev. Lett.* **62**(11), 1259–1262 (1989). DOI: [10.1103/PhysRevLett.62.1259](https://dx.doi.org/10.1103/PhysRevLett.62.1259)

[^Hu_1997]: B. Hu *et al.*, Plateau in Above-Threshold-Ionization Spectra and Chaotic Behavior in Rescattering Processes. *Phys. Lett. A* **236**, 533–542 (1997). DOI: [10.1016/S0375-9601(97)00811-6](https://dx.doi.org/10.1016/S0375-9601(97)00811-6)

---------------------------

## Features

- *Versatile* : *SemiclassicalSFI.jl* supports a wide range of functions. As for initial conditions (rate method), the library supports (for atoms) *ADK*, *SFA* and *SFA-AE*, (for molecules) *MOADK* and *WFAT*. As for the trajectory phase method, the library supports *CTMC*, *QTMC* and *SCTS*. Non-dipole effects can also be included during the trajectory simulation.
- *Out-of-box* : The usage of *SemiclassicalSFI.jl* is simple and straightforward.
- *Extensible* : *SemiclassicalSFI.jl* has a well-defined structure, which makes it easy to include new features.

---------------------------

## Installation

### Prerequisites

- *Minimum prerequisites* : Julia ≥1.7

- *GPU acceleration of traj. simulation* : a supported graphic card (NVIDIA)

- *MOADK and WFAT features* : Linux or macOS platform, Python 3 with the [PySCF](https://github.com/pyscf/pyscf) python package installed and the [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) package successfully built.

### Installing the package

This package is currently not in julia's general registry, but can be added through the repository URL:

```julia
using Pkg
Pkg.add(url="https://github.com/TheStarAlight/SemiclassicalSFI.jl.git")
# In pkg mode of REPL:
# (@v1.8) pkg> add https://github.com/TheStarAlight/SemiclassicalSFI.jl.git
```

It is suggested to test the package to check if the functions check if some special features (e.g., GPU acceleration and molecular calculation) work on your platform:

```julia
Pkg.test("SemiclassicalSFI")
# In pkg mode of REPL:
# (@v1.8) pkg> test SemiclassicalSFI
```

### Configuring Python and PySCF

Currently the MO-ADK and WFAT features related to molecules rely on the [PySCF](https://github.com/pyscf/pyscf) python package. *SemiclassicalSFI.jl* calls the PySCF using the [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) package. There are two ways to set up the Python environment used by PyCall, here we suggest using your local Python environment for convenience.

To correctly set up the configuration of PyCall, first, set the `PYTHON` environment variable to your Python executable, and build the PyCall package:

```julia
ENV["PYTHON"] = "path/to/python_exec"
using Pkg
Pkg.build("PyCall")
```

And don't forget to install PySCF in your Python via pip:

```
$ pip3 install pyscf
```


**Note**: Since the PySCF does not support the Windows, the molecular calculation must be performed on a Linux or macOS platform.
However, for Windows users, they may install the [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) (Windows Subsystem for Linux), which supports the PySCF.

---------------------------

## Usage

The usage of the program is simple and straightforward.
The main function of the program is wrapped in a method `performSFI`, with all simulation parameters provided, the program would automatically generate initial electron samples and perform trajectory simulations, collecting the electrons and mapping them to the momentum spectrum, and finally save the spectrum to the specified location.

The following are parameters of the method `performSFI` (for more detailed documentation, see [here](https://TheStarAlight.github.io/SemiclassicalSFI.jl/stable/manual3_main_method/)):

### Required params. for all methods:
- `init_cond_method = <:ADK|:SFA|:SFAAE|:WFAT|:MOADK>`  : Method of electrons' initial conditions. Currently supports `:ADK`, `:SFA`, `:SFAAE` for atoms and `:WFAT`, `:MOADK` for molecules.
- `laser::Laser`                                        : A `Lasers.Laser` object containing information of the laser field.
- `target::Target`                                      : A `Targets.Target` object containing information of the atom/molecule target.
- `sample_t_intv = (start,stop)`                        : Time interval in which the initial electrons are sampled.
- `sample_t_num`                                        : Number of time samples.
- `traj_t_final`                                        : Time when every trajectory simulation ends.
- `final_p_max = (pxMax,pyMax,pzMax)`                   : Boundaries of final momentum spectrum collected in three dimensions.
- `final_p_num = (pxNum,pyNum,pzNum)`                   : Numbers of final momentum spectrum collected in three dimensions.

### Required params. for step-sampling methods:
- `ss_kd_max`   : Boundary of kd (momentum's component along transverse direction (in xy plane)) samples.
- `ss_kd_num`   : Number of kd (momentum's component along transverse direction (in xy plane)) samples.
- `ss_kz_max`   : Boundary of kz (momentum's component along propagation direction (z ax.)) samples.
- `ss_kz_num`   : Number of kz (momentum's component along propagation direction (z ax.)) samples.

### Required params. for Monte-Carlo-sampling methods:
- `mc_kp_num`   : Number of kp (initial momentum which is perpendicular to field direction, two dimensional) samples in a single time sample.
- `mc_kp_max`   : Maximum value of momentum's transversal component (perpendicular to field direction).

### Optional params. for all methods:
- `save_path`                                       : Output HDF5 file path.
- `save_3D_spec = false`                            : Determines whether to save the 3D momentum spectrum (otherwise 2D) (default `false`).
- `traj_phase_method = <:CTMC|:QTMC|:SCTS>`         : Method of classical trajectories' phase (default `:CTMC`). Currently `:QTMC` and `:SCTS` only support atomic cases.
- `traj_rtol = 1e-6`                                : Relative error tolerance when solving classical trajectories using adaptive methods (default `1e-6`).
- `traj_nondipole = false`                          : Determines whether the non-dipole effect is taken account in the simulation (default `false`).
- `traj_GPU = false`                                : [Experimental] Determines whether to enable GPU acceleration in trajectory simulation (default `false`).
- `sample_monte_carlo = false`                      : Determines whether Monte-Carlo sampling is used when generating electron samples (default `false`). Currently only supports ADK.
- `final_ryd_collect = false`                       : Determines whether the rydberg final states are collected (default `false`).
- `final_ryd_n_max`                                 : Determines the maximum principle quantum number n for rydberg final states to be collected.

### Optional params. for atomic SFA, SFA-AE and ADK methods:
- `rate_prefix = <:ExpRate|:ExpPre|:ExpJac|:Full>`  : Prefix of the exponential term in the ionization rate (default `:ExpRate`).

### Optional params. for target `Molecule`:
- `mol_orbit_idx = 0`   : Index of the ionizing orbit relative to the HOMO (e.g., `0` indicates HOMO, and `-1` indicates HOMO-1) (default `0`).

### Optional params. for MO-ADK method:
- `moadk_orbit_m = 0`   : Magnetic quantum number m of the ionizing orbital along the z axis. m = 0,1,2 indicate σ, π and δ respectively (default `0`).

### Optional params. for ADK method:
- `adk_tun_exit = <:IpF|:FDM|:Para>` : Tunneling exit method for ADK methods (when `init_cond_method==:ADK`) (default `:IpF`).


---------------------------

## Example

This section presents a short minimal example.

```julia
using SemiclassicalSFI

filename = "SCSFI_HydLike_Ip_0.5662.h5"  # output file name

# The target is a Hydrogen-like atom with ionization potential of 0.5662 a.u. and single nuclear charge.
t = SemiclassicalSFI.Targets.HydrogenLikeAtom(Ip=0.5662, Z=1)

# The laser is a 800nm (NIR) 2-cycle circularly polarized pulse with a cos^4-shaped envelope (propagating in z direction).
l = SemiclassicalSFI.Lasers.Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)

# Invokes the main method.
SemiclassicalSFI.performSFI(
    # Using ADK as the rate method to give initial conditions of electron samples.
    init_cond_method = :ADK,
    laser = l,
    target = t,
    # Electrons ejected from the target between time range -80 and 80 a.u. would be sampled, and 5000 time samples would be picked at equal distances in the above time range.
    sample_t_intv = (-80,80),
    sample_t_num = 5000,
    # Trajectory simulation of the ejected electrons would end at t=120 a.u. (a short span after the laser ends)
    traj_t_final = 120,
    # Parameters of the final momentum spectrum. Electrons with momentum within (±2,±2,±2) a.u. would be collected, and mapped to a 800×800×1 cartesian grid (grid size in z direction is 1 and save_3D is false, which indicates that the simulation only saves 2D momentum spectrum).
    final_p_max = (2,2,2),
    final_p_num = (500,500,1),
    save_3D_spec = false,
    # Step-sampling parameters of electrons' initial momentum. `kd` refers to the momentum's component along transverse direction (in xy plane), while `kz` refers to the momentum's component along propagation direction (z ax.).
    ss_kd_max = 2.,
    ss_kd_num = 500,
    ss_kz_max = 2.,
    ss_kz_num = 150,
    # Classical Trajectory Monte-Carlo (CTMC), which indicates no account of phase.
    traj_phase_method = :CTMC,
    # GPU acceleration in trajectory simulation is enabled.
    traj_GPU = true,
    # output file name.
    save_path = filename
    )
```

---------------------------

## Troubleshooting

### Possible solution to precompilation failure

Sometimes the precompilation of the package and its dependencies fails, which usually happens on SciML's packages.
Under such circumstances, try to delete the compiled julia code (usually stored in `~/.julia/compiled/<julia_version>`) and precompile again.
If the problem still exists after precompiling from scratch, you may try switching the dependencies' versions in the julia, which is done by specifying the version when adding the packages:
```julia
using Pkg
Pkg.add(name="package_name", version="1.0")
# In pkg mode of REPL:
# (@v1.8) pkg> add package_name@1.0
```

It is shown that *OrdinaryDiffEq@6.53.3* and *DiffEqGPU@2.4.1* runs well on Windows 10 (10.0.19044) and WSL Ubuntu (22.04.1 LTS).

---------------------------

## Contributors

- [Mingyu Zhu](https://github.com/TheStarAlight) @ ECNU
- [Hongcheng Ni](https://faculty.ecnu.edu.cn/_s29/nhc_en/main.psp) @ ECNU

---------------------------

## License

This package is licensed under the Apache 2.0 license, and is copyrighted by Mingyu Zhu and the other contributors.