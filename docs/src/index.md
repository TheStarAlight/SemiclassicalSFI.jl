# 🎆SemiclassicalSFI.jl

*Implementation of classical/semiclassical methods in strong-field ionization of atoms and molecules.*

```@contents
Pages = ["index.md"]
```

## Background

The interaction between laser and matter has attracted widespread interest since the invention of laser technology decades ago.
To study the interaction between an ultrafast and intense laser pulse and atoms/molecules, where the electrons are ionized from the targets through multi-photon or tunneling/over-barrier processes, a time-dependent Schrödinger equation (TDSE) simulation is usually required to be carried out.
However, its high demand in computational resources and limited application scope (atoms and simple molecules) prevents it from its extensive application.

To overcome the shortcomings of TDSE, Corkum *et al.* [^Corkum_1989] proposed a scheme, where the electron is first ionized from the target through the tunneling mechanism, and then acts as a classical electron in the laser field.
This scheme was further developed by Hu *et al.* [^Hu_1997], in which the initial conditions of the classical electrons and the Coulomb potential of the parent ion are more appropriatedly taken account.
This scheme is named after the *Classical Trajectory Monte-Carlo (CTMC)* method, which has been widely adopted for research in interaction between high-intensity ultra-fast laser pulses and atoms/molecules.
Compared with TDSE, trajectory simulation schemes including CTMC and its variants, are less demanding in computational resources, which, in addition, provides a clear physical picture of strong-field ionization.

The essence of the trajectory simulation scheme lies in two aspects:
(1) The initial conditions of the classical electron samples at the beginning of the classical trajectories, which consists of initial position $\bm{r}_0$ (i.e., the tunneling exit position), initial momenta $\bm{p}_0$, and the corresponding ionization probability $W$ carried by the electron sample.
(2) The quantum phase property of classical trajectories, while the full classical trajectory (i.e., the CTMC) is widely adopted, there are schemes (e.g., QTMC and SCTS, which would be discussed further in the documentation) which introduce quantum phases in the electron trajectories and develop a semiclassical method for trajectory simulations.

After decades of accumulation of research and development, the trajectory simulation has grown to a complete solution of research on strong-field ionization of atoms and molecules. Developing a library with implementation of existing methods, efficiency of calculation, extensibility for future development and ease of maintenance would provide great convenience for theoretical research on strong-field ionization. With such aim, here we present *SemiclassicalSFI.jl*, a program package written in julia language, which provides a general, efficient and out-of-box solution of performing trajectory simulations.

[^Corkum_1989]: P. B. Corkum *et al.*, Above-Threshold Ionization in the Long-Wavelength Limit. *Phys. Rev. Lett.* **62**(11), 1259–1262 (1989). DOI: [10.1103/PhysRevLett.62.1259](https://dx.doi.org/10.1103/PhysRevLett.62.1259)

[^Hu_1997]: B. Hu *et al.*, Plateau in Above-Threshold-Ionization Spectra and Chaotic Behavior in Rescattering Processes. *Phys. Lett. A* **236**, 533–542 (1997). DOI: [10.1016/S0375-9601(97)00811-6](https://dx.doi.org/10.1016/S0375-9601(97)00811-6)

## Features

- *Versatile* :     *SemiclassicalSFI.jl* supports a wide range of functions. As for initial conditions (rate method), the library supports (for atoms) *ADK*, *SFA* and *SFA-AE*, (for molecules) *MOADK* and *WFAT*. As for the trajectory phase method, the library supports *CTMC*, *QTMC* and *SCTS*. Non-dipole effects can also be included during the trajectory simulation.
- *Out-of-box* :    The usage of *SemiclassicalSFI.jl* is simple and straightforward.
- *Extensible* :    *SemiclassicalSFI.jl* has a well-defined structure, which makes it easy to include new features.

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

!!! note "Possible solution to precompilation failure"

    Sometimes the precompilation of the package and its dependencies fails, which usually happens on SciML's packages.
    Under such circumstances, try to delete the compiled julia code (usually stored in `~/.julia/compiled/<julia_version>`) and precompile again.
    If the problem still exists after precompiling from scratch, you may try switching the dependencies' versions in the julia, which is done by specifying the version when adding the packages:
    ```julia
    using Pkg
    Pkg.add(name="package_name", version="1.0")
    # In pkg mode of REPL:
    # (@v1.8) pkg> add package_name@1.0
    ```

    It is shown that *OrdinaryDiffEq@6.53.3* and *DiffEqGPU@2.4.1* runs well on Windows 10 (10.0.19044), WSL Ubuntu (22.04.1 LTS) and Manjaro Linux (Uranos 23.0.0).

### Configuring Python and PySCF

Currently the [MO-ADK](@ref MOADK) and [WFAT](@ref WFAT) features related to molecules rely on the [PySCF](https://github.com/pyscf/pyscf) python package. *SemiclassicalSFI.jl* calls the PySCF using the [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) package. There are two ways to set up the Python environment used by PyCall, here we suggest using your local Python environment for convenience.

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

!!! tip "Note"
    Since the PySCF does not support the Windows, the molecular calculation must be performed on a Linux or macOS platform.
    However, for Windows users, they may install the [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) (Windows Subsystem for Linux), which supports the PySCF.


## Contributors

- [Mingyu Zhu](https://github.com/TheStarAlight) @ ECNU
- [Hongcheng Ni](https://faculty.ecnu.edu.cn/_s29/nhc_en/main.psp) @ ECNU

## License

This package is licensed under the Apache 2.0 license, and is copyrighted by Mingyu Zhu, Hongcheng Ni and the other contributors.
