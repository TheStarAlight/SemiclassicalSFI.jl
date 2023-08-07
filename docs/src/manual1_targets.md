# [Targets](@id targets_doc)

*This section provides information of available targets (atoms/molecules) in the library.*

A target interacts with the laser field and release the electron through multi-photon or tunneling processes.
Here we list available targets implemented in the [`Targets`](@ref) module of the library.

```@docs
Targets
```

```@contents
Pages = ["manual1_targets.md"]
Depth = 3
```

```@meta
CurrentModule = SemiclassicalSFI.Targets
```

```@setup manual_targets
using SemiclassicalSFI
using SemiclassicalSFI.Targets
```

## List of Available Properties and Methods

The available properties and methods of the targets are listed below.
To obtain a property of the target, invoke the property as a method and pass the target object as an argument. The following shows an example.

```@repl manual_targets
t = HeAtom()
IonPotential(t)
V = TargetPotential(t)
V(1.0,0.0,0.0)
```

### Atoms

|                   |[`HydrogenLikeAtom`](@ref) | [`SAEAtom`](@ref) |
|:------------------|:-:|:-:|
|`IonPotential`     | ✔ | ✔ |
|`AsympNuclCharge`  | ✔ | ✔ |
|`SoftCore`         | ✔ |   |
|`TargetName`       | ✔ | ✔ |
|`TargetPotential`  | ✔ | ✔ |
|`TargetForce`      | ✔ | ✔ |

### Molecule

#### Properties

| property name                 | params                        | example                           |
|:------------------------------|:------------------------------|:----------------------------------|
|`MolAtoms`                     | --                            |`["H","H"]`                        |
|`MolAtomCoords`                | --                            |`[0.0 0.0 -0.375; 0.0 0.0 0.375]`  |
|`MolCharge`                    | --                            |`0`                                |
|`MolEnergyDataAvailable`       | --                            |`true`                             |
|`MolEnergyLevels`              | --                            |`[-0.590975, 0.174188, ...]`       |
|`MolHOMOIndex`                 | --                            |`1`                                |
|`MolHOMOEnergy`                | --                            |`-0.590975`                        |
|`MolWFATAvailableIndices`      | --                            |`Set([0])`                         |
|`MolWFATData`                  |`orbitIdx_relHOMO=0`           | --                                |
|`MolWFATStructureFactor_G`     |`orbitIdx_relHOMO=0,nξ,m,β,γ`  |`1.91369 - 1.46476e-33im`          |
|`MolAsympCoeffAvailableIndices`| --                            |`Set([0])`                         |
|`MolAsympCoeff`                |`orbitIdx_relHOMO=0`           | --                                |
|`MolMOADKStructureFactor_B`    |`orbitIdx_relHOMO=0,m_,β,γ`    |`2.21496 + 0.0im`                  |
|[`MolRotation`](@ref)          | --                            |`(0.0,0.0,0.0)`                    |
|`IonPotential`                 |`orbitIdx_relHOMO=0`           |`0.590975`                         |
|`AsympNuclCharge`              | --                            |`1`                                |
|`TargetName`                   | --                            |`"Hydrogen"`                       |
|`TargetPotential`              | --                            | --                                |
|`TargetForce`                  | --                            | --                                |

#### Methods

| method name                   | params                                |
|:------------------------------|:--------------------------------------|
|[`SetMolRotation`](@ref)       |`α,β,γ` or `(α,β,γ)`                   |
|`MolCalcEnergyData!`           |`orbitIdx_relHOMO=0, MCType, kwargs...`|
|[`MolCalcWFATData!`](@ref)     |`orbitIdx_relHOMO=0, MCType, kwargs...`|
|[`MolCalcAsympCoeff!`](@ref)   |`orbitIdx_relHOMO=0, MCType, kwargs...`|

## Hydrogen-Like Atom

A hydrogen-like atom has a potential of the form
```math
V(r) = -\frac{Z}{\sqrt{r^2+a}},
```
where ``Z`` is the nuclear charge number;
``a`` denotes the soft-core parameter, which is applied to avoid singularity of the potential and can be adjusted to fit the actual ionization potential of the atom (obtained by TDSE).

The hydrogen-like atom is implemented in the library as [`HydrogenLikeAtom`](@ref).
```@docs
Targets.HydrogenLikeAtom
```


## Single-Active-Electron (SAE) Atom

The single-active-electron (SAE) atom is an implementation of the empirical atomic SAE model potential proposed by Tong *et al.* [^Tong_2005]
The model potential of the SAE atom has the form
```math
V(r) = -\frac{Z + a_1 \mathrm{e}^{-b_1 r} + a_2 r \mathrm{e}^{-b_2 r} + a_3 \mathrm{e}^{-b_3 r}}{r},
```
where the ``a_i`` and ``b_i`` are tunable model potential parameters [^note].

The SAE atom is implemented in the library as [`SAEAtom`](@ref).
```@docs
Targets.SAEAtom
```

[^Tong_2005]: X. M. Tong *et al.*, Empirical Formula for Static Field Ionization Rates of Atoms and Molecules by Lasers in the Barrier-Suppression Regime. *J. Phys. B: At. Mol. Opt. Phys.* **38**, 2593–2600. DOI: [10.1088/0953-4075/38/15/001](https://dx.doi.org/10.1088/0953-4075/38/15/001)
[^note]: The symbols of the parameters are different from that in the original article. ``a_1, b_1, a_2, b_2, a_3, b_3`` correspond to ``a_1, a_2, a_3, a_4, a_5, a_6`` in the original article respectively.


## Preset Atom Library

The library provides some preset commonly-used atoms or atomic ions for convenience.

- [`HydrogenLikeAtom`](@ref) : H, He⁺, Li²⁺

- [`SAEAtom`](@ref) :          He, Ne, Ne⁺, Ne²⁺, Ar, Ar⁺, Ar²⁺, V, Ni, Kr, Kr⁺, Rb, Nb, Pd, Xe, Xe⁺, Ta

These atoms/ions can be obtained by invoking `Targets.**Atom()` (for neutral atoms) or `Targets.**#pAtom()` (for positive atomic ions), where `**` denotes the symbol of the nucleus and `#` denotes the positive charge the ion carries.

Example:

```@setup manual_targets
using SemiclassicalSFI
```
```@repl manual_targets
t1 = Targets.HAtom()
t2 = Targets.Xe1pAtom()
```


## Molecule

The `Molecule` object represents a generic molecule, which is implemented in the library as [`Molecule`](@ref).
The structure of `Molecule` is much more complex than that of atoms because the [Molecular ADK (MO-ADK) theory](@ref MOADK) and [Weak-Field Asymptotic Theory (WFAT)](@ref WFAT) features for molecular strong-field ionization require a number of coefficients, which are saved to files for convenience.

### Initialization, saving and loading

The `Molecule` object can be initialized either by providing necessary information of the molecule (mainly atoms, coordinates of the atoms and the charge of the molecule) or from external data (stored in the HDF5 format), cf. the documentation of [`Molecule`](@ref):

```@docs
Targets.Molecule
```

The `Molecule` object, after modification, can be manually saved to a HDF5 file via [`MolSaveDataAs`](@ref).

```@docs
Targets.MolSaveDataAs
```

### Molecular-SFI Data Preparation

To use the molecular strong-field ionization theories such as the MO-SFA, MO-ADK and WFAT to provide the intitial conditions of the electrons, the structure coefficients of the `Molecule` have to be calculated beforehand and stored in the object.
Cf. the documentation of [`MolCalcAsympCoeff!`](@ref) and [`MolCalcWFATData!`](@ref).

Evaluation of the structure coefficients depends on the external quantum chemistry packages.
The [`Targets.MolecularCalculators`](@ref) module undertakes the task of communication with the external quantum packages.
Currently only the [`PySCFMolecularCalculator`](@ref) is implemented.

!!! tip "Customized calculation parameters"
    When invoking `MolCalcAsympCoeff!` and `MolCalcWFATData!` to perform calculation of structure coefficients, customized calculation parameters can be passed to the `kwargs` of these methods.
    These parameters would be passed to the constructor method of the `MolecularCalculator` (e.g., the `basis` parameter of the [`PySCFMolecularCalculator`](@ref)), as well as the [`MolecularCalculators.calc_asymp_coeff`](@ref), [`MolecularCalculators.calc_WFAT_data`](@ref) methods.
    Refer to their documentation below for more information.

```@docs
Targets.MolCalcAsympCoeff!
Targets.MolCalcWFATData!
```

```@docs
Targets.MolecularCalculators
Targets.MolecularCalculators.PySCFMolecularCalculator
```

```@docs
Targets.MolecularCalculators.calc_asymp_coeff
Targets.MolecularCalculators.calc_WFAT_data
```

### Molecule's Orientation

The molecule's orientation is described by a set of Euler angles (``z-y'-z''`` convention), which defines a rotational transformation from the molecular frame (MF) to the lab frame (LF).
This property of `Molecule` is NOT included in the saved file and thus needs to be specified each time upon initialization of the `Molecule` object from external files.

!!! note "Note"
    Here the three Euler angles `(α,β,γ)` that describe the `Molecule`'s orientation are completely different from that of the Euler angles `(β',γ')` in the [WFAT](@ref WFAT) and [MO-ADK](@ref MOADK) theory.
    These theories' "lab frame" is chosen for convenience of theoretical formulation, where the electric field is assumed to be static, pointing towards the ``+z`` direction,
    and has no relation with the lab frame mentioned above.

The orientation of the molecule can be obtained and set via the [`MolRotation`](@ref) and [`SetMolRotation`](@ref) methods.

```@docs
Targets.MolRotation
Targets.SetMolRotation
```
