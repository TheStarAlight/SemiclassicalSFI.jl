# [Theory: Initial Conditions](@id theory_init_cond)

*This section reviews commonly-used theories used to provide initial conditions in the trajectory simulations.*

A number of theories can be adapted to provide initial conditions of the classical electrons in the trajectory simulation scheme.
The initial condition usually consists of three properties:

- Initial position ``\bm{r}_0`` (i.e., the tunneling exit position);
- Initial momentum ``\bm{p}_0``, we note that in the trajectory simulation schemes, initial momenta are usually denoted using ``\bm{k}_0``;
- The corresponding ionization probability ``W`` carried by each electron sample, depending on the time-dependent laser field and the properties of the target atoms/molecules.

In the following we will give a brief review on the available theories we implemented in *SemiclassicalSFI.jl*.
Atomic units (a.u.) are used throughout unless stated otherwise.

```@contents
Pages = ["theory1_initial_conditions.md"]
```

## [Strong-Field Approximation (SFA)](@id SFA)

The Strong-Field Approximation (SFA) [^Popruzhenko_2014] is originated from the Keldysh theory of strong-field ionization.
Compared with the pertubative methods and adiabatic tunneling theories, the SFA is able to predict both the multi-photon and the tunneling process during the laser-atom interaction, as well as high-order non-pertubative phenomenona such as the above-threshold ionization (ATI) because it fully includes the non-adiabatic effect of the laser-atom interaction.
The broad scope of SFA has contributed to its widespread application in theoretical investigations of strong-field ionization.

Considering the electron evolving under a combined field of the Coulomb field ``V(\bm{r})`` and the laser field ``\bm{F}(t)=-\partial_t \bm{A}(t)``, under the length gauge (LG), its Hamiltonian reads
```math
H = \frac12 \bm{p}^2 + V(\bm{r}) + \bm{F}(t)\cdot\bm{r}.
```
Denoting ``\ket{\Psi_0} = \ket{\psi_0} \mathrm{e}^{\mathrm{i}I_{\mathrm{p}}t}`` as the unperturbed initial state with ionization potential of ``I_{\mathrm{p}}``, ``\ket{\Psi_{\bm{p}}}`` as the continuum state of momentum ``\bm{p}``, and
```math
U(t_{\mathrm{f}},t_0) = \exp \left[ -\mathrm{i} \int_{t_0}^{t_{\mathrm{f}}} H(\tau) \mathrm{d}\tau \right]
```
the time-evolution operator, the transition amplitude between the initial state (at ``t_0``) and the final state of momentum ``\bm{p}`` (at ``t_{\mathrm{f}}``) is written as
```math
M_{\bm{p}} = \braket{ \Psi_{\bm{p}} | U(t_{\mathrm{f}},t_0) | \Psi_0 }.
```

Here lies the key idea of SFA: when the influence of the Coulomb field to the ionized electrons is weak compared with that of the external laser field, we may neglect the influence of the Coulomb field in the expression of ``M_{\bm{p}}`` by replacing the time-evolution operator with a Coulomb-free one ``U_{\mathrm{f}}``, and meanwhile replacing the continuum state with the Volkov state ``\ket{\Psi^{\mathrm{V}}_{\bm{p}}}`` which represents a free electron evolving under the same laser field:
```math
M_{\bm{p}} = \braket{ \Psi^{\mathrm{V}}_{\bm{p}} | U_{\mathrm{f}}(t_{\mathrm{f}},t_{\mathrm{0}}) | \Psi_0 },
```
where the Volkov state under the LG is the product of a plane wave and a phase factor:
```math
\ket{ \Psi^{\mathrm{V}}_{\bm{p}} } = \ket{ \bm{p}+\bm{A}(t) } \mathrm{e}^{\mathrm{i}S_{\bm{p}}(t)},
```
and the phase has the expression:
```math
S_{\bm{p}}(t) = - \int^{t} \frac12 [\bm{p}+\bm{A}(\tau)]^2 \mathrm{d}\tau.
```
In this way the ``M_{\bm{p}}`` is expressed as
```math
M_{\bm{p}} = -\mathrm{i} \int_{t_{\mathrm{0}}}^{t_{\mathrm{f}}} \braket{ \bm{p}+\bm{A}(\tau) | \bm{F}(\tau)\cdot\bm{r} | \psi_0 } \mathrm{e}^{-\mathrm{i}\tilde{S}_{\bm{p}}(\tau)} \mathrm{d}\tau,
```
and we note that here we have extracted the phase factor of ``\ket{\Psi_0}`` and combined it with the former ``\mathrm{e}^{-\mathrm{i}S_{\bm{p}}(t)}``, giving
```math
\tilde{S}_{\bm{p}}(t) = - \int^{t} \left[ \frac12 [\bm{p}+\bm{A}(\tau)]^2 + I_{\mathrm{p}} \right] \mathrm{d}\tau.
```

An additional saddle-point approximation facilitates preparation of initial conditions of the electron trajectories.
The variation of ``\tilde{S}_{\bm{p}}(t)`` is much more sensitive than that of ``\braket{ \bm{p}+\bm{A}(t) | \bm{F}(t)\cdot\bm{r} | \psi_0 }`` as ``t`` varies, which leads to a fact that the whole integrand in our latest expression of ``M_{\bm{p}}`` oscillates in its complex phase and its values cancel out each other in most cases, except when the variation of the phase ``\tilde{S}_{\bm{p}}(t)`` becomes stable, i.e., at the saddle points of ``\tilde{S}_{\bm{p}}(t)``. The saddle points ``t_{\mathrm{s}}=t_{\mathrm{r}}+\mathrm{i}t_{\mathrm{i}}`` are the zeroes of the derivative of the complex function ``\tilde{S}_{\bm{p}}(t)``, which satisfy
```math
-\partial_t \tilde{S}_{\bm{p}}(t) |_{t=t_{\mathrm{s}}} = \frac12 [\bm{p}+\bm{A}(t_{\mathrm{s}})]^2 + I_{\mathrm{p}} = 0.
```
The integral can be approximated by a summation over the saddle points:
```math
M_{\bm{p}} \approx \sum_{t_{\mathrm{s}}} P_{\bm{p}}(t_{\mathrm{s}}) \mathrm{e}^{-\mathrm{i}\tilde{S}_{\bm{p}}(t_{\mathrm{s}})},
```
where ``P_{\bm{p}}(t_{\mathrm{s}})`` denotes the prefactor.
Here we use a modified version of SFA which takes account of the Coulomb tail [^Kjeldsen_2006] [^Milosevic_2006], which gives the prefactor
```math
P_{\bm{p}}(t_{\mathrm{s}}) = \{ [\bm{p}+\bm{A}(t_{\mathrm{s}})] \cdot \bm{F}(t_{\mathrm{s}}) \}^{-\alpha/2},
```
where ``\alpha = 1+Z/\sqrt{2I_{\mathrm{p}}}``, and ``Z`` is the asymptotic charge of the nucleus.
The phase ``\tilde{S}_{\bm{p}}(t_{\mathrm{s}})`` is obtained by solving the integral
```math
\begin{aligned}
    \tilde{S}_{\bm{p}}(t_{\mathrm{s}})
    &= \int_{t_{\mathrm{s}}}^{\infty} \left[ \frac12 [\bm{p}+\bm{A}(\tau)]^2 + I_{\mathrm{p}} \right] \mathrm{d}\tau \\
    &= \left( \int_{t_{\mathrm{s}}}^{t_{\mathrm{r}}} + \int_{t_{\mathrm{r}}}^{\infty} \right) \left[ \frac12 [\bm{p}+\bm{A}(\tau)]^2 + I_{\mathrm{p}} \right] \mathrm{d}\tau \\
    &= \Phi_{\mathrm{tun}} + \Phi_{\mathrm{traj}},
\end{aligned}
```
where ``\Phi_{\mathrm{tun}}`` represents the complex phase accumulation during the tunneling process, whose real part denotes the quantum phase, while its imaginary part, is related to the ionization probability; the ``\Phi_{\mathrm{traj}}``, is the phase accumulation during the electron trajectory motion in the continuum.

The SFA provides the final momentum distribution, while the trajectory simulation requires initial conditions of the eletrons.
To utilize the SFA to give initial conditions, we suppose that the classical electron is ejected at time ``t_{\mathrm{r}}`` at tunneling exit ``\bm{r}_0`` with momentum ``\bm{k}_0``.
The initial momentum ``\bm{k}_0``, neglecting the Coulomb potential, is related to the final momentum ``\bm{p}`` through
```math
\bm{p} = \bm{k}_0 - \int_{t_{\mathrm{r}}}^{\infty} \bm{F}(\tau) \mathrm{d}\tau = \bm{k}_0 - \bm{A}(t_{\mathrm{r}}).
```
The initial position ``\bm{r}_0``, i.e., the tunneling exit, is found by constructing a quantum tunneling trajectory.
The beginning of the trajectory, i.e., the tunneling entrance, has a real part of zero; the electron tunnels through the barrier during the time interval ``t_{\mathrm{s}}`` to ``t_{\mathrm{r}}`` and emerges as a classical electron at the tunneling exit with real position and momentum.
In this way we obtain the expression of the initial position:
```math
\bm{r}_0 = \mathrm{Re} \int_{t_{\mathrm{s}}}^{t_{\mathrm{r}}} \bm{A}(\tau) \mathrm{d}\tau = \mathrm{Im} \int_{0}^{t_{\mathrm{i}}} \bm{A}(t_{\mathrm{r}}+\mathrm{i}\tau) \mathrm{d}\tau.
```
The probablity density (in the final momentum space) carried by the electron sample is
```math
\mathrm{d}W/\mathrm{d}\bm{p} = \lvert P_{\bm{p}}(t_{\mathrm{s}}) \rvert^2 \exp(2\ \mathrm{Im}\ \Phi_{\mathrm{tun}}).
```

We note that the ionization probability is expressed in the coordinate of final momentum ``(p_x,p_y,p_z)``.
However, in the trajectory simulation, the initial electrons are sampled in the coordinate ``(t_{\mathrm{r}},k_d,k_z)``, where ``k_d`` denotes the initial momentum's component in the ``xy`` plane (which is perpendicular to the electric field).
Thus, adding a Jacobian in the prefix of the ionization probability is required if we sample the initial electrons within such coordinate, the transformed expression reads
```math
\mathrm{d}W/\mathrm{d}\bm{k}_\perp \mathrm{d}t_{\mathrm{r}} = \bm{J}(k_d,t_{\mathrm{r}}) \lvert P_{\bm{p}}(t_{\mathrm{s}}) \rvert^2 \exp(2\ \mathrm{Im}\ \Phi_{\mathrm{tun}}),
```
where the Jacobian is
```math
\bm{J}(k_d,t_{\mathrm{r}}) = \frac{\partial(p_x,p_y)}{\partial(k_d,t_{\mathrm{r}})} =
\begin{vmatrix}
    \partial p_x/\partial k_d & \partial p_x/\partial t_{\mathrm{r}} \\
    \partial p_y/\partial k_d & \partial p_y/\partial t_{\mathrm{r}}
\end{vmatrix}
.
```

[^Popruzhenko_2014]: S. V. Popruzhenko, Keldysh Theory of Strong Field Ionization: History, Applications, Difficulties and Perspectives. *J. Phys. B: At. Mol. Opt. Phys.* **47**, 204001 (2014). DOI:[10.1088/0953-4075/47/20/204001](https://dx.doi.org/10.1088/0953-4075/47/20/204001)
[^Kjeldsen_2006]: T. K. Kjeldsen *et al.*, Strong-Field Ionization of Atoms and Molecules: The Two-Term Saddle-Point Method. *Phys. Rev. A* **74**, 023407 (2006). DOI:[10.1103/PhysRevA.74.023407](https://dx.doi.org/10.1103/PhysRevA.74.023407)
[^Milosevic_2006]: D. B. Milošević *et al.*, Above-Threshold Ionization by Few-Cycle Pulses. *J. Phys. B: At. Mol. Opt. Phys.* **39**, R203–R262 (2006). DOI: [10.1088/0953-4075/39/14/R01](https://dx.doi.org/10.1088/0953-4075/39/14/R01)


## [SFA with Adiabatic Expansion (SFA-AE)](@id SFAAE)

For small Keldysh parameter ``\gamma``, the non-adiabatic effect is not significant, thus an adiabatic expansion scheme can be carried out to develop a modified theory based on the SFA, which is named after the SFA with adiabatic expansion (SFA-AE) [^Ni_2018]. It includes the non-adiabatic effect to a large extent and is capable of giving similar results compared with that given by the SFA under small Keldysh parameters.

The SFA-AE is applicable when the Keldysh parameter is small or the non-adiabatic effect is insignificant, and we recall that in the SFA there is a corresponding quantity ``t_{\mathrm{i}}`` which quantifies the non-adiabacity of tunneling.
For small ``t_{\mathrm{i}}``, we expand the vector potential ``\bm{A}(t_{\mathrm{r}} + \mathrm{i}t_{\mathrm{i}})`` at ``t_{\mathrm{r}}``, up to the second order of ``t_{\mathrm{i}}``:
```math
\bm{A}(t_{\mathrm{r}} + \mathrm{i}t_{\mathrm{i}}) = \bm{A}(t_{\mathrm{r}}) - \mathrm{i}t_{\mathrm{i}}\bm{F}(t_{\mathrm{r}}) + \frac12 t_{\mathrm{i}}^2 \bm{F}'(t_{\mathrm{r}}) + o(t_{\mathrm{i}}^2).
```
Inserting the above expression into the saddle-point equation in the SFA gives
```math
\bm{k}(t_{\mathrm{r}}) \cdot \bm{F}(t_{\mathrm{r}}) = 0,
```
and
```math
t_{\mathrm{i}} = \sqrt{\frac{k^2(t_{\mathrm{r}})+2I_{\mathrm{p}}}{F^2(t_{\mathrm{r}})-\bm{k}(t_{\mathrm{r}}) \cdot \bm{F}'(t_{\mathrm{r}})}}.
```

The ``\mathrm{Im}\ \Phi_{\mathrm{tun}}`` related to the ionization rate, in the SFA-AE, is
```math
\mathrm{Im}\ \Phi_{\mathrm{tun}} \approx -\frac13 \frac{[k^2(t_{\mathrm{r}})+2I_{\mathrm{p}}]^{3/2}}{\sqrt{F^2(t_{\mathrm{r}})-\bm{k}(t_{\mathrm{r}}) \cdot \bm{F}'(t_{\mathrm{r}})}},
```
and we obtain
```math
\begin{aligned}
    \mathrm{d}W/\mathrm{d}\bm{p}
    &= \lvert P_{\bm{p}}(t_{\mathrm{s}}) \rvert^2 \exp(2\ \mathrm{Im}\ \Phi_{\mathrm{tun}}) \\
    &= \lvert P_{\bm{p}}(t_{\mathrm{s}}) \rvert^2 \exp \left[ -\frac23 \frac{(k_\perp^2+2I_{\mathrm{p}})^{3/2}}{\sqrt{F^2-\bm{k}_\perp \cdot \bm{F}'}} \right],
\end{aligned}
```
where the ``\bm{k}_{\perp}`` denotes the transverse momentum at the tunneling exit, which is actually equivalent to ``\bm{k}(t_{\mathrm{r}})`` in the SFA-AE because the above saddle-point equation requires ``\bm{k}(t_{\mathrm{r}}) \cdot \bm{F}(t_{\mathrm{r}}) = 0``. We note that the initial momentum, ``\bm{k}_0``, is exactly ``\bm{k}_{\perp}``.
The prefactor ``P_{\bm{p}}(t_{\mathrm{s}})`` which included the Coulomb tail correction, in the SFA-AE, has the expression [^Frolov_2017]:
```math
P_{\bm{p}}(t_{\mathrm{s}}) = \left[ (k_\perp^2+2I_{\mathrm{p}})(F^2-\bm{k}_\perp \cdot \bm{F}') \right]^{-\alpha/4}.
```

The initial position has the expression
```math
\bm{r}_0 = \mathrm{Im} \int_{0}^{t_{\mathrm{i}}} \bm{A}(t_{\mathrm{r}}+\mathrm{i}\tau) \mathrm{d}\tau = -\frac{\bm{F}}{2} \frac{k_\perp^2+2I_{\mathrm{p}}}{F^2-\bm{k}_\perp \cdot \bm{F}'}.
```

[^Ni_2018]: H. Ni *et al.*, Tunneling Criteria and a Nonadiabatic Term for Strong-Field Ionization. *Phys. Rev. A* **98**, 013411 (2018). DOI:[10.1103/PhysRevA.98.013411](https://dx.doi.org/10.1103/PhysRevA.98.013411)
[^Frolov_2017]: Frolov *et al.*, Adiabatic-Limit Coulomb Factors for Photoelectron and High-Order-Harmonic Spectra. *Phys. Rev. A* **96**, 023406 (2017). DOI:[10.1103/PhysRevA.96.023406](https://dx.doi.org/10.1103/PhysRevA.96.023406)



## [Ammosov-Delone-Krainov (ADK)](@id ADK)

The Ammosov-Delone-Krainov (ADK) theory [^Ammosov_1986] [^Delone_1998] is used to study the adiabatic tunneling in the strong-field ionization, and is, in a sense, the adiabatic limit of the SFA.

In the adiabatic limit, the laser field can be treated as static, thus we have ``\bm{F}'(t)=\bm{0}`` (higher order derivatives of ``\bm{F}(t)`` remains zero as well).
Substuting it into the expressions of SFA-AE yields the ADK rate
```math
\mathrm{d}W/\mathrm{d}\bm{p} = \lvert P_{\bm{p}}(t_{\mathrm{s}}) \rvert^2 \exp \left[ -\frac23 \frac{[k_\perp^2+2I_{\mathrm{p}}]^{3/2}}{F} \right],
```
where the prefactor reads
```math
\lvert P_{\bm{p}}(t_{\mathrm{s}}) \rvert^2 = \left[ (k_\perp^2+2I_{\mathrm{p}})F^2\right]^{-\alpha/2}.
```
The tunneling exit position can be obtained in the same approach:
```math
\bm{r}_0 = \mathrm{Im} \int_{0}^{t_{\mathrm{i}}} \bm{A}(t_{\mathrm{r}}+\mathrm{i}\tau) \mathrm{d}\tau = -\frac{\bm{F}}{2} \frac{k_\perp^2+2I_{\mathrm{p}}}{F^2}.
```

[^Ammosov_1986]: M. V. Ammosov *et al.*, Tunnel Ionization of Complex Atoms and of Atomic Ions in an Alternating Electromagnetic Field. *Sov. Phys. JETP* **64**, 1191 (1986).
[^Delone_1998]: N. B. Delone *et al.*, Tunneling and Barrier-Suppression Ionization of Atoms and Ions in a Laser Radiation Field. *Phys.-Usp.* **41**, 469–485. DOI: [10.1070/PU1998v041n05ABEH000393](https://dx.doi.org/10.1070/PU1998v041n05ABEH000393)

### [Tunneling Exit Methods for Atomic ADK](@id tun_exit_atomic_adk)

The tunneling exit, i.e., the initial position ``\bm{r}_0`` of the classical trajectories, remains a controversial problem.
Below we briefly review three methods for determination of the tunneling exit, which apply for the ADK method.

(1) The ``I_{\mathrm{p}}/F`` model.

The ``I_{\mathrm{p}}/F`` model is initially derived from the SFA theory and is the adiabatic limit of it.
In this model the distance ``r_0`` between the tunneling exit and the nucleus is simply expressed as the quotient of the ionization potential ``I_{\mathrm{p}}`` and the electric field strength ``F``.
Taking account of the initial kinetic energy ``k_\perp^2/2``, ``r_0`` reads
```math
r_0 = \frac{I_{\mathrm{p}}+k_\perp^2/2}{F}.
```
However, the shortcoming of the ``I_{\mathrm{p}}/F`` model is obvious: derived from the SFA which assumes a short-range potential (i.e., a delta potential), this model neglects the complexity of the Coulomb characteristics of the potential and thus the conclusion lacks persuasiveness.

(2) The field-direction model (FDM).

The field-direction model (FDM) determines the tunneling exit by treating the problem as a one-dimensional tunneling problem: the tunneling exit is found in the one-dimensional cut of the combined potential of the laser and nucleus along the field direction.
We assume that the field points towards the ``-z`` direction, and the tunneling exit ``r_0`` is found from
```math
V(r_0 \bm{e}_z) - F r_0 = - I_{\mathrm{p}}.
```
For sufficiently large ``r_0``, the potential ``V(r)`` behaves as ``-Z/r``, and the tunneling exit is expressed as
```math
r_0 = \frac{I_{\mathrm{p}}+\sqrt{I_{\mathrm{p}}^2-4FZ}}{2F}.
```
The FDM approach takes the Coulomb field into account, however, not in such a scientific way because this problem is actually not independent of the transverse dimensions.

(3) The parabolic-coordinate model.

For a static electric field towards the ``-z`` direction, the parabolic-coordinate model introduces the parabolic coordinate ``\xi = r+z``, ``\eta = r-z``, ``\phi = \tan^{-1}{y/x}``, in which the Schrödinger equation becomes separable in the three dimensions.
The parabolic-coordinate approach gives expression of the tunneling exit:
```math
r_0 = \frac{I_{\mathrm{p}}+\sqrt{I_{\mathrm{p}}^2-4F \left[ Z-(1+|m|)\sqrt{I_{\mathrm{p}}/2} \right] }}{2F},
```
where ``m`` is the magnetic quantum number along the ``z`` axis (due to the isotropic angular distribution of the atomic wavefunction, we set ``m=0``).


## [Molecular ADK (MO-ADK)](@id MOADK)

The molecular ADK (MO-ADK) theory generalizes the original ADK theory by extending the application scope from atomic to simple linear molecules [^Tong_2002].

In the MO-ADK theory, the wavefunction of a linear molecule's ionizing orbital behaves asymptotically as
```math
\psi_0(\bm{r}) \sim \sum_{lm} F_{lm}(r) Y_{lm}(\theta,\phi)
```
in the molecular frame (MF) when ``r\rightarrow\infty``.
Assigning ``\kappa=\sqrt{2I_{\mathrm{p}}}``, the ``F_{lm}(r)`` has the following asymptotic behavior when ``r\rightarrow\infty``:
```math
F_{lm}(r) \sim C_{lm} r^{Z/\kappa-1} \mathrm{e}^{-\kappa r}.
```
In numerical implementation we obtain the parameters ``C_{lm}`` by fitting the above expression [^Zhang_2015], and the ``F_{lm}(r)`` is found by the spherical-harmonic expansion of the wavefunction:
```math
F_{lm}(r) = \int \mathrm{d}\bm{\Omega} Y_{lm}^{*}(\bm{\Omega}) \psi_0(\bm{r}).
```

[^Zhang_2015]: Zhang, B. *et al.*, SLIMP: Strong Laser Interaction Model Package for Atoms and Molecules. *Comp. Phys. Comm.* **192**, 330–341 (2015). DOI: [10.1016/j.cpc.2015.02.031](https://dx.doi.org/10.1016/j.cpc.2015.02.031)

We assume the electric field is pointing towards the ``z`` axis in the laboratory frame (LF).
The angle-dependent tunneling ionization rate in the MO-ADK theory reads
```math
\Gamma(\beta,\gamma) = \mathrm{d}W/\mathrm{d}t = \sum_{m'} \frac{|B_{m'}(\beta,\gamma)|^2}{2^{|m'|}|m'|!} \kappa^{-|m'|} \left(\frac{2\kappa^2}{F}\right)^{2Z/\kappa-|m'|-1} \mathrm{e}^{-2\kappa^3/3F},
```
where the molecule's orientation is described using a set of Euler angles ``\hat{\bm{R}} = (\alpha,\beta,\gamma)`` (``z-y'-z''`` convention), which represents the rotational transformation from the MF to the LF; ``B_{m'}(\beta,\gamma)`` are the structural parameters which depend on the molecule's orbital wavefunction (here we omitted the ``\alpha`` dependence because the structural parameters are independent of ``\alpha``).
The structural parameters ``B_{m'}(\beta,\gamma)`` have the following expression:
```math
B_{m'}(\beta,\gamma) = \sum_{lm} C_{lm} d_{m' m}^{l}(\beta) \mathrm{e}^{-\mathrm{i}m\gamma} Q_{l m'},
```
with ``d_{m' m}^{l}(\beta)`` being the Wigner-``d`` rotation matrix, and
```math
Q_{l m'} = (-1)^{(m'+|m'|)/2} \sqrt{\frac{(2l+1)(l+|m'|)!}{2(l-|m'|)!}}.
```

To ultilize the MO-ADK theory to provide the initial conditions in the trajectory simulation, we simply adopt the result of the atomic ADK theory:
```math
\bm{r}_0 = -\frac{\bm{F}}{2} \frac{k_\perp^2+2I_{\mathrm{p}}}{F^2}.
```
As for the ionization probability, we include the influence of the initial kinetic energy ``k_\perp^2/2`` by replacing the ``\kappa=\sqrt{2I_{\mathrm{p}}}`` with ``\kappa'(k_\perp)=\sqrt{2I_{\mathrm{p}}+k_\perp^2}`` in the exponential term of the ionization probability in the MO-ADK theory, giving
```math
\mathrm{d}W/\mathrm{d}\bm{k}_\perp \mathrm{d}t = \sum_{m'} \frac{|B_{m'}(\beta,\gamma)|^2}{2^{|m'|}|m'|!} \kappa^{-|m'|} \left(\frac{2\kappa^2}{F}\right)^{2Z/\kappa-|m'|-1} \mathrm{e}^{-2\kappa'^3(k_\perp)/3F}.
```

[^Tong_2002]: X. M. Tong *et al.*, Theory of Molecular Tunneling Ionization. *Phys. Rev. A* **66**, 033402 (2002). DOI: [10.1103/PhysRevA.66.033402](https://dx.doi.org/10.1103/PhysRevA.66.033402)



## [Weak-Field Asymptotic Theory (WFAT)](@id WFAT)

The weak-field asymptotic theory (WFAT) generalizes the tunneling ionization from isotropic atomic potentials to arbitrary molecular potentials [^Tolstikhin_2011]. Compared with the MO-ADK theory, the WFAT accounts for the influence of the molecules' permanent dipole moment, and is applicable for complex molecules other than simple linear molecules.

The formulation of the WFAT is based on the expansion in the parabolic coordinates.
The total ionization rate ``\Gamma(\beta,\gamma) = \mathrm{d}W/\mathrm{d}t``, is split into different parabolic channels:
```math
\Gamma(\beta,\gamma) = \sum_\nu \Gamma_\nu(\beta,\gamma),
```
where ``\Gamma_\nu(\beta,\gamma)`` are partial rates of parabolic quantum number indices ``\nu=(n_\xi,m)``, and ``n_\xi=0,1,2,\cdots``, ``m=0,\pm 1,\pm 2,\cdots``.
In the leading-order approximation of the WFAT, the partial rates can be separated into two factors, namely the structural part ``|G_\nu(\beta,\gamma)|^2`` and the field part ``W_\nu(F)``:
```math
\Gamma_\nu(\beta,\gamma) = |G_\nu(\beta,\gamma)|^2 W_\nu(F).
```
The field factor is expressed as
```math
W_\nu(F) = \frac{\kappa}{2} \left(\frac{4\kappa^2}{F}\right)^{2Z/\kappa-2n_\xi-|m|-1} \mathrm{e}^{-2\kappa^3/3F}.
```
The structure factor, in the integral representation of the WFAT [^Dnestryan_2018], is given as an integral:
```math
G_\nu (\beta,\gamma) = \mathrm{e}^{-\kappa\mu_z} \int \Omega_\nu^* \left(\hat{\bm{R}}^{-1} \bm{r}\right) \hat{V}_{\mathrm{c}}(\bm{r}) \psi_0(\bm{r}) \mathrm{d} \bm{r},
```
where ``\psi_0`` is the wavefunction of the ionizing orbital,
```math
\bm{\mu} = \int \psi_0^*(\bm{r}) \bm{r} \psi_0(\bm{r}) \mathrm{d} \bm{r}
```
denotes the orbital dipole moment in the LF, with ``\mu_z`` being its component along the field direction;
```math
\Omega_\nu(\bm{r}) = \sum_{l=|m|}^{\infty} \Omega^\nu_{lm}(\bm{r}) = \sum_{l=|m|}^{\infty} R_l^\nu(r) Y_{lm}(\theta, \phi)
```
is a reference function which can be expanded into spherical harmonics, its radial part is expressed as
```math
R_l^\nu(r)=\omega_l^\nu \ (\kappa r)^l \ \mathrm{e}^{-\kappa r} \ \mathrm{M}(l+1-Z/\kappa, 2l+2, 2 \kappa r),
```
with ``\mathrm{M}(a,b,x)`` being the confluent hyper-geometric function and
```math
\begin{aligned}
    \omega_l^\nu = & \      (-1)^{l+(|m|-m)/2+1}\ 2^{l+3/2}\ \kappa^{Z/\kappa-(|m|+1)/2-n_\xi} \\
                   & \times \sqrt{(2l+1)(l+m)!(l-m)!(|m|+n_\xi)!n_\xi!}\ \frac{l!}{(2l+1)!} \\
                   & \times \!\!\!\!\!\! \sum_{k=0}^{\min{(n_\xi,l-|m|)}} \!\!\!\!\!\!\!\!\!\! \frac{\Gamma(l+1-Z/\kappa+n_\xi-k)}{k!(l-k)!(|m|+k)!(l-|m|-k)!(n_\xi-k)!}
\end{aligned}
```
the normalization coefficient;
``\hat{V}_{\mathrm{c}}(\bm{r})=\hat{V}(\bm{r})+Z/r`` is the core potential with the Coulomb tail removed, where ``Z`` is the asymptotic charge of the residual ion.

The effective potential ``\hat{V}(\bm{r})`` describes the interaction between the ionizing electron and the residual parent ion.
Under the framework of the Hartree-Fock method, the effective potential consists of three parts, namely the nuclear Coulomb potential (``V_{\mathrm{nuc}}``), the direct (``V_{\mathrm{d}}``) and exchange (``V_{\mathrm{ex}}``) parts of inter-electron interactions:
```math
\hat{V}(\bm{r}) = V_{\mathrm{nuc}}(\bm{r}) + V_{\mathrm{d}}(\bm{r}) + \hat{V}_{\mathrm{ex}}(\bm{r}),
```
and
```math
\begin{aligned}
    V_{\mathrm{nuc}}(\bm{r}) &= -\sum_{A=1}^{N_{\mathrm{atm}}} \frac{Z_A}{\left|\bm{r}-\bm{R}_A\right|}, \\
    V_{\mathrm{d}}(\bm{r}) &= \sum_{i=1}^N \int \frac{\psi_i^*(\bm{r}') \psi_i(\bm{r}')}{|\bm{r}-\bm{r}'|} \mathrm{d} \bm{r}', \\
    \hat{V}_{\mathrm{ex}}(\bm{r}) \psi_0(\bm{r}) &= -\sum_{i=1}^N \psi_i(\bm{r}) \int \frac{\psi_i^*(\bm{r}') \psi_0(\bm{r}')}{|\bm{r}-\bm{r}'|} \braket{\sigma_i | \sigma_0} \mathrm{d} \bm{r}',
\end{aligned}
```
where ``N`` is the number of electrons, ``N_{\mathrm{atm}}`` is the number of nuclei, ``\psi_i(\bm{r})`` and ``\sigma_i`` denote the molecular orbital and the spin state of the electron of index ``i`` (``\braket{\sigma_i|\sigma_j}=1`` if electrons of index ``i`` and ``j`` have the same spin, and ``\braket{\sigma_i|\sigma_j}=0`` otherwise), ``Z_A`` and ``\bm{R}_A`` are the nuclear charge and position of atom of index ``A``.

As the WFAT provides only the ionization rate ``\Gamma = \mathrm{d}W/\mathrm{d}t`` as the MO-ADK does, we adopt the same procedure as we did in the MO-ADK theory to provide the initial conditions for the trajectory simulation.
The initial position ``\bm{r}_0`` is the same as that in the MO-ADK theory, and the ionization rate reads
```math
\mathrm{d}W/\mathrm{d}\bm{k}_\perp \mathrm{d}t = \sum_\nu |G_\nu(\beta,\gamma)|^2 \cdot \frac{\kappa}{2} \left(\frac{4\kappa^2}{F}\right)^{2Z/\kappa-2n_\xi-|m|-1} \mathrm{e}^{-2\kappa'^3(k_\perp)/3F}.
```

[^Tolstikhin_2011]: O. I. Tolstikhin *et al.*, Theory of Tunneling Ionization of Molecules: Weak-Field Asymptotics Including Dipole Effects. *Phys. Rev. A* **84**, 053423 (2011). DOI: [10.1103/PhysRevA.84.053423](https://dx.doi.org/10.1103/PhysRevA.84.053423)
[^Dnestryan_2018]: A. I. Dnestryan *et al.*, Structure Factors for Tunneling Ionization Rates of Molecules: General Grid-Based Methodology and Convergence Studies. *J. Chem. Phys.* **149**, 164107. DOI: [10.1063/1.5046902](https://dx.doi.org/10.1063/1.5046902)

