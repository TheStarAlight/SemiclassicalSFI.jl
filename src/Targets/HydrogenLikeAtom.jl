using StaticArrays

"""
```
struct HydrogenLikeAtom <: SAEAtomBase
```

Represents a hydrogen-like atom.

An instance of `HydrogenLikeAtom` can be initialized via the constructor method:
```julia
HydrogenLikeAtom(Ip, Z, soft_core=1.0, name="[NA]")
```

# Example:
```jldoctest
julia> t = Targets.HydrogenLikeAtom(Ip=0.5, Z=1.0, soft_core=1.0, name="H")
[HydrogenLikeAtom] Atom H, Ip=0.5, Z=1.0, SoftCore=1.0
```
"""
struct HydrogenLikeAtom <: SAEAtomBase
    "Ionization potential of the atom."
    Ip;
    "Asymptotic charge of the inner nucleus."
    nucl_charge;
    "Soft core parameter of the atom."
    soft_core;
    "Name of the atom."
    name::String;
    "Initializes a new instance of HydrogenLikeAtom."
    HydrogenLikeAtom( Ip, Z, soft_core=0.2, name="[NA]") = new(Ip, Z, soft_core, name)
    HydrogenLikeAtom(;Ip, Z, soft_core=0.2, name="[NA]") = new(Ip, Z, soft_core, name)
end

"Gets the ionization potential of the atom."
IonPotential(t::HydrogenLikeAtom) = t.Ip
"Gets the asymptotic nuclear charge of the atom."
AsympNuclCharge(t::HydrogenLikeAtom) = t.nucl_charge
"Gets the soft core parameter of the atom."
SoftCore(t::HydrogenLikeAtom) = t.soft_core
"Gets the name of the atom."
TargetName(t::HydrogenLikeAtom) = t.name
"Gets the potential function of the atom."
TargetPotential(t::HydrogenLikeAtom) = (x,y,z) -> -t.nucl_charge*(x^2+y^2+z^2+t.soft_core)^(-0.5)
"Gets the force exerted on the electron from the atom (which is the neg-grad of potential)."
TargetForce(t::HydrogenLikeAtom) = (x,y,z) -> -t.nucl_charge*(x^2+y^2+z^2+t.soft_core)^(-1.5) .* (x,y,z)
"Gets the trajectory function according to given parameter."
function TrajectoryFunction(t::HydrogenLikeAtom, laserFx::Function, laserFy::Function, phase_method::Symbol, non_dipole::Bool; kwargs...)
    Z  = t.nucl_charge
    Ip = t.Ip
    soft_core = t.soft_core
    # including external function call is infeasible in GPU, thus the external targetF & targetP are replaced by pure Coulomb ones.
    return if ! non_dipole
        if phase_method == :CTMC
            function traj_dipole_ctmc(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                tFx, tFy, tFz = -Z*(u[1]^2+u[2]^2+u[3]^2+soft_core)^(-1.5) .* (u[1],u[2],u[3])
                du1 = u[4]
                du2 = u[5]
                du3 = u[6]
                du4 = tFx-laserFx(t)
                du5 = tFy-laserFy(t)
                du6 = tFz
                @SVector [du1,du2,du3,du4,du5,du6]
            end
        elseif phase_method == :QTMC
            function traj_dipole_qtmc(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                tFx, tFy, tFz = -Z*(u[1]^2+u[2]^2+u[3]^2+soft_core)^(-1.5) .* (u[1],u[2],u[3])
                du1 = u[4]
                du2 = u[5]
                du3 = u[6]
                du4 = tFx-laserFx(t)
                du5 = tFy-laserFy(t)
                du6 = tFz
                # du7 = -(Ip + (du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]))
                du7 = -(Ip + (du1^2+du2^2+du3^2)/2 - Z*(u[1]^2+u[2]^2+u[3]^2+soft_core)^(-0.5))
                @SVector [du1,du2,du3,du4,du5,du6,du7]
            end
        elseif phase_method == :SCTS
            function traj_dipole_scts(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                tFx, tFy, tFz = -Z*(u[1]^2+u[2]^2+u[3]^2+soft_core)^(-1.5) .* (u[1],u[2],u[3])
                du1 = u[4]
                du2 = u[5]
                du3 = u[6]
                du4 = tFx-laserFx(t)
                du5 = tFy-laserFy(t)
                du6 = tFz
                # du7 = -(Ip + (du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]) + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                du7 = -(Ip + (du1^2+du2^2+du3^2)/2 - Z*(u[1]^2+u[2]^2+u[3]^2+soft_core)^(-0.5) + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                @SVector [du1,du2,du3,du4,du5,du6,du7]
            end
        end
    else
        if phase_method == :CTMC
            function traj_nondipole_ctmc(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                tFx, tFy, tFz = -Z*(u[1]^2+u[2]^2+u[3]^2+soft_core)^(-1.5) .* (u[1],u[2],u[3])
                c0 = 137.035999173
                du1 = u[4] + u[3]*laserFx(t)/c0
                du2 = u[5] + u[3]*laserFy(t)/c0
                du3 = u[6]
                du4 = tFx - laserFx(t)
                du5 = tFy - laserFy(t)
                du6 = tFz - (u[4]*laserFx(t)+u[5]*laserFy(t))/c0
                @SVector [du1,du2,du3,du4,du5,du6]
            end
        elseif phase_method == :QTMC
            function traj_nondipole_qtmc(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                tFx, tFy, tFz = -Z*(u[1]^2+u[2]^2+u[3]^2+soft_core)^(-1.5) .* (u[1],u[2],u[3])
                c0 = 137.035999173
                du1 = u[4] + u[3]*laserFx(t)/c0
                du2 = u[5] + u[3]*laserFy(t)/c0
                du3 = u[6]
                du4 = tFx - laserFx(t)
                du5 = tFy - laserFy(t)
                du6 = tFz - (u[4]*laserFx(t)+u[5]*laserFy(t))/c0
                # du7 = -(Ip + (du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]))
                du7 = -(Ip + (du1^2+du2^2+du3^2)/2 - Z*(u[1]^2+u[2]^2+u[3]^2+soft_core)^(-0.5))
                @SVector [du1,du2,du3,du4,du5,du6,du7]
            end
        elseif phase_method == :SCTS
            function traj_nondipole_scts(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                tFx, tFy, tFz = -Z*(u[1]^2+u[2]^2+u[3]^2+soft_core)^(-1.5) .* (u[1],u[2],u[3])
                c0 = 137.035999173
                du1 = u[4] + u[3]*laserFx(t)/c0
                du2 = u[5] + u[3]*laserFy(t)/c0
                du3 = u[6]
                du4 = tFx - laserFx(t)
                du5 = tFy - laserFy(t)
                du6 = tFz - (u[4]*laserFx(t)+u[5]*laserFy(t))/c0
                # du7 = -(Ip + (du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]) + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                du7 = -(Ip + (du1^2+du2^2+du3^2)/2 - Z*(u[1]^2+u[2]^2+u[3]^2+soft_core)^(-0.5) + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                @SVector [du1,du2,du3,du4,du5,du6,du7]
            end
        end
    end
end
"""
Gets the exponential term of ADK rate which depends on
Field strength `F`,
momentum's transverse component `kd` (in xy plane),
and propagation-direction (which is Z axis) component `kz`.
"""
ADKRateExp(t::HydrogenLikeAtom) = (F,kd,kz) -> exp(-2(kd^2+kz^2+2*t.Ip)^1.5/3F)

"Prints the information of the atom."
Base.show(io::IO, t::HydrogenLikeAtom) = print(io, "[HydrogenLikeAtom] Atom $(t.name), Ip=$(t.Ip), Z=$(t.nucl_charge), soft_core=$(t.soft_core)\n")

using Parameters, OrderedCollections
"Returns a `Dict{Symbol,Any}` containing properties of the object."
function Serialize(t::HydrogenLikeAtom)
    dict = OrderedDict{Symbol,Any}()
    type        = typeof(t)
    Ip          = t.Ip
    nucl_charge = t.nucl_charge
    soft_core   = t.soft_core
    name        = t.name
    @pack! dict = (type, Ip, nucl_charge, soft_core, name)
    return dict
end