using StaticArrays

"""
```
struct SAEAtom <: SAEAtomBase
```

Represents an atom under single-active-electron (SAE) approximation.

An instance of `SAEAtom` can be initialized via the constructor method:
```julia
SAEAtom(Ip, Z, a1=0.0, b1=0.0, a2=0.0, b2=0.0, a3=0.0, b3=0.0, name="[NA]")
```

# Example:
```jldoctest
julia> t = Targets.SAEAtom(Ip=0.9035698802, Z=1.0, a1=1.230723, b1=0.6620055, a2=-1.325040, b2=1.236224, a3=-0.2307230, b3=0.4804286, name="He")
[SAEAtom] Atom He, Ip=0.9035698802, Z=1.0
```
"""
struct SAEAtom <: SAEAtomBase
# should implement TargetPotential, TargetForce, TrajectoryFunction, ADKRateExp.
    "Ionization potential of the atom."
    Ip;
    "Asymptotic charge of the inner nucleus."
    nucl_charge;
    "Atomic parameters used to fit the atomic potential. See [J. Phys. B 38 2593 (2005)]"
    a1;
    b1;
    a2;
    b2;
    a3;
    b3;
    "Name of the atom."
    name::String;
    "Initializes a new instance of SAEAtom."
    SAEAtom( Ip, Z, a1=0., b1=0., a2=0., b2=0., a3=0., b3=0., name="[NA]") = new(Ip,Z,a1,b1,a2,b2,a3,b3,name)
    SAEAtom(;Ip, Z, a1=0., b1=0., a2=0., b2=0., a3=0., b3=0., name="[NA]") = new(Ip,Z,a1,b1,a2,b2,a3,b3,name)
end

"Gets the ionization potential of the atom."
IonPotential(t::SAEAtom) = t.Ip
"Gets the asymptotic nuclear charge of the atom."
AsympNuclCharge(t::SAEAtom) = t.nucl_charge
"Gets the name of the atom."
TargetName(t::SAEAtom) = t.name
"""
Gets the potential function of the atom.
Expression: V(r) = - [Z + a1*exp(-b1*r) + a2*r*exp(-b2*r) + a3*exp(-b3*r)] / r.
"""
function TargetPotential(t::SAEAtom)
    return function(x,y,z)
        r = sqrt(x^2+y^2+z^2)
        return - (t.nucl_charge + t.a1*exp(-t.b1*r) + t.a2*r*exp(-t.b2*r) + t.a3*exp(-t.b3*r)) / r
    end
end
"""
Gets the force exerted on the electron from the atom (which is the neg-grad of potential).
Expression: F(rvec) = - rvec / r^3 * [ Z + a1*(1+b1*r)*exp(-b1*r) + a3*(1+b3*r)*exp(-b3*r) ] - rvec * a2*b2/r * exp(-b2*r).
"""
function TargetForce(t::SAEAtom)
    return function(x,y,z)
        r = sqrt(x^2+y^2+z^2)
        return (-x,-y,-z) .* (r^(-3) * (t.nucl_charge + t.a1*(1+t.b1*r)*exp(-t.b1*r) + t.a3*(1+t.b3*r)*exp(-t.b3*r)) + t.a2*t.b2/r * exp(-t.b2*r))
    end
end
"Gets the trajectory function according to given parameter."
function TrajectoryFunction(t::SAEAtom, laserFx::Function, laserFy::Function, phase_method::Symbol, non_dipole::Bool; kwargs...)
    Z  = t.nucl_charge
    Ip = t.Ip
    a1,b1,a2,b2,a3,b3 = t.a1,t.b1,t.a2,t.b2,t.a3,t.b3
    return if ! non_dipole
        if phase_method == :CTMC
            function traj_dipole_ctmc(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                r = sqrt(u[1]^2+u[2]^2+u[3]^2)
                tFx, tFy, tFz = (u[1],u[2],u[3]) .* -(r^(-3)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/r*exp(-b2*r))
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
                r = sqrt(u[1]^2+u[2]^2+u[3]^2)
                tFx, tFy, tFz = (u[1],u[2],u[3]) .* -(r^(-3)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/r*exp(-b2*r))
                du1 = u[4]
                du2 = u[5]
                du3 = u[6]
                du4 = tFx-laserFx(t)
                du5 = tFy-laserFy(t)
                du6 = tFz
                # du7 = -(Ip + (du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]))
                du7 = -(Ip + (du1^2+du2^2+du3^2)/2 - (Z+a1*exp(-b1*r)+a2*r*exp(-b2*r)+a3*exp(-b3*r))/r)
                @SVector [du1,du2,du3,du4,du5,du6,du7]
            end
        elseif phase_method == :SCTS
            function traj_dipole_scts(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                r = sqrt(u[1]^2+u[2]^2+u[3]^2)
                tFx, tFy, tFz = (u[1],u[2],u[3]) .* -(r^(-3)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/r*exp(-b2*r))
                du1 = u[4]
                du2 = u[5]
                du3 = u[6]
                du4 = tFx-laserFx(t)
                du5 = tFy-laserFy(t)
                du6 = tFz
                # du7 = -(Ip + (du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]) + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                du7 = -(Ip + (du1^2+du2^2+du3^2)/2 - (Z+a1*exp(-b1*r)+a2*r*exp(-b2*r)+a3*exp(-b3*r))/r + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                @SVector [du1,du2,du3,du4,du5,du6,du7]
            end
        end
    else
        if phase_method == :CTMC
            function traj_nondipole_ctmc(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                r = sqrt(u[1]^2+u[2]^2+u[3]^2)
                tFx, tFy, tFz = (u[1],u[2],u[3]) .* -(r^(-3)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/r*exp(-b2*r))
                c0 = 137.035999173
                du1 = u[4] + u[3]*laserFx(t)/c0
                du2 = u[5] + u[3]*laserFy(t)/c0
                du3 = u[6]
                du4 = tFx-laserFx(t)
                du5 = tFy-laserFy(t)
                du6 = tFz - (u[4]*laserFx(t)+u[5]*laserFy(t))/c0
                @SVector [du1,du2,du3,du4,du5,du6]
            end
        elseif phase_method == :QTMC
            function traj_nondipole_qtmc(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                r = sqrt(u[1]^2+u[2]^2+u[3]^2)
                tFx, tFy, tFz = (u[1],u[2],u[3]) .* -(r^(-3)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/r*exp(-b2*r))
                c0 = 137.035999173
                du1 = u[4] + u[3]*laserFx(t)/c0
                du2 = u[5] + u[3]*laserFy(t)/c0
                du3 = u[6]
                du4 = tFx-laserFx(t)
                du5 = tFy-laserFy(t)
                du6 = tFz - (u[4]*laserFx(t)+u[5]*laserFy(t))/c0
                # du7 = -(Ip + (du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]))
                du7 = -(Ip + (du1^2+du2^2+du3^2)/2 - (Z+a1*exp(-b1*r)+a2*r*exp(-b2*r)+a3*exp(-b3*r))/r)
                @SVector [du1,du2,du3,du4,du5,du6,du7]
            end
        elseif phase_method == :SCTS
            function traj_nondipole_scts(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                r = sqrt(u[1]^2+u[2]^2+u[3]^2)
                tFx, tFy, tFz = (u[1],u[2],u[3]) .* -(r^(-3)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/r*exp(-b2*r))
                c0 = 137.035999173
                du1 = u[4] + u[3]*laserFx(t)/c0
                du2 = u[5] + u[3]*laserFy(t)/c0
                du3 = u[6]
                du4 = tFx - laserFx(t)
                du5 = tFy - laserFy(t)
                du6 = tFz - (u[4]*laserFx(t)+u[5]*laserFy(t))/c0
                # du7 = -(Ip + (du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]) + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                du7 = -(Ip + (du1^2+du2^2+du3^2)/2 - (Z+a1*exp(-b1*r)+a2*r*exp(-b2*r)+a3*exp(-b3*r))/r + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
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
ADKRateExp(t::SAEAtom) = (F,kd,kz) -> exp(-2(kd^2+kz^2+2*t.Ip)^1.5/3F)

"Prints the information of the atom."
Base.show(io::IO, t::SAEAtom) = print(io, "[SAEAtom] Atom $(t.name), Ip=$(t.Ip), Z=$(t.nucl_charge)\n")

using Parameters, OrderedCollections
"Returns a `Dict{Symbol,Any}` containing properties of the object."
function Serialize(t::SAEAtom)
    dict = OrderedDict{Symbol,Any}()
    type        = typeof(t)
    Ip          = t.Ip
    nucl_charge = t.nucl_charge
    name        = t.name
    a1 = t.a1
    b1 = t.b1
    a2 = t.a2
    b2 = t.b2
    a3 = t.a3
    b3 = t.b3
    @pack! dict = (type, Ip, nucl_charge, a1,b1,a2,b2,a3,b3, name)
    return dict
end