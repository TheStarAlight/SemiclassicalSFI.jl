# code to generate plots used in manual2_lasers.md

using SemiclassicalSFI
using SemiclassicalSFI.Lasers
using Plots
pyplot()

# =======================

l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0)
F0 = LaserF0(l)
Fx = LaserFx(l)
Fy = LaserFy(l)
t = -400:1:400
Fxt = Fx.(t)/F0
Fyt = Fy.(t)/F0
color = RGB(0/255,154/255,250/255)
color_trans = RGBA(0/255,154/255,250/255,0.3)
f1=plot(t,Fxt,Fyt, label="ε=0.0", color=color)
plot!(repeat([-400],length(t)),Fxt,Fyt, label=:none, color=color_trans)
plot!(t,Fxt,repeat([-1.0],length(t)), label=:none, color=color_trans)
plot!(t,repeat([1.0],length(t)),Fyt, label=:none, color=color_trans)
plot!(xlim=(-400,400),ylim=(-1,1),zlim=(-1,1))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$F_x/F_0$", zlabel=raw"$F_y/F_0$")

l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0.5)
F0 = LaserF0(l)
Fx = LaserFx(l)
Fy = LaserFy(l)
t = -400:1:400
Fxt = Fx.(t)/F0
Fyt = Fy.(t)/F0
color = RGB(255/255,127/255,24/255)
color_trans = RGBA(255/255,127/255,24/255,0.3)
f2=plot(t,Fxt,Fyt, label="ε=0.5", color=color)
plot!(repeat([-400],length(t)),Fxt,Fyt, label=:none, color=color_trans)
plot!(t,Fxt,repeat([-1.0],length(t)), label=:none, color=color_trans)
plot!(t,repeat([1.0],length(t)),Fyt, label=:none, color=color_trans)
plot!(xlim=(-400,400),ylim=(-1,1),zlim=(-1,1))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$F_x/F_0$", zlabel=raw"$F_y/F_0$")

l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=1.0)
F0 = LaserF0(l)
Fx = LaserFx(l)
Fy = LaserFy(l)
t = -400:1:400
Fxt = Fx.(t)/F0
Fyt = Fy.(t)/F0
color = RGB(44/255,160/255,44/255)
color_trans = RGBA(44/255,160/255,44/255,0.3)
f3=plot(t,Fxt,Fyt, label="ε=1.0", color=color)
plot!(repeat([-400],length(t)),Fxt,Fyt, label=:none, color=color_trans)
plot!(t,Fxt,repeat([-1.0],length(t)), label=:none, color=color_trans)
plot!(t,repeat([1.0],length(t)),Fyt, label=:none, color=color_trans)
plot!(xlim=(-400,400),ylim=(-1,1),zlim=(-1,1))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$F_x/F_0$", zlabel=raw"$F_y/F_0$")

l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=-1.0)
F0 = LaserF0(l)
Fx = LaserFx(l)
Fy = LaserFy(l)
t = -400:1:400
Fxt = Fx.(t)/F0
Fyt = Fy.(t)/F0
color = RGB(214/255,39/255,40/255)
color_trans = RGBA(214/255,39/255,40/255,0.3)
f4=plot(t,Fxt,Fyt, label="ε=-1.0", color=color)
plot!(repeat([-400],length(t)),Fxt,Fyt, label=:none, color=color_trans)
plot!(t,Fxt,repeat([-1.0],length(t)), label=:none, color=color_trans)
plot!(t,repeat([1.0],length(t)),Fyt, label=:none, color=color_trans)
plot!(xlim=(-400,400),ylim=(-1,1),zlim=(-1,1))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$F_x/F_0$", zlabel=raw"$F_y/F_0$")

f=plot(f1,f2,f3,f4, layout=(2,2), size=(1000,800), fmt=:svg)
savefig("./docs/src/assets/manual2_lasers_ellip.svg")

# =======================

l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0, azi=0π)
F0 = LaserF0(l)
Fx = LaserFx(l)
Fy = LaserFy(l)
t = -400:1:400
Fxt = Fx.(t)/F0
Fyt = Fy.(t)/F0
color = RGB(0/255,154/255,250/255)
color_trans = RGBA(0/255,154/255,250/255,0.3)
f1=plot(t,Fxt,Fyt, label="azi=0", color=color)
plot!(repeat([-400],length(t)),Fxt,Fyt, label=:none, color=color_trans)
plot!(t,Fxt,repeat([-1.0],length(t)), label=:none, color=color_trans)
plot!(t,repeat([1.0],length(t)),Fyt, label=:none, color=color_trans)
plot!(xlim=(-400,400),ylim=(-1,1),zlim=(-1,1))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$F_x/F_0$", zlabel=raw"$F_y/F_0$")

l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0, azi=π/6)
F0 = LaserF0(l)
Fx = LaserFx(l)
Fy = LaserFy(l)
t = -400:1:400
Fxt = Fx.(t)/F0
Fyt = Fy.(t)/F0
color = RGB(255/255,127/255,24/255)
color_trans = RGBA(255/255,127/255,24/255,0.3)
f2=plot(t,Fxt,Fyt, label="azi=π/6 (30°)", color=color)
plot!(repeat([-400],length(t)),Fxt,Fyt, label=:none, color=color_trans)
plot!(t,Fxt,repeat([-1.0],length(t)), label=:none, color=color_trans)
plot!(t,repeat([1.0],length(t)),Fyt, label=:none, color=color_trans)
plot!(xlim=(-400,400),ylim=(-1,1),zlim=(-1,1))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$F_x/F_0$", zlabel=raw"$F_y/F_0$")

l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0, azi=π/3)
F0 = LaserF0(l)
Fx = LaserFx(l)
Fy = LaserFy(l)
t = -400:1:400
Fxt = Fx.(t)/F0
Fyt = Fy.(t)/F0
color = RGB(44/255,160/255,44/255)
color_trans = RGBA(44/255,160/255,44/255,0.3)
f3=plot(t,Fxt,Fyt, label="azi=π/3 (60°)", color=color)
plot!(repeat([-400],length(t)),Fxt,Fyt, label=:none, color=color_trans)
plot!(t,Fxt,repeat([-1.0],length(t)), label=:none, color=color_trans)
plot!(t,repeat([1.0],length(t)),Fyt, label=:none, color=color_trans)
plot!(xlim=(-400,400),ylim=(-1,1),zlim=(-1,1))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$F_x/F_0$", zlabel=raw"$F_y/F_0$")

l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0, azi=π/2)
F0 = LaserF0(l)
Fx = LaserFx(l)
Fy = LaserFy(l)
t = -400:1:400
Fxt = Fx.(t)/F0
Fyt = Fy.(t)/F0
color = RGB(214/255,39/255,40/255)
color_trans = RGBA(214/255,39/255,40/255,0.3)
f4=plot(t,Fxt,Fyt, label="azi=π/2 (90°)", color=color)
plot!(repeat([-400],length(t)),Fxt,Fyt, label=:none, color=color_trans)
plot!(t,Fxt,repeat([-1.0],length(t)), label=:none, color=color_trans)
plot!(t,repeat([1.0],length(t)),Fyt, label=:none, color=color_trans)
plot!(xlim=(-400,400),ylim=(-1,1),zlim=(-1,1))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$F_x/F_0$", zlabel=raw"$F_y/F_0$")

f=plot(f1,f2,f3,f4, layout=(2,2), size=(1000,800), fmt=:svg)
savefig("./docs/src/assets/manual2_lasers_azimuth.svg")

# =======================

t = -400:1:400
l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0, cep=0)
plot(t, UnitEnvelope(l).(t), color=:gray, fill=0, α=0.1, label=:none)
plot!(t, -1 .* UnitEnvelope(l).(t), color=:gray, fill=0, α=0.1, label=:none)

l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0, cep=0)
plot!(t, LaserAx(l).(t) ./ LaserA0(l), label="CEP=0", color=RGB(64/255,183/255,173/255))
l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0, cep=π/3)
plot!(t, LaserAx(l).(t) ./ LaserA0(l), label="CEP=π/3", color=RGB(52/255,143/255,167/255))
l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0, cep=2π/3)
plot!(t, LaserAx(l).(t) ./ LaserA0(l), label="CEP=2π/3", color=RGB(55/255,101/255,158/255))
l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0, cep=π)
plot!(t, LaserAx(l).(t) ./ LaserA0(l), label="CEP=π", color=RGB(65/255,61/255,163/255))

f=plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$f_\mathrm{env}$", xlim=(-300,300), ylim=(-1.05,1.05), size=(500,400), framestyle=:box, fmt=:svg)
savefig("./docs/src/assets/manual2_lasers_cep.svg")

# =======================

t = -600:1:600
l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=10, ellip=0)
plot(t, UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)
plot!(t, -1 .* UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)
plot!(t, LaserAx(l).(t) ./ LaserA0(l), color=RGB(0/255,154/255,250/255))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$f_\mathrm{env}$", xlim=(-600,600), ylim=(-1.05,1.05), legend=:none, size=(500,400), framestyle=:box, fmt=:svg)
savefig("./docs/src/assets/manual2_lasers_env_cos4.svg")

# =======================

t = -600:1:600
l = Cos2Laser(peak_int=1e14, wave_len=800.0, cyc_num=10, ellip=0)
plot(t, UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)
plot!(t, -1 .* UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)
plot!(t, LaserAx(l).(t) ./ LaserA0(l), color=RGB(0/255,154/255,250/255))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$f_\mathrm{env}$", xlim=(-600,600), ylim=(-1.05,1.05), legend=:none, size=(500,400), framestyle=:box, fmt=:svg)
savefig("./docs/src/assets/manual2_lasers_env_cos2.svg")

# =======================

t = -1000:1:1000
l = GaussianLaser(peak_int=1e14, wave_len=800.0, FWHM_duration=1103.2, ellip=0)
plot(t, UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)
plot!(t, -1 .* UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)
plot!(t, LaserAx(l).(t) ./ LaserA0(l), color=RGB(0/255,154/255,250/255))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$f_\mathrm{env}$", xlim=(-1000,1000), ylim=(-1.05,1.05), legend=:none, size=(500,400), framestyle=:box, fmt=:svg)
savefig("./docs/src/assets/manual2_lasers_env_gaussian.svg")

# =======================

t = -600:1:600
l = TrapezoidalLaser(peak_int=1e14, wave_len=800.0, cyc_num_turn_on=3, cyc_num_turn_off=3, cyc_num_const=4, ellip=0, t_shift=-551.6)
plot(t, UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)
plot!(t, -1 .* UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)
plot!(t, LaserAx(l).(t) ./ LaserA0(l), color=RGB(0/255,154/255,250/255))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$f_\mathrm{env}$", xlim=(-600,600), ylim=(-1.05,1.05), legend=:none, size=(500,400), framestyle=:box, fmt=:svg)
savefig("./docs/src/assets/manual2_lasers_env_trapezoidal.svg")
