# === example_phase_methods_plot.jl ===
#
# To correctly generate the plots using the code in the example, the HDF5, Colors, ColorSchemes, Plots and PyPlot (which is successfully precompiled) packages are also required to be installed via the package manager of julia.


@info "Plotting..."
using HDF5
using Plots
pyplot() # PyPlot backend, need to install the PyPlot package manually. Switching to other backends is okay but the output might be slightly different.
include("ijet.jl") # modified jet colormap which has white color for the bottom value.
clim = (2,6)

function plot_heatmap(h5file; kwargs...)
    ion_prob = read(h5file, "ion_prob")
    px = read(h5file, "px")
    py = read(h5file, "py")
    pz = read(h5file, "pz")
    spec = read(h5file, "momentum_spec_3D")
    spec ./= ion_prob # normalize the total ionization yield
    @. spec = max(1e-100, spec) # remove zero value because the figure is displayed in logarithmic scale.
    # plots the slice on the pz-px plane near py=0.
    return heatmap(pz,px,log10.(mapreduce(i->spec[:,i,:],+,length(py)÷2-5:length(py)÷2+5))', xlim=(-1,1), ylim=(-1,1), c=:ijet, clim=clim,
                aspectratio=:equal, tick_direction=:out, bordercolor=:black, framestyle=:box; kwargs...)
end

figure_list = Vector{Any}(nothing,2)
figure_list[1] = plot_heatmap(h5open("SCSFI_H_0.9e14_800nm_8cyc_LP_ADKPara_ExpRate_QTMC.h5");   title="QTMC", xlabel=raw"$p_z$ (a.u.)", ylabel=raw"$p_x$ (a.u.)", legend=:none)
figure_list[2] = plot_heatmap(h5open("SCSFI_H_0.9e14_800nm_8cyc_LP_ADKPara_ExpRate_SCTS.h5");   title="SCTS", xlabel=raw"$p_z$ (a.u.)")

l = @layout [a{0.45w} b{0.55w}]
plot(figure_list..., size=(600,300), layout=l, fmt=:svg)
savefig("example_phase_methods.pdf")
@info "Plot saved as \"example_phase_methods.pdf\"."