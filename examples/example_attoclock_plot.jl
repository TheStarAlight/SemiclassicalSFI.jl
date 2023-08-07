# === example_attoclock_plot.jl ===
#
# To correctly generate the plots using the code in the example, the HDF5, Colors, ColorSchemes, Plots and PyPlot (which is successfully precompiled) packages are also required to be installed via the package manager of julia.


@info "Plotting..."
using HDF5
using Plots
pyplot()    # PyPlot backend, need to install the PyPlot package manually. Switching to other backends is okay but the output might be slightly different.
include("ijet.jl")  # modified jet colormap which has white color for the bottom value.
clim = (-5,0)

function plot_heatmap(h5file; kwargs...)
    ion_prob = read(h5file, "ion_prob")
    px = read(h5file, "px")
    py = read(h5file, "py")
    spec = read(h5file, "momentum_spec_2D")
    spec ./= ion_prob
    @. spec = max(1e-100, spec)     # remove zero value because the figure is displayed in logarithmic scale.
    return heatmap(px,py,log10.(spec'), xlim=(-2,2), ylim=(-2,2), c=:ijet, clim=clim,
                aspectratio=:equal, tick_direction=:out, bordercolor=:black, framestyle=:box; kwargs...)
end

figure_list = Vector{Any}(nothing,3)
figure_list[1] = plot_heatmap(h5open("SCSFI_He_4e14_800nm_2cyc_CP_ADK_CTMC.h5");    title="ADK",   xlabel=raw"$p_x$ (a.u.)", ylabel=raw"$p_y$ (a.u.)", legend=:none)
figure_list[2] = plot_heatmap(h5open("SCSFI_He_4e14_800nm_2cyc_CP_SFAAE_CTMC.h5");  title="SFAAE", xlabel=raw"$p_x$ (a.u.)", legend=:none)
figure_list[3] = plot_heatmap(h5open("SCSFI_He_4e14_800nm_2cyc_CP_SFA_CTMC.h5");    title="SFA",   xlabel=raw"$p_x$ (a.u.)")

l = @layout [a{0.32w} b{0.32w} c{0.36w}]
plot(figure_list..., size=(900,300), layout=l, fmt=:svg)
savefig("example_attoclock.pdf")
@info "Plot saved as \"example_attoclock.pdf\"."