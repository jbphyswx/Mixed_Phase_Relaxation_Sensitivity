
export_plot(z, all_data, ϕ_all, filename, ylabel; xlabel) = nothing
export_plot_snapshot(z, all_data, ϕ_all, filename, ylabel) = nothing

"""
    plot_friendly_name(ϕ)

Get plot-friendly string, since many Unicode
characters do not render in plot labels.
"""
function plot_friendly_name(ϕ)
    s = ϕ
    s = replace(s, "ρ" => "rho")
    s = replace(s, "α" => "alpha")
    s = replace(s, "∂" => "partial")
    s = replace(s, "∇" => "nabla")
    return s
end

"""
    export_plot(z, all_data, ϕ_all, filename, ylabel)

Export plot of all variables, or all
available time-steps in `all_data`.
"""
function export_plot(z, all_data, ϕ_all, filename, ylabel; xlabel = nothing)
    ϕ_all isa Tuple || (ϕ_all = (ϕ_all,))
    p = plot()
    for n in 0:(length(keys(all_data)) - 1)
        for ϕ in ϕ_all
            ϕ_string = String(ϕ)
            ϕ_name = plot_friendly_name(ϕ_string)
            ϕ_data = all_data[n][ϕ_string][:]
            if !isnothing(xlabel)
                plot!(ϕ_data, z, xlabel = xlabel, ylabel = ylabel)
            else
                plot!(ϕ_data, z, xlabel = ϕ_name, ylabel = ylabel)
            end
        end
    end
    savefig(filename)
end

"""
    export_plot(z, all_data, ϕ_all, filename, ylabel)

Export plot of all variables, or all
available time-steps in `all_data`.
"""
function export_plot(
    z,
    all_data::Array,
    ϕ_all,
    filename,
    ylabel;
    xlabel = nothing,
)
    ϕ_all isa Tuple || (ϕ_all = (ϕ_all,))
    p = plot()
    for data in all_data
        for ϕ in ϕ_all
            ϕ_string = String(ϕ)
            ϕ_name = plot_friendly_name(ϕ_string)
            ϕ_data = data[ϕ_string][:]
            if !isnothing(xlabel)
                plot!(ϕ_data, z, xlabel = xlabel, ylabel = ylabel)
            else
                plot!(ϕ_data, z, xlabel = ϕ_name, ylabel = ylabel)
            end
        end
    end
    savefig(filename)
end

"""
    export_plot_snapshot(z, all_data, ϕ_all, filename, ylabel)

Export plot of all variables in `all_data`
"""
function export_plot_snapshot(z, all_data, ϕ_all, filename, ylabel)
    ϕ_all isa Tuple || (ϕ_all = (ϕ_all,))
    p = plot()
    for ϕ in ϕ_all
        ϕ_string = String(ϕ)
        ϕ_name = plot_friendly_name(ϕ_string)
        ϕ_data = all_data[ϕ_string][:]
        plot!(ϕ_data, z, xlabel = ϕ_name, ylabel = ylabel)
    end
    savefig(filename)
end

function save_binned_surface_plots(
    x,
    y,
    z,
    title,
    filename,
    n_plots = (3, 3),
    z_label_prefix = "z",
    n_digits = 5,
)
    n_z_partitions = prod(n_plots)
    z_min_global = min(z...)
    z_max_global = max(z...)
    Δz = (z_max_global - z_min_global) / n_z_partitions
    z_min = ntuple(i -> z_min_global + (i - 1) * Δz, n_z_partitions)
    z_max = ntuple(i -> z_min_global + (i) * Δz, n_z_partitions)
    p = []
    for i in 1:n_z_partitions
        mask = z_min[i] .<= z .<= z_max[i]
        x_i = x[mask]
        y_i = y[mask]
        sz_min = string(z_min[i])[1:min(n_digits, length(string(z_min[i])))]
        sz_max = string(z_max[i])[1:min(n_digits, length(string(z_max[i])))]
        p_i = plot(
            x_i,
            y_i,
            title = "$(title), in ($sz_min, $sz_max)",
            seriestype = :scatter,
            markersize = 5,
        )
        push!(p, p_i)
    end
    plot(p..., layout = n_plots, legend = false)
    savefig(filename)
end;
