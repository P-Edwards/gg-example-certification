using Plots, HomotopyContinuation, DynamicPolynomials, CSV, DataFrames

# High-precision grid: pass the DynamicPolynomials form of the curve (exact
# rational coefficients, via produce_curve_from_data(...; dympoly=true)) and
# evaluate in BigFloat. Use this for large denominator bases where Float64
# cancellation can distort the zero level set.
function curve_grid_exact(curve_dympoly, dyvars, xlims, ylims; resolution=400)
    xs = range(xlims[1], xlims[2]; length=resolution)
    ys = range(ylims[1], ylims[2]; length=resolution)
    Z = Matrix{Float64}(undef, length(ys), length(xs))
    for (j, xv) in enumerate(xs), (i, yv) in enumerate(ys)
        Z[i, j] = Float64(subs(curve_dympoly[1], dyvars => [BigFloat(xv), BigFloat(yv)]))
    end
    return xs, ys, Z
end

# Evaluates the curve on a grid through a compiled HomotopyContinuation system.
function curve_grid_fast(curve_modelkit, vars, xlims, ylims; resolution=800)
    sys = System(curve_modelkit, variables=vars)
    xs = range(xlims[1], xlims[2]; length=resolution)
    ys = range(ylims[1], ylims[2]; length=resolution)
    Z = Matrix{Float64}(undef, length(ys), length(xs))
    for (j, xv) in enumerate(xs), (i, yv) in enumerate(ys)
        Z[i, j] = real(sys([xv, yv])[1])
    end
    return xs, ys, Z
end

# Draws a view larger than the given box by `margin`, a fraction of each side's
# width. Affects the drawn view only, not the box.
function plot_curve(curve_modelkit, vars, xlims, ylims;
                    resolution=800, size=(900, 500), margin=0.15, kwargs...)
    xpad, ypad = margin * (xlims[2] - xlims[1]), margin * (ylims[2] - ylims[1])
    view_x = [xlims[1] - xpad, xlims[2] + xpad]
    view_y = [ylims[1] - ypad, ylims[2] + ypad]
    xs, ys, Z = curve_grid_fast(curve_modelkit, vars, view_x, view_y; resolution=resolution)
    contour(xs, ys, Z; levels=[0.0], color=:black, colorbar=false,
        linewidth=1.5, size=size, xlabel="x", ylabel="y",
        xlims=Tuple(view_x), ylims=Tuple(view_y), kwargs...)
end

# Overlays both kinds of maxima: distance-function maxima as hollow rings and
# maxima of f as smaller filled markers, so coincident points remain readable.
function plot_maxima!(plt, distance_maxima, polynomial_maxima;
                      distance_color=:red, polynomial_color=:green, marker_scale=1.0)
    if !isempty(distance_maxima)
        scatter!(plt, [p[1] for p in distance_maxima], [p[2] for p in distance_maxima];
            markershape=:circle, markersize=10 * marker_scale, markercolor=RGBA(0, 0, 0, 0),
            markerstrokecolor=distance_color, markerstrokewidth=2,
            label="distance-function maxima")
    end
    if !isempty(polynomial_maxima)
        scatter!(plt, [p[1] for p in polynomial_maxima], [p[2] for p in polynomial_maxima];
            markershape=:utriangle, markersize=5 * marker_scale, markercolor=polynomial_color,
            markerstrokecolor=polynomial_color, label="maxima of u")
    end
    return plt
end

_center(b) = haskey(b, "circumcenter_candidate") ? real.(b["circumcenter_candidate"]) : real.(b["circumcenter"])
_endpoints(b) = haskey(b, "endpoint_candidates") ? [real.(e) for e in b["endpoint_candidates"]] : [real.(e) for e in b["endpoints"]]

# Overlays bottlenecks (certified standard form or floating-point dicts) on an
# existing curve plot: circumcenters, endpoints, and connecting segments.
function plot_bottlenecks!(plt, bottlenecks; center_color=:red, endpoint_color=:blue,
                           segment_color=:gray, markersize=4)
    for b in bottlenecks
        c = _center(b)
        for e in _endpoints(b)
            plot!(plt, [c[1], e[1]], [c[2], e[2]]; color=segment_color, label="", linewidth=1)
            scatter!(plt, [e[1]], [e[2]]; color=endpoint_color, label="", markersize=markersize)
        end
        scatter!(plt, [c[1]], [c[2]]; color=center_color, label="", markersize=markersize, markershape=:diamond)
    end
    return plt
end

# Draws the curve at true geometric scale: one unit in x is one unit in y, so
# the shape of the component is undistorted. The canvas is sized from the data
# ratio so that a long thin component still fills the frame. Type is set in
# points and does not scale with the canvas, so on a figure this wide the
# defaults read as tiny; `fontsize` sets the x tick labels and is the base for
# the x axis label, title and legend. The y axis carries few ticks, since its
# labels stack along the short side of the frame.
function plot_curve_to_scale(curve_modelkit, vars, xlims, ylims;
                             resolution=1600, width=1800, min_height=220,
                             margin=0.05, fontsize=18, yticks=[-1, 0, 1],
                             furniture=nothing, kwargs...)
    xpad, ypad = margin * (xlims[2] - xlims[1]), margin * (ylims[2] - ylims[1])
    view_x = [xlims[1] - xpad, xlims[2] + xpad]
    view_y = [ylims[1] - ypad, ylims[2] + ypad]
    data_ratio = (view_y[2] - view_y[1]) / (view_x[2] - view_x[1])
    furniture === nothing && (furniture = round(Int, 12 * fontsize))
    height = max(min_height, round(Int, width * data_ratio) + furniture)
    xs, ys, Z = curve_grid_fast(curve_modelkit, vars, view_x, view_y; resolution=resolution)
    contour(xs, ys, Z; levels=[0.0], color=:black, colorbar=false,
        linewidth=1.5, size=(width, height), xlabel="x", ylabel="y",
        xlims=Tuple(view_x), ylims=Tuple(view_y), aspect_ratio=:equal,
        legend=:outertopright,
        titlefontsize=fontsize + 6, legendfontsize=fontsize,
        tickfontsize=fontsize, guidefontsize=fontsize + 2, yticks=yticks,
        left_margin=10Plots.mm, bottom_margin=10Plots.mm, top_margin=5Plots.mm,
        kwargs...)
end

function plot_curve_with_bottlenecks(curve_modelkit, vars, bottlenecks;
                                     xlims=nothing, ylims=nothing, resolution=800, kwargs...)
    if xlims === nothing || ylims === nothing
        xrational, yrational = certification_box(curve_modelkit, vars)
        xlims, ylims = float_box(xrational), float_box(yrational)
    end
    plt = plot_curve(curve_modelkit, vars, xlims, ylims; resolution=resolution, kwargs...)
    plot_bottlenecks!(plt, bottlenecks)
    return plt
end

# Bottlenecks from the CSV format written by the exploration scripts:
# columns bneck-x, bneck-y, endpt-j-x, endpt-j-y
function bottlenecks_from_csv(path)
    using_frame = CSV.read(path, DataFrame)
    k = (ncol(using_frame) - 2) ÷ 2
    bottlenecks = []
    for row in eachrow(using_frame)
        push!(bottlenecks, Dict(
            "circumcenter" => [row[1], row[2]],
            "endpoints" => [[row[2j + 1], row[2j + 2]] for j in 1:k]))
    end
    return bottlenecks
end
