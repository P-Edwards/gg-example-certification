using LinearAlgebra, HomotopyContinuation
import Contour

# Points of {u = 0} inside the box, sampled by marching squares.
function curve_point_cloud(curve, vars, xlims, ylims; resolution=1500)
    xs, ys, Z = curve_grid_fast(curve, vars, xlims, ylims; resolution=resolution)
    points = Vector{Vector{Float64}}()
    for line in Contour.lines(Contour.contour(collect(xs), collect(ys), Matrix(transpose(Z)), 0.0))
        lx, ly = Contour.coordinates(line)
        append!(points, [[lx[i], ly[i]] for i in eachindex(lx)])
    end
    @assert !isempty(points) "contour extraction found no curve points"
    return points
end

# The nearest points of the cloud to z, one per branch: after each contact is
# taken, all points within cluster_radius of it are dropped.
function nearest_contacts(z, points, order; cluster_radius=0.5)
    remaining = collect(eachindex(points))
    contacts, distances = Vector{Vector{Float64}}(), Float64[]
    for _ in 1:order
        isempty(remaining) && return nothing
        candidate_distances = [norm(points[i] - z) for i in remaining]
        nearest = argmin(candidate_distances)
        contact = points[remaining[nearest]]
        push!(contacts, contact)
        push!(distances, candidate_distances[nearest])
        remaining = [i for i in remaining if norm(points[i] - contact) > cluster_radius]
    end
    return contacts, distances
end

# The point of C = R cap {u > 0} at which the nearest contacts are most nearly
# equidistant. The search is confined to C and the spread is measured relative
# to the contact distance: at a point far from the curve all contacts are
# nearly equidistant, so an absolute spread would be minimized by running away
# from the curve. Returns (z, contacts, distances), or nothing.
function approximate_junction(curve, vars, xlims, ylims;
                              order=3, resolution=1500, cluster_radius=0.5,
                              coarse=200, refine=60)
    points = curve_point_cloud(curve, vars, xlims, ylims; resolution=resolution)
    curve_system = System(curve, variables=vars)
    inside_component(z) =
        xlims[1] <= z[1] <= xlims[2] && ylims[1] <= z[2] <= ylims[2] &&
        real(curve_system(z)[1]) > 0

    function relative_spread(z)
        inside_component(z) || return Inf
        found = nearest_contacts(z, points, order; cluster_radius=cluster_radius)
        found === nothing && return Inf
        distances = found[2]
        return (maximum(distances) - minimum(distances)) / minimum(distances)
    end

    function grid_argmin(xl, yl, n)
        best, best_z = Inf, nothing
        for xv in range(xl[1], xl[2]; length=n), yv in range(yl[1], yl[2]; length=n)
            value = relative_spread([xv, yv])
            value < best && (best = value; best_z = [xv, yv])
        end
        return best_z
    end

    z = grid_argmin(xlims, ylims, coarse)
    z === nothing && return nothing
    for scale in (0.2, 0.02)
        dx, dy = scale * (xlims[2] - xlims[1]), scale * (ylims[2] - ylims[1])
        candidate = grid_argmin([z[1] - dx, z[1] + dx], [z[2] - dy, z[2] + dy], refine)
        candidate === nothing || (z = candidate)
    end
    contacts, distances = nearest_contacts(z, points, order; cluster_radius=cluster_radius)
    return z, contacts, distances
end
