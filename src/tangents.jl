using HomotopyContinuation

# Points of {u = 0} in the closed rectangle R at which curvature vanishes. 
# Characterized by the equation u_y^2*u_{xx} - 2u_x*u_y*u_{xy} + u_x^2*u_{yy} = 0. 
function curvature_vanishing(curve, vars, xlims, ylims)
    u = curve[1]
    ux = differentiate(u, vars[1])
    uy = differentiate(u, vars[2])
    N = uy^2 * differentiate(ux, vars[1]) - 2 * ux * uy * differentiate(ux, vars[2]) +
        ux^2 * differentiate(uy, vars[2])
    system = System([u, N], variables=vars)
    result = solve(system; start_system=:polyhedral,
        tracker_options=TrackerOptions(; parameters=:conservative),
        endgame_options=EndgameOptions())
    in_box(p) = xlims[1] <= p[1] <= xlims[2] && ylims[1] <= p[2] <= ylims[2]
    sign_changes = [p for p in real_solutions(result; only_nonsingular=true) if in_box(p)]    
    return sign_changes
end
