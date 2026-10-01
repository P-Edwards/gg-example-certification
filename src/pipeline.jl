using HomotopyContinuation, Arblib, LinearAlgebra

# Solutions of the bottleneck correspondence occur in orbits of size
# factorial(order) under the symmetric group permuting the endpoints, all
# representing the same geometric bottleneck. Counts are reported per orbit.
function geometric_count(n, order)
    orbit = factorial(order)
    # Occasionally some of the solutions in a particular
    # orbit are misclassified as singular solutions when 
    # doing a numerical solve. This ceiling gives a correct
    # count on the minimum number of orbits found. 
    return ceil(n / orbit)
end

function distinct_points(points; tolerance)
    representatives = Vector{eltype(points)}()
    for p in points
        any(q -> maximum(abs.(p - q)) <= tolerance, representatives) && continue
        push!(representatives, p)
    end
    return representatives
end


function certified_distance_maxima(k, bottleneck_order; denominator=nothing, anchors=nothing)
    denom_base = denominator_base(k, denominator)
    curve, vars = curve_for(k, denom_base; anchors=anchors)

    if bottleneck_order == 2
        results, correspondence, variable_map, alg_correspondence, alg_var_map =
            solve_unparameterized_correspondence(curve, vars, bottleneck_order)
    else
        results, _, _, _, _ = solve_correspondence_with_monodromy(curve, vars, bottleneck_order)
        correspondence, variable_map, alg_correspondence, alg_var_map, _ =
            bottleneck_correspondence_equations_with_first_deriv(curve, vars, bottleneck_order; multi_system=false)
    end

    harvested = solutions(results)
    certified = certify(correspondence, harvested)
    standard_form = (c -> convert_correspondence_solution_point_to_standard_form(c, variable_map)).(filter(is_certified, certificates(certified)))

    multiplier_filter = bottleneck_order == 2 ? check_bottleneck_multipliers : check_bottleneck_multipliers_strict
    not_on_quadric = filter(multiplier_filter, standard_form)

    diee_expected, diee_satisfied = nothing, nothing
    if bottleneck_order == 2
        d = anchors === nothing ? 2 * k : max(3, length(anchors))
        diee_expected = d^4 - 5 * d^2 + 4 * d
        diee_satisfied = length(not_on_quadric) == diee_expected
    end

    real_certs = filter(b -> b["solution_is_real"], not_on_quadric)
    xrational, yrational = certification_box(curve, vars)
    xlims, ylims = float_box(xrational), float_box(yrational)
    in_box = filter(b -> check_the_bounds_cert(b, xlims, ylims), real_certs)
    on_component = filter(check_circumcenter_on_component, in_box)

    distance_system = minimum_distance_start_system(curve, vars)
    start = randn(ComplexF64, length(vars) + length(curve))
    initial_solution = solutions(solve(distance_system, target_parameters=start,
        start_system=:polyhedral, tracker_options=ESPECIALLY_CONSERVATIVE,
        endgame_options=PERMISSIVE_ENDGAME))
    for b in on_component
        b["minimum_distance_interval"] = compute_minimum_distance_interval(distance_system, initial_solution, start, b)
    end
    candidates = filter(b -> Arblib.overlaps(real(b["distance_squared"]), b["minimum_distance_interval"]), on_component)

    is_max = Bool[]
    if !isempty(candidates)
        with_mults, _ = unzip((b -> find_multipliers_for_critical_on_algebraic_medial_axis(
            alg_correspondence, alg_var_map, variable_map, b)).(candidates))
        certified_mults = certify(alg_correspondence, with_mults)
        solution_intervals = (c -> certified_solution_interval(c)[:]).(certificates(certified_mults))
        # The dympoly versions here are to consult a technical issue. HomotopyContinuation.jl supplies
        # an expression manipulation library called ModelKit which, amongst other things, does not
        # handle some arblib types at all and does not perform interval arithmetic after the substitution 
        # when it does.

        # Another reasonable question may be why these variables and necessary function values 
        # aren't just included directly in the original correspondence so that they can be certified that way. 
        # The issue is that the bordered Hessian determinants are already too large to do a cofactor expansion
        # with polynomial expressions efficiently, and are also fairly high degree. This order is the only one that
        # works, then.
        dpoly, dvars = curve_for(k, denom_base; anchors=anchors, dympoly=true)
        sd_matrix, sd_vars, num_vars, num_constraints = dympoly_second_derivative(dpoly, dvars, bottleneck_order)
        computed_tests = (si -> substituted_second_deriv_test(sd_matrix, si, sd_vars;
            number_of_variables=num_vars, number_of_constraints=num_constraints)).(solution_intervals)
        is_max = [is_certified(c) && certified_constrained_maximum(t, num_constraints) for (c,t) in zip(certificates(certified_mults),computed_tests)]
    end

    box_scale = max(xlims[2] - xlims[1], ylims[2] - ylims[1])
    maxima_points = distinct_points(
        [real.(candidates[i]["circumcenter_candidate"]) for i in eachindex(candidates) if is_max[i]];
        tolerance=1e-8 * box_scale)

    return Dict(
        "curve" => curve, "vars" => vars, "k" => k, "denom_base" => denom_base,
        "bottleneck_order" => bottleneck_order, "xlims" => xlims, "ylims" => ylims,
        "xlims_rational" => xrational, "ylims_rational" => yrational,
        "n_harvested" => geometric_count(length(harvested), bottleneck_order),
        "n_certified" => geometric_count(length(standard_form), bottleneck_order),
        "n_not_on_quadric" => geometric_count(length(not_on_quadric), bottleneck_order),
        "n_real" => geometric_count(length(real_certs), bottleneck_order),
        "n_in_box" => geometric_count(length(in_box), bottleneck_order),
        "n_on_component" => geometric_count(length(on_component), bottleneck_order),
        "n_minimizer_candidates" => geometric_count(length(candidates), bottleneck_order),
        "maxima_lower_bound" => geometric_count(count(is_max), bottleneck_order),
        "maxima_points" => maxima_points,
        "diee_expected" => diee_expected,
        "diee_satisfied" => diee_satisfied)
end

function box_status(interval, xlims, ylims)
    x, y = real(interval[1]), real(interval[2])
    (x > xlims[1] && x < xlims[2] && y > ylims[1] && y < ylims[2]) && return :inside
    (x < xlims[1] || x > xlims[2] || y < ylims[1] || y > ylims[2]) && return :outside
    return :boundary_ambiguous
end

# This function provides a partially certified alternative to the CAD-based
# local maxima check for u on C = R\cap {u > 0}. Its count is not a fully certified
# uppper bound because the numerical solve ("critical_results = ...") 
# may in principle improperly filter some solutions as singular which are not actually in practice. 
# This does not seem to affect the outputs on our particular examples in practice.
function certified_polynomial_maxima(k; xlims=nothing, ylims=nothing, denominator=nothing, anchors=nothing)
    denom_base = denominator_base(k, denominator)
    curve, vars = curve_for(k, denom_base; anchors=anchors)
    dpoly, dvars = curve_for(k, denom_base; anchors=anchors, dympoly=true)
    if xlims === nothing || ylims === nothing
        xrational, yrational = certification_box(curve, vars)
        xlims, ylims = float_box(xrational), float_box(yrational)
    end

    gradient_system = System(differentiate(curve[1], vars), variables=vars)
    critical_results = solve(gradient_system; start_system=:polyhedral,
        tracker_options=ESPECIALLY_CONSERVATIVE, endgame_options=PERMISSIVE_ENDGAME)
    certification = certify(gradient_system, solutions(critical_results))
    real_certs = filter(c -> is_certified(c) && c.real, certificates(certification))

    dhessian = differentiate(differentiate(dpoly[1], dvars), dvars)
    counts = Dict(:maximum => 0, :minimum => 0, :saddle => 0, :inconclusive => 0,
                  :boundary_ambiguous => 0, :outside => 0,
                  :off_component => 0, :component_ambiguous => 0)
    maxima = []
    for cert in real_certs
        interval = certified_solution_interval(cert)[:]
        status = box_status(interval, xlims, ylims)
        if status != :inside
            counts[status] += 1
            continue
        end
        value = real(convert_dympoly_expression_to_acb(subs(dpoly[1], dvars => interval)))
        if Arblib.isnegative(value)
            counts[:off_component] += 1
            continue
        elseif !Arblib.ispositive(value)
            counts[:component_ambiguous] += 1
            continue
        end
        H = convert_dympoly_expression_to_acb.(subs(dhessian, dvars => interval))
        h11 = real(H[1, 1])
        detH = real(H[1, 1] * H[2, 2] - H[1, 2] * H[2, 1])
        classification = if Arblib.ispositive(detH) && Arblib.isnegative(h11)
            :maximum
        elseif Arblib.ispositive(detH) && Arblib.ispositive(h11)
            :minimum
        elseif Arblib.isnegative(detH)
            :saddle
        else
            :inconclusive
        end
        counts[classification] += 1
        classification == :maximum && push!(maxima, Dict(
            "point" => real.(solution_candidate(cert)),
            "interval" => interval,
            "value_interval" => value))
    end

    return Dict(
        "counts" => counts,
        "maxima" => maxima,
        "maxima_points" => [m["point"] for m in maxima],
        "n_distinct_certified" => ndistinct_certified(certification),
        "n_real_certified" => length(real_certs),
        "xlims" => xlims, "ylims" => ylims)
end
