using HomotopyContinuation, Arblib, DynamicPolynomials, LinearAlgebra

function convert_correspondence_solution_point_to_standard_form(solution_certificate, variable_map)
    k = length(variable_map["barycentric_coordinates"])
    interval = certified_solution_interval(solution_certificate)
    candidate = solution_candidate(solution_certificate)
    endpoint_intervals = [interval[variable_map["endpoint_coordinates" * string(i)]] for i in 1:k]
    endpoint_candidates = [candidate[variable_map["endpoint_coordinates" * string(i)]] for i in 1:k]
    mult_intervals = interval[variable_map["barycentric_coordinates"]]
    mult_candidates = candidate[variable_map["barycentric_coordinates"]]
    return Dict{String,Any}(
        "endpoints" => endpoint_intervals,
        "endpoint_candidates" => endpoint_candidates,
        "multipliers" => mult_intervals,
        "multiplier_candidates" => mult_candidates,
        "circumcenter_interval" => sum(mult_intervals[i] * endpoint_intervals[i] for i in 1:k),
        "circumcenter_candidate" => sum(mult_candidates[i] * endpoint_candidates[i] for i in 1:k),
        "lambda_coordinates" => candidate[variable_map["lambda_coordinates"]],
        "distance_squared" => interval[variable_map["distance_function"]][1],
        "function_value" => real(interval[variable_map["function_value"]][1]),
        "function_value_candidate" => real(candidate[variable_map["function_value"]][1]),
        "distance_squared_candidate" => candidate[variable_map["distance_function"]][1],
        "full_solution" => candidate,
        "solution_is_real" => solution_certificate.real)
end

check_two_bottleneck_multipliers(b) =
    all(m -> Arblib.overlaps(Acb(1 / 2 + 0im), m), b["multipliers"])

check_bottleneck_multipliers(b) =
    all(m -> !Arblib.overlaps(Acb(1 + 0im), m), b["multipliers"])

check_bottleneck_multipliers_strict(b) =
    all(m -> !Arblib.overlaps(Acb(1 + 0im), m) && !Arblib.overlaps(Acb(0im), m), b["multipliers"])

function minimum_distance_start_system(F, x)
    N, m = length(x), length(F)
    @unique_var lambdas[0:m]
    @var epsilon[1:m]
    @var p[1:N]
    @var D
    perturbed_equations = [F[i] - epsilon[i] for i in 1:m]
    gradient_of_original = differentiate(F, x)
    normal_equations = vcat([[lambdas[1] * (x[j] - p[j]) for j in 1:N] +
                             sum([[lambdas[i + 1] * gradient_of_original[i, j] for j in 1:N] for i in 1:m])]...)
    patch_equation = [dot(rand(m + 1), lambdas) - 1]
    distance_equation = [D - dot(p - x, p - x)]
    return System([perturbed_equations; patch_equation; normal_equations; distance_equation],
        variables=[x; lambdas; D], parameters=[p; epsilon])
end

# Enclosure of the minimum of a collection of real balls. The minimum lies
# between the smallest of the lower bounds and the smallest of the upper
# bounds. 
function minimum_enclosure(values)
    lower = minimum(Arblib.lbound(Arblib.Arf, v) for v in values)
    upper = minimum(Arblib.ubound(Arblib.Arf, v) for v in values)
    return Arblib.union(Arb(lower), Arb(upper))
end

# Upper bound for the distance between the certified circumcenter and the
# floating point candidate that the homotopy was tracked to.
function displacement_bound(point_interval, point_candidate)
    magnitudes = [Arblib.abs_ubound(Arb, coordinate - candidate)
                  for (coordinate, candidate) in zip(point_interval, point_candidate)]
    return sqrt(sum(m * m for m in magnitudes))
end


function compute_minimum_distance_interval(distance_system, initial_solution, start, b)
    point_candidate = real(b["circumcenter_candidate"])
    point_interval = real(b["circumcenter_interval"])
    target = [point_candidate; 0.0]
    distance_results = solve(distance_system, initial_solution;
        start_parameters=start, target_parameters=target,
        tracker_options=TrackerOptions(; parameters=:conservative))
    subbed = System(subs(distance_system.expressions, parameters(distance_system) => target),
        variables=variables(distance_system))
    certs = filter(c -> c.real, certificates(certify(subbed, distance_results)))
    certified_values = [real(certified_solution_interval(c))[end] for c in certs]
    minimum_squared = minimum_enclosure(certified_values)
    Arblib.ispositive(minimum_squared) ||
        error("enclosure of the minimum squared distance does not exclude zero")


    shift = displacement_bound(point_interval, point_candidate)
    # distance is an interval enclosing the 
    # minimum distance between the numerical candidate
    # circumcenter and the curve
    distance = sqrt(minimum_squared)
    # The distance-to-curve function is 1-Lipschitz, 
    # so the following encloses the minimum distance
    # between the true circumcenter and the curve
    transported = Arblib.union(distance - shift, distance + shift)
    # We square this for convenience in comparing to 
    # squared distances later
    return transported * transported
end


check_circumcenter_on_component(b) = Arblib.ispositive(b["function_value"])

flatten_bottleneck_intervals(b) = [b["endpoints"]; [b["circumcenter_interval"]]]

function check_the_bounds_cert(b, xlims, ylims)
    for point in flatten_bottleneck_intervals(b)
        (real(point[1]) > xlims[2] || real(point[1]) < xlims[1]) && return false
        (real(point[2]) > ylims[2] || real(point[2]) < ylims[1]) && return false
    end
    return true
end

function find_multipliers_for_critical_on_algebraic_medial_axis(
        alg_correspondence, alg_variable_map, bottleneck_variable_map, b)
    shared = intersect(keys(alg_variable_map), keys(bottleneck_variable_map))
    alg_vars = variables(alg_correspondence)
    subbed_variables = alg_vars
    for label in shared
        subbed_variables = subs(subbed_variables,
            alg_vars[alg_variable_map[label]] => b["full_solution"][bottleneck_variable_map[label]])
    end
    subbed_variables = subs(subbed_variables,
        alg_vars[alg_variable_map["circumcenter_coordinates"]] => b["circumcenter_candidate"])
    substituted = subs(expressions(alg_correspondence), alg_vars => subbed_variables)
    lagrange_system = System(substituted, variables=alg_vars[alg_variable_map["lagrange_multipliers"]])
    results = solve(lagrange_system; tracker_options=ESPECIALLY_CONSERVATIVE, endgame_options=PERMISSIVE_ENDGAME)
    return subs(subbed_variables, variables(lagrange_system) => solutions(results)[1]), results
end

# Bordered Hessian of the Lagrangian for the interval second-derivative test.
# Built with DynamicPolynomials so that it can be evaluated over Acb intervals.
function dympoly_second_derivative(F, x, k; multi_system=false)
    @assert k > 1
    N = length(x)
    systems = multi_system ? F : [F for _ in 1:k]

    @polyvar y[1:k, 1:N]
    correspondence_original_functions = vcat([subs(systems[i], x => y[i, :]) for i in 1:k]...)
    @polyvar z[1:N]

    max_codimension = maximum(length(system) for system in systems)
    @polyvar lambdas[0:max_codimension, 1:k]
    final_lambda_variables = vcat([lambdas[2:(length(systems[i]) + 1), i] for i in 1:k]...)

    flattened_normals = []
    for (i, system) in enumerate(systems)
        n = length(system)
        gradients = subs(differentiate(system, x), x => y[i, :])
        normals = (y[i, :] - z) - sum(lambdas[j + 1, i] * gradients[j, :] for j in 1:n)
        push!(flattened_normals, vcat(normals...))
    end
    flattened_normals = vcat(flattened_normals...)

    distance_equations = [dot(y[i, :] - z, y[i, :] - z) - dot(y[i + 1, :] - z, y[i + 1, :] - z) for i in 1:(k - 1)]
    constraints = [correspondence_original_functions; distance_equations; flattened_normals]
    constraint_vars = [y'[:]; z; final_lambda_variables]

    objective = sum(dot(y[i, :] - z, y[i, :] - z) for i in 1:k)
    @polyvar active_lambdas[1:length(constraints)]
    base_lagrangian = objective - dot(active_lambdas, constraints)
    bordered_hessian = differentiate(differentiate(base_lagrangian,
        [active_lambdas; constraint_vars]), [active_lambdas; constraint_vars])

    return bordered_hessian,
           [constraint_vars; active_lambdas],
           length(constraint_vars),
           length(constraints)
end

# Constant value of a fully-substituted DynamicPolynomials expression.
function convert_dympoly_expression_to_acb(e)
    cs = coefficients(e)
    return isempty(cs) ? Acb(0.0) : Acb(cs[1])
end

function substituted_second_deriv_test(symbolic_matrix, solution_interval, variable_order;
                                       number_of_variables, number_of_constraints)
    subbed = convert_dympoly_expression_to_acb.(subs(symbolic_matrix, variable_order => solution_interval))
    return [det(subbed[1:j, 1:j]) for j in (2 * number_of_constraints + 1):(number_of_variables + number_of_constraints)]
end

# Sufficient second-order condition for a constrained local maximum: with n
# variables and m constraints, the leading principal minors of the bordered
# Hessian of orders 2m+1, ..., n+m satisfy sign(minor of order m+j) = (-1)^j.
# Returns true only when every required sign is certified over the whole
# interval. An empty minor sequence (n == m) returns false.
function certified_constrained_maximum(minors, number_of_constraints)
    isempty(minors) && return false
    for (index, minor) in enumerate(minors)
        j = number_of_constraints + index
        certified = iseven(j) ? Arblib.ispositive(real(minor)) : Arblib.isnegative(real(minor))
        certified || return false
    end
    return true
end
