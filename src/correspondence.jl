using HomotopyContinuation, LinearAlgebra

const ESPECIALLY_CONSERVATIVE = TrackerOptions(; automatic_differentiation=3, parameters=:conservative, max_steps = 100_000, extended_precision = true)
const PERMISSIVE_ENDGAME = EndgameOptions()

function lagrange_multiplier_system(objective_function, constraints; x=variables(constraints))
    constraints_differential = differentiate(constraints, x)
    objective_differential = transpose(differentiate(objective_function, x))
    @unique_var lambdas[1:length(constraints)]
    base_system = objective_differential - transpose(lambdas) * constraints_differential
    return (base_system[:], lambdas)
end


function bottleneck_correspondence_equations_with_first_deriv(F, x, k;
        parameters=nothing, extra_equations=[], multi_system=true, target_parameters=nothing)
    @assert k > 1
    N = length(x)
    systems = multi_system ? F : [F for _ in 1:k]
    multi_system && @assert length(F) == k

    @unique_var y[1:k, 1:N]
    correspondence_original_functions = vcat([subs(systems[i], x => y[i, :]) for i in 1:k]...)

    @var t[1:k]
    barycenter_equation = [sum(t) - 1]
    @unique_var z[1:N]

    max_codimension = maximum(length(system) for system in systems)
    @unique_var lambdas[0:max_codimension, 1:k]
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
    @var D
    squared_distance_to_first_endpoint = [D - dot(y[1, :] - z, y[1, :] - z)]
  
    @var U
    function_value_at_circumcenter = [U - subs(systems[1][1], x => z)]

    lagrangian_eqns, lagrange_multiplier_variables = lagrange_multiplier_system(
        sum(dot(y[i, :] - z, y[i, :] - z) for i in 1:k),
        [correspondence_original_functions; distance_equations; flattened_normals];
        x=[y'[:]; z; final_lambda_variables])

    variable_map = build_variable_map(
        vcat(["endpoint_coordinates" * string(i) for i in 1:k],
             ["barycentric_coordinates", "lambda_coordinates", "distance_function", "function_value"]),
        vcat([y[i, :] for i in 1:k], [t, final_lambda_variables, D, U]))

    core = [correspondence_original_functions; barycenter_equation; flattened_normals;
            distance_equations; squared_distance_to_first_endpoint;
            function_value_at_circumcenter]
    substitute_z = eqn -> subs(eqn, z => sum(t[i] * y[i, :] for i in 1:k))
    bottleneck_equations = substitute_z.([core; extra_equations])
    bottleneck_equations_no_extras = substitute_z.(core)
    bottleneck_system = System(bottleneck_equations,
        variables=[y'[:]; t; final_lambda_variables; D; U], parameters=parameters)
    bottleneck_system_at_target = bottleneck_system
    if target_parameters !== nothing
        bottleneck_system_at_target = System(
            subs(bottleneck_equations_no_extras, parameters => target_parameters),
            variables=variables(bottleneck_system))
    end

    algebraic_variable_map = build_variable_map(
        vcat(["endpoint_coordinates" * string(i) for i in 1:k],
             ["circumcenter_coordinates", "lambda_coordinates", "lagrange_multipliers"]),
        vcat([y[i, :] for i in 1:k], [z, final_lambda_variables, lagrange_multiplier_variables]))
    algebraic_medial_axis_equations = [correspondence_original_functions; flattened_normals;
                                       distance_equations; lagrangian_eqns]
    if target_parameters !== nothing
        algebraic_medial_axis_equations = subs(algebraic_medial_axis_equations, parameters => target_parameters)
    end
    algebraic_medial_axis_system = System(algebraic_medial_axis_equations,
        variables=[y'[:]; z; final_lambda_variables; lagrange_multiplier_variables])

    return (bottleneck_system, variable_map, algebraic_medial_axis_system,
            algebraic_variable_map, bottleneck_system_at_target)
end

function build_variable_map(names, variable_lists)
    counter = 1
    mapping = Dict()
    for (label, vs) in zip(names, variable_lists)
        mapping[label] = counter:(counter + length(vs) - 1)
        counter += length(vs)
    end
    return mapping
end

function dense_parameterization(systems, x)
    split_polys, params, target = [], [], []
    for system in systems
        current = []
        for f in system
            current_dense_poly, current_dense_params = dense_poly(x, degree(f))
            push!(current, current_dense_poly)
            append!(params, current_dense_params)
            append!(target, coeffs_as_dense_poly(f, x, degree(f)))
        end
        push!(split_polys, current)
    end
    return split_polys, [params...], [target...]
end

function solve_unparameterized_correspondence(F, x, k)
    correspondence, var_map, alg_correspondence, alg_var_map, _ =
        bottleneck_correspondence_equations_with_first_deriv(F, x, k; multi_system=false)
    solved = solve(correspondence; tracker_options=ESPECIALLY_CONSERVATIVE, endgame_options=PERMISSIVE_ENDGAME)
    return solved, correspondence, var_map, alg_correspondence, alg_var_map
end

function solve_correspondence_with_monodromy(F, x, k; monodromy_options...)
    split_polys, params, target = dense_parameterization([F for _ in 1:k], x)
    correspondence, var_map, alg_correspondence, alg_var_map, unparameterized =
        bottleneck_correspondence_equations_with_first_deriv(split_polys, x, k;
            parameters=params, target_parameters=target)
    initial = monodromy_solve(correspondence; monodromy_options...)
    solved = solve(correspondence, solutions(initial);
        start_parameters=parameters(initial), target_parameters=target,
        tracker_options=ESPECIALLY_CONSERVATIVE)
    return solved, unparameterized, var_map, alg_correspondence, alg_var_map
end
