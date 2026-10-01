include(joinpath(@__DIR__, "..", "src", "WFSCertification.jl"))
using Dates

# Usage:
#   julia --project scripts/run_full_pipeline.jl --k 5 --order 3
#   julia --project scripts/run_full_pipeline.jl --points -3, -1, 1, 3 --order 2 --denominator 300
# Options:
#   --k N                 number of torsion points / run index (required unless
#                         --points is given)
#   --order N             bottleneck order. Default 2.
#   --points P1, P2, ...  the full list of real anchors a_j of f = -prod(x + iy - a_j),
#                         as exact rationals (e.g. 1/2, 1), replacing the tabulated
#                         anchors; excludes --k and requires --denominator
#   --denominator B       override the denominator base tabulated in
#                         src/curves.jl, giving eps = 1/B^2, eps2 = 1/B^3
#   --skip-mathematica    skip the CAD step (metadata is still written)
#   --skip-smoothness     skip the affine smoothness certificate (needs Oscar)
#   --curvanishing        count the points of {u = 0} in R where curvature is 0
#   --wolframscript PATH  Wolfram executable to use; accepts wolframscript OR a
#                         kernel binary (WolframKernel/math), auto-detected by name.
#                         If omitted, standard install locations are searched.

function parse_point(token)
    parts = split(token, "/")
    if length(parts) == 1
        return Rational{Int128}(parse(Int128, parts[1]))
    elseif length(parts) == 2
        return parse(Int128, parts[1]) // parse(Int128, parts[2])
    end
    error("Could not parse point: " * token)
end


function parse_arguments(args)
    options = Dict{String,Any}(
        "skip_mathematica" => false,
        "skip_smoothness" => false,
        "curvanishing" => false)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--k"
            options["k"] = parse(Int, args[i + 1])
            i += 2
        elseif a == "--order"
            options["order"] = parse(Int, args[i + 1])
            i += 2
        elseif a == "--denominator"
            options["denominator"] = parse(Int128, args[i + 1])
            i += 2
        elseif a == "--points"
            j = i + 1
            while j <= length(args) && !startswith(args[j], "--")
                j += 1
            end
            tokens = split(join(args[(i + 1):(j - 1)], " "), r"[,\s]+"; keepempty=false)
            if isempty(tokens)
                error("--points needs at least one point.")
            end
            options["points"] = parse_point.(tokens)
            i = j
        elseif a == "--skip-mathematica"
            options["skip_mathematica"] = true
            i += 1
        elseif a == "--skip-smoothness"
            options["skip_smoothness"] = true
            i += 1
        elseif a == "--curvanishing"
            options["curvanishing"] = true
            i += 1
        elseif a == "--wolframscript"
            options["wolframscript"] = args[i + 1]
            i += 2
        else
            error("Unknown option: " * a)
        end
    end

    if !haskey(options, "order")
        options["order"] = 2
    end
    if haskey(options, "points")
        if haskey(options, "k")
            error("--k and --points are mutually exclusive.")
        end
        if !haskey(options, "denominator")
            error("--points requires --denominator.")
        end
    elseif !haskey(options, "k")
        error("--k is required unless --points is given.")
    end
    return options
end

function mathematica_bool(value)
    if value
        return "True"
    end
    return "False"
end

options = parse_arguments(ARGS)


order = options["order"]
manual_points = haskey(options, "points")
if manual_points
    anchors = options["points"]
    k = div(length(anchors), 2)
    base = options["denominator"]
    anchor_source = "--points"
    curve_tag = "manpts"
else
    k = options["k"]
    anchors = Rational{Int128}.(default_anchors(TORSION_NUMBERS[k]))
    if haskey(options, "denominator")
        base = options["denominator"]
    else
        base = Int128(DENOM_BASES[k])
    end
    anchor_source = "tabulated for k = $k"
    curve_tag = "k$(k)"
end
run_mathematica = !options["skip_mathematica"]
run_smoothness = !options["skip_smoothness"]
run_curvature = options["curvanishing"]
preferred_wolfram_binaries = String[]
if haskey(options, "wolframscript")
    push!(preferred_wolfram_binaries, options["wolframscript"])
end
curve_label = "k = $k"
run_tag = "$(curve_tag)_order$(order)_denom$(base)"
results_dir = joinpath(@__DIR__, "..", "results")
exports_dir = joinpath(@__DIR__, "..", "mathematica", "exports")
certificates_dir = joinpath(@__DIR__, "..", "certificates")

println("Anchors (", anchor_source, "): ", anchor_string(anchors))
if manual_points && sort(anchors) != sort(-anchors)
    println("  warning: the anchors are not symmetric about 0; the heuristic box R assumes they are")
end

if run_smoothness
    include(joinpath(@__DIR__, "..", "src", "smoothness.jl"))
end
mkpath(results_dir)
mkpath(exports_dir)

# 1 & 2: certified distance-function maxima and the heuristic box R
println("Running certified distance-function computation ($curve_label, order = $order, denominator = $base)...")
result = certified_distance_maxima(k, order; denominator=base, anchors=anchors)
curve, vars = result["curve"], result["vars"]
xlims, ylims = result["xlims"], result["ylims"]

curvature_vanishing_points = Vector{Vector{Float64}}()
if run_curvature
    println("Computing curvature vanishing points of {u = 0} in R...")
    curvature_vanishing_points = curvature_vanishing(curve, vars, xlims, ylims)
end

# 3: local maxima of u inside R (Julia cross-check on the CAD count)
println("Computing local maxima of u (Julia cross-check)...")
polynomial_result = certified_polynomial_maxima(k; xlims=xlims, ylims=ylims,
    denominator=base, anchors=anchors)
f_maxima = polynomial_result["maxima_points"]
f_counts = polynomial_result["counts"]

# 4: plot curve + both sets of maxima
plt = plot_curve(curve, vars, xlims, ylims; resolution=1200,
    title="$curve_label, order = $order, eps = 1/$base^2")
plot_maxima!(plt, result["maxima_points"], f_maxima)
plot_path = joinpath(results_dir, "plot_$(run_tag).png")
savefig(plt, plot_path)

scaled = plot_curve_to_scale(curve, vars, xlims, ylims;
    title="$curve_label, order = $order, eps = 1/$base^2", margin=0.01)
plot_maxima!(scaled, result["maxima_points"], f_maxima; marker_scale=0.6)
scaled_path = joinpath(results_dir, "plot_to_scale_$(run_tag).png")
savefig(scaled, scaled_path)

# 5: metadata export for Mathematica (exact rational curve and box)
dpoly, dvars = curve_for(k, base; anchors=anchors, dympoly=true)
metadata_path = joinpath(exports_dir, "run_$(run_tag).wl")
write_mathematica_metadata(metadata_path, [
    ("K", string(k)),
    ("ManualPoints", mathematica_bool(manual_points)),
    ("Anchors", "{" * anchor_string(anchors) * "}"),
    ("DenominatorBase", string(base)),
    ("BottleneckOrder", string(order)),
    ("GeneratedAt", "\"" * string(now()) * "\""),
    ("DistanceMaximaLowerBound", string(result["maxima_lower_bound"])),
    ("CertifiedPolynomialMaximaInBox", string(f_counts[:maximum])),
    ("SmoothnessChecked", mathematica_bool(run_smoothness)),
    ("Box", mathematica_box_string(result["xlims_rational"], result["ylims_rational"])),
    ("Curve", mathematica_polynomial_string(dpoly[1], dvars))])

# 6: CAD checks. wolframscript takes `-file`, a raw kernel takes `-script`.
# Returns (found, path) for the first usable Wolfram executable.
function find_wolfram_binary(preferred)
    candidates = copy(preferred)
    append!(candidates, ["wolframscript", "wolfram", "math", "WolframKernel", "MathKernel"])
    if Sys.isapple()
        for app in ("Mathematica", "Wolfram", "Wolfram Desktop"), bin in ("wolframscript", "WolframKernel")
            push!(candidates, "/Applications/$app.app/Contents/MacOS/$bin")
        end
    elseif Sys.islinux()
        for root in ("/usr/local/Wolfram/Mathematica", "/usr/local/Wolfram/WolframEngine", "/opt/Wolfram/Mathematica")
            if !isdir(root)
                continue
            end
            for version in sort(readdir(root); rev=true)
                push!(candidates, joinpath(root, version, "Executables", "wolframscript"))
                push!(candidates, joinpath(root, version, "Executables", "math"))
            end
        end
    elseif Sys.iswindows()
        for root in ("C:\\Program Files\\Wolfram Research\\Mathematica",)
            if !isdir(root)
                continue
            end
            for version in sort(readdir(root); rev=true)
                push!(candidates, joinpath(root, version, "wolframscript.exe"))
                push!(candidates, joinpath(root, version, "math.exe"))
            end
        end
    end
    for candidate in candidates
        resolved = Sys.which(candidate)
        if resolved !== nothing
            return true, resolved
        end
        if isfile(candidate)
            return true, candidate
        end
    end
    return false, ""
end

# CAD_* values must be integers or True/False; anything else means the
# Mathematica side did not evaluate and the numbers must not be reported.
const CAD_INTEGER_KEYS = ["CAD_COMPONENTS", "CAD_MAXIMA_NONDEGENERATE",
                          "CAD_MAXIMA_DEGENERATE", "CAD_MAXIMA_TOTAL"]
const CAD_BOOLEAN_KEYS = ["CAD_EDGES_CLEAR", "CAD_CONTAINS_TEST_POINT",
                          "CAD_BOX_CERTIFIED", "CAD_MAXIMA_ALL_ISOLATED"]

function parse_cad_output(output)
    values = Dict{String,String}()
    for line in split(output, '\n')
        if !startswith(line, "CAD_")
            continue
        end
        parts = split(line, ' '; limit=2)
        if length(parts) == 2
            values[parts[1]] = strip(parts[2])
        end
    end
    return values
end

function cad_problems(values)
    problems = String[]
    if haskey(values, "CAD_ERROR")
        push!(problems, values["CAD_ERROR"])
    end
    for key in CAD_INTEGER_KEYS
        if !haskey(values, key)
            push!(problems, "missing $key")
        elseif tryparse(Int, values[key]) === nothing
            push!(problems, "$key is not an integer: $(values[key])")
        end
    end
    for key in CAD_BOOLEAN_KEYS
        if !haskey(values, key)
            push!(problems, "missing $key")
        elseif !(values[key] in ("True", "False"))
            push!(problems, "$key is not True/False: $(values[key])")
        end
    end
    return problems
end

# Returns (ran, values, status, details). Output is captured even when the
# kernel exits nonzero: the script reports failures on CAD_ERROR lines, which
# would otherwise be discarded.
function run_cad_checks(preferred_binaries, script, package, metadata)
    found, binary = find_wolfram_binary(preferred_binaries)
    if !found
        return false, Dict{String,String}(), "no Wolfram executable found",
            ["searched PATH and the standard install locations",
             "give one with --wolframscript PATH (wolframscript or a kernel binary)",
             "run manually: wolframscript -file $script $package $metadata"]
    end

    if occursin("wolframscript", lowercase(basename(binary)))
        flag = "-file"
    else
        flag = "-script"
    end
    # The paths are also supplied in the environment: a raw kernel invoked with
    # -script does not populate $ScriptCommandLine.
    command = addenv(`$binary $flag $script $package $metadata`,
                     "WFS_CAD_PACKAGE" => abspath(package),
                     "WFS_CAD_METADATA" => abspath(metadata))
    manual = "run manually: WFS_CAD_PACKAGE=$(abspath(package)) " *
             "WFS_CAD_METADATA=$(abspath(metadata)) $binary $flag $script $package $metadata"

    println("Running CAD checks via ", binary, "...")
    output = ""
    try
        output = read(ignorestatus(command), String)
    catch e
        return false, Dict{String,String}(), "Wolfram could not be run", [string(e), manual]
    end

    parsed = parse_cad_output(output)
    problems = cad_problems(parsed)
    if isempty(problems)
        return true, parsed, "", String[]
    end

    details = copy(problems)
    for line in split(strip(output), '\n')
        if !isempty(strip(line))
            push!(details, "output: " * line)
        end
    end
    push!(details, manual)
    return false, Dict{String,String}(), "Wolfram ran but the checks did not complete", details
end

cad_script = joinpath(@__DIR__, "..", "mathematica", "cad_checks.wls")
cad_package = joinpath(@__DIR__, "..", "mathematica", "CADChecks.wl")
if run_mathematica
    cad_ran, cad, cad_status, cad_details =
        run_cad_checks(preferred_wolfram_binaries, cad_script, cad_package, metadata_path)
else
    cad_ran = false
    cad = Dict{String,String}()
    cad_status = "skipped (--skip-mathematica)"
    cad_details = ["run manually: wolframscript -file $cad_script $cad_package $metadata_path"]
end

# 7: affine smoothness of {u = 0}
# Returns (certified, status, certificate path).
function certify_affine_smoothness(k, curve_tag, base, dympoly, dyvars, output_directory)
    println("Certifying affine smoothness of {u = 0}...")
    certified, cofactors, generators, _ = smoothness_certificate(dympoly, dyvars)
    if !certified
        return false, "NOT certified: <u, u_x, u_y> is not the unit ideal", ""
    end
    if !verify_certificate(generators, cofactors)
        return false, "certificate identity failed to verify", ""
    end
    mkpath(output_directory)
    path = write_smoothness_certificate(
        joinpath(output_directory, "smoothness_$(curve_tag)_denom$(base).txt"), k, generators, cofactors)
    return true, "certified", path
end

if run_smoothness
    smoothness_certified, smoothness_status, smoothness_path =
        certify_affine_smoothness(k, curve_tag, base, dpoly[1], dvars, certificates_dir)
else
    smoothness_certified = false
    smoothness_status = "skipped (--skip-smoothness)"
    smoothness_path = ""
end

# 8: console summary
println("\n", "="^70)
println("Run: $curve_label, bottleneck order = $order")
println("Anchors: ", anchor_string(anchors))
println("Bounding box R: x in $(xlims), y in $(ylims)")
println("-"^70)
println("Section 1: Certified distance-function pipeline (counts are geometric bottlenecks)")
println("  harvested path solutions:            ", result["n_harvested"])
println("  certified:                           ", result["n_certified"])
println("  after degenerate-multiplier filter:  ", result["n_not_on_quadric"])
println("  certified real:                      ", result["n_real"])
println("  inside R:                            ", result["n_in_box"])
println("  circumcenter in C (u > 0):           ", result["n_on_component"])
println("  confirmed local maxima of d_X (lower bound): ", result["maxima_lower_bound"])
println("-"^70)
println("Section 2: Local maxima of u inside R (Julia cross-check): ", f_counts[:maximum])
println("-"^70)
println("Section 3: Affine smoothness of {u = 0}")
println("  ", smoothness_status)
if smoothness_certified
    println("  certificate: ", smoothness_path)
end
println("-"^70)
println("Section 4: CAD certification (Mathematica)")
if cad_ran
    blocked_edges = get(cad, "CAD_BLOCKED_EDGES", "none")
    if blocked_edges == "none"
        blocked_note = ""
    else
        blocked_note = "  (blocked: " * blocked_edges * ")"
    end
    if string(f_counts[:maximum]) == cad["CAD_MAXIMA_TOTAL"]
        agreement = "yes"
    else
        agreement = "NO"
    end
    println("  components of {u >= 0} in R:  ", cad["CAD_COMPONENTS"])
    println("  edges of R disjoint from {u = 0}:", cad["CAD_EDGES_CLEAR"], blocked_note)
    println("  box certified:                ", cad["CAD_BOX_CERTIFIED"])
    println("  local maxima of u in component:     ", cad["CAD_MAXIMA_TOTAL"],
        "  (nondegenerate ", cad["CAD_MAXIMA_NONDEGENERATE"],
        " + degenerate ", cad["CAD_MAXIMA_DEGENERATE"], ")")
    if cad["CAD_MAXIMA_ALL_ISOLATED"] != "True"
        println("  WARNING: a solution component is positive-dimensional;")
        println("  the figures above count components, not points.")
    end
    println("  agrees with Julia count:      ", agreement)
else
    println("  NOT RUN: ", cad_status)
    for detail in cad_details
        println("    ", detail)
    end
end
if run_curvature
    println("-"^70)
    println("Section 5: Curvature vanishing points of {u = 0} in R (numerical): ", length(curvature_vanishing_points))
    for p in curvature_vanishing_points
        println("  ", p)
    end
end
println("="^70)
println("Plot:     ", plot_path)
println("To scale: ", scaled_path)
println("Metadata: ", metadata_path)