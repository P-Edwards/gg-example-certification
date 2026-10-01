using HomotopyContinuation, DynamicPolynomials, LinearAlgebra

# Run tables: run index == number of torsion points k; eps = 1/b^2, eps2 = 1/b^3
const TORSION_NUMBERS = [1, 2, 3, 4, 5, 6, 7, 8, 9]
const DENOM_BASES = [2^1, Int128(3800), Int128(3800), Int128(3800), Int128(50000),
                     Int128(20000), Int128(200000), Int128(2e6), Int128(2^26)]

# Exact rational division for integer inputs, float division otherwise
op(x, y) = (x isa Integer || x isa Int128) && (y isa Integer || y isa Int128) ? x // y : x / y

unzip(a) = map(x -> getfield.(a, x), fieldnames(eltype(a)))

function default_anchors(k; range=1)
    h = 1 // (2 * (2 * Int128(k) - 3))
    positive = [1 - (2j - 1) * h for j in (k - 1):-1:1]
    push!(positive, Rational{Int128}(1))
    return range * [-reverse(positive); positive]
end

anchor_string(anchors) = join([denominator(a) == 1 ? string(numerator(a)) :
                               string(numerator(a), "/", denominator(a)) for a in anchors], ", ")


function potential_counterexample_curve(k, eps; range=1, eps2=false, dympoly=false, anchors=nothing)
    @polyvar x y
    anchor_numbers = anchors === nothing ? default_anchors(k; range=range) : anchors
    factors_for_fk = [x + y * im - a for a in anchor_numbers]
    fk = -1 * prod(factors_for_fk)
    cs, ms = coefficients(fk), monomials(fk)
    real_part_fk = sum(real(cs[i]) * ms[i] for i in eachindex(cs) if abs(imag(cs[i])) < 1e-5)
    if eps2 == false
        eps2 = eps^(3 / 2)
    end
    curve = 1 - y^2 + eps * y * (y^2 - 3 * x^2) + eps2 * real_part_fk
    dympoly && return ([curve], [x; y])
    @var a b
    modelkit_eqn = Expression(subs(curve, [x; y] => [a; b]))
    @var x y
    return ([subs(modelkit_eqn, [a; b] => [x; y])], [x; y])
end

function produce_curve_from_data(denom_base_for_frac, number_of_torsion_points;
                                 range_of_torsion_points=1, dympoly=false, anchors=nothing)
    return potential_counterexample_curve(number_of_torsion_points,
        op(1, denom_base_for_frac^2);
        range=op(1, range_of_torsion_points),
        eps2=op(1, denom_base_for_frac^3),
        dympoly=dympoly,
        anchors=anchors)
end

function curve_for(k, base; anchors=nothing, dympoly=false)
    if anchors === nothing
        return produce_curve_from_data(base, TORSION_NUMBERS[k]; dympoly=dympoly)
    end
    return produce_curve_from_data(base, nothing; anchors=anchors, dympoly=dympoly)
end

function compact_component_bounds(curve, vars; expand=0.05)
    hsols = real_solutions(solve(System(subs(curve, vars[2] => 0), variables=[vars[1]])))
    vsols = real_solutions(solve(System(subs(curve, vars[1] => 0), variables=[vars[2]])))
    h = minimum(abs(s[1]) for s in hsols)
    v = minimum(abs(s[1]) for s in vsols)
    return [-(1 + expand) * h, (1 + expand) * h], [-(1 + expand) * v, (1 + expand) * v]
end

# The denominator base for a run, or an override supplied by the caller.
denominator_base(k, override=nothing) =
    override === nothing ? DENOM_BASES[k] : Int128(override)

# The rectangle R is rational with a power-of-two denominator, so that its
# endpoints convert to Float64 exactly and the numerical filtering and the CAD
# check use the same rectangle.
const BOX_DENOMINATOR = Int128(2)^12

function rationalize_box(xlims, ylims; denominator=BOX_DENOMINATOR)
    function round_outward(v, mode)
        n = mode === :down ? floor(Int128, v * denominator) : ceil(Int128, v * denominator)
        @assert abs(n) < Int128(2)^53 "box endpoint not exactly representable in Float64"
        return Rational{Int128}(n, denominator)
    end
    return [round_outward(xlims[1], :down), round_outward(xlims[2], :up)],
           [round_outward(ylims[1], :down), round_outward(ylims[2], :up)]
end

# The rectangle R: heuristic bounds rounded outward to rationals.
function certification_box(curve, vars; expand=0.25, denominator=BOX_DENOMINATOR)
    xlims, ylims = compact_component_bounds(curve, vars; expand=expand)
    return rationalize_box(xlims, ylims; denominator=denominator)
end

float_box(lims) = [Float64(lims[1]), Float64(lims[2])]
