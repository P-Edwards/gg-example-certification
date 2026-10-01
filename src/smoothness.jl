import Oscar

# Affine smoothness of {u = 0} in C^2, established by showing that the ideal
# I = <u, u_x, u_y> in Q[x,y] is the unit ideal. By the weak Nullstellensatz
# this is equivalent to the complex variety V(I) being empty.
#
# A positive result is accompanied by a Nullstellensatz certificate: cofactors
# a, b, c with a*u + b*u_x + c*u_y = 1. Verifying the certificate requires only
# polynomial arithmetic.
#
# Oscar is imported rather than brought into the namespace: it exports names
# such as coefficients and degree that also come from MultivariatePolynomials,
# which the rest of the pipeline uses.

const SMOOTHNESS_VARIABLE_NAMES = [:x, :y]

# Clears denominators, giving a generator with coprime integer coefficients.
# Scaling a generator by a nonzero rational leaves the ideal unchanged; the
# scale factor is returned so certificates can be related to the original
# generators.
function integralize(p)
    cs = collect(Oscar.coefficients(p))
    isempty(cs) && return p, Oscar.QQ(1)
    scale = Oscar.QQ(lcm([Oscar.denominator(c) for c in cs]),
                     gcd([Oscar.numerator(c) for c in cs]))
    return scale * p, scale
end

# Converts the curve, given in DynamicPolynomials form with exact rational
# coefficients, into an Oscar polynomial over QQ.
function oscar_curve(dympoly, dyvars)
    R, (x, y) = Oscar.polynomial_ring(Oscar.QQ, SMOOTHNESS_VARIABLE_NAMES)
    oscar_variables = Dict(dyvars[1] => x, dyvars[2] => y)
    result = zero(R)
    for (c, m) in zip(coefficients(dympoly), monomials(dympoly))
        monomial = one(R)
        for v in dyvars
            monomial *= oscar_variables[v]^degree(m, v)
        end
        result += Oscar.QQ(numerator(c), denominator(c)) * monomial
    end
    return R, [x, y], result
end

# Generators of the singular locus ideal, scaled to integer coefficients.
function singular_locus_generators(R, vars, u)
    scaled = integralize.([u, Oscar.derivative(u, vars[1]), Oscar.derivative(u, vars[2])])
    return [p for (p, _) in scaled], [s for (_, s) in scaled]
end

# Expresses 1 as a combination of the ideal generators, or returns nothing if
# the ideal is not the unit ideal. The transformation matrix T of a Groebner
# basis satisfies gens(I) * T == gens(G); a constant element of G therefore
# gives the cofactors as the corresponding column of T.
function unit_cofactors(I)
    basis, transformation =
        Oscar.groebner_basis_with_transformation_matrix(I; complete_reduction=true)
    basis_elements = Oscar.gens(basis)
    index = findfirst(p -> !iszero(p) && Oscar.total_degree(p) == 0, basis_elements)
    index === nothing && return nothing
    scale = inv(Oscar.leading_coefficient(basis_elements[index]))
    return [scale * transformation[i, index] for i in 1:Oscar.number_of_generators(I)]
end

# Returns (certified, cofactors, generators, scales).
function smoothness_certificate(dympoly, dyvars)
    R, vars, u = oscar_curve(dympoly, dyvars)
    generators, scales = singular_locus_generators(R, vars, u)
    cofactors = unit_cofactors(Oscar.ideal(R, generators))
    cofactors === nothing && return false, nothing, generators, scales
    return true, cofactors, generators, scales
end

# Checks a*u + b*u_x + c*u_y == 1 by direct polynomial arithmetic.
function verify_certificate(generators, cofactors)
    length(generators) == length(cofactors) || return false
    combination = sum(cofactors[i] * generators[i] for i in eachindex(generators))
    return combination == one(parent(generators[1]))
end

function polynomial_to_terms(p)
    return [(exponents[1], exponents[2],
             BigInt(Oscar.numerator(c)), BigInt(Oscar.denominator(c)))
            for (c, exponents) in Oscar.coefficients_and_exponents(p)]
end

function terms_to_polynomial(R, vars, rows)
    result = zero(R)
    for (ex, ey, num, den) in rows
        result += Oscar.QQ(num, den) * vars[1]^ex * vars[2]^ey
    end
    return result
end

# Certificate file format: a header line per polynomial giving its role and
# term count, followed by one line per term as
#   <x exponent> <y exponent> <coefficient numerator> <coefficient denominator>
function write_smoothness_certificate(path, k, generators, cofactors)
    open(path, "w") do io
        println(io, "# affine smoothness certificate for run k = ", k)
        println(io, "# sum_i cofactor_i * generator_i = 1 in Q[x,y]")
        println(io, "k ", k)
        for (label, polynomials) in (("generator", generators), ("cofactor", cofactors))
            for (index, p) in enumerate(polynomials)
                rows = polynomial_to_terms(p)
                println(io, label, " ", index, " ", length(rows))
                for (ex, ey, num, den) in rows
                    println(io, ex, " ", ey, " ", num, " ", den)
                end
            end
        end
    end
    return path
end

function read_smoothness_certificate(path)
    R, (x, y) = Oscar.polynomial_ring(Oscar.QQ, SMOOTHNESS_VARIABLE_NAMES)
    vars = [x, y]
    lines = filter(line -> !isempty(line) && !startswith(line, "#"), readlines(path))
    k, generators, cofactors = nothing, [], []
    index = 1
    while index <= length(lines)
        fields = split(lines[index])
        if fields[1] == "k"
            k = parse(Int, fields[2])
            index += 1
            continue
        end
        count = parse(Int, fields[3])
        rows = [(parse(Int, f[1]), parse(Int, f[2]), parse(BigInt, f[3]), parse(BigInt, f[4]))
                for f in (split(lines[index + offset]) for offset in 1:count)]
        push!(fields[1] == "generator" ? generators : cofactors,
              terms_to_polynomial(R, vars, rows))
        index += count + 1
    end
    return k, R, vars, generators, cofactors
end

# Recomputes the generators for run k, checks that they agree with those stored
# in the certificate, then checks the certificate identity.
function check_smoothness_certificate(path)
    k, R, vars, stored_generators, cofactors = read_smoothness_certificate(path)
    dympoly, dyvars = produce_curve_from_data(DENOM_BASES[k], TORSION_NUMBERS[k]; dympoly=true)
    _, curve_vars, u = oscar_curve(dympoly[1], dyvars)
    recomputed_generators, _ = singular_locus_generators(parent(u), curve_vars, u)
    generators_match = length(recomputed_generators) == length(stored_generators) &&
        all(recomputed_generators[i] == stored_generators[i] for i in eachindex(stored_generators))
    return k, generators_match, verify_certificate(stored_generators, cofactors)
end
