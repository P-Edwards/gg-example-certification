Certified Magnanini's conjecture computations for Gladiali-Grossi curves
===================================================================

Copyright (C) 2026 [Parker Edwards](https://parkeredw.com)

Overview
--------

Certified computation checking Magnanini's conjecture for the Gladiali-Grossi family of plane curves
defined by

    1 - y^2 + (1/b^2) * y * (y^2 - 3x^2) + (1/b^3) * Re(f_k) = 0

where f_k is a product of 2k complex-linear factors with equally spaced real
anchors and b is an exact integer denominator base. See the manuscript *Counterexamples to Magnanini's Conjecture Concerning Critical Points of the Torsion Function* [here](https://doi.org/10.48550/arXiv.2608.17081) for details. 


Requirements
------------

- Julia 1.9 or later.
- Mathematica 11.3 or later (optional, for cylindrical algebraic decomposition).

Installation
------------

From the repository root, start Julia in interactive mode and run:

    julia

    julia> ]                                    
    (@v1.x) pkg> activate .
    (@v1.x) pkg> instantiate .

This will download and install the Julia dependencies indicated in the `Project.toml` file.


Usage
-----
If you would like to use CPU parallelization, run `export JULIA_NUM_THREADS=<number>` for whatever `<number>` of cores you would like before running other commands.

To run one of the examples reported in the above manuscript for `k = 2, 3, ..., 6`, 
use the following command template:

    julia --project=. scripts/run_full_pipeline.jl --k 5

which will use pre-chosen settings. The 3-bottleneck related example from the manuscript can be run using

    julia --project=. scripts/run_full_pipeline.jl --points -9,-7,-5,-3,-1,1,3,5,7,9  --denominator 1300

The following additional options are available:

- `--k N` - number of torsion points (required; 2-8).
- `--points P1, P2, ...` - This argument can be used to replace the `2k` "anchor" points
in the G-G construction with different list of `2k` points. They must be rational, specified
as a comma separated list like `-1, -1/2, 1/2, 1`. 
- `--denominator B` - controls epsilon used in the construction as `eps = 1/B^2`.
- `--skip-mathematica` - skip the CAD step.
- `--wolframscript PATH` - Wolfram executable to use. If omitted, script searches common locations.
- `--skip-smoothness` - skips certifying that the curve is algebraically smooth.
- `--curvanishing` - Runs a utility which numerically computes the number of places on the curve that the curvature vanishes.

Outputs
-------
By default, 2 plots (scaled and unscaled) similar to those in the manuscript are placed in `results/*`. If the smoothness
test was run, a machine-readable smoothness certificate is placed in `certificates/*`. Additionally, the console
prints results from each certification step. The primary outputs of interest are: 

- Section 1 (S1): confirmed local maxima of d_X (lower bound). 
- Section 2: Local maxima of $u$ inside $R$. This count is done numerically and is not certified to be an upper bound, 
but is useful if Mathematica is unavailable. 
- Section 4 (S4): local maxima of u in component. This is a symbolically certified count. 

If S4 <= S1, then Magnanini's conjecture is satisfied for the G-G example corresponding to the given inputs.