(* ::Package:: *)

(* CAD-based checks for the WFS pipeline.
   Load with Get["CADChecks.wl"] from a notebook or wolframscript. *)

ClearAll[certifyBox, countMaximaInBox];

(* Certify that R ∩ {f > 0} has exactly one connected component, that it
   contains a given test point, and that f is strictly negative on the whole
   boundary of the box. Together these give C = R ∩ {f > 0} compact.

   The component count and the test point use f > 0. The edge tests use
   f >= 0: an edge is clear exactly when {f >= 0} does not meet it, which is
   the statement that f < 0 everywhere on that edge.

   Returns an Association; "Certified" is the conjunction of all three. *)
certifyBox[f_, {{xmin_, xmax_}, {ymin_, ymax_}}, {x_, y_}, testPoint_ : {0, 0}] :=
  Module[{box, containsTestPoint, edgeMeetsRegion, edgeStatus, blockedEdges, components},

    box = xmin <= x <= xmax && ymin <= y <= ymax;

    containsTestPoint =
      TrueQ[Simplify[(f /. {x -> testPoint[[1]], y -> testPoint[[2]]}) > 0]] &&
      TrueQ[xmin <= testPoint[[1]] <= xmax && ymin <= testPoint[[2]] <= ymax];

    (* Each edge test is a univariate real quantifier elimination:
       Reduce returns False exactly when the region misses the edge. *)
    edgeMeetsRegion[condition_, variable_] :=
      Reduce[condition, variable, Reals] =!= False;

    edgeStatus = {
      edgeMeetsRegion[(f /. x -> xmin) >= 0 && ymin <= y <= ymax, y],
      edgeMeetsRegion[(f /. x -> xmax) >= 0 && ymin <= y <= ymax, y],
      edgeMeetsRegion[(f /. y -> ymin) >= 0 && xmin <= x <= xmax, x],
      edgeMeetsRegion[(f /. y -> ymax) >= 0 && xmin <= x <= xmax, x]};

    blockedEdges = Pick[{"left", "right", "bottom", "top"}, edgeStatus];

    components = CylindricalDecomposition[f > 0 && box, {x, y}, "Components"];

    <|
      "Components" -> Length[components],
      "EdgesClear" -> blockedEdges === {},
      "BlockedEdges" -> blockedEdges,
      "ContainsTestPoint" -> containsTestPoint,
      "Certified" -> Length[components] == 1 && blockedEdges === {} && containsTestPoint
    |>];

(* Counts the local maxima of f inside the closed box, restricted to f > 0 so
   that only maxima belonging to the connected component of interest are
   counted. The count is split into two semialgebraic systems:

     nondegenerate:  grad f == 0, det Hessian > 0, f_xx < 0, f > 0, in box
     degenerate:     grad f == 0, det Hessian == 0,           f > 0, in box

   Degenerate critical points cannot be classified by the second-derivative
   test, so they are counted separately and added, making the total an upper
   bound on the number of local maxima.

   Each CAD component of a zero-dimensional solution set is a conjunction of
   equations. "AllIsolated" reports whether that holds for every component; if
   it is False some component is positive-dimensional and the counts are counts
   of components, not of points. *)
countMaximaInBox[f_, {{xmin_, xmax_}, {ymin_, ymax_}}, {x_, y_}] :=
  Module[{fx, fy, fxx, fxy, fyy, detHessian, box, criticalOnComponent,
          nondegenerate, degenerate, isolatedQ, components},

    {fx, fy} = {D[f, x], D[f, y]};
    {fxx, fxy, fyy} = {D[f, {x, 2}], D[D[f, x], y], D[f, {y, 2}]};
    detHessian = fxx fyy - fxy^2;

    box = xmin <= x <= xmax && ymin <= y <= ymax;
    criticalOnComponent = fx == 0 && fy == 0 && f > 0 && box;

    nondegenerate = CylindricalDecomposition[
      criticalOnComponent && detHessian > 0 && fxx < 0, {x, y}, "Components"];
    degenerate = CylindricalDecomposition[
      criticalOnComponent && detHessian == 0, {x, y}, "Components"];

    isolatedQ[component_] :=
      AllTrue[If[Head[component] === And, List @@ component, {component}],
              Head[#] === Equal &];
    components = Join[nondegenerate, degenerate];

    <|
      "Nondegenerate" -> Length[nondegenerate],
      "Degenerate" -> Length[degenerate],
      "Total" -> Length[nondegenerate] + Length[degenerate],
      "AllIsolated" -> AllTrue[components, isolatedQ]
    |>];
