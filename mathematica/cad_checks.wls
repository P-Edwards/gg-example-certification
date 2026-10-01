#!/usr/bin/env wolframscript
(* Runs the CAD checks for one pipeline run.
   Usage: wolframscript -file cad_checks.wls <CADChecks.wl> <metadata.wl>
   With a raw kernel, pass the paths in the environment instead:
     WFS_CAD_PACKAGE=<CADChecks.wl> WFS_CAD_METADATA=<metadata.wl> \
       wolfram -script cad_checks.wls
   The metadata file defines runData with keys "Curve" (polynomial in x, y with
   exact rational coefficients) and "Box" ({{xmin, xmax}, {ymin, ymax}} with
   exact rational entries). Results are printed as CAD_* key-value lines parsed
   by the Julia driver; failures are reported as CAD_ERROR lines. *)

abort[message_] := (Print["CAD_ERROR ", message]; Exit[1]);

(* A raw kernel invoked with -script leaves $ScriptCommandLine empty, so the
   paths are taken from the environment first and from the script arguments
   only as a fallback. *)
packagePath = Environment["WFS_CAD_PACKAGE"];
metadataPath = Environment["WFS_CAD_METADATA"];
If[!StringQ[packagePath] || !StringQ[metadataPath],
  If[Length[$ScriptCommandLine] >= 3,
    {packagePath, metadataPath} = Take[$ScriptCommandLine, -2],
    abort["no input paths: set WFS_CAD_PACKAGE and WFS_CAD_METADATA, " <>
          "or pass them as the last two script arguments"]]];

If[!FileExistsQ[packagePath], abort["package not found: " <> packagePath]];
If[!FileExistsQ[metadataPath], abort["metadata not found: " <> metadataPath]];
If[Get[packagePath] === $Failed, abort["failed to load " <> packagePath]];
If[Get[metadataPath] === $Failed, abort["failed to load " <> metadataPath]];

If[Head[runData] =!= Association, abort["metadata did not define runData"]];
If[!KeyExistsQ[runData, "Curve"], abort["metadata has no Curve entry"]];
If[!KeyExistsQ[runData, "Box"], abort["metadata has no Box entry"]];
If[Length[DownValues[certifyBox]] == 0, abort["certifyBox is not defined"]];
If[Length[DownValues[countMaximaInBox]] == 0, abort["countMaximaInBox is not defined"]];

f = runData["Curve"];
box = runData["Box"];
If[!MatchQ[box, {{_, _}, {_, _}}], abort["Box is not of the form {{a,b},{c,d}}"]];

boxResult = certifyBox[f, box, {x, y}];
maximaResult = countMaximaInBox[f, box, {x, y}];
If[Head[boxResult] =!= Association, abort["certifyBox did not evaluate"]];
If[Head[maximaResult] =!= Association, abort["countMaximaInBox did not evaluate"]];

Print["CAD_COMPONENTS ", boxResult["Components"]];
Print["CAD_EDGES_CLEAR ", boxResult["EdgesClear"]];
Print["CAD_BLOCKED_EDGES ",
  If[boxResult["BlockedEdges"] === {}, "none",
    StringRiffle[boxResult["BlockedEdges"], ","]]];
Print["CAD_CONTAINS_TEST_POINT ", boxResult["ContainsTestPoint"]];
Print["CAD_BOX_CERTIFIED ", boxResult["Certified"]];
Print["CAD_MAXIMA_NONDEGENERATE ", maximaResult["Nondegenerate"]];
Print["CAD_MAXIMA_DEGENERATE ", maximaResult["Degenerate"]];
Print["CAD_MAXIMA_TOTAL ", maximaResult["Total"]];
Print["CAD_MAXIMA_ALL_ISOLATED ", maximaResult["AllIsolated"]];
