# CLAUDE.md

HCIT/coronagraph instrument software (FALCO-based EFC testbed control + reporting).

## Class/file layout conventions
- `C`-prefixed classdef files = handle classes: `CRunData` (generic base,
  in "HCIT utils toolbox/CRunData.m") is the base for run/iteration data;
  `CfalcoRunData < CRunData` (same folder) is the FALCO-specific concrete
  subclass. If EFC software ever replaces falco, a new CRunData subclass
  replaces CfalcoRunData — don't hardcode CfalcoRunData-specific behavior
  into generic code paths.
- `CConstants` (classes toolbox/CConstants.m) holds constants like `NM`;
  `CRunData < handle & CConstants` so `S.NM` etc. work by inheritance.
- `Cppt` (PPT ActiveX toolbox/Cppt.m) wraps PowerPoint ActiveX automation:
  methods `NewSlide(slide_num)` (append if `[]`), `CopyFigSlide(slide,hfig)`,
  `CopyFigNewSlide(hfig)`, `AddPictureNewSlide(fn)`. Windows/ActiveX only
  (`ispc`). Note: `classdef Cppt` is NOT `< handle` (value class), but its
  properties are COM object handles so mutation-by-reference still works in
  practice. Constructor supports `'open'` option (default `false`) to reopen
  an existing .pptx via `invoke(S.Presentations, 'Open', fn)` instead of
  creating a new blank presentation.
- `CheckOption(varstring, defaultval, varargin{:})` is a shared utility
  function in "functions toolbox/CheckOption.m" (on the Matlab path,
  not a local function) — the standard varargin option-parsing idiom used
  everywhere in this codebase.
- `CEFCReport` (HCIT utils toolbox/CEFCReport.m) is the class-based EFC
  report generator (replaced the old free-function GenerateEFCReport_falco.m,
  which is now a thin backward-compat wrapper). It's a handle class: a
  constructed `R` persists in the workspace, so later work (e.g. plotting a
  specific iteration subset, or calling `R.Update()` after the testbed
  produces more iterations) is just further method calls on the same `R`,
  not re-passing `listS`/`Sppt` the way the old free-function `GenReport`
  required.
- OMC_MSWC repo (separate from this Matlab repo) is at
  `C:\Users\dmarx\Documents\OMC_MSWC`; `run_testbed/GenReport.m` there is a
  thin preset layer on top of GenerateEFCReport_falco.m.

## MATLAB gotchas
- A function handle `@obj.method` checks method access (public/private) at
  handle-CREATION time, not call time. You cannot create `@obj.method` from
  outside the classdef if `method` is private — make it public if external
  code needs a handle to it.
- MATLAB classdef hot-reloads an edited .m file for an already-existing
  handle instance within the same session, as long as properties/attributes
  are unchanged (adding/modifying methods is fine). No need to `clear
  classes` for that kind of change. One case where `clear classes` would be
  destructive: it drops a live Cppt COM handle inside a CEFCReport's R.Sppt.
- `CfalcoRunData`/`CRunData` objects cache `S.falcoData = out` (the run's
  cumulative `_snippet.mat`) fresh at construction time — it's a frozen
  snapshot, not live. Methods that need the latest cumulative run history
  (e.g. normIntHist-type arrays) must read from `listS(end).falcoData`
  (the most-recently-constructed/loaded iteration), NOT `listS(1)` (the
  first-loaded, likely stale). `CEFCReport.PlotBeta` uses `listS(end)`
  correctly — use it as the reference pattern. (Bug of this shape was found
  and fixed in `CEFCReport.PlotNormIntensity`, which caused summary plots to
  silently truncate at a stale iteration count after `R.Update()`.)

## Workflow notes
- User sometimes captures a plan-mode plan into a manual .txt snapshot if a
  session gets interrupted mid-execution. If asked to "continue where we
  left off" and nothing obvious is in progress, ask the user if they saved
  a plan snapshot before re-deriving one from scratch.
