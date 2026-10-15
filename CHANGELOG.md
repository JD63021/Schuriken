# Changelog

## Unreleased

- Fix spectral recompression: OpenFOAM pressure matrices are negative
  definite, so `W^T A W` failed the positive-definite check and recompression
  never ran. `definiteInvSqrt` handles either sign. Recompression is now opt-in
  (`recompress yes`); the default keeps the full basis at `maxRank`, which is
  what the published results actually used. If the spectral selection is
  impossible the newest vectors are kept instead of freezing the basis.
- A singular or non-finite `W^T A W` no longer aborts the run: the basis is
  reduced to its numerically independent A-orthogonal part, or cleared so that
  the solve falls back to stock PCG + GAMG. `G^-1` is symmetrized.
- Cheaper iterations: `z = g + W (a - b)` without a separate coarse vector,
  cache-blocked multi-vector kernels (deflation kernels 1.6-2x faster; Y24
  serial pressure time about -5% on the same solve sequence), the trainer
  reuses the solve's `A W`/`G` and grows it one column at a time, and the
  eigen-solver uses cyclic Jacobi sweeps.
- Three collectives per recycled PCG iteration instead of five; batches use a
  single `MPI_Allreduce` when built against MPI.
- Explicit `learn yes|no` per solver dictionary. One adaptive state per field
  and mesh region, held in the mesh object registry (two `foamMultiRun` regions
  of equal size previously shared one basis).
- `persist yes` writes the basis to the time directories and restores it on
  restart.
- Per-solve and cumulative timing (`seconds`, `setupSeconds`, `iterSeconds`,
  `trainSeconds`) and region/field names in the log markers.
- Tests: linear-algebra unit tests, forced recompression, preconditioner
  symmetry (serial and MPI), learn switch, multi-region, restart. Test scripts
  move old case directories aside instead of deleting them.

## 0.3.0-dev

- Refactor the single large `PCGNinja.C` implementation into solver, AMGsensei,
  recycled-preconditioner, and linear-algebra modules.
- Preserve the v0.2 PCG/GAMG/POD algorithm and OpenFOAM-facing controls.
- Document native OpenFOAM tolerance and MPI behavior.
- Add Y24 16-rank validation results.

## 0.2.0

- Add distributed-memory MPI support for the adaptive POD basis.
- Batch reduced-vector and reduced-matrix global communications.
- Add 2-rank MPI cavity smoke test.

## 0.1.0

- Initial serial `pcgninja` + `AMGsensei` package.
