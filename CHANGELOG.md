# Changelog

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
