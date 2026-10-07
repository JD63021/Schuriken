# Schuriken architecture

Schuriken is split so the OpenFOAM runtime solver, adaptive learning, recycled
preconditioning, and low-level algebra can evolve independently.

## Solve path

`PCGNinja` owns the normal PCG iteration and OpenFOAM solver interface.

For every pressure solve it:

1. reads ordinary OpenFOAM solver controls and the `AMGsensei` configuration;
2. constructs the configured OpenFOAM preconditioner (normally GAMG);
3. constructs `RecycledGAMG` from the current learned basis;
4. executes PCG using the balanced recycled preconditioner;
5. passes the completed solve to `AMGsenseiTrainer` for optional learning.

The OpenFOAM `solverPerformance::checkConvergence` path remains responsible for
`tolerance`, `relTol`, and `maxIter` behavior.

## Recycled preconditioner

For basis `W`, Schuriken forms

```text
AW = A W
G  = W^T A W
```

and applies the balanced correction around stock GAMG. `G` and its inverse are
small replicated dense objects; `W` and `AW` remain distributed.

## AMGsensei training

The trainer learns only on configured solves (`learnOddOnly=1` by default).
Initial training uses full solve corrections. After a basis exists, the trainer
can use the tail correction after `tailIter` PCG iterations.

A small snapshot POD proposes new modes. Candidate modes are A-orthogonalized
against the current basis. At `maxRank`, the current implementation attempts a
GAMG-informed spectral recompression to `recompressRank` and then allows the
basis to grow again.

## MPI model

Every rank holds the local cells for every basis vector. Distributed operations
live in `linalg/DistributedOps`.

The implementation batches:

- `W^T r` into one reduced-vector collective;
- `W^T A W` and related cross-Gram products into one reduced-matrix collective.

Thus the learned vectors are never gathered to a master rank.

## Extension points

The current separation is intended to support future solvers without copying
AMGsensei:

```text
PCGNinja          -> symmetric PCG backend
BiCGStabNinja     -> future nonsymmetric backend
FGMRESNinja       -> future flexible/matrix-free backend
        |                 |
        +------ AMGsensei-+
                 |
          DistributedOps
          DenseOps / PODBasis
```

A future matrix-free Krylov backend should reuse `DistributedOps` for batched
Arnoldi reductions and can reuse the existing learned basis/trainer.
