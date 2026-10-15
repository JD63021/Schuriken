# Schuriken architecture

Schuriken is split so the OpenFOAM runtime solver, adaptive learning, recycled
preconditioning, and low-level algebra can evolve independently.

## Solve path

`PCGNinja` owns the normal PCG iteration and OpenFOAM solver interface.

For every pressure solve it:

1. reads ordinary OpenFOAM solver controls and the `AMGsensei` configuration;
2. looks up the adaptive state of this field in its mesh's object registry
   (created on first use, read back from the start time with `persist`);
3. constructs the configured OpenFOAM preconditioner (normally GAMG);
4. builds the `ReducedSpace` (`A W`, `W^T A W` and its inverse) for the
   current matrix, repairing the basis if `W^T A W` is singular;
5. executes PCG using the balanced recycled preconditioner `RecycledGAMG`;
6. passes the completed solve and its reduced space to `AMGsenseiTrainer` for
   optional learning.

The OpenFOAM `solverPerformance::checkConvergence` path remains responsible for
`tolerance`, `relTol`, and `maxIter` behavior.

## Recycled preconditioner

For basis `W`, Schuriken forms

```text
AW = A W
G  = W^T A W
```

and applies the balanced correction around stock GAMG (`M`):

```text
a = G^-1 W^T r,   t = r - AW a,   g = M t,   b = G^-1 (AW)^T g
z = g + W (a - b)          = Q r + (I - Q A) M (I - A Q) r,  Q = W G^-1 W^T
```

`G` and its inverse are small replicated dense objects; `W` and `AW` remain
distributed. `G^-1` is symmetrized so `z` is a symmetric function of `r` to
round-off, as PCG requires (`tests/symmetry_cavity.sh`).

OpenFOAM pressure matrices are negative definite, so `G` normally is too. Code
that needs a definite square root (`definiteInvSqrt`) works on `sign*G`.

## AMGsensei training

The trainer learns on the solves selected by `learn` or, without it, by
`learnOddOnly`. Initial training uses full solve corrections. After a basis
exists, the trainer uses the tail correction after `tailIter` PCG iterations.

A small snapshot POD proposes new modes. Candidate modes are A-orthogonalized
against the current basis using the solve's reduced space; each accepted mode
borders `G` with one new row and column (`appendToReducedSpace`), so a mode
costs two matrix products regardless of the basis size.

At `maxRank`, with `recompress yes`, the trainer recompresses to
`recompressRank` by generalized Ritz values `mu` of `M A` on `span(W)` and keeps
the smallest (the modes GAMG reduces worst), then lets the basis grow again.
With `recompress no` (default) the full basis is kept and learning stops.

## State

`AMGsenseiState` holds one field's basis, pending snapshots and counters. It
lives in a `regIOobject` named `AMGsensei_<field>` in the mesh registry, so each
region and field is independent and the state is destroyed with its mesh. With
`persist` the object is `AUTO_WRITE` and written in binary at write times; a
restart reads it from the start time if the cell count matches.

## MPI model

Every rank holds the local cells for every basis vector. Distributed operations
live in `linalg/DistributedOps` and `linalg/AllReduce`.

A recycled PCG iteration needs three collectives, as stock PCG:

- `[sum|r|, W^T r]` after the residual update;
- `[(AW)^T g, r.g]` inside the preconditioner, from which `r.z` follows
  without another reduction;
- `p.Ap`.

Setup batches `W^T A W` and related cross-Gram products into one reduced-matrix
collective. Built against MPI, each batch is one `MPI_Allreduce` on OpenFOAM's
communicator; otherwise a gather/scatter is used. The learned vectors are never
gathered to a master rank.

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
          DistributedOps / AllReduce
          DenseOps / PODBasis / ReducedSpace
```

A future matrix-free Krylov backend should reuse `DistributedOps` for batched
Arnoldi reductions and can reuse the existing learned basis/trainer.
