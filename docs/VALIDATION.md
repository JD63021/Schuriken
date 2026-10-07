# Validation notes

## Smoke tests

The bundled smoke sequence validates:

- OpenFOAM runtime registration of `pcgninja`;
- library loading through `controlDict`;
- stock GAMG construction beneath `pcgninja`;
- AMGsensei snapshot collection and rank growth;
- recycled pressure solves after training;
- decomposed MPI operation on 2 ranks.

## Y24 MPI-16 regression

The 22,848-cell VMFL003 Y24 pipe case was run on 16 MPI ranks with the same
normal OpenFOAM solver controls used in serial.

Result:

```text
outer iterations          126
pressure solves           252
total PCG iterations      5739
stock/learning PCG        806
POD PCG                   4933
mean PCG / solve          22.7738
max PCG / solve           318
DeltaP                    18834.856899
mean y+                   24.0449682316
rank                      24
accepted updates          12
wall time                 6.278 s
```

Representative serial reference:

```text
outer iterations          126
pressure solves           252
total PCG iterations      4115
DeltaP                    18835.8283963
mean y+                   24.0450293652
rank                      18
accepted updates          9
wall time                 14.637 s
```

The parallel solution preserved the nonlinear iteration count and final physics.
The parallel GAMG hierarchy required more Krylov work and a larger learned
basis, but the total wall time still decreased by about 2.33x.
