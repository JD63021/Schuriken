# Y24 pipe MPI-16 example

This is the validated 22,848-cell turbulent pipe case used for the Schuriken MPI regression.

The example is an ordinary OpenFOAM Foundation v12 case. There are no Schuriken-specific command-line arguments.

## Configuration

- pressure solver: `pcgninja`
- preconditioner: stock OpenFOAM GAMG
- AU/AP/AK/AO = 0.70 / 0.70 / 0.70 / 0.70
- pressure tolerance = 1e-5
- pressure relTol = 0.10
- outer residual gate: U, k, omega = 1e-4
- validated decomposition: 16 Scotch subdomains

## Run

```bash
cd ~/Schuriken
./Allwmake
cd examples/Y24_MPI16
./Allrun
```

A different rank count can be tried with `NP=8 ./Allrun`, but the published validation target is MPI-16.

## Validated result

```text
MPI ranks                 16
outer iterations          126
pressure solves           252
total PCG iterations      5739
stock/learning PCG        806
POD PCG                   4933
mean PCG / solve          22.7738095
max PCG / solve           318
final recycled rank       24
accepted updates          12
recompress count          0
DeltaP                    18834.856899
mean y+                   24.0449682316
wall time                 6.146348116 s
```

Wall time is machine-dependent. Convergence history and final physics are the stronger regression checks.

## Clean

```bash
./Allclean
```
