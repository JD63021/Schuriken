# Schuriken

**Adaptive recycled pressure solves for OpenFOAM Foundation v12.**

Schuriken provides a runtime-selectable symmetric pressure solver:

```text
solver pcgninja;
```

`pcgninja` keeps the standard OpenFOAM **PCG + GAMG** solve path and adds an
adaptive recycled POD coarse correction. The learning/recycling component is
called **AMGsensei**.

The current implementation supports both serial and decomposed MPI runs. No
Schuriken-specific MPI launch mode or nonlinear convergence mechanism is
required.

## Why Schuriken

Repeated pressure systems in steady SIMPLE/PIMPLE calculations are closely
related from one outer iteration to the next. AMGsensei learns a small
subspace from earlier pressure corrections and inserts a balanced recycled
correction around the normal GAMG preconditioner.

The large POD vectors stay distributed in MPI. Only small reduced vectors and
matrices are globally reduced.

## Requirements

- OpenFOAM Foundation **v12**
- GCC/OpenFOAM build environment matching the local installation
- `mpirun` for parallel use

The development system used `/opt/openfoam12`.

## Build

```bash
source /opt/openfoam12/etc/bashrc
cd ~/Schuriken
./Allwmake
```

The library is installed as:

```text
$FOAM_USER_LIBBIN/libSchuriken.so
```

Run the bundled serial and MPI smoke tests with:

```bash
./install_and_test.sh
```

## Normal OpenFOAM use

Load the library in `system/controlDict`:

```text
libs
(
    "libSchuriken.so"
);
```

Use `pcgninja` in `system/fvSolution`:

```text
p
{
    solver          pcgninja;

    preconditioner
    {
        preconditioner          GAMG;
        smoother                symGaussSeidel;
        nVcycles                1;

        cacheAgglomeration      true;
        agglomerator            faceAreaPair;
        nCellsInCoarsestLevel   10;
        mergeLevels             1;

        nPreSweeps              0;
        nPostSweeps             2;
        nFinestSweeps           2;
    }

    tolerance       1e-5;
    relTol          0.1;
    maxIter         500;

    AMGsensei
    {
        enabled             true;
        trainSnapshots      3;
        addRank             2;
        maxRank             24;
        recompressRank      16;
        tailIter            10;
        learnOddOnly        1;
    }
}

pFinal
{
    $p;
    relTol          0;
}
```

### AMGsensei controls

| Entry | Default | Meaning |
|---|---|---|
| `enabled` | `true` | `false` gives plain PCG + GAMG |
| `trainSnapshots` | 3 | snapshots per POD batch |
| `addRank` | 2 | modes added per batch |
| `maxRank` | 24 | basis size limit |
| `recompress` | `no` | at `maxRank`: `no` keeps the full basis and stops learning; `yes` recompresses to `recompressRank` and keeps learning |
| `recompressRank` | 16 | size after a spectral recompression |
| `tailIter` | 10 | PCG depth of the tail snapshot |
| `learn` | (unset) | `yes`/`no` per solver dictionary; overrides `learnOddOnly` |
| `learnOddOnly` | 1 | without `learn`: learn on odd solves of the field (1) or every solve (0) |
| `persist` | `no` | write the basis to the time directories and read it back on restart |
| `checkSymmetry` | `false` | test diagnostic: print the asymmetry of the recycled preconditioner |

`learn` replaces the odd/even rule when the number of pressure solves per outer
iteration is not two. OpenFOAM merges sub-dictionaries, so `pFinal` can inherit
everything and switch learning off alone:

```text
pFinal
{
    $p;
    relTol          0;
    AMGsensei { learn no; }
}
```

Each field keeps its own adaptive state in its mesh region, so several fields or
regions (e.g. `foamMultiRun`) learn independently.

With `persist yes` the basis of field `p` is written in binary as
`<time>/AMGsensei_p` (per processor in parallel) at every write time. A restart
reads it from the start time when the cell count matches; otherwise learning
starts afresh. Set it in the first-solved dictionary (`p`; `$p` passes it on).

`recompress` defaults to `no` because that is how the published results were
obtained: before the sign fix in this version recompression never ran, so
the basis froze at `maxRank`. On the Y24 pipe, enabling it cost about 13% more
PCG iterations on 2 ranks. It may suit long transients, where an old basis can
go stale.

### Tolerances

Schuriken does **not** add a second tolerance system.

The pressure linear solve uses the standard OpenFOAM entries:

```text
tolerance
relTol
maxIter
```

Outer CFD convergence remains the normal OpenFOAM mechanism, for example:

```text
PIMPLE
{
    residualControl
    {
        U       1e-4;
        k       1e-4;
        omega   1e-4;
    }
}
```

Benchmark wrapper variables such as `SIMPLE_TOL` or `P_RELTOL` are only script
conveniences used to write these ordinary OpenFOAM dictionaries. They are not
required by Schuriken.

## MPI use

There are no MPI-specific Schuriken controls.

Use the case exactly as a normal decomposed OpenFOAM case:

```bash
decomposePar -force
mpirun -np 16 foamRun -parallel
```

Each rank stores its local part of every recycled vector. Global POD inner
products and reduced matrices are formed with batched collectives on the mesh
communicator. The full pressure or POD vectors are never gathered to one rank.

A recycled PCG iteration needs three collectives, the same as stock PCG:
`W^T r` travels with the residual norm and `(AW)^T M t` with `r.z`. When
OpenFOAM is built with MPI each batch is one `MPI_Allreduce` (OpenFOAM-12
itself only all-reduces single scalars); otherwise a gather/scatter is used.


## Reproducible Y24 example

A complete 22,848-cell Y24 pipe case, including the OpenFOAM mesh, is provided under:

```text
examples/Y24_MPI16/
```

After building Schuriken, reproduce the validated 16-rank run with:

```bash
cd ~/Schuriken/examples/Y24_MPI16
./Allrun
```

The example is a normal OpenFOAM case. `pcgninja` is selected in `fvSolution`, `libSchuriken.so` is loaded from `controlDict`, and MPI is launched normally. No Schuriken-specific tolerance or parallel arguments are required.

## Repository layout

```text
src/
├── solver/
│   ├── PCGNinja.H
│   └── PCGNinja.C
├── sensei/
│   ├── AMGsensei.H
│   ├── AMGsensei.C
│   ├── AMGsenseiState.H
│   ├── AMGsenseiState.C
│   ├── AMGsenseiTrainer.H
│   └── AMGsenseiTrainer.C
├── preconditioner/
│   ├── RecycledGAMG.H
│   └── RecycledGAMG.C
├── linalg/
│   ├── AllReduce.H
│   ├── AllReduce.C
│   ├── DistributedOps.H
│   ├── DistributedOps.C
│   ├── DenseOps.H
│   ├── DenseOps.C
│   ├── PODBasis.H
│   ├── PODBasis.C
│   ├── ReducedSpace.H
│   └── ReducedSpace.C
└── Make/
```

The responsibilities are intentionally separated:

- `solver/`: OpenFOAM runtime registration and the PCG iteration.
- `sensei/`: user configuration, per-field adaptive state (registry-held,
  optionally persisted), snapshot training, basis augmentation and spectral
  recompression.
- `preconditioner/`: balanced recycled correction wrapped around stock GAMG.
- `linalg/`: batched distributed reductions, cache-blocked multi-vector kernels,
  small dense algebra, snapshot POD and the reduced space `A W`, `W^T A W`.

See [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) for details.

## Validation

### Serial/MPI smoke

The stock OpenFOAM cavity test has passed in serial and on 2 MPI ranks. The MPI
smoke confirms that AMGsensei forms a common distributed basis and then uses it
inside subsequent recycled pressure solves.

### Y24 pipe, 16 MPI ranks

Representative settings:

```text
AU/AP/AK/AO = 0.70 / 0.70 / 0.70 / 0.70
p tolerance = 1e-5
p relTol    = 0.10
outer gate  = U,k,omega < 1e-4
```

Observed MPI-16 result:

```text
outer iterations          126
pressure solves           252
total PCG iterations      5739
final recycled rank       24
accepted updates          12
DeltaP                    18834.856899
mean y+                   24.0449682316
wall time                 6.15 s
```

Representative serial reference from the development campaign:

```text
outer iterations          126
pressure solves           252
total PCG iterations      4115
final recycled rank       18
accepted updates          9
DeltaP                    18835.8283963
mean y+                   24.0450293652
wall time                 14.64 s
```

The MPI run reproduced the same outer iteration count and essentially the same
physics while reducing wall time by about **2.38x** on this relatively small
22,848-cell case. MPI required a richer recycled space and more Krylov work,
which is expected because decomposition changes the parallel GAMG hierarchy and
floating-point reduction order.

Additional development data are in `benchmarks/results.csv`.

## Current scope

- OpenFOAM Foundation v12
- symmetric LDU systems
- PCG + stock GAMG + adaptive recycled POD
- serial and distributed-memory MPI
- one adaptive state per field and mesh region, optionally persisted for
  restarts (same decomposition only)

The experimental polynomial/matrix-free FGMRES work is intentionally not part
of the default solver yet. The modular layout is designed so future BiCGStab or
FGMRES backends can reuse AMGsensei and the distributed linear-algebra layer.

## Tests

Static repository audit:

```bash
./tests/static_selfcheck.sh
```

Serial smoke:

```bash
./tests/smoke_cavity.sh
```

MPI smoke:

```bash
./tests/smoke_cavity_mpi.sh 2
```

Algorithm and feature tests (cavity based, a few seconds each; `NP` runs the
case on `NP` MPI ranks):

```bash
./tests/unit_linalg.sh             # dense/distributed algebra unit tests
./tests/recompress_cavity.sh       # spectral recompression runs; recompress no keeps the basis
./tests/symmetry_cavity.sh [NP]    # recycled preconditioner symmetric to round-off
./tests/learn_cavity.sh            # explicit learn switch
./tests/multiregion_cavity.sh      # one state per foamMultiRun region
./tests/restart_cavity.sh [NP]     # persist: basis written and restored
```

Test cases are created under `/tmp` (override with `SCHURIKEN_TEST_CASE`); an
existing case directory is moved aside to `<case>.oldN`, never deleted.

Full build + test sequence:

```bash
./install_and_test.sh
```

## License

GPL-3.0-or-later. See `LICENSE`.
