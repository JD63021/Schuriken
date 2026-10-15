/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "AllReduce.H"
#include "List.H"
#include "Pstream.H"
#include "ops.H"

#ifdef SCHURIKEN_HAVE_MPI
    #include "PstreamGlobals.H"
    #include <mpi.h>

    #if defined(WM_SP)
        #define SCHURIKEN_MPI_SCALAR MPI_FLOAT
    #elif defined(WM_DP)
        #define SCHURIKEN_MPI_SCALAR MPI_DOUBLE
    #elif defined(WM_LP)
        #define SCHURIKEN_MPI_SCALAR MPI_LONG_DOUBLE
    #endif
#endif

void Foam::schuriken::sumAllReduce
(
    scalar* values,
    const label n,
    const label comm
)
{
    if (!Pstream::parRun() || n <= 0)
    {
        return;
    }

#ifdef SCHURIKEN_HAVE_MPI
    // The communicator table is the one OpenFOAM's own reductions use
    if
    (
        MPI_Allreduce
        (
            MPI_IN_PLACE,
            values,
            n,
            SCHURIKEN_MPI_SCALAR,
            MPI_SUM,
            PstreamGlobals::MPICommunicators_[comm]
        )
    )
    {
        FatalErrorInFunction
            << "MPI_Allreduce failed"
            << Foam::abort(FatalError);
    }
#else
    List<scalar> list(n);
    forAll(list, i)
    {
        list[i] = values[i];
    }
    Pstream::listCombineGather
    (
        list,
        plusEqOp<scalar>(),
        UPstream::msgType(),
        comm
    );
    Pstream::listCombineScatter(list, UPstream::msgType(), comm);
    forAll(list, i)
    {
        values[i] = list[i];
    }
#endif
}
