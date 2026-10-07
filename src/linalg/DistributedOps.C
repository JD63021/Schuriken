/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "DistributedOps.H"
#include "List.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "ops.H"

#include <cmath>

Foam::scalar Foam::schuriken::dotLocal
(
    const scalarField& a,
    const scalarField& b
)
{
    scalar s = 0.0;
    forAll(a, i)
    {
        s += a[i]*b[i];
    }
    return s;
}

Foam::scalar Foam::schuriken::dotGlobal
(
    const scalarField& a,
    const scalarField& b,
    const label comm
)
{
    scalar s = dotLocal(a, b);
    reduce(s, sumOp<scalar>(), UPstream::msgType(), comm);
    return s;
}

Foam::scalar Foam::schuriken::normGlobal
(
    const scalarField& a,
    const label comm
)
{
    const scalar s = dotGlobal(a, a, comm);
    return std::sqrt(max(s, scalar(0)));
}

std::vector<Foam::scalar> Foam::schuriken::dotBatchGlobal
(
    const std::vector<scalarField>& left,
    const scalarField& right,
    const label comm
)
{
    List<scalar> values(static_cast<label>(left.size()), scalar(0));
    for (label i=0; i<values.size(); ++i)
    {
        values[i] = dotLocal(left[static_cast<std::size_t>(i)], right);
    }

    if (Pstream::parRun())
    {
        Pstream::listCombineGather
        (
            values,
            plusEqOp<scalar>(),
            UPstream::msgType(),
            comm
        );
        Pstream::listCombineScatter(values, UPstream::msgType(), comm);
    }

    std::vector<scalar> result(static_cast<std::size_t>(values.size()), 0.0);
    forAll(values, i)
    {
        result[static_cast<std::size_t>(i)] = values[i];
    }
    return result;
}

std::vector<Foam::scalar> Foam::schuriken::crossGramGlobal
(
    const std::vector<scalarField>& left,
    const std::vector<scalarField>& right,
    const label comm
)
{
    const label nr = static_cast<label>(left.size());
    const label nc = static_cast<label>(right.size());
    List<scalar> values(nr*nc, scalar(0));

    for (label i=0; i<nr; ++i)
    {
        for (label j=0; j<nc; ++j)
        {
            values[i*nc+j] = dotLocal
            (
                left[static_cast<std::size_t>(i)],
                right[static_cast<std::size_t>(j)]
            );
        }
    }

    if (Pstream::parRun())
    {
        Pstream::listCombineGather
        (
            values,
            plusEqOp<scalar>(),
            UPstream::msgType(),
            comm
        );
        Pstream::listCombineScatter(values, UPstream::msgType(), comm);
    }

    std::vector<scalar> result(static_cast<std::size_t>(nr*nc), 0.0);
    forAll(values, i)
    {
        result[static_cast<std::size_t>(i)] = values[i];
    }
    return result;
}

bool Foam::schuriken::normalizeEuclid
(
    scalarField& a,
    const label comm,
    const scalar tol
)
{
    const scalar n = normGlobal(a, comm);
    if (!(n > tol) || !std::isfinite(n))
    {
        return false;
    }

    forAll(a, i)
    {
        a[i] /= n;
    }
    return true;
}
