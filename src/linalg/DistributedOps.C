/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "DistributedOps.H"
#include "AllReduce.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"

#include <cmath>

namespace
{

using Foam::label;
using Foam::scalar;

// Cells per block in the fused multi-vector loops: one block of the shared
// vector stays in L1 while every basis vector streams past it once.
const label blockSize = 512;

// Four partial sums break the floating-point add dependency chain
inline scalar dotRange
(
    const scalar* __restrict__ a,
    const scalar* __restrict__ b,
    const label start,
    const label end
)
{
    scalar s0 = 0.0;
    scalar s1 = 0.0;
    scalar s2 = 0.0;
    scalar s3 = 0.0;
    label i = start;
    for (; i + 3 < end; i += 4)
    {
        s0 += a[i]*b[i];
        s1 += a[i+1]*b[i+1];
        s2 += a[i+2]*b[i+2];
        s3 += a[i+3]*b[i+3];
    }
    for (; i < end; ++i)
    {
        s0 += a[i]*b[i];
    }
    return (s0 + s1) + (s2 + s3);
}

}


Foam::scalar Foam::schuriken::dotLocal
(
    const scalarField& a,
    const scalarField& b
)
{
    return dotRange(a.cdata(), b.cdata(), 0, a.size());
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

std::vector<Foam::scalar> Foam::schuriken::dotBatchLocal
(
    const std::vector<scalarField>& V,
    const scalarField& x
)
{
    const label n = x.size();
    std::vector<scalar> values(V.size(), 0.0);

    for (label start=0; start<n; start+=blockSize)
    {
        const label end = min(start + blockSize, n);
        for (std::size_t k=0; k<V.size(); ++k)
        {
            values[k] += dotRange(V[k].cdata(), x.cdata(), start, end);
        }
    }
    return values;
}

std::vector<Foam::scalar> Foam::schuriken::dotBatchGlobal
(
    const std::vector<scalarField>& left,
    const scalarField& right,
    const label comm
)
{
    std::vector<scalar> values = dotBatchLocal(left, right);
    sumAllReduce(values.data(), static_cast<label>(values.size()), comm);
    return values;
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
    const label n = nr ? left[0].size() : 0;
    std::vector<scalar> values(static_cast<std::size_t>(nr*nc), 0.0);

    for (label start=0; start<n; start+=blockSize)
    {
        const label end = min(start + blockSize, n);
        for (label i=0; i<nr; ++i)
        {
            const scalar* li = left[static_cast<std::size_t>(i)].cdata();
            for (label j=0; j<nc; ++j)
            {
                values[static_cast<std::size_t>(i*nc+j)] += dotRange
                (
                    li,
                    right[static_cast<std::size_t>(j)].cdata(),
                    start,
                    end
                );
            }
        }
    }

    sumAllReduce(values.data(), nr*nc, comm);
    return values;
}

void Foam::schuriken::addCombination
(
    scalarField& y,
    const std::vector<scalarField>& V,
    const std::vector<scalar>& c
)
{
    const label n = y.size();
    scalar* __restrict__ yp = y.data();

    for (label start=0; start<n; start+=blockSize)
    {
        const label end = min(start + blockSize, n);
        for (std::size_t k=0; k<V.size(); ++k)
        {
            const scalar ck = c[k];
            const scalar* __restrict__ vk = V[k].cdata();
            for (label i=start; i<end; ++i)
            {
                yp[i] += ck*vk[i];
            }
        }
    }
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
