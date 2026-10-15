/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "DenseOps.H"

#include <algorithm>
#include <cmath>
#include <limits>

namespace Foam
{
namespace schuriken
{

std::vector<scalar> identityDense(const label n)
{
    std::vector<scalar> I(static_cast<std::size_t>(n*n), 0.0);
    for (label i=0; i<n; ++i)
    {
        I[static_cast<std::size_t>(i*n+i)] = 1.0;
    }
    return I;
}

void symmetrize(std::vector<scalar>& A, const label n)
{
    for (label i=0; i<n; ++i)
    {
        for (label j=i+1; j<n; ++j)
        {
            const scalar s = 0.5*
            (
                A[static_cast<std::size_t>(i*n+j)]
              + A[static_cast<std::size_t>(j*n+i)]
            );
            A[static_cast<std::size_t>(i*n+j)] = s;
            A[static_cast<std::size_t>(j*n+i)] = s;
        }
    }
}

bool jacobiEigenSymmetric
(
    std::vector<scalar> A,
    const label n,
    std::vector<scalar>& eval,
    std::vector<scalar>& V
)
{
    V = identityDense(n);
    if (n <= 0)
    {
        eval.clear();
        return false;
    }

    scalar scale = 0.0;
    for (label i=0; i<n; ++i)
    {
        scale = Foam::max(scale, Foam::mag(A[static_cast<std::size_t>(i*n+i)]));
    }
    scale = Foam::max(scale, scalar(1e-300));

    // Zero A(p,q) by a plane rotation and accumulate it into V
    const auto rotate = [&](const label p, const label q)
    {
        const scalar app = A[static_cast<std::size_t>(p*n+p)];
        const scalar aqq = A[static_cast<std::size_t>(q*n+q)];
        const scalar apq = A[static_cast<std::size_t>(p*n+q)];

        const scalar phi = 0.5*std::atan2(2.0*apq, aqq-app);
        const scalar c = std::cos(phi);
        const scalar s = std::sin(phi);

        for (label k=0; k<n; ++k)
        {
            if (k == p || k == q) continue;

            const scalar akp = A[static_cast<std::size_t>(k*n+p)];
            const scalar akq = A[static_cast<std::size_t>(k*n+q)];

            const scalar nkp = c*akp - s*akq;
            const scalar nkq = s*akp + c*akq;

            A[static_cast<std::size_t>(k*n+p)] = nkp;
            A[static_cast<std::size_t>(p*n+k)] = nkp;
            A[static_cast<std::size_t>(k*n+q)] = nkq;
            A[static_cast<std::size_t>(q*n+k)] = nkq;
        }

        A[static_cast<std::size_t>(p*n+p)] =
            c*c*app - 2.0*s*c*apq + s*s*aqq;
        A[static_cast<std::size_t>(q*n+q)] =
            s*s*app + 2.0*s*c*apq + c*c*aqq;
        A[static_cast<std::size_t>(p*n+q)] = 0.0;
        A[static_cast<std::size_t>(q*n+p)] = 0.0;

        for (label k=0; k<n; ++k)
        {
            const scalar vkp = V[static_cast<std::size_t>(k*n+p)];
            const scalar vkq = V[static_cast<std::size_t>(k*n+q)];
            V[static_cast<std::size_t>(k*n+p)] = c*vkp - s*vkq;
            V[static_cast<std::size_t>(k*n+q)] = s*vkp + c*vkq;
        }
    };

    // Cyclic Jacobi: rotate every off-diagonal pair once per sweep instead of
    // searching the whole matrix for the largest element before each rotation
    const label maxSweeps = 100;

    for (label sweep=0; sweep<maxSweeps; ++sweep)
    {
        scalar mx = 0.0;
        for (label i=0; i<n; ++i)
        {
            for (label j=i+1; j<n; ++j)
            {
                mx = Foam::max(mx, Foam::mag(A[static_cast<std::size_t>(i*n+j)]));
            }
        }

        if (!(mx > 1e-13*scale))
        {
            break;
        }

        for (label p=0; p<n-1; ++p)
        {
            for (label q=p+1; q<n; ++q)
            {
                if (Foam::mag(A[static_cast<std::size_t>(p*n+q)]) > 1e-15*scale)
                {
                    rotate(p, q);
                }
            }
        }
    }

    eval.resize(static_cast<std::size_t>(n));
    for (label i=0; i<n; ++i)
    {
        eval[static_cast<std::size_t>(i)] = A[static_cast<std::size_t>(i*n+i)];
        if (!std::isfinite(eval[static_cast<std::size_t>(i)]))
        {
            return false;
        }
    }
    return true;
}

bool invertDense
(
    const std::vector<scalar>& Ain,
    const label n,
    std::vector<scalar>& inv
)
{
    std::vector<scalar> A = Ain;
    inv = identityDense(n);

    scalar scale = 0.0;
    for (label i=0; i<n*n; ++i)
    {
        scale = Foam::max(scale, Foam::mag(A[static_cast<std::size_t>(i)]));
    }
    scale = Foam::max(scale, scalar(1e-300));

    for (label col=0; col<n; ++col)
    {
        label piv = col;
        scalar best = Foam::mag(A[static_cast<std::size_t>(col*n+col)]);
        for (label r=col+1; r<n; ++r)
        {
            const scalar v = Foam::mag(A[static_cast<std::size_t>(r*n+col)]);
            if (v > best)
            {
                best = v;
                piv = r;
            }
        }

        if (!(best > 1e-14*scale))
        {
            return false;
        }

        if (piv != col)
        {
            for (label j=0; j<n; ++j)
            {
                std::swap
                (
                    A[static_cast<std::size_t>(col*n+j)],
                    A[static_cast<std::size_t>(piv*n+j)]
                );
                std::swap
                (
                    inv[static_cast<std::size_t>(col*n+j)],
                    inv[static_cast<std::size_t>(piv*n+j)]
                );
            }
        }

        const scalar d = A[static_cast<std::size_t>(col*n+col)];
        for (label j=0; j<n; ++j)
        {
            A[static_cast<std::size_t>(col*n+j)] /= d;
            inv[static_cast<std::size_t>(col*n+j)] /= d;
        }

        for (label r=0; r<n; ++r)
        {
            if (r == col) continue;
            const scalar f = A[static_cast<std::size_t>(r*n+col)];
            if (f == 0.0) continue;
            for (label j=0; j<n; ++j)
            {
                A[static_cast<std::size_t>(r*n+j)] -=
                    f*A[static_cast<std::size_t>(col*n+j)];
                inv[static_cast<std::size_t>(r*n+j)] -=
                    f*inv[static_cast<std::size_t>(col*n+j)];
            }
        }
    }

    return true;
}

std::vector<scalar> denseMatVec
(
    const std::vector<scalar>& A,
    const std::vector<scalar>& x,
    const label n
)
{
    std::vector<scalar> y(static_cast<std::size_t>(n), 0.0);
    for (label i=0; i<n; ++i)
    {
        scalar s = 0.0;
        for (label j=0; j<n; ++j)
        {
            s += A[static_cast<std::size_t>(i*n+j)]*x[static_cast<std::size_t>(j)];
        }
        y[static_cast<std::size_t>(i)] = s;
    }
    return y;
}

std::vector<scalar> denseMul
(
    const std::vector<scalar>& A,
    const std::vector<scalar>& B,
    const label n
)
{
    std::vector<scalar> C(static_cast<std::size_t>(n*n), 0.0);
    for (label i=0; i<n; ++i)
    {
        for (label k=0; k<n; ++k)
        {
            const scalar aik = A[static_cast<std::size_t>(i*n+k)];
            for (label j=0; j<n; ++j)
            {
                C[static_cast<std::size_t>(i*n+j)] +=
                    aik*B[static_cast<std::size_t>(k*n+j)];
            }
        }
    }
    return C;
}

bool definiteInvSqrt
(
    const std::vector<scalar>& G,
    const label n,
    std::vector<scalar>& Gmhalf,
    scalar& sign,
    scalar& gmin,
    scalar& gmax,
    scalar& gcond
)
{
    std::vector<scalar> eval;
    std::vector<scalar> V;
    if (!jacobiEigenSymmetric(G, n, eval, V))
    {
        return false;
    }

    scalar emin = std::numeric_limits<scalar>::max();
    scalar emax = -std::numeric_limits<scalar>::max();
    for (label i=0; i<n; ++i)
    {
        emin = Foam::min(emin, eval[static_cast<std::size_t>(i)]);
        emax = Foam::max(emax, eval[static_cast<std::size_t>(i)]);
    }

    // Both extremes must share one sign; flip a negative definite G
    sign = emax > 0.0 ? 1.0 : -1.0;
    gmin = Foam::min(sign*emin, sign*emax);
    gmax = Foam::max(sign*emin, sign*emax);

    if (!(gmax > 0.0) || !(gmin > 1e-13*gmax))
    {
        return false;
    }

    gcond = gmax/gmin;
    Gmhalf.assign(static_cast<std::size_t>(n*n), 0.0);

    for (label k=0; k<n; ++k)
    {
        const scalar d =
            1.0/std::sqrt(sign*eval[static_cast<std::size_t>(k)]);
        for (label i=0; i<n; ++i)
        {
            const scalar vik = V[static_cast<std::size_t>(i*n+k)];
            for (label j=0; j<n; ++j)
            {
                Gmhalf[static_cast<std::size_t>(i*n+j)] +=
                    d*vik*V[static_cast<std::size_t>(j*n+k)];
            }
        }
    }

    return true;
}

}
}
