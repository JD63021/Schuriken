/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Serial unit tests for the dense and distributed linear algebra used by
 * AMGsensei.  Needs no case:  Test-schurikenLinalg
 */

#include "DenseOps.H"
#include "DistributedOps.H"
#include "PODBasis.H"
#include "ReducedSpace.H"
#include "UPstream.H"
#include "IOstreams.H"

#include <cmath>
#include <limits>
#include <vector>

using namespace Foam;
using namespace Foam::schuriken;

namespace
{

label nFailed = 0;
const label comm = UPstream::worldComm;

void check(const char* name, const scalar err, const scalar tol)
{
    const bool ok = std::isfinite(err) && err <= tol;
    if (!ok)
    {
        ++nFailed;
    }
    Info<< "SCHURIKEN_UNIT test=" << name << " err=" << err << " tol=" << tol
        << " status=" << (ok ? "PASS" : "FAIL") << endl;
}

void checkTrue(const char* name, const bool ok)
{
    check(name, ok ? 0 : 1, 0);
}

// Deterministic pseudo-random numbers in [-0.5, 0.5)
scalar pseudoRandom(const label i, const label seed)
{
    const scalar s =
        std::sin(12.9898*scalar(i + 1) + 78.233*scalar(seed))*43758.5453;
    return s - std::floor(s) - 0.5;
}

// Symmetric positive definite B B^T + n I
std::vector<scalar> randomSPD(const label n, const label seed)
{
    const std::size_t N = n;
    std::vector<scalar> B(N*N), A(N*N, 0.0);
    for (std::size_t i=0; i<N*N; ++i)
    {
        B[i] = pseudoRandom(i, seed);
    }
    for (std::size_t i=0; i<N; ++i)
    {
        for (std::size_t j=0; j<N; ++j)
        {
            for (std::size_t k=0; k<N; ++k)
            {
                A[i*N+j] += B[i*N+k]*B[j*N+k];
            }
        }
        A[i*N+i] += n;
    }
    return A;
}

std::vector<scalar> negated(std::vector<scalar> A)
{
    for (scalar& a : A)
    {
        a = -a;
    }
    return A;
}

scalar maxAbsDiff(const std::vector<scalar>& A, const std::vector<scalar>& B)
{
    scalar e = 0;
    for (std::size_t i=0; i<A.size(); ++i)
    {
        e = max(e, mag(A[i] - B[i]));
    }
    return A.size() == B.size() ? e : GREAT;
}

scalarField randomField(const label n, const label seed)
{
    scalarField f(n);
    forAll(f, i)
    {
        f[i] = pseudoRandom(i, seed);
    }
    return f;
}

void testInvertDense()
{
    const label n = 6;
    const std::vector<scalar> I = identityDense(n);
    for (const scalar sign : {1.0, -1.0})
    {
        std::vector<scalar> A = randomSPD(n, 3);
        if (sign < 0)
        {
            A = negated(A);
        }
        std::vector<scalar> Ainv;
        checkTrue("invertDense_ok", invertDense(A, n, Ainv));
        check("invertDense_AAinv", maxAbsDiff(denseMul(A, Ainv, n), I), 1e-12);
    }

    std::vector<scalar> S(n*n, 1.0);
    std::vector<scalar> Sinv;
    checkTrue("invertDense_rejects_singular", !invertDense(S, n, Sinv));
}

void testJacobi(const label n)
{
    const std::size_t N = n;
    const std::vector<scalar> A = randomSPD(n, 5);
    std::vector<scalar> eval, V;
    checkTrue("jacobi_ok", jacobiEigenSymmetric(A, n, eval, V));

    std::vector<scalar> R(N*N, 0.0), VtV(N*N, 0.0);
    for (std::size_t i=0; i<N; ++i)
    {
        for (std::size_t j=0; j<N; ++j)
        {
            for (std::size_t k=0; k<N; ++k)
            {
                R[i*N+j] += V[i*N+k]*eval[k]*V[j*N+k];
                VtV[i*N+j] += V[k*N+i]*V[k*N+j];
            }
        }
    }
    check("jacobi_reconstruct", maxAbsDiff(R, A)/n, 1e-12);
    check("jacobi_orthogonal", maxAbsDiff(VtV, identityDense(n)), 1e-12);
}

void testDefiniteInvSqrt()
{
    const label n = 5;
    for (const scalar expected : {1.0, -1.0})
    {
        std::vector<scalar> G = randomSPD(n, 7);
        if (expected < 0)
        {
            G = negated(G);
        }
        std::vector<scalar> S;
        scalar sign = 0, gmin = 0, gmax = 0, gcond = 0;
        checkTrue
        (
            "definiteInvSqrt_ok",
            definiteInvSqrt(G, n, S, sign, gmin, gmax, gcond)
        );
        check("definiteInvSqrt_sign", mag(sign - expected), 0);
        check("definiteInvSqrt_gmin_positive", gmin > 0 ? 0 : 1, 0);

        // S (sign G) S = I
        std::vector<scalar> sG(G);
        for (scalar& g : sG)
        {
            g *= sign;
        }
        check
        (
            "definiteInvSqrt_SGS",
            maxAbsDiff(denseMul(denseMul(S, sG, n), S, n), identityDense(n)),
            1e-11
        );
    }

    std::vector<scalar> D(identityDense(n));
    D[(n - 1)*(n + 1)] = -1.0;
    std::vector<scalar> S;
    scalar sign = 0, gmin = 0, gmax = 0, gcond = 0;
    checkTrue
    (
        "definiteInvSqrt_rejects_indefinite",
        !definiteInvSqrt(D, n, S, sign, gmin, gmax, gcond)
    );
}

void testDistributedOps(const label n)
{
    std::vector<scalarField> L, R;
    for (label k=0; k<3; ++k)
    {
        L.push_back(randomField(n, 10 + k));
        R.push_back(randomField(n, 20 + k));
    }

    const std::vector<scalar> G = crossGramGlobal(L, R, comm);
    const std::vector<scalar> d = dotBatchGlobal(L, R[1], comm);
    scalar eG = 0, ed = 0, edot = 0;
    for (std::size_t i=0; i<L.size(); ++i)
    {
        for (std::size_t j=0; j<R.size(); ++j)
        {
            eG = max(eG, mag(G[i*R.size()+j] - sum(L[i]*R[j])));
        }
        ed = max(ed, mag(d[i] - sum(L[i]*R[1])));
        edot = max(edot, mag(dotLocal(L[i], R[i]) - sum(L[i]*R[i])));
    }
    check("crossGramGlobal", eG/n, 1e-15);
    check("dotBatchGlobal", ed/n, 1e-15);
    check("dotLocal", edot/n, 1e-15);

    // y += sum_k c_k L_k
    const std::vector<scalar> c{0.5, -2.0, 3.0};
    scalarField y(R[0]);
    addCombination(y, L, c);
    scalarField yRef(R[0] + 0.5*L[0] - 2.0*L[1] + 3.0*L[2]);
    check("addCombination", max(mag(y - yRef)), 1e-14);

    scalarField v(randomField(n, 30));
    checkTrue("normalizeEuclid_ok", normalizeEuclid(v, comm));
    check("normalizeEuclid_unit", mag(normGlobal(v, comm) - 1), 1e-14);
    scalarField z(n, 0.0);
    checkTrue("normalizeEuclid_rejects_zero", !normalizeEuclid(z, comm));
}

void testSnapshotPOD()
{
    // Three snapshots in a two-dimensional subspace
    const label n = 50;
    const scalarField a(randomField(n, 40));
    const scalarField b(randomField(n, 41));
    std::vector<scalarField> snaps;
    snaps.push_back(scalarField(a + 0.5*b));
    snaps.push_back(scalarField(2.0*a - b));
    snaps.push_back(scalarField(-a + 3.0*b));

    std::vector<scalarField> modes;
    scalar energy = 0;
    checkTrue("snapshotPOD_ok", snapshotPOD(snaps, 3, modes, energy, comm));
    checkTrue("snapshotPOD_rank2", modes.size() == 2);
    check("snapshotPOD_energy", mag(energy - 1), 1e-10);

    if (modes.size() == 2)
    {
        check("snapshotPOD_orthogonal", mag(sum(modes[0]*modes[1])), 1e-10);

        // Every snapshot lies in span(modes)
        scalar e = 0;
        for (const scalarField& s : snaps)
        {
            scalarField r(s);
            for (const scalarField& m : modes)
            {
                r -= sum(m*s)*m;
            }
            e = max(e, Foam::sqrt(sum(r*r)/sum(s*s)));
        }
        check("snapshotPOD_span", e, 1e-10);
    }
}

// AW = d*W for a diagonal "matrix" d
ReducedSpace diagonalSpace
(
    const std::vector<scalarField>& W,
    const scalarField& d
)
{
    ReducedSpace space;
    for (const scalarField& w : W)
    {
        space.AW.push_back(scalarField(d*w));
    }
    return space;
}

void testRepairBasis()
{
    const label n = 30;

    // Negative definite diagonal, like an OpenFOAM pressure matrix
    scalarField d(n);
    forAll(d, i)
    {
        d[i] = -(1.0 + i);
    }

    // w2 = w0 + 2 w1 makes W^T A W singular
    std::vector<scalarField> W;
    W.push_back(randomField(n, 50));
    W.push_back(randomField(n, 51));
    W.push_back(scalarField(W[0] + 2.0*W[1]));
    const std::vector<scalarField> W0(W);

    ReducedSpace space = diagonalSpace(W, d);
    checkTrue
    (
        "factorReducedSpace_rejects_singular",
        !factorReducedSpace(W, space, comm)
    );

    const label rank = repairBasis(W, space, comm);
    checkTrue("repairBasis_rank2", rank == 2 && space.rank() == 2);

    if (rank == 2)
    {
        // AW is still A W
        scalar eAW = 0;
        for (std::size_t k=0; k<W.size(); ++k)
        {
            eAW = max(eAW, max(mag(space.AW[k] - d*W[k])));
        }
        check("repairBasis_AW_consistent", eAW, 1e-12);

        // G Ginv = I
        check
        (
            "repairBasis_Ginv",
            maxAbsDiff(denseMul(space.G, space.Ginv, 2), identityDense(2)),
            1e-10
        );

        // The original vectors still lie in span(W): A-projection residual
        scalar eSpan = 0;
        for (const scalarField& w : W0)
        {
            const std::vector<scalar> h = dotBatchGlobal(space.AW, w, comm);
            const std::vector<scalar> c = denseMatVec(space.Ginv, h, 2);
            scalarField r(w);
            for (std::size_t k=0; k<W.size(); ++k)
            {
                r -= c[k]*W[k];
            }
            eSpan = max(eSpan, Foam::sqrt(sum(r*r)/sum(w*w)));
        }
        check("repairBasis_span", eSpan, 1e-8);
    }

    // Growing the space one vector at a time matches a full factorisation
    {
        std::vector<scalarField> Wg(1, randomField(n, 70));
        ReducedSpace grown = diagonalSpace(Wg, d);
        checkTrue
        (
            "appendToReducedSpace_seed",
            factorReducedSpace(Wg, grown, comm)
        );
        for (label k=1; k<4; ++k)
        {
            const scalarField v(randomField(n, 70 + k));
            scalar vAv = 0;
            checkTrue
            (
                "appendToReducedSpace_ok",
                appendToReducedSpace(Wg, grown, v, scalarField(d*v), vAv, comm)
            );
        }

        ReducedSpace full = diagonalSpace(Wg, d);
        checkTrue("appendToReducedSpace_full", factorReducedSpace(Wg, full, comm));
        check("appendToReducedSpace_G", maxAbsDiff(grown.G, full.G), 1e-13);
        check
        (
            "appendToReducedSpace_Ginv",
            maxAbsDiff(grown.Ginv, full.Ginv)/max(mag(full.Ginv[0]), 1.0),
            1e-10
        );

        // A dependent vector is refused and leaves W unchanged
        const scalarField dep(Wg[0] - Wg[1]);
        scalar vAv = 0;
        checkTrue
        (
            "appendToReducedSpace_rejects_dependent",
            !appendToReducedSpace(Wg, grown, dep, scalarField(d*dep), vAv, comm)
         && Wg.size() == 4 && grown.rank() == 4
        );
    }

    // A non-finite basis is dropped entirely
    std::vector<scalarField> Wnan(1, randomField(n, 60));
    Wnan[0][3] = std::numeric_limits<scalar>::quiet_NaN();
    ReducedSpace nanSpace = diagonalSpace(Wnan, d);
    checkTrue
    (
        "factorReducedSpace_rejects_nan",
        !factorReducedSpace(Wnan, nanSpace, comm)
    );
    checkTrue
    (
        "repairBasis_drops_nan",
        repairBasis(Wnan, nanSpace, comm) == 0 && Wnan.empty()
    );
}

} // End anonymous namespace


int main(int argc, char* argv[])
{
    testInvertDense();
    testJacobi(7);
    testJacobi(24);
    testDefiniteInvSqrt();
    testDistributedOps(40);
    testDistributedOps(1300);   // several partial cache blocks
    testSnapshotPOD();
    testRepairBasis();

    Info<< nl << "SCHURIKEN_UNIT_TEST="
        << (nFailed ? "FAIL" : "PASS") << " failed=" << nFailed << endl;

    return nFailed ? 1 : 0;
}
