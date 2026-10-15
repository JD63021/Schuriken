/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "AMGsenseiState.H"
#include "List.H"
#include "objectRegistry.H"
#include "regIOobject.H"
#include "Time.H"

#include <map>

namespace Foam
{
namespace schuriken
{

// Owner of one field's state inside an object registry.  With AUTO_WRITE
// the basis is written to the time directories (always binary); with
// READ_IF_PRESENT a written basis is read back on construction.
class AMGsenseiStateHolder
:
    public regIOobject
{
public:

    TypeName("AMGsenseiState");

    AMGsenseiState data;

    explicit AMGsenseiStateHolder(const IOobject& io)
    :
        regIOobject(io),
        data()
    {
        if (readOpt() == IOobject::READ_IF_PRESENT && headerOk())
        {
            Istream& is = readStream(typeName);
            List<scalarField> basis;
            is >> data.fieldName >> data.nCells >> basis;
            close();

            data.W.assign(basis.begin(), basis.end());
            data.initialized = true;
            data.restoredRank = basis.size();
        }
    }

    virtual bool writeData(Ostream& os) const
    {
        List<scalarField> basis(static_cast<label>(data.W.size()));
        forAll(basis, k)
        {
            basis[k] = data.W[static_cast<std::size_t>(k)];
        }
        os  << data.fieldName << token::SPACE << data.nCells << nl
            << basis << nl;
        return os.good();
    }

    // The basis is large, so it is written in binary whatever writeFormat is
    virtual bool writeObject
    (
        IOstream::streamFormat,
        IOstream::versionNumber ver,
        IOstream::compressionType cmp,
        const bool write
    ) const
    {
        return regIOobject::writeObject(IOstream::BINARY, ver, cmp, write);
    }
};

defineTypeNameAndDebug(AMGsenseiStateHolder, 0);

}
}


Foam::schuriken::AMGsenseiState::AMGsenseiState()
:
    initialized(false),
    nCells(0),
    fieldName(),
    region("none"),
    solveCount(0),
    stockSolves(0),
    stockPCGIts(0),
    podSolves(0),
    podPCGIts(0),
    acceptedUpdates(0),
    recompressCount(0),
    setupSeconds(0),
    iterSeconds(0),
    trainSeconds(0),
    restoredRank(-1),
    W(),
    batch()
{}

void Foam::schuriken::AMGsenseiState::reset
(
    const label n,
    const word& field
)
{
    initialized = true;
    nCells = n;
    fieldName = field;
    solveCount = 0;
    stockSolves = 0;
    stockPCGIts = 0;
    podSolves = 0;
    podPCGIts = 0;
    acceptedUpdates = 0;
    recompressCount = 0;
    setupSeconds = 0;
    iterSeconds = 0;
    trainSeconds = 0;
    restoredRank = -1;
    W.clear();
    batch.clear();
}

Foam::schuriken::AMGsenseiState& Foam::schuriken::state
(
    const lduMesh& mesh,
    const word& fieldName,
    const bool persist
)
{
    const objectRegistry* db = dynamic_cast<const objectRegistry*>(&mesh);

    if (!db)
    {
        static std::map<word, AMGsenseiState> table;
        return table[fieldName];
    }

    const word name("AMGsensei_" + fieldName);

    if (!db->foundObject<AMGsenseiStateHolder>(name))
    {
        // A persisted basis is read from the start time of the run
        const Time& runTime = db->time();
        AMGsenseiStateHolder* holder = new AMGsenseiStateHolder
        (
            IOobject
            (
                name,
                Time::timeName(runTime.startTime().value()),
                *db,
                persist ? IOobject::READ_IF_PRESENT : IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        holder->store();
        holder->data.region = db->name();
    }

    AMGsenseiStateHolder& holder =
        db->lookupObjectRef<AMGsenseiStateHolder>(name);

    // Once any solver dictionary of the field asks for it, keep writing
    if (persist)
    {
        holder.writeOpt() = IOobject::AUTO_WRITE;
    }

    return holder.data;
}
