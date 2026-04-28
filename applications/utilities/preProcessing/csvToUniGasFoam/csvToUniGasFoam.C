/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2024 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    csvToUniGasFoam

Description
    Reads a CSV file containing gas state variables and maps them to
    uniGasFoam initial fields using nearest-neighbour interpolation.

    Column names and the gas species name are specified in the
    system/csvToUniGasFoam dictionary.

    Writes the following fields to the 0 directory:
      - numberDensity_<species>  [m^-3]
      - transT_<species>         [K]
      - rotT_<species>           [K]
      - vibT_<species>           [K]
      - U_<species>              [m/s]

Usage
    \b csvToUniGasFoam \<csv-file\>

    Example dictionary (system/csvToUniGasFoam):
    \verbatim
    gasSpecies      Ar;

    csvColumns
    {
        x           x;
        y           y;
        z           z;
        numberDensity  rhoN;
        transT      T_trans;
        rotT        T_rot;
        vibT        T_vib;
        elecT       T_elec;
        Ux          Ux;
        Uy          Uy;
        Uz          Uz;
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOdictionary.H"
#include "indexedOctree.H"
#include "treeDataPoint.H"

#include <fstream>
#include <sstream>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Split a CSV line into tokens, stripping surrounding whitespace
static List<std::string> splitCSVLine(const std::string& line)
{
    List<std::string> tokens;
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, ','))
    {
        const size_t s = tok.find_first_not_of(" \t\r\n");
        const size_t e = tok.find_last_not_of(" \t\r\n");
        tokens.append(s == std::string::npos ? "" : tok.substr(s, e - s + 1));
    }
    return tokens;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Reads a CSV file of gas state variables and maps them to "
        "uniGasFoam initial fields using nearest-neighbour interpolation."
    );

    argList::noParallel();
    argList::addArgument("csv-file", "The input .csv file");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOstream::minPrecision(10);

    // ------------------------------------------------------------------ //
    //  Read control dictionary                                            //
    // ------------------------------------------------------------------ //

    IOdictionary csvDict
    (
        IOobject
        (
            "csvToUniGasFoamDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word gasSpecies(csvDict.get<word>("gasSpecies"));

    const dictionary& colDict = csvDict.subDict("csvColumns");

    const word xCol(colDict.getOrDefault<word>("x", "x"));
    const word yCol(colDict.getOrDefault<word>("y", "y"));
    const word zCol(colDict.getOrDefault<word>("z", "z"));

    Info<< "Gas species : " << gasSpecies << nl
        << "Coordinates : x=" << xCol << "  y=" << yCol << "  z=" << zCol
        << nl << endl;

    // ------------------------------------------------------------------ //
    //  Read and parse CSV                                                 //
    // ------------------------------------------------------------------ //

    const fileName csvPath = args.get<fileName>(1);
    if (!isFile(csvPath, false))
    {
        FatalErrorInFunction
            << "Cannot read file " << csvPath
            << exit(FatalError);
    }

    std::ifstream csvFile(csvPath);
    if (!csvFile.good())
    {
        FatalErrorInFunction
            << "Cannot open file " << csvPath
            << exit(FatalError);
    }

    // Parse header line and build column-name-to-index map
    std::string headerLine;
    std::getline(csvFile, headerLine);
    const List<std::string> headers = splitCSVLine(headerLine);

    HashTable<label, word> colId;
    forAll(headers, i)
    {
        colId.insert(word(headers[i]), i);
    }

    auto requireCol = [&](const word& name) -> label
    {
        const auto it = colId.cfind(name);
        if (!it.good())
        {
            FatalErrorInFunction
                << "Column '" << name << "' not found in CSV header." << nl
                << "Available columns: " << headers << nl
                << exit(FatalError);
        }
        return *it;
    };

    // Returns -1 if the key is absent from csvColumns; the field stays zero
    auto optCol = [&](const word& dictKey) -> label
    {
        if (!colDict.found(dictKey)) return -1;
        return requireCol(colDict.get<word>(dictKey));
    };

    const label idX    = requireCol(xCol);
    const label idY    = requireCol(yCol);
    const label idZ    = requireCol(zCol);
    const label idRhoN = optCol("numberDensity");
    const label idTransT  = optCol("transT");
    const label idRotT = optCol("rotT");
    const label idVibT  = optCol("vibT");
    const label idElecT = optCol("elecT");
    const label idUx    = optCol("Ux");
    const label idUy   = optCol("Uy");
    const label idUz   = optCol("Uz");

    Info<< "Column mapping:" << nl
        << "  numberDensity : " << (idRhoN  >= 0 ? word(colDict.get<word>("numberDensity")) : word("(zero)")) << nl
        << "  transT        : " << (idTransT >= 0 ? word(colDict.get<word>("transT"))        : word("(zero)")) << nl
        << "  rotT          : " << (idRotT  >= 0 ? word(colDict.get<word>("rotT"))          : word("(zero)")) << nl
        << "  vibT          : " << (idVibT  >= 0 ? word(colDict.get<word>("vibT"))  : word("(zero)")) << nl
        << "  elecT         : " << (idElecT >= 0 ? word(colDict.get<word>("elecT")) : word("(zero)")) << nl
        << "  Ux            : " << (idUx    >= 0 ? word(colDict.get<word>("Ux"))    : word("(zero)")) << nl
        << "  Uy            : " << (idUy    >= 0 ? word(colDict.get<word>("Uy"))            : word("(zero)")) << nl
        << "  Uz            : " << (idUz    >= 0 ? word(colDict.get<word>("Uz"))            : word("(zero)")) << nl
        << endl;

    // Read data rows
    DynamicList<point>  csvPts;
    DynamicList<scalar> csvRhoN;
    DynamicList<scalar> csvTransT;
    DynamicList<scalar> csvRotT;
    DynamicList<scalar> csvVibT;
    DynamicList<scalar> csvElecT;
    DynamicList<vector> csvU;

    std::string dataLine;
    while (std::getline(csvFile, dataLine))
    {
        if (dataLine.empty() || dataLine.find_first_not_of(" \t\r\n,") == std::string::npos)
            continue;

        const List<std::string> vals = splitCSVLine(dataLine);

        auto val = [&](label id) -> scalar
        {
            return id >= 0 ? std::stod(vals[id]) : 0.0;
        };

        csvPts.append
        (
            point
            (
                std::stod(vals[idX]),
                std::stod(vals[idY]),
                std::stod(vals[idZ])
            )
        );
        csvRhoN.append(val(idRhoN));
        csvTransT .append(val(idTransT));
        csvRotT.append(val(idRotT));
        csvVibT .append(val(idVibT));
        csvElecT.append(val(idElecT));
        csvU.append(vector(val(idUx), val(idUy), val(idUz)));
    }

    const label nPts = csvPts.size();
    Info<< "Read " << nPts << " data points from " << csvPath << nl << endl;

    if (nPts == 0)
    {
        FatalErrorInFunction
            << "No data points read from " << csvPath
            << exit(FatalError);
    }

    // ------------------------------------------------------------------ //
    //  Build nearest-neighbour octree over CSV points                    //
    // ------------------------------------------------------------------ //

    const pointField ptField(csvPts);

    treeBoundBox bb(ptField);
    bb.inflate(0.01);

    indexedOctree<treeDataPoint> octree
    (
        treeDataPoint(ptField),
        bb,
        8,    // maxLevel
        10,   // leafSize
        3.0   // duplicateTol
    );

    // ------------------------------------------------------------------ //
    //  Declare output fields                                              //
    // ------------------------------------------------------------------ //

    runTime.setTime(0, 0);

    volScalarField rhoN
    (
        IOobject
        (
            "numberDensity_" + gasSpecies,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(pow(dimLength, -3), Zero)
    );

    volScalarField transT
    (
        IOobject
        (
            "transT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimTemperature, Zero)
    );

    volScalarField rotT
    (
        IOobject
        (
            "rotT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimTemperature, Zero)
    );

    volScalarField vibT
    (
        IOobject
        (
            "vibT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimTemperature, Zero)
    );

    volScalarField elecT
    (
        IOobject
        (
            "elecT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimTemperature, Zero)
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(dimVelocity, vector::zero)
    );

    // ------------------------------------------------------------------ //
    //  Map CSV data to internal cell centres                             //
    // ------------------------------------------------------------------ //

    Info<< "Mapping to " << mesh.nCells() << " internal cells ..." << nl;

    const vectorField& cellCentres = mesh.cellCentres();

    forAll(cellCentres, cellI)
    {
        const pointIndexHit hit = octree.findNearest(cellCentres[cellI], GREAT);
        const label id = hit.index();

        rhoN [cellI] = csvRhoN[id];
        transT[cellI] = csvTransT [id];
        rotT [cellI] = csvRotT[id];
        vibT [cellI] = csvVibT [id];
        elecT[cellI] = csvElecT[id];
        U    [cellI] = csvU    [id];
    }

    // ------------------------------------------------------------------ //
    //  Map CSV data to boundary face centres                             //
    // ------------------------------------------------------------------ //

    Info<< "Mapping to boundary faces ..." << nl;

    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();

    forAll(boundaryMesh, patchi)
    {
        const polyPatch& patch = boundaryMesh[patchi];

        if (!polyPatch::constraintType(patch.type()))
        {

            const vectorField fC = patch.faceCentres();

            fvPatchScalarField& rhoNPatch   = rhoN.boundaryFieldRef()[patchi];
            fvPatchScalarField& transTPatch = transT.boundaryFieldRef()[patchi];
            fvPatchScalarField& rotTPatch   = rotT.boundaryFieldRef()[patchi];
            fvPatchScalarField& vibTPatch   = vibT.boundaryFieldRef()[patchi];
            fvPatchScalarField& elecTPatch  = elecT.boundaryFieldRef()[patchi];
            fvPatchVectorField& UPatch      = U.boundaryFieldRef()[patchi];

            forAll(fC, faceI)
            {
                const pointIndexHit hit = octree.findNearest(fC[faceI], GREAT);
                const label id = hit.index();

                rhoNPatch [faceI] = csvRhoN[id];
                transTPatch[faceI] = csvTransT [id];
                rotTPatch  [faceI] = csvRotT[id];
                vibTPatch  [faceI] = csvVibT [id];
                elecTPatch [faceI] = csvElecT[id];
                UPatch    [faceI] = csvU    [id];
            }

        }

    }

    // ------------------------------------------------------------------ //
    //  Write fields                                                       //
    // ------------------------------------------------------------------ //

    Info<< "\nWriting fields to " << runTime.timeName() << " ..." << nl << endl;

    rhoN.write();
    transT.write();
    rotT.write();
    vibT.write();
    elecT.write();
    U.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
