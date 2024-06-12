IOdictionary uniGasInitialiseDict
(
    IOobject
    (
        "uniGasInitialiseDict",
        mesh.time().system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

IOstream::defaultPrecision(15);

Info<< "Initialising uniGasFoam for Time = " << runTime.timeName() << nl << endl;

uniGasCloud uniGas("uniGas", uniGasInitialiseDict, runTime, mesh);

Info<< nl << "ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;
