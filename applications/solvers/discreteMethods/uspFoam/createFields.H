IOdictionary uspInitialiseDict
(
    IOobject
    (
        "uspInitialiseDict",
        mesh.time().system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

IOstream::defaultPrecision(15);

Info<< "Initialising usp for Time = " << runTime.timeName() << nl << endl;

uspCloud usp(runTime, "usp", mesh, uspInitialiseDict);

Info<< nl << "ClockTime = " << runTime.elapsedClockTime() << " s"
    << nl << endl;
