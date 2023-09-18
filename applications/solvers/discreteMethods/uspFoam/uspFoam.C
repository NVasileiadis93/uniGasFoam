/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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
    uspFoam

Group
    grpDiscreteMethodsSolvers

Description
    Direct simulation Monte Carlo (USP) solver for 3D, transient, multi-
    species flows

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvCFD.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"

    argList::addBoolOption
    (
        "no-lagrangian",  // noLagrangian
        "Skip reconstructing lagrangian positions and fields"
    );
    argList::addOptionCompat("no-lagrangian", {"noLagrangian", 2106});

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    
    label infoCounter = 0;

    const bool writeLagrangian = !args.found("no-lagrangian");

    while (runTime.loop())
    {          
        infoCounter++;
        
        if (infoCounter >= usp.nTerminalOutputs())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;
        }

        usp.evolve();

        if (infoCounter >= usp.nTerminalOutputs())
        {
            usp.info();   
        }

        if (writeLagrangian)
        {
            runTime.write();
        }

        if (infoCounter >= usp.nTerminalOutputs())
        {
            usp.loadBalanceCheck();
        }
        
        if (infoCounter >= usp.nTerminalOutputs())
        {
            Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
                
            infoCounter = 0;
        }
    
    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
