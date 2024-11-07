/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    uniGasFoam

Group
    grpDiscreteMethodsSolvers

Description
    Unified particle-based solver for 3D multiscale rarefied gas flows

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvCFD.H"
#include "uniGasCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #define NO_CONTROL
    #include "postProcess.H"
 
     argList::addBoolOption
    (
        "keep-lagrangian", 
        "Keeping lagrangian data in all write directories"
    );
 
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    const bool cleanLagrangian = !args.found("keep-lagrangian");
    
    if (cleanLagrangian)
    {
        uniGas.cleanLagrangian();
    }

    label infoCounter = 0;

    while (runTime.loop())
    {          
        infoCounter++;
        
        if (infoCounter >= uniGas.nTerminalOutputs())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;
        }

        uniGas.evolve();

        if (infoCounter >= uniGas.nTerminalOutputs())
        {
            uniGas.info();   
        }

        runTime.write();

        if (infoCounter >= uniGas.nTerminalOutputs())
        {
            uniGas.loadBalanceCheck();
        }
        
        if (infoCounter >= uniGas.nTerminalOutputs())
        {
            Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
                
            infoCounter = 0;
        }
    
        if (cleanLagrangian && runTime.write())
        {
            uniGas.cleanLagrangian();
        }

    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
