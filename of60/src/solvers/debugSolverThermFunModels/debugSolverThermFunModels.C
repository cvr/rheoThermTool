/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    debugFoam

Description
    Testing application to be used with library lconstitutiveEquations. This 
    application returns the principal components of the extra-stress tensor,
    given the velocity gradient tensor, which is defined by the user.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IFstream.H"
#include "OFstream.H"
#include "simpleControl.H"

#include "constitutiveModel.H"

#include "thermFunModel.H"
#include "thermFunScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    Info<< "mesh.schemesDict() = " << mesh.schemesDict() << endl;
    #include "createControl.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    
    Info<< "cp.F() =" << cp.F() << endl;
    cp.calcThermFun(T, T0);
    Info<< "cp.F() =" << cp.F() << endl;
    
    Info<< "kt.F() =" << kt.F() << endl;
    kt.calcThermFun(T, T0);
    Info<< "kt.F() =" << kt.F() << endl;

    Info<< "eta.F() =" << eta.F() << endl;
    eta.calcThermFun(T, T0);
    Info<< "eta.F() =" << eta.F() << endl;
    
    Info<< "lambda.F() =" << lambda.F() << endl;
    lambda.calcThermFun(T, T0);
    Info<< "lambda.F() =" << lambda.F() << endl;

    Info<< "eta.F0() = " << eta.F0() << endl;
    //Info<< "eta.c1() = " << eta.c1() << endl;  //error, as c1() is undefined in thermFunScheme.H


    int nIter = 0;
    //while (simple.loop(runTime))
    while (nIter < 2)
    {
        nIter += 1;
        
        simple.loop(runTime);

        Info<< "Time = " << runTime.timeName() << nl << endl;

        for (int i=0; i<1; i++)
        {
            Info<< "Iteration  " << i << nl << endl; 
        }
        
        //runTime.write(); // Uncomment if needed
      
        Info<< runTime.timeOutputValue() << tab 
            << endl;  
         
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    runTime.write();
    runTime.writeNow();

    return 0;
}


// ************************************************************************* //
