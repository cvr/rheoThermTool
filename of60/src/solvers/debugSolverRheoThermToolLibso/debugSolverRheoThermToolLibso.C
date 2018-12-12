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

#include "ppUtilInterface.H"
#include "constitutiveModel.H"

//#include "thermFunModel.H"
//#include "thermFunScheme.H"

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

    
/*    
    Info<< "schemesDict() = " << mesh.schemesDict() << endl;
    Info<< "schemesDict().dictName() = " << mesh.schemesDict().dictName() << endl;
    Info<< "schemesDict().name() = " << mesh.schemesDict().name() << endl;
    Info<< "schemesDict().digest() = " << mesh.schemesDict().digest() << endl;
    Info<< "schemesDict().tokens() = " << mesh.schemesDict().tokens() << endl;
    Info<< "schemesDict().found(\"gradSchemes\") = " << mesh.schemesDict().found("gradSchemes") << endl;
    Info<< "schemesDict().isDict(\"gradSchemes\") = " << mesh.schemesDict().isDict("gradSchemes") << endl;
    Info<< "schemesDict().subDict(\"gradSchemes\") = " << mesh.schemesDict().subDict("gradSchemes") << endl;
    Info<< "schemesDict().subDict(\"gradSchemes\").toc() = " << mesh.schemesDict().subDict("gradSchemes").toc() << endl;
    Info<< "schemesDict().subDict(\"gradSchemes\")[\"default\"] = " << mesh.schemesDict().subDict("gradSchemes")["default"] << endl;
    Info<< "schemesDict().subDict(\"gradSchemes\").lookup(\"default\") = " << mesh.schemesDict().subDict("gradSchemes").lookup("default") << endl;
    Info<< "schemesDict().subDict(\"gradSchemes\")[\"default\"][0] = " << mesh.schemesDict().subDict("gradSchemes")["default"][0] << endl;
    Info<< "schemesDict().subDict(\"gradSchemes\")[\"default\"][1] = " << mesh.schemesDict().subDict("gradSchemes")["default"][1] << endl;
    Info<< "schemesDict().subDict(\"gradSchemes\").parent() = " << mesh.schemesDict().subDict("gradSchemes").parent() << endl;
    Info<< "schemesDict().subDict(\"gradSchemes\").parent().dictName() = " << mesh.schemesDict().subDict("gradSchemes").parent().dictName() << endl;



    //dictionary d0 = mesh.schemesDict();
    dictionary d0(mesh.schemesDict());
    Info<< "d0 = " << d0 << endl;
    Info<< "d0.dictName() = " << d0.dictName() << endl;
    Info<< "d0.name() = " << d0.name() << endl;
    Info<< "d0.digest() = " << d0.digest() << endl;
    Info<< "d0.found(\"gradSchemes\") = " << d0.found("gradSchemes") << endl;
    Info<< "d0.isDict(\"gradSchemes\") = " << d0.isDict("gradSchemes") << endl;
    
    //dictionary d1 = mesh.schemesDict().subDict("gradSchemes");
    dictionary d1(d0.subDict("gradSchemes"));
    Info<< "d1 = " << d1 << endl;
    Info<< "d1.toc() = " << d1.toc() << endl;
    Info<< "d1[\"default\"] = " << d1["default"] << endl;
    Info<< "d1.lookup(\"default\") = " << d1.lookup("default") << endl;
    Info<< "d1[\"default\"][0] = " << d1["default"][0] << endl;
    Info<< "d1[\"default\"][1] = " << d1["default"][1] << endl;
    Info<< "d1.parent() = " << d1.parent() << endl;
    Info<< "d1.parent().dictName() = " << d1.parent().dictName() << endl;



    //Info<< mesh.schemesDict().lookup("gradSchemes") << endl;

//    Info<< "ddtScheme(U) = " << mesh.schemesDict().ddtSchemes("U") << endl;
//    Info<< "thermFunScheme(p) = " << mesh.schemesDict().thermFunScheme("p") << endl;
//    Info<< "thermFunScheme(phi) = " << mesh.schemesDict().thermFunScheme("phi") << endl;

*/

    /*
    Info<< "mu.F() =" << mu.F() << endl;
    mu.calcThermFun(T, T0);
    Info<< "mu.F() =" << mu.F() << endl;
    
    Info<< "ert.F() =" << ert.F() << endl;
    ert.calcThermFun(T, T0);
    Info<< "ert.F() =" << ert.F() << endl;

    Info<< "cp.F() =" << cp.F() << endl;
    cp.calcThermFun(T, T0);
    Info<< "cp.F() =" << cp.F() << endl;
    
    Info<< "k.F() =" << k.F() << endl;
    k.calcThermFun(T, T0);
    Info<< "k.F() =" << k.F() << endl;

    Info<< "mu.F0() = " << mu.F0() << endl;
    //Info<< "mu.c1() = " << mu.c1() << endl;  //error, as c1() is undefined in thermFunScheme.H
    */

    int nIter = 0;
    //while (simple.loop(runTime))
    while (nIter < 5)
    {
        nIter += 1;
        
        simple.loop(runTime);

        Info<< "Time = " << runTime.timeName() << nl << endl;

        for (int i=0; i<1; i++)
        {
            Info<< "Iteration  " << i << nl << endl; 
            
            // --- Solve only the constitutive equation
            constEq.correct();

            //increase etaS, etaP and lambda by a factor
            if (constEq.isGNF()) {
                constEq.etaRef() *= Xadim;
            } else {
                constEq.etaSRef() *= Xadim;
                constEq.etaPRef() *= Yadim;
                constEq.lambdaRef() *= XYadim;
            }
        }
        
        //extraStress = constEq.tauTotal();

        if(runTime.outputTime()) {
            if (constEq.isGNF()) {
                constEq.etaRef().write();
            } else {
                constEq.etaSRef().write();
                constEq.etaPRef().write();
                constEq.lambdaRef().write();
            }
        }

        runTime.write(); // Uncomment if needed
      
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
