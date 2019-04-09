Info<< "Updating thermal functions\n" << endl;

//Info<< "Update thermFun for cp (isobaric specific heat)" << endl;
cp.calcThermFun(T, T0);

//Info<< "Update thermFun for kt (thermal conductivity)" << endl;
kt.calcThermFun(T, T0);

//- Comment: dirty way to do it but if I create
//           thermFunModel k("k", T) inside the
//           if statement, it gets within that
//           scope only and I cannot get it out!
//           This works as a kind of hack.

if ( (&constEq.etaRef() != nullptr) && (constEq.eqTypeName() == "Newtonian") )
{
    //Info<< "Update thermFun for eta (GNF Newtonian)" << endl;
    (*eta).calcThermFun(T, T0);
    constEq.etaRef() = (*eta).F();
}

if (&constEq.kRef() != nullptr)
{
    //Info<< "Update thermFun for k (GNF)" << endl;
    (*k).calcThermFun(T, T0);
    constEq.kRef() = (*k).F();
}

if (&constEq.etaMinRef() != nullptr)
{
    //Info<< "Update thermFun for etaMin (GNF)" << endl;
    (*etaMin).calcThermFun(T, T0);
    constEq.etaMinRef() = (*etaMin).F();
}

if (&constEq.etaMaxRef() != nullptr)
{
    //Info<< "Update thermFun for etaMax (GNF)" << endl;
    (*etaMax).calcThermFun(T, T0);
    constEq.etaMaxRef() = (*etaMax).F();
}

if (&constEq.etaSRef() != nullptr)
{
    //Info<< "Update thermFun for etaS (VE)" << endl;
    (*etaS).calcThermFun(T, T0);
    constEq.etaSRef() = (*etaS).F();
}

if (&constEq.etaPRef() != nullptr)
{
    //Info<< "Update thermFun for etaP (VE)" << endl;
    (*etaP).calcThermFun(T, T0);
    constEq.etaPRef() = (*etaP).F();
}

if (&constEq.lambdaRef() != nullptr) {
    //Info<< "Update thermFun for lambda (VE)" << endl;
    (*lambda).calcThermFun(T, T0);
    constEq.lambdaRef() = (*lambda).F();
}

