Info<< "Creating thermal functions\n" << endl;

Info<< "Setting thermFun for cp (isobaric specific heat)" << endl;
thermFunModel cp("cp", T);
cp.calcThermFun(T, T0);

Info<< "Setting thermFun for kt (thermal conductivity)" << endl;
thermFunModel kt("kt", T);
kt.calcThermFun(T, T0);

//- Comment: dirty way to do it but if I create
//           thermFunModel k("k", T) inside the
//           if statement, it gets within that
//           scope only and I cannot get it out!
//           This works as a kind of hack.

thermFunModel* k = nullptr;
if (&constEq.kRef() != nullptr)
{
    Info<< "Setting thermFun for k (GNF)" << endl;
    k = new thermFunModel("k", T);
    (*k).calcThermFun(T, T0);
}

thermFunModel* etaMin = nullptr;
if (&constEq.etaMinRef() != nullptr)
{
    Info<< "Setting thermFun for etaMin (GNF)" << endl;
    etaMin = new thermFunModel("etaMin", T);
    (*etaMin).calcThermFun(T, T0);
}

thermFunModel* etaMax = nullptr;
if (&constEq.etaMaxRef() != nullptr)
{
    Info<< "Setting thermFun for etaMax (GNF)" << endl;
    etaMax = new thermFunModel("etaMax", T);
    (*etaMax).calcThermFun(T, T0);
}

thermFunModel* etaS = nullptr;
if (&constEq.etaSRef() != nullptr)
{
    Info<< "Setting thermFun for etaS (VE)" << endl;
    etaS = new thermFunModel("etaS", T);
    (*etaS).calcThermFun(T, T0);
}

thermFunModel* etaP = nullptr;
if (&constEq.etaPRef() != nullptr)
{
    Info<< "Setting thermFun for etaP (VE)" << endl;
    etaP = new thermFunModel("etaP", T);
    (*etaP).calcThermFun(T, T0);
}

thermFunModel* lambda = nullptr;
if (&constEq.lambdaRef() != nullptr) {
    Info<< "Setting thermFun for lambda (VE)" << endl;
    lambda = new thermFunModel("lambda", T);
    (*lambda).calcThermFun(T, T0);
}



