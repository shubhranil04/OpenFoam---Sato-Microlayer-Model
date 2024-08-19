{
// Calculate mass source due to microlayer evaporation

//const scalarField& mlFaceAreas = mesh.magSf().boundaryField()[mlPatchID];

const volScalarField& TSat = phaseChange.satProp().TSat();

const volScalarField& L = phaseChange.satProp().L();

const scalar kappal = max(mixture.thermo1().kappa()).value();

volScalarField sharpmdotml(mdotml * 0.0);

const fvMesh& mesh = alpha2.mesh();


forAll(mlCellLabels,i){

	const label cellI = mlCellLabels[i];

	if ((alpha1[cellI] < 0.5) && (del[cellI]>1e-8)){
		
		mfluxml[cellI] = (kappal*(TL.boundaryField()[mlPatchID][i] - TSat[cellI])/del[cellI])/L[cellI];
		
		sharpmdotml[cellI] = mfluxml[cellI]*mlFaceAreas[i]/mesh.V()[cellI];

	}
	else {
		sharpmdotml[cellI] = 0;
		
		mfluxml[cellI] = 0;
	}


}

// Spread mass source similar to Hardt-and-Wandra

    dimensionedScalar DPsi_
    (
        "DPsi_",
        dimensionSet(0,2,0,0,0,0,0),
        3/sqr(gAverage(mesh.nonOrthDeltaCoeffs())
    ));

    dimensionedScalar intPsi0_ = fvc::domainIntegrate(sharpmdotml);

    volScalarField psi_
    (
        IOobject
        (
            "psi_",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimDensity/dimTime, 0),
        "zeroGradient"
    );

    //- Smearing of source term field
    fvScalarMatrix psi_Eqn
    (
        fvm::Sp(scalar(1),psi_) - fvm::laplacian(DPsi_,psi_) == sharpmdotml
    );

    psi_Eqn.solve();

    // Cut cells with cutoff < alpha1 < 1-cutoff
    // and rescale remaining source term field
   
    dimensionedScalar intPsiVapor_
    (
        "intPsiVapor_",
        dimensionSet(1,0,-1,0,0,0,0),
        0.0
    );

    forAll(mesh.C(),celli)
    {
        if (alpha1[celli] < 0.001)
        {
            intPsiVapor_.value() +=
                (1.0-alpha1[celli])*psi_[celli]*mesh.V()[celli];
        }
    }

    //- Calculate Nl and Nv
    dimensionedScalar Nv_ ("Nv_", dimless, 2.0);

    reduce(intPsiVapor_.value(),sumOp<scalar>());

    if (intPsiVapor_.value() > 1e-99)
    {
        Nv_ = intPsi0_/intPsiVapor_;
    }


    //- Set source terms in cells with alpha1 < cutoff or alpha1 > 1-cutoff
    forAll(mesh.C(),celli)
    {
        if (alpha1[celli] < 0.001)
        {
            mdotml[celli] = Nv_.value()*(1.0-alpha1[celli])*psi_[celli];
        }
        else
        {
            mdotml[celli] = 0.0;
        }
    }
}
