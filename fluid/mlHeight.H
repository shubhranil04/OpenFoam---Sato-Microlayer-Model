{

// Update microlayer height

// Get time step
const scalar dt = mesh.time().deltaTValue();

// Store previous time step heights and face areas
const volScalarField delOld = del;

const volScalarField mlFaceAreasOld = mlFaceAreas;

forAll(mlCellLabels,i){
	
	const label cellI = mlCellLabels[i];	
	
	mlFaceAreas[i] = alpha2[cellI]*mesh.magSf().boundaryField()[mlPatchID][i];
	
	// Pure vapor cells
	if (alpha1[cellI] < 0.001){
		del[cellI] = delOld[cellI] - dt*mfluxml[cellI]/max(rho1).value();		
	}
	
	// Contact Line cells
	else if (alpha1[cellI] < 0.999){
	
		if (mlFaceAreas[i] - mlFaceAreasOld[i] > 0)
		{
			scalar delStar = delOld[cellI]*mlFaceAreasOld[i]/mlFaceAreas[i] + del0[cellI]*(mlFaceAreas[i]-mlFaceAreasOld[i])/mlFaceAreas[i];

			del[cellI] = delStar - dt*mfluxml[cellI]/max(rho1).value();
		}
		
		else
		{
			del[cellI] = delOld[cellI] - dt*mfluxml[cellI]/max(rho1).value();
		}
	}
	
	// Liquid Cells
	else 
	{
		del[cellI] = 0;
	}
	
}

}


