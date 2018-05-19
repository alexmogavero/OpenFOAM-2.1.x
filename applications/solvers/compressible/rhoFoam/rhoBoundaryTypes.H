const volScalarField::GeometricBoundaryField& pbf = p.boundaryField();
wordList rhoBoundaryTypes = pbf.types();

forAll(rhoBoundaryTypes, patchi)
{
	//rhoBoundaryTypes[patchi] = calculatedFvPatchScalarField::typeName;
    /*if (rhoBoundaryTypes[patchi] == "waveTransmissive")
    {
        rhoBoundaryTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
    }
    else if (pbf[patchi].fixesValue())
    {
        rhoBoundaryTypes[patchi] = fixedRhoFvPatchScalarField::typeName;
    }*/
	rhoBoundaryTypes[patchi] = calculatedFvPatchScalarField::typeName;
}
