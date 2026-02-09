// Reconstruct cell-centered velocity from face-normal velocity Un
// using cached least-squares reconstruction tensor.
//
// Algorithm: U_cell = reconTensor & surfaceSum(SfHat * Un)
//            phi    = Un * magSf
//            Un_bnd = U_bnd & nhat
{
    // Face-normal weighted by unit normal -> vector contribution per face
    surfaceVectorField SfHatUn("SfHatUn", SfHat * Un);

    // Least-squares reconstruction: sum face contributions, apply inverse tensor
    U.primitiveFieldRef() =
        (reconTensor & fvc::surfaceSum(SfHatUn))().primitiveField();

    // Apply boundary conditions from the volVectorField U
    U.correctBoundaryConditions();

    // Update face flux from face-normal velocity
    phi = Un * mesh.magSf();

    // Project U boundary values back onto face normals for consistency
    const surfaceVectorField::Boundary& SfBf = mesh.Sf().boundaryField();
    const surfaceScalarField::Boundary& magSfBf = mesh.magSf().boundaryField();
    surfaceScalarField::Boundary& UnBf = Un.boundaryFieldRef();
    const volVectorField::Boundary& UBf = U.boundaryField();

    forAll(mesh.boundary(), patchi)
    {
        UnBf[patchi] = UBf[patchi] & (SfBf[patchi]/magSfBf[patchi]);
    }
}
