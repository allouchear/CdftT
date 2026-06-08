// Backbone integration workflows across multiple jobs.

void test_density_workflow_backbone() {
    // P1: MakeDensityCube -> ComputePartialCharges -> ComputeDescriptors.
}

void test_orbital_conversion_workflow_backbone() {
    // P1: ConvertOrbitals -> MakeOrbitalsCube.
}

void test_grid_arithmetic_workflow_backbone() {
    // P1: Density cubes -> ComputeGridDifference -> ComputeIntegrals.
}

void test_excited_state_workflow_backbone() {
    // P1: ComputeEnergyWithPointCharges + LambdaDiagnostic.
}
