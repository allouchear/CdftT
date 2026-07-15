Computing energies with added point charges: ``ComputeEnergyWithPointCharges``
==============================================================================

This job computes the new energy levels of a system when one or many point charges are added. This computation can be performed analytically, on a cartesian grid or on a Becke grid.

The number of energy levels computed is equal to the number of states given in the input file that describes the electronic transitions in the unperturbed system.


Input file parameters
---------------------

Below are the available parameters for this job.

| **Mandatory** parameters are indicated with an **asterisk: "\*"**.
| An asterisk enclosed in parentheses "(\*)" indicates that the parameter can be mandatory or optional, depending on the conditions mentioned in its description below.

Optional parameters that are not specified in the input file will take their default value (this will be announced with a "note" during the run).


Mandatory parameters
^^^^^^^^^^^^^^^^^^^^

``AnalyticFiles*``
""""""""""""""""""

Specifies the name of the input file containing the information about the system. Supported formats are:

- ``.fchk``
- ``.gab``
- ``.log``
- ``.molden``
- ``.wfx``.

| Only one file must be provided when using a ``.fchk`` or a ``.log`` file, as these store the Ground State (GS) energy of the system.
| For other formats (that do not store the GS energy), you must either:

- add a second file (``.fchk`` or ``.log``) in this parameter value, that will only be used to read the GS energy
- define the GS energy using the ``GroundStateEnergy`` parameter (see below)
- use a file format that store both electronic transitions and the GS energy in the ``TransitionsFileName`` parameter (see below).


``RunType*``
""""""""""""

Specifies the CdftT module to run. To run this job, set this parameter to ``ComputeEnergyWithPointCharges``.


Mandatory/optional parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``CustomGridSize(*)``
"""""""""""""""""""""

This parameter is mandatory to define the grid parameters when the ``Size`` parameter is set to ``Custom`` (see below).

The format of this parameter value is: ``Nx, Ny, Nz, Ox, Oy, Oz, T11,T12, T13, T21, T22, T23, T31, T32, T33``, where:

- ``Nx``, ``Ny``, ``Nz`` are the number of grid points in each direction
- ``Ox``, ``Oy``, ``Oz`` define the origin of the grid
- ``Tij`` define the transformation matrix elements.


``GroundStateEnergy(*)``
""""""""""""""""""""""""

Used to specify the Ground State (GS) energy of the system, **in Hartree**.

This parameter is either mandatory or optional depending on the file format used in the ``AnalyticFiles`` parameter. GS energy can also be read from an analytic file (if provided) or from a file describing the electronic transitions in the unperturbed system (if provided). If none of these files are provided, then this parameter must be specified in the input file.


``TransitionsFileName(*)``
""""""""""""""""""""""""""

Used to specify the name of the input file describing the electronic transitions in the unperturbed system. Supported formats are:

- ``.log``
- ``.out``
- ``.txt`` (see below for the format).

This parameter is either mandatory or optional depending on the file format used in the ``AnalyticFiles`` parameter. If the file specified in the ``AnalyticFiles`` parameter does not contain information about the electronic transitions in the unperturbed system, then this parameter must be specified in the input file.


Optional parameters
^^^^^^^^^^^^^^^^^^^

``Becke``
"""""""""

Used to define the parameters of the Becke grid.

The format of this parameter is: ``kmax, lebedev_order, radial_grid_factor``, where:

- ``kmax`` (default ``3``): Controls the smoothness of the transition zone between two atoms. It is recommended to keep the default value, as modifying it has a minor impact on the results compared to the other two parameters.
- ``lebedev_order`` (default ``41``): Controls the number of points on the angular part of the grid (i.e., on the surface of concentric spheres around each atom). A higher value increases the density of points and thus the precision.
- ``radial_grid_factor`` (default ``5``): Adjusts the number of points along the radial direction (i.e., the number of concentric spheres) for each atom. This number is determined automatically based on the atom's identity, so this parameter acts as a multiplier (i.e., a higher value increases the density of points).

The program will only perform computations on a Becke grid if this parameter is specified in the input file.
Therefore, to use default values without having to learn them by heart, you can simply write ``Becke = default`` in the input file.


``Charges``
"""""""""""

Used to specify the charge(s) to add, **in electron charge units**.

The default value is -1.0.


``ChargesPositionsBijections``
""""""""""""""""""""""""""""""

Used to specify whether the relationship between the charges and their positions is a list of bijections (i.e. each charge has one and only one position, and each position is associated with one and only one charge).

Possible values are ``True`` and ``False`` (default).

| In the case of bijective relationships, if positions are specified for the charges, the number of positions must be a multiple of the number of charges. The i-th charge will be put at the (m + i)-th positions, where m is a integer between 0 and (number of positions / number of charges) - 1.
| Example with 2 charges and 6 positions: charge 1 will be put at positions 1, 3 and 5, while charge 2 will be put at positions 2, 4 and 6.

| In the case of bijective relationships, if positions are not specified for the charges (charges are put on atoms, see ``Positions`` below), the number of charges must be a multiple of the number of atoms. The i-th charge will be put on the (m + i)-th atoms, where m is a integer between 0 and (number of charges / number of atoms) - 1.
| Example with 2 charges and 6 atoms: charge 1 will be put on atoms 1, 3 and 5, while charge 2 will be put on atoms 2, 4 and 6.

In the case of non-bijective relationships, all charges will be put at/on all positions/atoms.


``EnergyPointChargeMethods``
""""""""""""""""""""""""""""

Used to specify the method(s) used to compute the new energy levels of the system when one or many point charges are added. Possible values are:

- ``LinearResponse``: the new energy levels are computed using linear response theory, within the Independant Particle Approximation (IPA)
- ``Perturbative``: the new energy levels are computed using perturbation theory
- ``Variational`` (default): the new energy levels are computed using a variational method.


``MaxNumberOfExcitedStates``
""""""""""""""""""""""""""""

Used to specify the maximum number of excited states to take into account when loading the electronic transitions.

The default value is -1 (all available excited states are considered).


``NuclearCutoff``
"""""""""""""""""

Used to specify the cutoff distance for the charge-nuclei interactions, **in angstroms**. This avoids numerical instabilities when a charge is very close to a nucleus.

The default value is 0.1 Å (≈ 0.18897 Bohr radius).


``OutputPrefix``
""""""""""""""""

Used to specify the prefix of the output files generated by this job. An underscore ("_") is automatically added at the end of the prefix if it is not already present.

The default value is an empty string.

Below are some examples of filenames obtained for the file that contains the new energy levels (``energies.cdftt``) with different values of the ``OutputPrefix`` parameter:

- ``OutputPrefix`` is not given in the input file (so it is set to its default value: an empty string): the filename will not be modified (``energies.cdftt``)
- ``OutputPrefix = output``: the filename will become ``output_energies.cdftt``
- ``OutputPrefix = results_``: the filename will become ``results_energies.cdftt``
- ``OutputPrefix = output_results``, the filename will become ``output_results_energies.cdftt``, etc.


``Positions``
"""""""""""""

Used to specify the coordinates of the charge(s), **in angstroms**.

A coordinate is defined by a triplet of numbers ``(x, y, z)`` (with the space dimensions in this order) enclosed in parentheses and separated by commas.

If this parameter is not provided in the input file, the coordinates of the atoms will be used.


``SavePseudoOrbitals``
""""""""""""""""""""""

Used to specify whether the pseudo-orbitals should be saved in a ``.cube`` file.
Two files are created: one for the alpha spin and one for the beta spin.
The resulting filenames will be ``[outputPrefix_]lrf_pseudoOrbitals_alpha.cube`` and ``[outputPrefix_]lrf_pseudoOrbitals_beta.cube``, with the prefix defined by the ``OutputPrefix`` parameter (see above), if any.

Possible values are ``True`` and ``False`` (default).


``ShowProgress``
""""""""""""""""

Used to specify whether the progress of the computation should be shown in the terminal, during an interactive session.

Possible values are ``True`` and ``False`` (default).


``SingleCharge``
""""""""""""""""

Used to specify whether the computation should include only one point charge.

Possible values are ``True`` (default) and ``False``.

In single charge mode, if multiple charges are given in the input file, the program will generate multiple runs and loop on the given charges, using only one charge at a time.


``Size``
""""""""

Used to specify the density of points for a regular cartesian grid.

| Possible values are ``coarse`` (3 points per Bohr radius), ``medium`` (6 points per Bohr radius, default) or ``fine`` (12 points per Bohr radius).
| It is also possible to define a custom grid using the ``custom`` value. In this case, the grid parameters must be defined with the ``CustomGridSize`` parameter (see above).

The program will only perform grid-based computations if this parameter is specified in the input file.


``Verbose``
"""""""""""

Used to specify the level of details in the log file.

Possible values are:

- ``0`` (default): No log file is generated.
- ``1``: A log file is generated containing:

    - the description of the excited states (read from the file describing the electronic transitions in the unperturbed system)
    - the matrix elements :math:`\langle \, i \, \middle| \, \hat{H} \, \middle| \, j \, \rangle` and :math:`\langle \, i \, \middle| \, \hat{H} - \hat{H}_0 \, \middle| \, j \, \rangle` (triangular matrix: only elements with :math:`j \leq i` are written)

- ``2``: In addition with the information given for verbose level 1, the log file includes:

    - the Slater determinants that contribute to each excited state with their coefficients
    - the detail of the computation of the matrix elements :math:`\langle \, i \, \middle| \, \hat{H} \, \middle| \, j \, \rangle` and :math:`\langle \, i \, \middle| \, \hat{H} - \hat{H}_0 \, \middle| \, j \, \rangle`: :math:`\langle \, i \, \middle| \, \hat{H} \, \middle| \, j \, \rangle`, :math:`\langle \, i \, \middle| \, V_{\mathrm{ions/nuclei}} \, \middle| \, j \, \rangle` and :math:`\langle \, i \, \middle| \, V_{\mathrm{ions/electrons}} \, \middle| \, j \, \rangle`. In the case of multiple charges, the contributions of each individual charge are summed.


- ``3``: In addition with the information given for verbose levels 1 and 2, the log file includes:

    - the information about the orbitals read from the analytic file
    - the contribution of each individual charge in the case of multiple charges (i.e. :math:`\langle \, i \, \middle| \, V_{1/\mathrm{nuclei}} \, | \, j \, \rangle`, :math:`\langle \, i \, \middle| \, V_{2/\mathrm{nuclei}} \, \middle| \, j \, \rangle`, ..., :math:`\langle \, i \, \middle| \, V_{N/\mathrm{nuclei}} \, \middle| \, j \, \rangle` and :math:`\langle i \, \middle| \, V_{1/\mathrm{electrons}} \, \middle| \, j \, \rangle`, :math:`\langle i \, \middle| \, V_{2/\mathrm{electrons}} \, \middle| \, j \, \rangle`, ..., :math:`\langle \, i \, \middle| \, V_{N/\mathrm{electrons}} \, \middle| \, j \, \rangle`, where :math:`N` is the number of charges.)
    - the eigenvalues (energies) and their associated eigenvectors before their sorting (the final results are sorted by increasing energy)
    - the projection of the perturbed states onto the unperturbed basis (with the contribution of each state). 



Transitions file format
-----------------------

The electronic transitions can be read from a quantum chemistry software output file (such as ``.log`` from Gaussian) or from a ``.cdftt`` file produced by ``ComputeEnergyWithPointCharges``.
They can also be described in a ``.txt`` file written manually. Below is an example of the expected format for this file, that includes two excited states:

.. code-block:: none
    :linenos:

    # This is a comment line, it will be ignored by the program. Blank lines are also ignored.
    
    # The line giving the GS energy below is optional. It can be used to specify the GS energy of the system if it is not given in an analytic file (see "AnalyticFiles" parameter above).
    Ground State Energy -191.91 H

    Energy 8.5844 eV
    7 A 9 A 0.69993
    7 B 9 B 0.69993

    Energy 9.7919 eV
    5 9 -0.49831
    6 8 0.49831

Each excited state starts with a line indicating its difference in energy with respect to the GS energy. Accepted units are ``eV`` (electronvolt) or ``H`` (Hartree). This line is followed by one or more lines describing the electronic transitions that contribute to this excited state. Each transition line contains three or five space-separated values. In the first excited state above, five values are given: the number of the occupied orbital and its spin (``A`` for alpha, ``B`` for beta), the number of the virtual orbital and its spin, and the coefficient of the transition. In the second excited state, only three values are given: the number of the occupied orbital, the number of the virtual orbital and the coefficient. In this case, it is assumed that the line describes both an alpha electron and a beta electron (i.e. the first line for this state could be expanded into two lines: ``5 A 9 A -0.49831`` and ``5 B 9 B -0.49831``).


Example input file
------------------

Here is an example input file for this job:

.. code-block:: none
    :linenos:

    # This is a comment line, it will be ignored by the program. Blank lines are also ignored.
    # RunType = Help
    RunType = ComputeEnergyWithPointCharges

    # Analytic file name
    AnalyticFiles = filename.fchk

    # File describing the electronic transitions in the unperturbed system
    TransitionsFile = transitions.txt
    MaxNumberOfExcitedStates = 49

    # Output control 
    OutputPrefix = myMolecule
    Verbose = 1
    ShowProgress = True

    # Charge(s) to add, in electron charge units
    Charges = -0.1, -0.3
    SingleCharge = True     # The -0.1 and -0.3 charges will be put in separate runs

    # Coordinates of the charge(s), in Angstrom
    Positions = (2.0, 2.0, 2.0), (-1.0, -1.5, -3.0)
    ChargesPositionsBijections = False   # The -0.1 and -0.3 charges will be put at both positions (2.0, 2.0, 2.0) and (-1.0, -1.5, -3.0) so there will be four runs in total:
                                         # Run 1: charge -0.1 at position (2.0, 2.0, 2.0)
                                         # Run 2: charge -0.1 at position (-1.0, -1.5, -3.0)
                                         # Run 3: charge -0.3 at position (2.0, 2.0, 2.0)
                                         # Run 4: charge -0.3 at position (-1.0, -1.5, -3.0)
    
    # Method to use for the computation of the new energy levels
    EnergyPointChargeMethods = Perturbative, Variational

    # Cutoff distance for the charge-nuclei interactions, in angstrom
    NuclearCutoff = 0.106

    # Compute using a default Becke grid (in addition with the analytical computation that is always performed)
    Becke = default

    # Compute using a fine regular grid (in addition with the Becke grid defined above, and the analytical computation that is always performed)
    Size = Fine
