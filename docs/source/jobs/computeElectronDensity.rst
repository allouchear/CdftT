<<<<<<< HEAD
Compute electronic densities and transition densities : ``ComputeElectronDensity``
=========================================================================================
=======
Compute electronic densities / electronic transition densities : ``ComputeElectronDensity``
===========================================================================================

This job analytically computes the first order reduced density matrix (RDM-1) of a system. It then uses it to compute the electronic density (or electronic transition density) on a cartesian grid and stores the results in a ``.cube`` file.
>>>>>>> 0937009a0c9004f4ab4da59ca61196de68eacffe

This Job computes (analytically) the first order reduced density matrix RDM-1 of a system and use it to compute the density (electronic density or transition density) on a cube grid. The program needs in input the molecular orbitals from a quantum chemistry software such as ``Gaussian`` or ``Orca``, and a file describing the electronic transitions.

Input file parameters
---------------------


Below are the available parameters for this job.

| **Mandatory** parameters are indicated with an **asterisk: "\*"**.
| An asterisk enclosed in parentheses "(\*)" indicates that the parameter can be mandatory or optional, depending on the conditions mentioned in its description below.

Optional parameters that are not specified in the input file will take their default value (this will be announced with a "note" during the run).


Mandatory parameters
^^^^^^^^^^^^^^^^^^^^

``RunType*``
""""""""""""""""""
Specifies wich of the modules available in CdftT you want to run, in this case it is ``ComputeElectronDensity``.


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
- use a file format that store both electronic transitions and the GS energy in the ``TransitionsFileName`` parameter (see below).
<<<<<<< HEAD
=======


``RunType*``
""""""""""""

Specifies the CdftT module to run. To run this job, set this parameter to ``ComputeElectronDensity``.
>>>>>>> 0937009a0c9004f4ab4da59ca61196de68eacffe


``TransitionsFile*``
""""""""""""""""""""""""""

Used to specify the name of the input file describing the electronic transitions. Supported formats are:

- ``.log``
- ``.out``
- ``.cdftt`` (created by ``ComputeEnergyWithPointCharges``)
- ``.txt`` (see below for the format).

This parameter is either mandatory or optional depending on the file format used in the ``AnalyticFiles`` parameter. If the file specified in the ``AnalyticFiles`` parameter does not contain information about the electronic transitions, then this parameter must be specified in the input file.


``CustomGridSize(*)``
"""""""""""""""""""""

This parameter is mandatory to define the grid parameters when the ``Size`` parameter is set to ``Custom`` (see below).

The format of this parameter value is: ``Nx, Ny, Nz, Ox, Oy, Oz, T11,T12, T13, T21, T22, T23, T31, T32, T33``, where:

- ``Nx``, ``Ny``, ``Nz`` are the number of grid points in each direction
- ``Ox``, ``Oy``, ``Oz`` define the origin of the grid
- ``Tij`` define the translation matrix elements.
<<<<<<< HEAD
=======


``TransitionsFile(*)``
""""""""""""""""""""""

Used to specify the name of the input file describing the electronic transitions. Supported formats are:

- ``.log``
- ``.out``
- ``.cdftt`` (created by ``ComputeEnergyWithPointCharges``)
- ``.txt`` (see below for the format).

This parameter is either mandatory or optional depending on the file format used in the ``AnalyticFiles`` parameter. If the file specified in the ``AnalyticFiles`` parameter does not contain information about the electronic transitions, then this parameter must be specified in the input file.
>>>>>>> 0937009a0c9004f4ab4da59ca61196de68eacffe


Optional parameters
^^^^^^^^^^^^^^^^^^^

``Size``
""""""""

Used to specify the density of points for a regular cartesian grid.

| Possible values are ``coarse`` (3 points per Bohr radius), ``medium`` (6 points per Bohr radius, default) or ``fine`` (12 points per Bohr radius).
| It is also possible to define a custom grid using the ``custom`` value. In this case, the grid parameters must be defined with the ``CustomGridSize`` parameter (see above).


``RDMMethod``
"""""""""""""
<<<<<<< HEAD
Choose which method to use to compute RDM-1, ``Gamma`` (default) or ``X``. Must be used is the case of mono-excitations and ``X`` should be use only for TDDFT-TDA calculation.
Computation time depends on the number of molecular orbitals nMO and number of transitions/Slater Determinants nSD. ``X`` should be faster (O(nMO^3) complexity vs O(nMO^2*nSD^2) for ``Gamma``). But Only ``Gamma`` is parallelized and nSD can be reduced by the ``SDCutoff`` parameter (see below). 
=======

Used to specify the method to use to compute the RDM-1.

Possible values are ``Gamma`` (default) or ``X``.
``Gamma`` is only valid in the case of mono-excitations. ``X`` must be used only for TDDFT calculations in the Tamm-Dancoff Approximation.

Computation time depends on the number N of molecular orbitals and number M of transitions/Slater Determinants. ``X`` should be faster (O(N^3) complexity vs O(N^2 * M^2) for ``Gamma``). However, ``Gamma`` only is parallelized and M can be reduced using the ``SDCutoff`` parameter (see below). 

>>>>>>> 0937009a0c9004f4ab4da59ca61196de68eacffe

``ExcitedStatesNumbers``
""""""""""""""""""""""""

<<<<<<< HEAD
List of states for which the density will be computed, all by default, starting at 0 for ground state.
=======
Used to specify the numbers of the states for which the density will be computed. Note that the Ground State corresponds to the state number 0.

By default, all available excited states are considered.

>>>>>>> 0937009a0c9004f4ab4da59ca61196de68eacffe

``ExcludedOrbitals``
""""""""""""""""""""

<<<<<<< HEAD
List of molecular orbitals to not take into account in calculations, 1-based index. ``none`` by default.
=======
Used to specify the numbers of the excluded molecular orbitals (not taken into account in calculations). Note that the first orbital is numbered 1.

By default, all available orbitals are considered.

>>>>>>> 0937009a0c9004f4ab4da59ca61196de68eacffe

``TransitionDensities``
"""""""""""""""""""""""

<<<<<<< HEAD
List of transitions Densities to compute in the following format : (a,b),(c,d)...
By default none are computed, and if choosing ``X`` as ``RDMMethod``, only (0,i) transitions can be computed.
=======
Used to specify the transitions densities to compute. The transition densities to compute are specified as a list of tuples, each tuple containing the numbers of the two states between which the transition density is computed. Note that the Ground State corresponds to the state number 0.

For example, to compute the transition densities between the Ground State and the first excited state, and between the Ground State and the tenth excited state, one would use: ``TransitionDensities = (0, 1) ; (0, 10)``.

By default, transition densities are not computed.

Note that if ``RDMMethod`` is set to ``X``, only the transitions from the Ground State to excited states can be computed: (0, i) where i is the number of the excited state.

>>>>>>> 0937009a0c9004f4ab4da59ca61196de68eacffe

``SDCutoff``
""""""""""""

<<<<<<< HEAD
Cutoff applied on transitions coefficients. Calculations will be carried using only those with a coefficient C such that C> ``SDCutoff``/max(C). Default is ``0``
=======
Used to set a cutoff on the transition coefficients.

Calculations will then be carried using only the transitions whose coefficient C satisfies C > ``SDCutoff``/max(C).

There is no cutoff by default.

>>>>>>> 0937009a0c9004f4ab4da59ca61196de68eacffe

``OutputPrefix``
""""""""""""""""

<<<<<<< HEAD
Prefix added to all files generated by ``ComputeElectronDensity``. None by default.
=======
Used to specify the prefix of the output files generated by this job. An underscore ("_") is automatically added at the end of the prefix if it is not already present.

The default value is an empty string.

Below are some examples of filenames obtained for the file that contains the new energy levels (``energies.cdftt``) with different values of the ``OutputPrefix`` parameter:

- ``OutputPrefix`` is not given in the input file (so it is set to its default value: an empty string): the filename will not be modified (``energies.cdftt``)
- ``OutputPrefix = output``: the filename will become ``output_energies.cdftt``
- ``OutputPrefix = results_``: the filename will become ``results_energies.cdftt``
- ``OutputPrefix = output_results``, the filename will become ``output_results_energies.cdftt``, etc.

>>>>>>> 0937009a0c9004f4ab4da59ca61196de68eacffe

``Precision``
"""""""""""""

<<<<<<< HEAD
number of significant digits used for output data, default is ``10``.
=======
Used to set the number of significant digits in output data.

The default value is ``10``.

>>>>>>> 0937009a0c9004f4ab4da59ca61196de68eacffe

``ShowProgress``
""""""""""""""""

<<<<<<< HEAD
If ``True`` (``False`` by default), adds a progress bar to the RDM and grid computation.
=======
Used to specify whether the progress of the computation should be shown in the terminal, during an interactive session.

Possible values are ``True`` and ``False`` (default).

>>>>>>> 0937009a0c9004f4ab4da59ca61196de68eacffe

``Verbose``
"""""""""""

By default is ``0``, if equals ``1``, generates a log files with additional data, such as for each state it's energy and the Slater determinant used.

``SaveReducedDensityMatrix``
""""""""""""""""""""""""""""

<<<<<<< HEAD
If ``True`` (``False`` by default), saves the RDM-1 matrix in a file.

=======
- ``0`` (default): No log file is generated.
- ``1``: A log file is generated containing:

    - the description of the excited states (read from the file describing the electronic transitions in the unperturbed system)
    - the matrix elements :math:`\langle \, i \, \middle| \, \hat{H} \, \middle| \, j \, \rangle` and :math:`\langle \, i \, \middle| \, \hat{H} - \hat{H}_0 \, \middle| \, j \, \rangle` (triangular matrix: only elements with :math:`j \leq i` are written)


``SaveReducedDensityMatrix``
""""""""""""""""""""""""""""

Used to specify whether the reduced density matrix should be saved in a file.

Possible values are ``True`` and ``False`` (default).
>>>>>>> 0937009a0c9004f4ab4da59ca61196de68eacffe


Transitions file format
-----------------------

<<<<<<< HEAD
Transitions can be obtained from quantum chemistry software (such as ``.log`` from Gaussian) or as a ``.cdftt`` file produced by ``ComputeEnergyWithPointCharges``.
The electronic transitions can also be described in a ``.txt`` file. Below is an example of the expected format for this file, that includes to excited states:
=======
The electronic transitions can be read from a quantum chemistry software output file (such as ``.log`` from Gaussian) or from a ``.cdftt`` file produced by ``ComputeEnergyWithPointCharges``.
They can also be described in a ``.txt`` file written manually. Below is an example of the expected format for this file, that includes two excited states:
>>>>>>> 0937009a0c9004f4ab4da59ca61196de68eacffe

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
<<<<<<< HEAD
    RunType=ComputeElectronDensity

    # input data
    AnalyticFiles=acrolein.fchk
    TransitionsFile=acrolein.log
    
    # Grid configuration
    Size=custom
    CustomSizeData=140,100,60,7.558904,5.669178,2.834589,0.094486,0,0,0,0.094486,0,0,0,0.094486
    
    # computation parameters
    RDMMethod=Gamma
    ExcitedStatesNumbers=0,1,10
    ExcludedOrbitals=1,2,3,4
    TransitionDensities=(0,1),(0,10)
    SDCutoff=0.0001
    
    # output parameters
    OutputPrefix=data/Gamma_SD0001
    Precision=6
    
    # verbosity
    ShowProgress=True
    verbose=0
    SaveReducedDensityMatrix=True
=======
    RunType = ComputeElectronDensity

    # Analytic file name
    AnalyticFiles = filename.fchk

    # File describing the electronic transitions in the unperturbed system
    TransitionsFile = transitions.txt

    # Output control 
    OutputPrefix = myMolecule
    Precision = 6
    SaveReducedDensityMatrix = True
    ShowProgress = True
    Verbose = 1
    
    # Custom cartesian grid configuration
    Size = Custom
    CustomSizeData = 140, 100, 60 ; 7.558904, 5.669178, 2.834589 ; 0.094486, 0, 0, 0, 0.094486, 0, 0, 0, 0.094486
    
    # Computation parameters
    RDMMethod = Gamma
    ExcitedStatesNumbers = 0, 1, 10
    ExcludedOrbitals = 1, 2, 3, 4
    TransitionDensities = (0, 1) ; (0, 10)
    SDCutoff = 0.0001
    
>>>>>>> 0937009a0c9004f4ab4da59ca61196de68eacffe
