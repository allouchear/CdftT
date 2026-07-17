Computing Molecular Orbital values:  ``MakeOrbitalsCube``
=========================================================

This job analytically computes Molecular Orbital (MO) values on a cartesian grid and saves it in a ``.cube`` file.


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


``OrbitalType*``
""""""""""""""""

Used to specify which MOs to compute. Possible values are:

- ``All`` (default): all MOs
- ``Occupied``: Occupied MOs
- ``Virtual``: Virtual MOs
- ``Homo``: Highest Occupied MO
- ``Lumo``: Lowest Unoccupied MO
- ``Homo-Lumo``: HOMO and LUMO
- ``Custom``: user-selected MOs.

When using ``Custom``, two parameters must be provided: ``OrbitalsList`` and ``OrbitalsSpins`` (see below).
Otherwise, the desired spins can be set using the ``SpinType`` parameter (see below).


``RunType*``
""""""""""""

Specifies the CdftT module to run. To run this job, set this parameter to ``MakeOrbitalsCube``.


Mandatory/optional parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``CustomGridSize(*)``
"""""""""""""""""""""

This parameter is mandatory to define the grid parameters when the ``Size`` parameter is set to ``Custom`` (see below).

The format of this parameter value is: ``Nx, Ny, Nz, Ox, Oy, Oz, T11,T12, T13, T21, T22, T23, T31, T32, T33``, where:

- ``Nx``, ``Ny``, ``Nz`` are the number of grid points in each direction
- ``Ox``, ``Oy``, ``Oz`` define the origin of the grid
- ``Tij`` define the transformation matrix elements.


``OrbitalsList(*)``
"""""""""""""""""""

This parameter is mandatory when the ``OrbitalType`` parameter is set to ``Custom``. It specifies the list of Molecular Orbitals (MOs) to compute, as a comma-separated list of MO indices (starting from 1).


``OrbitalsSpins(*)``
""""""""""""""""""""

This parameter is mandatory when the ``OrbitalType`` parameter is set to ``Custom``. It specifies the list of spin types for the selected Molecular Orbitals (MOs), as a comma-separated list of spin types (``Alpha`` or ``Beta``).


Optional parameters
^^^^^^^^^^^^^^^^^^^

``ShowProgress``
""""""""""""""""

Used to specify whether the progress of the computation should be shown in the terminal, during an interactive session.

Possible values are ``True`` and ``False`` (default).


``Size``
""""""""

Used to specify the density of points for a regular cartesian grid.

| Possible values are ``coarse`` (3 points per Bohr radius), ``medium`` (6 points per Bohr radius, default) or ``fine`` (12 points per Bohr radius).
| It is also possible to define a custom grid using the ``custom`` value. In this case, the grid parameters must be defined with the ``CustomGridSize`` parameter (see above).


``SpinType``
""""""""""""

Used to specify which spin type to compute. Possible values are ``Alpha``, ``Beta`` or ``Alpha-Beta`` (default).


Example input file
------------------

Here is an example input file for this job:

.. code-block:: none
    :linenos:

    # This is a comment line, it will be ignored by the program. Blank lines are also ignored.
    # RunType = Help
    RunType = MakeOrbitalsCube

    # Analytic file name
    AnalyticFiles = filename.fchk

    # Custom cartesian grid configuration and output file name
    GridFilesNames = output.cube
    Size = Custom
    CustomSizeData = 140, 100, 60 ; 7.558904, 5.669178, 2.834589 ; 0.094486, 0, 0, 0, 0.094486, 0, 0, 0, 0.094486

    # Show a progress bar in the terminal
    ShowProgress = True

    # Orbitals selection: example 1 - Select all spinorbitals
    OrbitalType = All
    SpinType = Alpha-Beta

    # OR

    # Orbitals selection, example 2 - Select specific spinorbitals
    # OrbitalType = Custom
    # OrbitalsNumbers = 1, 1, 2, 2, 3, 3, 10
    # SpinList = Alpha, Beta, Alpha, Beta, Alpha, Beta, Alpha
