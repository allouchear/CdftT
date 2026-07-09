Creating a density grid: ``MakeDensityCube``
============================================

This job creates a density grid and saves it in a ``.cube`` file.


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


``GridFilesNames*``
"""""""""""""""""""

Specifies the name of the output file.

Note that there is no check on the file extension, but it is recommended to use the ``.cube`` extension for this file.


``RunType*``
""""""""""""

Specifies the CdftT module to run. To run this job, set this parameter to ``MakeDensityCube``.


Mandatory/optional parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``CustomGridSize(*)``
"""""""""""""""""""""

This parameter is mandatory to define the grid parameters when the ``Size`` parameter is set to ``Custom`` (see below).

The format of this parameter value is: ``Nx, Ny, Nz, Ox, Oy, Oz, T11,T12, T13, T21, T22, T23, T31, T32, T33``, where:

- ``Nx``, ``Ny``, ``Nz`` are the number of grid points in each direction
- ``Ox``, ``Oy``, ``Oz`` define the origin of the grid
- ``Tij`` define the transformation matrix elements.


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


Example input file
------------------

Here is an example input file for this job:

.. code-block:: none
    :linenos:

    # This is a comment line, it will be ignored by the program. Blank lines are also ignored.
    # RunType = Help
    RunType = MakeDensityCube

    # Analytic file name
    AnalyticFiles = filename.fchk

    # Custom cartesian grid configuration and output file name
    GridFilesNames = output.cube
    Size = Custom
    CustomSizeData = 80, 80, 80 ; 5, 5, 5 ; 0.15, 0, 0, 0, 0.15, 0, 0, 0, 0.15

    # Show a progress bar in the terminal
    ShowProgress = True
