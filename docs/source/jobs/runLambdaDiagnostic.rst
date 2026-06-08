Run the Lambda Diagnostic: ``RunLambdaDiagnostic``
==================================================

This job prints the result of the Lambda diagnostic test, as described by Peach et al., that judges the reliability of TDDFT excited states calculations. It also allows to validate the grid size configuration by computing overlap integrals between the orbitals involved in the excited states.


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

Only one file must be provided.


``TransitionsFileName*``
""""""""""""""""""""""""

Used to specify the name of the input file describing the electronic transitions in the unperturbed system. Supported formats are:

- ``.log``
- ``.out``
- ``.txt`` (see below for the format).

This parameter must be specified in the input file.


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

``Size``
""""""""

Used to specify the density of points for the regular cartesian grid.

| Possible values are ``coarse`` (3 points per Bohr radius), ``medium`` (6 points per Bohr radius, default) or ``fine`` (12 points per Bohr radius).
| It is also possible to define a custom grid using the ``custom`` value. In this case, the grid parameters must be defined with the ``CustomGridSize`` parameter (see above).

The default value is ``medium``.


Transitions file format
-----------------------

The electronic transitions can be described in a ``.txt`` file. Below is an example of the expected format for this file, that includes to excited states:

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
    RunType = RunLambdaDiagnostic

    # Analytic file name
    AnalyticFiles = filename.wfx

    # Grid configuration
    Size = Custom
    CustomGridSize = 128, 128, 128 ; 5, 5, 5 ; 7.87401574803148e-02, 0, 0 ; 0, 7.87401574803148e-02, 0 ; 0, 0, 7.87401574803148e-02

    # File describing the electronic transitions in the unperturbed system
    TransitionsFileName = transitions.txt