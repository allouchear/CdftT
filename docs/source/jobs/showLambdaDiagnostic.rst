Show the Lambda Diagnostic
==========================

This job prints the result of the Lambda diagnostic test, as described by Peach et al., that judges the reliability of TDDFT excited states calculations. It also allows to validate the grid size configuration by computing overlap integrals between the orbitals involved in the excited states.


Transitions file format
-----------------------

This job requires a file containing the list of electronic transitions to analyse. Below is an example of the expected format for this file, that includes to excited states:

.. code-block:: none
    :linenos:

    Energy 8.5844 eV
    7 A 9 A 0.69993
    7 B 9 B 0.69993

    Energy 9.7919 eV
    5 A 9 A -0.49831
    5 B 9 B -0.49831
    6 A 8 A 0.49831
    6 B 8 B 0.49831

Each excited state starts with a line indicating its energy. Accepted units are ``eV`` (electronvolt) or ``H`` (Hartree). This line is followed by one or more lines describing the electronic transitions that contribute to this excited state. Each transition line contains five space-separated values: the number of the occupied orbital and its spin (``A`` for alpha, ``B`` for beta), the number of the virtual orbital and its spin, and the coefficient of the transition.


Available Grid sizes
--------------------

Three standard grid sizes are available:

- ``coarse``: 3 points per Bohr radius
- ``medium``: 6 points per Bohr radius
- ``fine``: 12 points per Bohr radius.

It is also possible to define a custom grid size using the ``custom`` option. In this case, the ``CustomSizeData`` parameter must be provided in the input file, with the following format: ``Nx,Ny,Nz,Ox,Oy,Oz,T11,T12,T13,T21,T22,T23,T31,T32,T33``, where ``Ni`` are the number of points in the ``i``-th direction, ``Oi`` are the coordinates of the bottom left corner of the cube and ``Tij`` are the coefficients of the translation vector.


Supported input formats
-----------------------

Below are the supported input file formats:

- ``.fchk``
- ``.gab``
- ``.log``
- ``.molden``
- ``.wfx``.


Example input file
------------------

Here is an example input file for this job:

.. code-block:: none
    :linenos:

    #Runtype=Help
    RunType=LambdaDiagnostic
    # Analytic file name
    AnalyticFile=filename.wfx
    # Grid configuration
    Size=Custom
    CustomSizeData=128,128,128,5,5,5,7.87401574803148e-02,0,0,0,7.87401574803148e-02,0,0,0,7.87401574803148e-02
    TransitionsFile=transitions.txt