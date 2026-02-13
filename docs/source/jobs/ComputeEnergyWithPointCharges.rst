Computing energies with added point charges
===========================================

This job computes the new energy levels of a system when one or many point charges are added. This computation can be performed analytically, on a cube grid or on a Becke grid.
The number of energy levels computed is equal to the number of states given in the input file that describes the electronic transitions in the unperturbed system.

If the positions of the charges are not given in the input file, they will be consecutively placed on the atoms of the system.



Supported input formats for the analytic file
---------------------------------------------

Below are the supported input file formats:

- ``.fchk``
- ``.gab``
- ``.log``
- ``.molden``
- ``.wfx``.


Supported input formats for the transitions file
------------------------------------------------

Below are the supported input file formats:

- ``.log``
- ``.molden``
- ``.txt`` (see below for the format).


Transitions file format
-----------------------

The electronic transitions can be described in a ``.txt`` file. Below is an example of the expected format for this file, that includes to excited states:

.. code-block:: none
    :linenos:

    Energy 8.5844 eV
    7 A 9 A 0.69993
    7 B 9 B 0.69993

    Energy 9.7919 eV
    5 9 -0.49831
    6 8 0.49831

Each excited state starts with a line indicating its energy. Accepted units are ``eV`` (electronvolt) or ``H`` (Hartree). This line is followed by one or more lines describing the electronic transitions that contribute to this excited state. Each transition line contains three or five space-separated values. In the first excited state above, five values are given: the number of the occupied orbital and its spin (``A`` for alpha, ``B`` for beta), the number of the virtual orbital and its spin, and the coefficient of the transition. In the second excited state, only three values are given: the number of the occupied orbital, the number of the virtual orbital and the coefficient. In this case, it is assumed that the line describes both an alpha and a beta electrons (i.e. the first line for this state could be expanded into two lines: ``5 A 9 A -0.49831`` and ``5 B 9 B -0.49831``).


Available Verbose levels
------------------------

There are four available Verbose levels for this job:

- ``0`` (default): No log file is generated.
- ``1``: A log file is generated containing:

  - the description of the excited states (read from the file describing the electronic transitions in the unperturbed system)
  - the matrix elements < i | H | j > and < i | H - H0 | j > (triangular matrix: only elements with j <= i are written)

- ``2``: In addition with the information given for verbose level 1, the log file includes:

  - the Slater determinants that contribute to each excited state with their coefficients
  - the detail of the computation of the matrix elements < i | H | j > and < i | H - H0 | j >: < i | H | j >, < i | V_ions/nuclei | j > and < i | V_ions/electrons | j >. In the case of multiple charges, the contributions of each individual charge are summed.

- ``3``: In addition with the information given for verbose levels 1 and 2, the log file includes:

  - the information about the orbitals read from the analytic file
  - the contribution of each individual charge in the case of multiple charges (i.e. < i | V_1/nuclei | j >, < i | V_2/nuclei | j >, ..., < i | V_N/nuclei | j > and < i | V_1/electrons | j >, < i | V_2/electrons | j >, ..., < i | V_N/electrons | j >)
  - the eigenvalues (energies) and their associated eigenvectors before their sorting (the final results are sorted by increasing energy)
  - the projection of the perturbed states onto the unperturbed basis (with the contribution of each state). 


Available Grid sizes
--------------------

Three standard grid sizes are available:

- ``coarse``: 3 points per Bohr radius
- ``medium``: 6 points per Bohr radius
- ``fine``: 12 points per Bohr radius.

It is also possible to define a custom grid size using the ``custom`` option. In this case, the ``CustomSizeData`` parameter must be provided in the input file, with the following format: ``Nx,Ny,Nz,Ox,Oy,Oz,T11,T12,T13,T21,T22,T23,T31,T32,T33``, where ``Ni`` are the number of points in the ``i``-th direction, ``Oi`` are the coordinates of the bottom left corner of the cube and ``Tij`` are the coefficients of the translation vector.


Example input file
------------------

Here is an example input file for this job:

.. code-block:: none
    :linenos:
    
    #RunType=Help
    RunType=ComputeEnergyWithPointCharges
    Verbose=1
    # Analytic file name
    AnalyticFiles=filename.fchk
    # File describing the electronic transitions in the unperturbed system
    TransitionsFile=transitions.txt
    # Output file prefix
    OutputPrefix=output
    # Charge(s) to add, in electron charge units
    Charges=-0.1,-0.3
    # Coordinates of the charge(s), in Angstrom
    Positions=2.0,2.0,2.0,-1.0,-1.5,-3.0
    # Cutoff distance for the charge-nuclei interactions, in angstrom
    NuclearCutoff=0.106
    
