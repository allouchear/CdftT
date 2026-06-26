Computing a chemical descriptor: ``ComputeDescriptors``
=======================================================

This job allows the computation of chemical descriptors from analytic or ``cube`` files using on-grid, near-grid, near-grid-refinement and Becke.
Frontier Molecular Orbitals (FMO) and finite difference (FD) are methods also provided for the computation. FMO requires one analytic file (``.log``, ``.wfx``, ``.molden``, ...). FD requires three analytic files.
The other methods require ``cube`` files of nucleophilic, electrophilic and radical attacks for the molecule.

Energies must also be given by the user: if two are given, they are assumed to be the ionization energy and the electron affinity. If three are given they are assumed to be the total energies of each file.


Input file parameters
---------------------

Below are the available parameters for this job.

| **Mandatory** parameters are indicated with an **asterisk: "\*"**.
| An asterisk enclosed in parentheses "(\*)" indicates that the parameter can be mandatory or optional, depending on the conditions mentioned in its description below.

Optional parameters that are not specified in the input file will take their default value (this will be announced with a "note" during the run).


Mandatory parameter
^^^^^^^^^^^^^^^^^^^

``PartitionMethod*``
""""""""""""""""""""

Specifies the method used to partition the volume. Possible values are:

- ``on-grid`` (default)
- ``near-grid``
- ``near-grid-refinement``
- ``Becke``
- ``FD`` (Finite Difference)
- ``FMO`` (Frontier Molecular Orbitals).


Mandatory/optional parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AnalyticFiles(*)``
""""""""""""""""""""

Specifies the name of the input file(s) containing the information about the system. This parameter is **mandatory** when ``PartitionMethod`` is set to ``FD`` or ``FMO``.

Supported formats are:

- ``.fchk``
- ``.gab``
- ``.log``
- ``.molden``
- ``.wfx``.

| For the ``FD`` method, three files must be provided: one for the nucleophilic attack, one for the electrophilic attack, and one for the radical attack on the molecule.
| For the ``FMO`` method, only one file must be provided, corresponding to the unperturbed system.


``Energies(*)``
"""""""""""""""

Specifies the energies required for the calculation, **in Hartree**. This parameter is **mandatory** when ``PartitionMethod`` is one of ``on-grid``, ``near-grid``, ``near-grid-refinement``, or ``Becke``.

Two or three energy values must be provided:

- If two values are given, they are interpreted as the Ionization Potential (I) and the Electron Affinity (A).
- If three values are given, they are interpreted as the total energies of the three systems provided in the ``Grids`` parameter.


``GridFilesNames(*)``
"""""""""""""""""""""

Specifies the name of the ``.cube`` files to use. This parameter is **mandatory** when ``PartitionMethod`` is one of ``on-grid``, ``near-grid``, ``near-grid-refinement``, or ``Becke``.

Three files must be provided, corresponding to the nucleophilic, electrophilic, and radical attacks on the molecule.


Example input file
------------------

Here is an example input file for this job:

.. code-block:: none
    :linenos:

    # This is a comment line, it will be ignored by the program. Blank lines are also ignored.
    # RunType = Help
    RunType = ComputeDescriptors
    
    # Partition method to use for the computation of the descriptors
    PartitionMethod = On-Grid
    
    # Grid files (mandatory for on-grid method)
    GridFilesNames = grid1.cube, grid2.cube, grid3.cube
    
    # Energies (mandatory for on-grid method)
    # Can be Ionization Potential and Electron Affinity (2 values) or total energies of the 3 systems (3 values)
    Energies = 0.25, 0.1
    # Energies= -190.0, -190.25, -189.9