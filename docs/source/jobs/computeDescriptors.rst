Computing a chemical descriptor
===============================

This job allows the computation of chemical descriptors from analytic or ``cube`` files using on-grid, near-grid, near-grid-refinement and Becke.
Frontier Molecular Orbitals (FMO) and finite difference (FD) are methods also provided for the computation. FMO requires one analytic file (``.log``, ``.wfx``, ``.molden``, ...). FD requires three analytic files.
The other methods require ``cube`` files of nucleophilic, electrophilic and radical attacks for the molecule.

Energies must also be given by the user: if two are given, they are assumed to be the ionisation potential and the electronic affinity. If three are given they are assumed to be the total energies of each file.


Example input file
------------------

Here is an example input file for this job:

.. code-block:: none
    :linenos:

    #RunType=Help
    RunType=ComputeDescriptors
    #GridFileName
    Grids=grid1.cube, grid2.cube, grid3.cube
    PartitionMethod=on-grid
    Energies=I, A or E1,E2,E3