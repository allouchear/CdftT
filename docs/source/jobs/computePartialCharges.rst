Computing partial charges
=========================

This job performs grid based computations of partial charges of a molecule.


Available atomic volumes computation methods
--------------------------------------------

The following methods are available to compute atomic volumes. The first three are based on Bader's Atoms In Molecule (AIM) theory:

- ``on-grid``: follows Tang's algorithm to find Bader volumes.
- ``near-grid``: more precise version of ``on-grid``.
- ``near-grid-refinement``: even more precise version, but requires more time.
- ``VDD``: uses a topological method. Assigns points to volumes by distance to the closest atom.
- ``Becke``: uses a regular density grid to interpolate Becke's atomic variable grids.


Example input file
------------------

Here is an example input file for this job:

.. code-block:: none
    :linenos:
    
    #RunType
    #RunType=Help
    RunType=ComputePartialCharges
    #GridFileName
    Grids=h2o_80_0.gcube 
    PartitionMethod=on-grid


References
----------

W. Tang, E. Sanville, G. Henkelman. A grid-based bader analysis algorithm without lattice bias. *Journal o Physics: Condensed Matter* 2009, 21(8):084204. DOI: `10.1088/0953-8984/21/8/084204 <https://doi.org/10.1088/0953-8984/21/8/084204>`_.