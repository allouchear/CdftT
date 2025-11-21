Computing local integrals
=========================

This job computes local integrals of grids on volumes defined by a method of choice. A grid is required to define the volumes.
The additional grids provided by the user should contain the quantities to be integrated.


Available partition methods
---------------------------

The following partition methods are available to define the volumes:

- ``on-grid``: defines volumes using on-grid AIM. This method requires an electronic density grid.
- ``near-grid``: defines volumes using near-grid AIM. This method requires electronic density grid.
- ``near-grid-refinement``: defines volumes using near-grid-refinement AIM. This method requires electronic density grid.
- ``VDD``: defines volumes by distance to atoms. This method can use any type of density.
- ``BBS``: builds Basins By SIGN. This method requires a grid of density difference (see :ref:`jobs/computeGridDifference`). An additional ``Cutoff`` parameter is also required in the input file for this method, to set a threshold for insignificant values.
- ``B2S``: builds two basins by SIGN. Same as ``BBS`` but only constructs two volumes.


Example input file
------------------

Here is an example input file for this job:

.. code-block:: none
    :linenos:
    
    #RunType=Help
    RunType=ComputeIntegrals
    #GridFileName
    Grids=gridDefiningVolumes.cube, grid1ToBeIntegrated.cube, grid2ToBeIntegrated.cube
    PartitionMethod=BBS
    Cutoff=1e-10