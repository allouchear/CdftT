Computing local integrals: ``ComputeIntegrals``
===============================================

This job computes local integrals of grids on volumes defined by a method of choice. A grid is required to define the volumes.
The additional grids provided by the user should contain the quantities to be integrated.


Input file parameters
---------------------

Below are the available parameters for this job.

| **Mandatory** parameters are indicated with an **asterisk: "\*"**.
| An asterisk enclosed in parentheses "(\*)" indicates that the parameter can be mandatory or optional, depending on the conditions mentioned in its description below.


Mandatory parameters
^^^^^^^^^^^^^^^^^^^^

``GridFilesNames*``
"""""""""""""""""""

Specifies the names of the ``.cube`` files to use. At least one file name must be provided.

- The first file in the list is used to define the integration volumes based on the chosen ``PartitionMethod``.
- The subsequent files (if any) are the grids whose quantities will be integrated over these volumes. If more than one file is provided, the integrals will be computed for each of them.


``PartitionMethod*``
""""""""""""""""""""

Specifies the method used to partition the volume. Possible values are:

- ``on-grid`` (default): Defines volumes using on-grid AIM. Requires an electronic density grid.
- ``near-grid``: Defines volumes using near-grid AIM. Requires an electronic density grid.
- ``near-grid-refinement``: Defines volumes using near-grid-refinement AIM. Requires an electronic density grid.
- ``VDD``: Defines volumes by distance to atoms. Can use any type of density.
- ``BBS``: Builds Basins By SIGN. Requires a grid of density difference (see :ref:`computeGridDifference`).
- ``B2S``: Builds two basins by SIGN. Same as ``BBS`` but only constructs two volumes.


Mandatory/optional parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``Cutoff(*)``
"""""""""""""

This parameter is **mandatory** when the ``PartitionMethod`` is set to ``BBS``. It defines a threshold to discard insignificant values during the basin construction. It is ignored for all other methods.


Example input file
------------------

Here is an example input file for this job:

.. code-block:: none
    :linenos:
    
    # This is a comment line, it will be ignored by the program. Blank lines are also ignored.
    # RunType = Help
    RunType = ComputeIntegrals
    
    # The first grid defines the volumes. The others are integrated over these volumes.
    GridFileNames = grid_defining_volumes.cube, grid_to_be_integrated_1.cube, grid_to_be_integrated_2.cube
    
    # Partition method to define the volumes
    PartitionMethod = BBS
    
    # Cutoff is mandatory for the BBS method.
    Cutoff = 1e-10