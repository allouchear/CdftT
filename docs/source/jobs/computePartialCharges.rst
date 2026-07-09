Computing partial charges: ``ComputePartialCharges``
====================================================

This job performs grid based computations of partial charges of a molecule.


Input file parameters
---------------------

Below are the available parameters for this job.

| **Mandatory** parameters are indicated with an **asterisk: "\*"**.
| An asterisk enclosed in parentheses "(\*)" indicates that the parameter can be mandatory or optional, depending on the conditions mentioned in its description below.


Mandatory parameters
^^^^^^^^^^^^^^^^^^^^

``GridFileNames*``
""""""""""""""""""

Specifies the name of the ``.cube`` file to use for the computation. Only one file should be provided.


``RunType*``
""""""""""""

Specifies the CdftT module to run. To run this job, set this parameter to ``ComputePartialCharges``.


Optional parameter
^^^^^^^^^^^^^^^^^^

``PartitionMethod``
"""""""""""""""""""

Specifies the method to be used for computing the atomic volumes. The available methods are:

- ``on-grid`` (default): Follows Tang's algorithm to find Bader volumes (based on AIM theory).
- ``near-grid``: A more precise version of ``on-grid``.
- ``near-grid-refinement``: An even more precise version, but requires more computation time.
- ``VDD``: Uses a topological method, assigning points to volumes by distance to the closest atom.
- ``Becke``: Uses a regular density grid to interpolate Becke's atomic variable grids.



Example input file
------------------

Here is an example input file for this job:

.. code-block:: none
    :linenos:
    
    # This is a comment line, it will be ignored by the program. Blank lines are also ignored.
    # RunType = Help
    RunType = ComputePartialCharges
    
    # A single grid file is required
    GridFileNames = grid.cube 
    
    # Use default partition method (on-grid): the line below is commented out, but it can be uncommented to use an other method.
    # PartitionMethod = On-Grid


References
----------

\W. Tang, E. Sanville, G. Henkelman. A grid-based bader analysis algorithm without lattice bias. *Journal of Physics: Condensed Matter* 2009, 21(8):084204. DOI: `10.1088/0953-8984/21/8/084204 <https://doi.org/10.1088/0953-8984/21/8/084204>`_.
