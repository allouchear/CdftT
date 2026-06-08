Computing the difference between two grids: ``ComputeGridDifference``
=====================================================================

This job computes the differences between the values of the first and the second provided grids.
It then assigns these to a third grid.


Input file parameters
---------------------

Below are the available parameters for this job.

| **Mandatory** parameters are indicated with an **asterisk: "\*"**.

``GridFilesNames*``
"""""""""""""""""""

Specifies the names of the ``.cube`` files to use. Three filenames must be provided:

- The first file is the grid from which the second will be subtracted (the minuend).
- The second file is the grid to subtract (the subtrahend).
- The third file is the output grid where the result of the difference will be saved.


Example input file
------------------

Here is an example input file for this job:

.. code-block:: none
    :linenos:
    
    # This is a comment line, it will be ignored by the program. Blank lines are also ignored.
    # RunType = Help
    RunType = ComputeGridDifference
    
    # Grid filenames
    GridFilesNames = minuend.cube, subtrahend.cube, difference_output.cube