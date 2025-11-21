Computing the difference between two grids
==========================================

This job computes the differences between the values of the first and the second provided grids.
It then assigns these to a third grid.


Example input file
------------------

Here is an example input file for this job:

.. code-block:: none
    :linenos:
    
    #Runtype=Help
    RunType=ComputeDifference
    #GridFileName
    Grids=in1.cube, in2.cube, out.cube