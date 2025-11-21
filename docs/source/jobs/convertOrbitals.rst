Converting analytical files
===========================

This job converts an analytical file into an other format.


Supported input formats
-----------------------

Below are the supported input formats for the analytical files:

- ``.fchk``
- ``.gab``
- ``.log``
- ``.molden``
- ``.wfx``.


Supported output formats
------------------------

Below are the supported output formats for the analytical files:

- ``.gab``
- ``.molden``
- ``.wfx``.


Example input file
------------------

Here is an example input file for this job:

.. code-block:: none
    :linenos:
    
    #RunType=Help
    RunType=ConvertOrbitals
    AnalyticFiles=input.wfx, output.molden