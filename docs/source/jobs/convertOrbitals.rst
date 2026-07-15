Converting analytical files: ``ConvertOrbitals``
================================================

This job converts an analytical file into an other format.


Input file parameters
---------------------

Below are the available parameters for this job.

| **Mandatory** parameters are indicated with an **asterisk: "\*"**.
| An asterisk enclosed in parentheses "(\*)" indicates that the parameter can be mandatory or optional, depending on the conditions mentioned in its description below.

Optional parameters that are not specified in the input file will take their default value (this will be announced with a "note" during the run).


Mandatory parameters
^^^^^^^^^^^^^^^^^^^^

``AnalyticFiles*``
""""""""""""""""""

Specifies the name of the input and output files. The first file in the list is the input file to convert, and the second one is the output file.

Supported input formats are:

- ``.fchk``
- ``.gab``
- ``.log``
- ``.molden``
- ``.wfx``.

Supported output formats are:

- ``.gab``
- ``.molden``
- ``.wfx``.


``RunType*``
""""""""""""

Specifies the CdftT module to run. To run this job, set this parameter to ``ConvertOrbitals``.


Example input file
------------------

Here is an example input file for this job:

.. code-block:: none
    :linenos:

    # This is a comment line, it will be ignored by the program. Blank lines are also ignored.
    # RunType = Help
    RunType = ConvertOrbitals

    # Analytic file names: first the file to convert, then the desired output file name
    AnalyticFiles = input.wfx, output.molden
