Creating an Electron Localization Function grid
===============================================

This job creates a grid, computes the Electron Localisation Function (ELF) using either Savin or Becke method, and saves it in a ``.cube`` file.


Available methods:
------------------

Below are the available methods to compute the ELF:

- ``Becke``
- ``Savin`` (default).


Available Grid sizes
--------------------

Three standard grid sizes are available:

- ``coarse``: 3 points per Bohr radius
- ``medium``: 6 points per Bohr radius
- ``fine``: 12 points per Bohr radius.

It is also possible to define a custom grid size using the ``custom`` option. In this case, the ``CustomSizeData`` parameter must be provided in the input file, with the following format: ``Nx,Ny,Nz,Ox,Oy,Oz,T11,T12,T13,T21,T22,T23,T31,T32,T33``, where ``Ni`` are the number of points in the ``i``-th direction, ``Oi`` are the coordinates of the bottom left corner of the cube and ``Tij`` are the coefficients of the translation vector.


Supported input formats
-----------------------

Below are the supported input file formats:

- ``.fchk``
- ``.gab``
- ``.log``
- ``.molden``
- ``.wfx``.


Example input file
------------------

Here is an example input file for this job:

.. code-block:: none
    :linenos:
    
    #RunType=Help
    RunType=MakeELFCube
    #GridFileName
    AnalyticFile=filename.wfx
    Size=Medium
    ELFmethod=Becke
    Grid=save.cube
