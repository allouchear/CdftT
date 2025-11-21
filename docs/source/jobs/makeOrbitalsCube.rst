Computing molecular orbital values
==================================

This job computes a grid of molecular orbital (MO) values and save it in a ``.cube`` file.


Available orbital selections
----------------------------

Below are the available options to select which MO values to compute:

- ``all`` (default): all MOs
- ``Occ``: Occupied MOs
- ``Virtual``: Virtual MOs
- ``Homo``: Highest Occupied MO
- ``Lumo``: Lowest Unoccupied MO
- ``Homo,Lumo``: HOMO and LUMO
- ``Custom``: user-selected MOs. In this case, the ``OrbitalsList`` parameter must be provided in the input file, with a comma-separated list of MO indices (starting from 1).


Available spin types
--------------------

Below are the available options to select which spin type to compute:

- ``Alpha``
- ``Beta``
- ``Alpha,Beta``.


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
    
    #Runtype=Help
    RunType=MakeOrbitalsCube
    ????? (To Be Done)
