CdftT documentation
===================

CdftT (CDFT Tools) is set of a library and a program to compute CDFT (Conceptual Density Functional Theory) descriptors. It is written in C++. It is still under developpement, however several tools are already implemented.

New tools will be implemented soon.


Quickstart
----------

To get started with CdftT, simply run ``cdftt jobInputFile.txt``, where ``jobInputFile.txt`` is an input file describing the job you want to perform. See the :ref:`Jobs` section for more details on available jobs and their input files.



.. _Jobs:

.. toctree::
    :maxdepth: 1
    :caption: Jobs

    jobs/computeDescriptors
    jobs/computeGridDifference
    jobs/computeIntegrals
    jobs/computePartialCharges
    jobs/convertOrbitals
    jobs/makeDensityCube
    jobs/makeELFCube
    jobs/makeOrbitalsCube
    