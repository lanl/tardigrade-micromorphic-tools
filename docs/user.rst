.. _user_manual:

###########
User Manual
###########

***********
Quick Start
***********

WORK IN PROGRESS!!!!
THE TEXT BELOW NEEDS TO BE UPDATED TO REFLECT  MATERIAL DEFINITIONS AND ASSOCIATED DATABASES
CONSISTENT WITH THE TARDEGRADE FRAMEWORK AND ASSOCIATED PROGRAMS (i.e. MOOSE)

This project is built and deployed to the `W-13 Python Environments`_ with continuous integration (CI) and continuous
deployment (CD). Most users will not need to build and install this project from source. Outside of the `W-13 Python
Environments`_, users may need to build and install directly from source. In that case, users are directed to the
:ref:`build` instructions.

With the `W-13 Python Environments`_, this project is installed in the Conda environment ``lib64`` and ``include``
directories, e.g. ``/path/to/my/conda/environment/{lib64,include}``

.. code:: bash

   $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path/to/conda/environment/lib64

Where the appropriate path can be found with

.. code:: bash

   $ find path/to/conda/environment -name "libmicromorphic_tools.so"

For instance, with the W-13 "release" environment on ``sstelmo``

.. code:: bash

   $ find /projects/python/release -name "libmicromorphic_tools.so"
   /projects/python/release/lib64/libmicromorphic_tools.so
   $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/projects/python/release/lib64

As a convenience, the following code may be used to determine the correct, active Conda environment at execution.
The following bash code is provided as an example for end users and not supported by this project. End users who wish to
learn more about bash scripting are directed to the online Bash documentation.

.. code:: bash

   # Get current conda environment information
   conda_env_path=$(conda info | grep "active env location" | cut -f 2 -d :)
   # Export the conda environment library path
   $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${conda_env_path}/lib64

***************************
Use after build from source
***************************

.. code:: bash


It is strongly recommended that anyone building from source make use of the CMake ``--install`` options in a local Conda
environment. It is also possible to install to more traditional system paths, but this may require significantly more
background reading in relevant system administration.

Unless the template repository and all upstream c++ libraries are built and installed to a common system path it is
recommended that the subroutines are left in the project build directory. However, it is possible to copy the shared
library files to any other directory provided the upstream projects ``{error,vector,constitutive}_tools``
are present in the build directory, e.g.
``micromorphic_tools/build/_deps/{error,vector,constitutive}_tools-build/``.

.. code:: bash

   $ pwd
   $ cp /path/to/micromorphic_tools/build/src/cpp/{libmicromorphic_tools.so} .

******************************
Input File Material Definition
******************************


.. warning::

   Constitutive modeler health warning! The integration tests use a ``STATEV`` and ``PROPS`` length of one as the
   "incorrect" lengths to check the thrown exceptions. If your real constitutive model actually using a length of one
   for either vector, the integration test expectation must be updated.

micromorphic_tools requires 2 material constants and 2 state variables. The c++ micromorphic_tools interface, material constants, and state
variables are described in the :ref:`sphinx_api`. The fixed expectations for the mooose interface are defined in the
"Variables" section of the :ref:`sphinx_api` for :ref:`micromorphic_tools_source`. A complete discussion about the constants and their
meaning is not included here. Instead users are directed to calibrated material parameters found in micromorphic_tools entries in the
`Granta/MIMS`_ `Material Database`_ :cite:`MIMS`. Material parameter calibration sets should be availble for download
with the correct mooose input file formatting from MIMS.

The micromorphic_tools project contains moose integration tests for the micromorphic_tools moose interface. These tests perform actual simulations using the same dummy parameters used for unit and integration testing of the micromorphic_tools c++ code. The micromorphic_tools
moose input files used for integration testing can be found in the micromorphic_tools source code repository with the following bash
command

.. code:: bash

    This section needs to be updated

The material definition from an integration test input file is included below for reference

.. warning::

   The material constants used in this example material definition are *NOT* calibrated for any real material data.
   For calibrated material parameters, see the micromorphic_tools entry for materials found in the `Granta/MIMS`_ `Material
   Database`_ :cite:`MIMS`.

.. literalinclude:: ../src/abaqus/single_element_c3d8.inp
   :linenos:
   :lines: 42-50
