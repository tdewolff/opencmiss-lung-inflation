==================
uniaxial_extension
==================

Uniaxial extension of a cube 

Building the example
====================

Instructions on how to configure and build with CMake::

  git clone https://github.com/OpenCMISS-Examples/uniaxial_extension.git
  mkdir build
  cmake -DOpenCMISSLibs_DIR=/path/to/opencmisslib/install ../uniaxial_extension
  make  # cmake --build . will also work here and is much more platform agnostic.

Running the example
===================

Explain how the example is run::

  cd build
  ./src/fortran/uniaxial_extension.F90

or maybe it is a Python only example::

  source /path/to/opencmisslibs/install/virtaul_environments/oclibs_venv_pyXY_release/bin/activate
  python src/python/uniaxial_extension.py

where the XY in the path are the Python major and minor versions respectively.

Prerequisites
=============

None

License
=======

Apache 2.0 License
