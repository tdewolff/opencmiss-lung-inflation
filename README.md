# OpenCMISS lung inflation

## Building the example

Instructions on how to configure and build with CMake:

  git clone https://github.com/tdewolff/opencmiss-lung-inflation.git
  mkdir build
  cmake -DOpenCMISSLibs_DIR=/path/to/opencmisslib/install ../opencmiss-lung-inflation
  make  # cmake --build . will also work here and is much more platform agnostic.

## Running the example

Explain how the example is run:

  python src/python/lung_inflation.py

## License

Apache 2.0
