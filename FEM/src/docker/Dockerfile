################################################################################
#                                                                              #
#  This Dockerfile is expected to be ran after it has been copied to the final #
#  PoisSolver structure by SiQAD, it will not run correctly if built in this   #
#  directory.                                                                  #
#                                                                              #
#  Sample build command:                                                       #
#      docker build --tag=poissolver .                                         #
#                                                                              #
#  Sample run command:                                                         #
#      docker run -v /path/to/jobtmp:/poissolver/runtmp poissolver (entrypoint)#
#                                                                              #
################################################################################

# Use a supported distribution with classic DOLFIN packages.
FROM ubuntu:22.04

# Keep classic DOLFIN. DOLFINx is a different API.
# Use signed distribution packages, so no Python package downloader runs here.
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates make gcc g++ libboost-dev python3-dev python3-dolfin python3-tk python3-distutils \
    gmsh swig python3-pil python3-networkx python3-matplotlib python3-numpy python3-shapely \
    python3-scipy python3-h5py python3-meshio \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory to /poissolver
WORKDIR /poissolver

# Copy the poissolver files to the working directory in the container
COPY . /poissolver

# Compile SWIG
RUN cd /poissolver/swig_siqadconn \
    && swig -python -c++ siqadconn.i \
    && python3 setup.py build_ext --inplace \
    && mkdir -p ../siqadtools \
    && mv siqadconn.py _siqadconn.*.so ../siqadtools/ \
    && touch ../siqadtools/__init__.py
RUN python3 /poissolver/poisson3D.py --help
RUN python3 /poissolver/docker_smoke_test.py

ENTRYPOINT ["python3", "/poissolver/poisson3D.py"]
