################################################################################
#                                                                              #
#  This Dockerfile is expected to be ran after it has been copied to the final #
#  PoisSolver structure by SiQAD, it will not run correctly if built in this   #
#  directory.                                                                  #
#                                                                              #
#  Sample build command:                                                       #
#      docker build --tag=poissolver                                           #
#                                                                              #
#  Sample run command:                                                         #
#      docker run -v /path/to/jobtmp:/poissolver/runtmp poissolver             #
#                                                                              #
################################################################################

# Base off the Ubuntu 18.04 LTS version
#FROM ubuntu:18.04
FROM quay.io/fenicsproject/stable:latest

# Set the working directory to /poissolver
WORKDIR /poissolver

# Copy the poissolver files to the working directory in the container
COPY . /poissolver

# Install dependencies
# TODO get pip packages from a requirements.txt instead
RUN apt-get update && apt-get install -y \
    python3-pip \
    python3-dolfin \
    gmsh \
    swig
RUN pip3 install \
    pillow \
    networkx \
    matplotlib \
    numpy \
    shapely

ENTRYPOINT ["python3", "/poissolver/poisson3D.py"]
