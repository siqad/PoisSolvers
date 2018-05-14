#!/bin/bash
swig -v -python -c++ phys_connector.i
python3 setup.py build_ext --inplace
cp _phys_connector.cpython-36m-x86_64-linux-gnu.so ../
cp phys_connector.py ../
