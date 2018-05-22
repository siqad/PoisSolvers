#!/bin/bash

swig -v -python -c++ phys_connector.i

if [ "$FOR_OS" == "win64" ]; then
    export PATH="$PATH:/usr/lib/mxe/usr/bin"
    MAKE_COMMAND=x86_64-w64-mingw32.static-g++

    $MAKE_COMMAND -O2 -fPIC -Wall -Wextra -std=c++11 -c phys_connector.cc
    $MAKE_COMMAND -O2 -fPIC -Wall -Wextra -std=c++11 -c phys_connector_wrap.cxx -I/home/samuelngsh/Python36-64/include -L/home/samuelngsh/Python36-64/libs -lpython36
    $MAKE_COMMAND -shared -o _phys_connector.cpython-36-x86_64-w64-mingw32.so phys_connector.o phys_connector_wrap.o -static-libstdc++ -I/home/samuelngsh/Python36-64/include -L/home/samuelngsh/Python36-64/libs -lpython36

    cp _phys_connector.cpython-36-x86_64-w64-mingw32.so ..
else
    python3 setup.py build_ext --inplace
    cp _phys_connector.cpython-36m-x86_64-linux-gnu.so ..
fi

cp phys_connector.py ..
