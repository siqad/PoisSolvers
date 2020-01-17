cmake_minimum_required(VERSION 3.10)

# for release builds, add flag "-DCMAKE_BUILD_TYPE=Release" (without quotes) when running cmake.

if(DEFINED SIQAD_INSTALL_ROOT)
    set(POIS_INSTALL_DIR "${SIQAD_INSTALL_ROOT}/phys/poissolver")
elseif(CMAKE_BUILD_TYPE MATCHES Release)
    set(POIS_INSTALL_DIR "release")
else()
    set(POIS_INSTALL_DIR "debug")
ENDIF()

set(POIS_FILES
    LICENSE
    poissolver.physeng
    FEM/src/python/ac.py
    FEM/src/python/capacitance.py
    FEM/src/python/charge_density.py
    FEM/src/python/dopant.py
    FEM/src/python/exporter.py
    FEM/src/python/genexp.py
    FEM/src/python/helpers.py
    FEM/src/python/mesher.py
    FEM/src/python/mesh_writer_3D.py
    FEM/src/python/plotter.py
    FEM/src/python/poisson3D.py
    FEM/src/python/poisson_class.py
    FEM/src/python/res_graph.py
    FEM/src/python/resistance.py
    FEM/src/python/subdomains.py
)

install(FILES ${POIS_FILES} DESTINATION ${POIS_INSTALL_DIR})
