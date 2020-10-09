message(STATUS "==== vtk_add_module CvtNeuronReader ====")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../")

vtk_module_add_module(
    CvtNeuronReader
    SOURCES vtkCvtNeuronReader.cxx
    HEADERS vtkCvtNeuronReader.h
    )

vtk_module_link(CvtNeuronReader PUBLIC neuron)

paraview_add_server_manager_xmls(XMLS PVCvtNeuronReader.xml)
