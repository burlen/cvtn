message(STATUS "==== paraview_add_plugin PVCvtNeuronReader ====")

paraview_add_plugin(PVCvtNeuronReader
  VERSION "1.0"
  MODULES CvtNeuronReader
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/CvtnReader/vtk.module"
  REQUIRED_ON_SERVER
  REQUIRED_ON_CLIENT
  )
