# Tools for visualization of BMTK output using ParaView

[cvtn](#cvtn) **A command line aplication to convert BMTK/neuron simulation data**

[PVCvtNeuronReader](#plugin) **A ParaView plugin for loading data directly into ParaView**

## First
The dependencies need to be installed.

### HDF5
Any recent version of HDF5 installed by a package manager should suffice.

### Intel TBB
A recent version installed from a package manager should suffice.

### VTK
If building the plugin skip this step, as the converter can get VTK from ParaView.
If only building the command line converter, then VTK is required.

A typical VTK build might be acheived as follows:
```
wget https://www.vtk.org/files/release/8.2/VTK-8.2.0.tar.gz
tar xzfv VTK-9.0.0.tar.gz
mkdir VTK-9.0.0-build
cd VTK-9.0.0-build/
cmake \
    -DCMAKE_INSTALL_PREFIX=../VTK-8.2.0-install \
    -DVTK_MODULE_USE_EXTERNAL_VTK_hdf5=ON \
    ../VTK-9.0.0
make -j
make -j install
```
### ParaView
If you are *not* building the plugin skip this step. ParaView is required for the plugin and can also satisfy the VTK dependency of the command line converter.

A typical ParaView build might be acheived as follows:
```bash
git clone https://gitlab.kitware.com/paraview/paraview.git
cd paraview
./Utilities/SetupForDevelopment.sh
git checkout v5.8.1
git submodule update --recursive
cd ..
mkdir paraview-5.8.1-build
cd paraview-5.8.1-build
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=`pwd`/../paraview-5.8.1-install  \
    -DVTK_MODULE_USE_EXTERNAL_VTK_hdf5=ON \
    -DPARAVIEW_USE_PYTHON=ON \
    -DPARAVIEW_USE_MPI=OFF \
    ../paraview
make -j
make -j install
```

## cvtn
**A command line aplication to convert BMTK/neuron simulation data**

The tool can convert neruons with geometry stored in a `seg_coords` folder and
interpolate time series data stored in a `im.h5` file onto the geometry.
After conversion the VTK file format is used to store data back on disk for
injestion into ParaView for visualization.

In the
`seg_coords` folder it is expected to have one file for each neuron amd the
files named by the neuron number. For instance if there are 32 neurons there
would be 32 files named `0.h5 ... 31.h5`.

Time series data is stored in a file named `im.h5` and there is a value for
each coordinate in each neuron at each time point.

The converter can generate 2 types of geometry, the neurons themselves or
a 3D Cartesian mesh with scalar fields interpolated onto it.

This project has a CMake build and depnds on HDF5, TBB, and VTK 8.2.0. Newer
versions of VTK are known to be problematic and will require a rewrite.



### Installing the converter
Once VTK, HDF5, and Intel TBB are installed one can install the converter.

```
git clone git@github.com:burlen/cvtn.git
mkdir cvtn-build
cmake -DVTK_DIR=../VTK-8.2.0-build/  ../cvtn
make
```

### Using the converter

```
cvt [input dir] [first neuron] [last neuron] [first step] [last step]
    [num cells] [output dir] [out file] [write geom] [write mesh] [n threads]
```
**input dir** : path where seg_coords and im.h live

**first neuron** : decimal number of the first neuron to load

**last neuron** : decimal number of the last neuron to load

**first step** : index of the first time step to load

**last step** : index of the last time step to load

**num cells** : number of cells in the longest side of the 3D Cartesian mesh.

**output dir** : path where to write the data as VTK files

**out file** : a name that's used when writeing output files

**write geometry** : 0/1 if 1 compute the neuron geometry

**write mesh** : 0/1 if 1 sample the scalar field onto a regular Cratesian mesh

**n threads** : number of threads to use (optional)


## plugin
**PVCvtNeuronReader - A reader plugin for ParaView**


### Building the plugin
```bash
git clone https://github.com/burlen/cvtn.git
mkdir cvtn-build
cd cvtn-build
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=`pwd`/../paraview-5.8.1-install  \
    -DParaView_DIR=`pwd`/../paraview-5.8.1-install/lib64/cmake/paraview-5.8/
    -DENABLE_PLUGIN=ON \
    ../cvtn

make -j
make -j install
```

### Loading the plugin in ParaView
Start ParaView open the *Tools->ManagePlugins* dialog box. Click *Load New*
button, and navigate to the plugin install and loacte and select the
*PVCvtnNeuronReader.so* file. Once back in the dialog expand the plugin's entry
and check *Autoload* so that the plugin is loaded when ParaView starts.

### Loadiing data
The plugin operates on a directory but must be given a file name due to how
ParaView works. In the top directory of the dataset make an empty file with the
*.bmtk* file extension.
```bash
touch master.bmtk
```
When this file is opened in ParaView by the *File->Open* dialog the
PVCvtnReaderPlugin will be used.
