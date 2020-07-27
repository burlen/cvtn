# Convert neuron simulation data into a format readable by ParaView

The tool can convert neruons with geometry stored in a `seg_coords` folder and
interpolate time series data stored in a `im.h5` file onto the geometry. In the
`seg_coords` folder it is expected to have one file for each neuron amd the
files named by the neuron number. For instance if there are 32 neurons there
would be 32 files named `0.h5 ... 31.h5`.

Time series data is stored in a file named `im.h5` and there is a value for
each coordinate in each neuron at each time point.

The converter can generate 2 types of geometry, the neurons themselves or
a 3D Cartesian mesh with scalar fields interpolated onto it.

This project has a CMake build and depnds on HDF5, TBB, and VTK 8.2.0. Newer
versions of VTK are known to be problematic and will require a rewrite.

## First
The dependencies need to be installed.

### VTK
After conversion the VTK file format is used to store data back on disk for
injestion into ParaView for visualization.

A typical VTK build might be acheived as follows:
```
wget https://www.vtk.org/files/release/8.2/VTK-8.2.0.tar.gz
tar xzfv VTK-8.2.0.tar.gz
mkdir VTK-8.2.0-build
cd VTK-8.2.0-build/
cmake -DCMAKE_INSTALL_PREFIX=../VTK-8.2.0-install ../VTK-8.2.0
make -j
make -j install
```

### HDF5
Any recent version of HDF5 should suffice. One might simply use a package manager.
for example: `sudo apt-get install gdf5-devel` or `brew install hdf5` etc etc.

### Intel TBB


## Install
Once VTK, HDF5, and Intel TBB are installed one can install the converter.

```
git clone git@github.com:burlen/cvtn.git
mkdir cvtn-build
cmake -DVTK_DIR=../VTK-8.2.0-build/  ../cvtn
make
```

## Use

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

