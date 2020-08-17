#include "neuron.h"

#include <chrono>
#include <iostream>

using index_t = long;
using coord_t = float;
using data_t = float;


int main(int argc, char **argv)
{
    std::chrono::high_resolution_clock::time_point t0;
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point st;
    std::chrono::high_resolution_clock::time_point et;
    st = std::chrono::high_resolution_clock::now();

    if (argc < 12)
    {
        std::cerr << "Usage: cvt [input dir] [first neuron] "
            << "[last neuron] [first step] [last step] [num cells] "
            << "[output dir] [out file] [write geom] [write mesh] "
            << "[n threads]" << std::endl
            << std::endl
            << "  input dir - path where seg_coords and im.h live" << std::endl
            << "  first neuron - decimal number of the first neuron to load" << std::endl
            << "  last neuron - decimal number of the last neuron to load" << std::endl
            << "  first step - index of the first time step to load" << std::endl
            << "  last step - index of the last time step to load" << std::endl
            << "  num cells - number of cells in the longest side of the Cartesian mesh" << std::endl
            << "  output dir - path where to write the data as VTK files" << std::endl
            << "  out file - name of the data set" << std::endl
            << "  write geometry - 0/1 if 1 compute the neuron geometry" << std::endl
            << "  write mesh - 0/1 if 1 sample the scalar field onto a " << std::endl
            << "               regular Cratesian mesh" << std::endl
            << "  n threads - number of threads to use (optional)" << std::endl
            << std::endl;
        return -1;
    }

    const char *inputDir = argv[1];
    int neuron0 = atoi(argv[2]);
    int neuron1 = atoi(argv[3]);
    int step0 = atoi(argv[4]);
    int step1 = atoi(argv[5]);
    index_t nCells = atoi(argv[6]);
    const char *outputDir = argv[7];
    const char *outputFile = argv[8];
    int writeGeom = atoi(argv[9]);
    int writeMesh = atoi(argv[10]);
    int nThreads = argc > 11 ? atoi(argv[11]) : -1;

    tbb::task_scheduler_init init(nThreads);
    std::cerr << "initializing with " << nThreads << " threads" << std::endl;

    // import the nerons
    int nNeuron = neuron1 - neuron0 + 1;

    std::cerr << "importing " << nNeuron << " neurons " << neuron0
        << " to " << neuron1 << " from " << inputDir << " write geometry"
        << EYESNO(!writeGeom) << " write mesh" << EYESNO(!writeMesh)
        << " data will be stored in " << outputDir << " at " << outputFile << "...";
    t0 = std::chrono::high_resolution_clock::now();

    std::vector<neuron::Neuron<index_t,coord_t,data_t>> neurons(nNeuron);

    tbb::parallel_for(tbb::blocked_range<int>(0,nNeuron),
        neuron::Importer<index_t,coord_t,data_t>(neuron0, inputDir, &neurons));

    t1 = std::chrono::high_resolution_clock::now();
    std::cerr << "done! ("
        << 1.e-6*std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()
        << "s)" << std::endl
        << "intializing time series, mesher, and exporter...";

    // disable hdf5 error spew.
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);

    t0 = std::chrono::high_resolution_clock::now();

    // initailze time series handlers
    neuron::TimeSeries<index_t, coord_t, data_t> timeSeries(inputDir, &neurons);
    if (timeSeries.initialize())
        return -1;

    neuron::Mesher<index_t, coord_t, data_t> mesher(&neurons);
    if (mesher.initialize(nCells))
        return -1;

    neuron::Exporter<index_t, coord_t, data_t> exporter;
    if (exporter.initialize(outputDir, outputFile, mesher))
        return -1;

    int nScalars = timeSeries.nScalars;
    std::vector<neuron::ReduceAverage<index_t, data_t>*> avg(nScalars);
    for (int i = 0; i < nScalars; ++i)
        avg[i] = new neuron::ReduceAverage<index_t, data_t>(timeSeries.scalarName[i]);

    std::vector<neuron::ReduceMaxAbs<index_t, data_t>*> mxa(nScalars);
    for (int i = 0; i < nScalars; ++i)
        mxa[i] = new neuron::ReduceMaxAbs<index_t, data_t>(timeSeries.scalarName[i]);

    int nReductions = 2*nScalars;
    std::vector<const char *>names(nReductions);
    std::vector<data_t *>arrays(nReductions);

    t1 = std::chrono::high_resolution_clock::now();
    std::cerr << "done! ("
        << 1.e-6*std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()
        << "s)" << std::endl;
    timeSeries.print();
    mesher.print();

    // process each time step
    step1 = (step1 > 0 ? step1 :timeSeries.nSteps - 1);
    int nSteps = step1 - step0 + 1;
    std::cerr << "processing " << nSteps << "..." << std::endl;
    for (int i = step0; i <= step1; ++i)
    {
        std::cerr << "processing step " << i << "...";
        t0 = std::chrono::high_resolution_clock::now();

        // update to this time step and interpolate new values to the neurons
        timeSeries.importTimeStep(i);

        if (writeMesh)
        {
            // interpolate neuron values onto the mesh
            for (int j = 0; j < nScalars; ++j)
            {
                // interpolate (average)
                mesher.mesh(*(avg[j]), j);
                // marshall results
                names[j] = avg[j]->getName();
                arrays[j] = avg[j]->getResult();

                // interpolate (maximum)
                mesher.mesh(*(mxa[j]), j);
                // marshall results
                int jj = 2*j;
                names[jj] = mxa[j]->getName();
                arrays[jj] = mxa[j]->getResult();
            }

            // export the mesh
            if (exporter.writeMesh(i, timeSeries.nScalars, arrays.data(), names.data()))
                return -1;

            for (int j = 0; j < nScalars; ++j)
            {
                avg[j]->freeMem();
                mxa[j]->freeMem();
            }
        }

        if (writeGeom)
        {
            // export the neurons
            if (exporter.writeCells(i, &neurons))
                return -1;
        }

        t1 = std::chrono::high_resolution_clock::now();
        std::cerr << "done! ("
            << 1.e-6*std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()
            << "s)" << std::endl;
    }

    for (int i = 0; i < nScalars; ++i)
    {
        delete avg[i];
        delete mxa[i];
    }


    std::cerr << "cleaning up...";
    t0 = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < nNeuron; ++i)
        neurons[i].freeMem();

    timeSeries.freeMem();
    mesher.freeMem();
    exporter.freeMem();

    t1 = std::chrono::high_resolution_clock::now();
    std::cerr << "done! ("
        << 1.e-6*std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()
        << "s)" << std::endl;

    et = std::chrono::high_resolution_clock::now();
    std::cerr << "total run time "
        << 1.e-6*std::chrono::duration_cast<std::chrono::microseconds>(et - st).count()
        << "s" << std::endl;

    return 0;
}
