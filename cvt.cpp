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
            << "[n threads]" << std::endl;
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

    // import the nerons
    int nNeuron = neuron1 - neuron0 + 1;

    std::cerr << "importing " << nNeuron << " neurons " << neuron0
        << " to " << neuron1 << " from " << inputDir << "...";
    t0 = std::chrono::high_resolution_clock::now();

    std::vector<neuron::Neuron<index_t,coord_t,data_t>> neurons(nNeuron);

    tbb::parallel_for(tbb::blocked_range<int>(0,nNeuron),
        neuron::Importer<index_t,coord_t,data_t>(neuron0, inputDir, &neurons));

    t1 = std::chrono::high_resolution_clock::now();
    std::cerr << "done! ("
        << 1.e-6*std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()
        << "s)" << std::endl
        << "intializing time series, mesher, and exporter...";
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

    neuron::ReduceAverage<index_t, data_t> avg(timeSeries.scalarName);
    neuron::ReduceMaxAbs<index_t, data_t> mxa(timeSeries.scalarName);

    t1 = std::chrono::high_resolution_clock::now();
    std::cerr << "done! ("
        << 1.e-6*std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()
        << "s)" << std::endl;
    timeSeries.print();
    mesher.print();

    // process each time step
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
            mesher.mesh(avg);
            mesher.mesh(mxa);

            // export the mesh
            const char *names[] = {avg.getName(), mxa.getName()};
            data_t *arrays[] = {avg.getResult(), mxa.getResult()};

            if (exporter.writeMesh(i, 2, arrays, names))
                return -1;

            avg.freeMem();
            mxa.freeMem();
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
