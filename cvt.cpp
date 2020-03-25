#include "import.h"



using index_t = long;
using coord_t = float;
using data_t = float;


int main(int argc, char **argv)
{
    if (argc < 5)
    {
        std::cerr << "Usage: cvt [input dir] [first neuron] "
            << "[last neuron] [first step] [last step] [num cells] "
            << "[output dir] [out file] [n threads]" << std::endl;
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
    int nThreads = argc > 6 ? atoi(argv[9]) : -1;

    tbb::task_scheduler_init init(nThreads);

    // import the nerons
    int nNeuron = neuron1 - neuron0 + 1;
    std::cerr << "importing " << nNeuron << " neurons " << neuron0
        << " to " << neuron1 << " from " << inputDir << "...";

    std::vector<import::Neuron<index_t,coord_t,data_t>> neurons(nNeuron);

    tbb::parallel_for(tbb::blocked_range<int>(0,nNeuron),
        import::Importer<index_t,coord_t,data_t>(neuron0, inputDir, &neurons));

    std::cerr << "done!" << std::endl
        << "intializing time series, mesher, and exporter...";

    // initailze time series handlers
    import::TimeSeries<index_t, coord_t, data_t> timeSeries(inputDir, &neurons);
    if (timeSeries.initialize())
        return -1;

    import::Mesher<index_t, coord_t, data_t> mesher(&neurons);
    if (mesher.initialize(nCells))
        return -1;

    import::Exporter<index_t, coord_t, data_t> exporter;
    if (exporter.initialize(outputDir, outputFile, mesher))
        return -1;

    import::ReduceAverage<index_t, data_t> avg;
    import::ReduceMaxAbs<index_t, data_t> mxa;

    std::cerr << "done!" << std::endl;
    timeSeries.print();
    mesher.print();

    // process each time step
    int nSteps = step1 - step0 + 1;
    for (int i = step0; i <= step1; ++i)
    {
        std::cerr << "processing step " << i << "...";

        // update to this time step and interpolate new values to the neurons
        timeSeries.importTimeStep(i);

        // interpolate neuron values onto the mesh
        mesher.mesh(avg);
        mesher.mesh(mxa);

        // export the mesh
        const char *names[] = {avg.getName(), mxa.getName()};
        data_t *arrays[] = {avg.getResult(), mxa.getResult()};

        if (exporter.write(i, 2, arrays, names))
            return -1;

        avg.freeMem();
        mxa.freeMem();

        // export the neurons
        if (exporter.write(i, &neurons))
            return -1;

        std::cerr << "done!" << std::endl;
    }


    std::cerr << "cleaning up...";

    for (int i = 0; i < nNeuron; ++i)
        neurons[i].freeMem();

    timeSeries.freeMem();
    mesher.freeMem();
    exporter.freeMem();

    return 0;
}
