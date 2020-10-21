#include <tbb/queuing_mutex.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <chrono>
#include <iostream>

namespace neuron
{

// this is meant to prevent multiple concurrent calls into HDF5
// and it is also used to serialize debug output when MEM_CHECK
// is defined
tbb::queuing_mutex ioMutex;


// --------------------------------------------------------------------------
int neuronFileExists(const char *dirName, int nid)
{
    char fn[512];
    snprintf(fn, 511, "%s/seg_coords/%d.h5", dirName, nid);
    int fd = open(fn, O_RDONLY);
    if (fd == -1)
    {
        return 0;
    }
    close(fd);
    return 1;
}

// --------------------------------------------------------------------------
int imFileExists(const char *dirName)
{
    char fn[512];
    snprintf(fn, 511, "%s/im.h5", dirName);
    int fd = open(fn, O_RDONLY);
    if (fd == -1)
    {
        return 0;
    }
    close(fd);
    return 1;
}

// --------------------------------------------------------------------------
int scanForNeurons(const char *dirName, int &lastId)
{
    // scan for an id that fails. this locates the next power of 2 above the
    // the last valid file
    int bot = 1;
    int top = 1;
    while (neuronFileExists(dirName, top - 1))
    {
        bot = top;
        top *= 2;
    }

    // use bisection to locate the last valid file.
    --bot;
    --top;
    while ((top - bot) > 1)
    {
        int next = bot + (top - bot) / 2;
        if (neuronFileExists(dirName, next))
        {
            bot = next;
        }
        else
        {
            top = next;
        }
    }

    // verify
    if (!neuronFileExists(dirName, bot))
        return -1;

    lastId = bot;

    return 0;
}

}
