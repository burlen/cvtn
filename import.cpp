#include <tbb/queuing_mutex.h>

namespace import
{

// this is meant to prevent multiple concurrent calls into HDF5
// and it is also used to serialize debug output when MEM_CHECK
// is defined
tbb::queuing_mutex ioMutex;

}
