#include <cstring>
#include <iostream>

#include <hdf5.h>
#include <hdf5_hl.h>
#include <tbb/tbb.h>

#include <vtkType.h>
#include <vtkSOADataArrayTemplate.h>
#include <vtkAOSDataArrayTemplate.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkIdTypeArray.h>
#include <vtkPolyData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkPolyDataWriter.h>

// defining MEM_CHECK will result in a walk of the dataset
// every point of every cell is sent to stderr stream used
// in conjunction with valgrind this validates the cells
//#define MEM_CHECK

// this is meant to prevent multiple concurrent calls into HDF5
// and it is also used to serialize debug output when MEM_CHECK
// is defined
tbb::mutex ioMutex;

// --------------------------------------------------------------------------
#define ERROR(_msg)                                                 \
{                                                                   \
    std::cerr << "Error: ["                                         \
        << __FILE__ << ":" << __LINE__ << "] " _msg << std::endl;   \
}

// --------------------------------------------------------------------------
template<typename cppT>
struct cppH5Tt {};

#define declareCppH5Tt(_cppT, _h5Name)                                  \
template<> struct cppH5Tt<_cppT>                                        \
{                                                                       \
    static                                                              \
    int read(hid_t fh, const char *dsn, _cppT *buf)                     \
    {                                                                   \
        herr_t ierr = H5LTread_dataset_ ## _h5Name (fh, dsn, buf);      \
        if (ierr < 0)                                                   \
            return -1;                                                  \
        return 0;                                                       \
    }                                                                   \
};

declareCppH5Tt(char, char)
declareCppH5Tt(short, short)
declareCppH5Tt(int, int)
declareCppH5Tt(long, long)
declareCppH5Tt(float, float)
declareCppH5Tt(double, double)

// --------------------------------------------------------------------------
template <typename T>
bool equal(T a, T b, T tol)
{
    T diff = std::abs(a - b);
    a = std::abs(a);
    b = std::abs(b);
    b = (b > a) ? b : a;
    if (diff <= (b*tol))
        return true;
    return false;
}

// --------------------------------------------------------------------------
template<typename index_t, typename coord_t>
int read(const std::string &baseDir, int nId, coord_t *&p0,
    coord_t *&p5, coord_t *&p1, index_t &nSimpleCells)
{
    tbb::mutex::scoped_lock lock(ioMutex);

    char fn[256];
    snprintf(fn, 256, "%s/%d.h5", baseDir.c_str(), nId);

    hid_t fh = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fh < 0)
    {
        ERROR("Failed to open " << fn)
        return -1;
    }

    int nDims = 0;
    if (H5LTget_dataset_ndims(fh, "/p0", &nDims) < 0)
    {
        H5Fclose(fh);
        ERROR("Failed to get coordinate dimensions")
        return -1;
    }

    if (nDims != 2)
    {
        H5Fclose(fh);
        ERROR("Coordinates have " << nDims << " but we require 2")
        return -1;
    }

    hsize_t dims[2] = {0};
    H5T_class_t elemClass;
    size_t elemSize = 0;
    if (H5LTget_dataset_info(fh, "/p0", dims, &elemClass, &elemSize) < 0)
    {
        H5Fclose(fh);
        ERROR("Failed to get dataset info")
        return -1;
    }

    index_t ptSize = dims[0];
    nSimpleCells = dims[1];
    index_t bufSize = ptSize*nSimpleCells;

    p0 = (coord_t*)malloc(bufSize*sizeof(coord_t));
    if (cppH5Tt<coord_t>::read(fh, "/p0", p0))
    {
        free(p0);
        H5Fclose(fh);
        ERROR("Failed to read p0")
        return -1;
    }

    p5 = (coord_t*)malloc(bufSize*sizeof(coord_t));
    if (cppH5Tt<coord_t>::read(fh, "/p05", p5))
    {
        free(p0);
        free(p5);
        H5Fclose(fh);
        ERROR("Failed to read p05")
        return -1;
    }

    p1 = (coord_t*)malloc(bufSize*sizeof(coord_t));
    if (cppH5Tt<coord_t>::read(fh, "/p1", p1))
    {
        free(p0);
        free(p5);
        free(p1);
        H5Fclose(fh);
        ERROR("Failed to read p1")
        return -1;
    }

    H5Fclose(fh);
    return 0;
}

// --------------------------------------------------------------------------
template<typename index_t, typename coord_t>
int thickness(index_t nComplexCells, index_t *complexCells,
    index_t *complexCellLens,  index_t *complexCellLocs,
    index_t nPts, coord_t *x, coord_t *y, coord_t *z,
    coord_t *dist, coord_t *&thick)
{
    thick = (coord_t*)malloc(nPts*sizeof(coord_t));

    for (index_t i = 0; i < nPts; ++i)
        thick[i] = 1.0 - dist[i];

    for (index_t i = 0; i < nComplexCells; ++i)
    {
        // last point in the cell has a 0 thickness
        index_t q = complexCellLocs[i] + complexCellLens[i] - 1;
        index_t ii = complexCells[q];
        thick[ii] = 0.0;
    }

    return 0;
}

// --------------------------------------------------------------------------
template<typename index_t, typename coord_t>
int distance(index_t nComplexCells, index_t *complexCells,
    index_t *complexCellLens,  index_t *complexCellLocs,
    index_t nPts, coord_t *x, coord_t *y, coord_t *z,
    coord_t *&dist)
{
    dist = (coord_t*)malloc(nPts*sizeof(coord_t));
    coord_t max_dist = 0.0;

    // compute the distance from root
    index_t p5c0 = complexCells[1];
    coord_t x0 = x[p5c0];
    coord_t y0 = y[p5c0];
    coord_t z0 = z[p5c0];

    for (index_t i = 0; i < nPts; ++i)
    {
        coord_t xi = x[i];
        coord_t yi = y[i];
        coord_t zi = z[i];

        coord_t r = xi - x0;
        r *= r;

        coord_t dy2 = yi - y0;
        dy2 *= dy2;
        r += dy2;

        coord_t dz2 = zi - z0;
        dz2 *= dz2;
        r += dz2;

        r = sqrt(r);

        max_dist = r > max_dist ? r : max_dist;
        dist[i] = r;
    }

    // normalize
    for (index_t i = 0; i < nPts; ++i)
        dist[i] /= max_dist;

    return 0;
}

#if defined(MEM_CHECK)
// --------------------------------------------------------------------------
template<typename index_t, typename coord_t>
void print(index_t nSimpleCells, index_t *simpleCells,
    index_t nComplexCells, index_t *complexCells,
    index_t *complexCellLens, index_t *complexCellLocs,
    index_t complexCellArraySize, index_t nPts,
    coord_t *x, coord_t *y, coord_t *z)
{
    tbb::mutex::scoped_lock lock(ioMutex);

    std::cerr << "simpleCells(" << nSimpleCells << ")=";
    for (index_t i = 0; i < 3*nSimpleCells; ++i)
        std::cerr << simpleCells[i] << ", ";
    std::cerr << std::endl;

    std::cerr << "complexCells(" << nComplexCells << ")=";
    for (index_t i = 0; i < complexCellArraySize; ++i)
        std::cerr << complexCells[i] << ", ";
    std::cerr << std::endl;

    std::cerr << "complexCellLens=";
    for (index_t i = 0; i < nComplexCells; ++i)
        std::cerr << complexCellLens[i] << ", ";
    std::cerr << std::endl;

    std::cerr << "complexCellLocs=";
    for (index_t i = 0; i < nComplexCells; ++i)
        std::cerr << complexCellLocs[i] << ", ";
    std::cerr << std::endl;

    std::cerr << "pts(" <<  nPts << ")=";
    /*for (index_t i = 0; i < nPts; ++i)
        std::cerr << "(" << x[i] << ", " << y[i] << ", " << z[i] << "), ";
    std::cerr << std::endl;*/

    index_t tot = 0;
    index_t *cc = complexCells;
    for (index_t i = 0; i < nComplexCells; ++i)
    {
        index_t celLen = complexCellLens[i];
        std::cerr << i << " " << celLen << std::endl;
        for (index_t j = 0; j < celLen; ++j)
            std::cerr << "    " << cc[j] << " -> " << x[cc[j]] << ", "
                << y[cc[j]] << ", " << z[cc[j]] << std::endl;
        cc += celLen;
        tot += celLen;
    }
    std::cerr << "tot=" << tot + nComplexCells << std::endl;
}
#endif

// --------------------------------------------------------------------------
template<typename index_t, typename coord_t>
int clean(index_t &nPts, coord_t *&x, coord_t *&y,
    coord_t *&z, index_t &nSimpleCells, index_t *&simpleCells,
    index_t &nComplexCells, index_t *&complexCells,
    index_t *&complexCellLens, index_t *&complexCellLocs,
    index_t &complexCellArraySize, coord_t tol=1.e-3)
{
    // allocate space for the complex cells.
    nComplexCells = 0;
    complexCellArraySize = 3*3*nSimpleCells;
    index_t nbytes = complexCellArraySize*sizeof(index_t);
    complexCells = (index_t*)malloc(nbytes);
    memset(complexCells, 0, nbytes);
    complexCellArraySize = 0;

    // random access structures
    nbytes = nSimpleCells*sizeof(index_t);
    complexCellLens = (index_t*)malloc(nbytes);
    memset(complexCellLens, 0, nbytes);

    // random access structure
    complexCellLocs = (index_t*)malloc(nbytes);
    memset(complexCellLocs, 0, nbytes);

    // this mask indicates if we deleted a given point already
    index_t nnPts = 0;
    index_t nDel = 0;
    nbytes = nPts*sizeof(int);
    int *del = (int*)malloc(nbytes);
    memset(del, 0, nbytes);

    // this mask indicates if the cell has been checked
    nbytes = nSimpleCells*sizeof(int);
    int *visit = (int*)malloc(nbytes);
    memset(visit, 0, nbytes);

    // check each point for duplicates
    index_t *cc = complexCells;
    for (index_t i = 0; i < nSimpleCells; ++i)
    {
        // skip cells we already checked
        if (visit[i])
            continue;

        // start the new complex cell
        complexCellLocs[nComplexCells] = complexCellArraySize;

        index_t *ccl = complexCellLens + nComplexCells;
        *ccl = 3;

        index_t ii = 3*i;
        cc[0] = simpleCells[ii];
        cc[1] = simpleCells[ii + 1];
        cc[2] = simpleCells[ii + 2];
        cc += 3;

        complexCellArraySize += 3;
        nComplexCells += 1;
        nnPts += 3;

        // check the last point in this cell
        ii = 3*i + 2;
        coord_t xi = x[ii];
        coord_t yi = y[ii];
        coord_t zi = z[ii];

        // mark this as visited
        visit[i] = 1;

        index_t j = 0;
        while (j < nSimpleCells)
        {
            // skip cells we already checked
            if (visit[j])
            {
                j += 1;
                continue;
            }

            // check first point in other cells
            index_t jj = 3*j;

            // this point was already deleted, don't bother testing it again.
            if (del[jj])
            {
                j += 1;
                continue;
            }

            coord_t xj = x[jj];
            coord_t yj = y[jj];
            coord_t zj = z[jj];

            // check if they are the same
            if (equal(xi,xj,tol) && equal(yi,yj,tol) && equal(zi,zj,tol))
            {
                // mark cell as visited
                visit[j] = 1;

                // mark point as deleted
                del[jj] = 1;
                nDel += 1;

                // update the map
                simpleCells[jj] = ii;

                // before merging the cell, finish the scan because we only can
                // have one merger per duplicate, but there may be more than
                // one deuplicate in the set and we need to remove all
                for (index_t k = j + 1; k < nSimpleCells; ++k)
                {
                    // skip cells we already checked
                    if (visit[k])
                        continue;

                    // first point in other cells
                    index_t kk = 3*k;

                    // this point was already deleted, don't bother testing it again.
                    if (del[kk])
                        continue;

                    coord_t xk = x[kk];
                    coord_t yk = y[kk];
                    coord_t zk = z[kk];

                    // check if they are the same
                    if (equal(xi,xk,tol) && equal(yi,yk,tol) && equal(zi,zk,tol))
                    {
                        // don't mark this cell as visited because it is not
                        // being merged here

                        // mark as deleted
                        del[kk] = 1;
                        nDel += 1;

                        // update the map
                        simpleCells[kk] = ii;
                    }
                }

                // merge the cell.
                *ccl += 2;

                cc[0] = 3*j + 1;
                cc[1] = 3*j + 2;
                cc += 2;

                complexCellArraySize += 2;
                nnPts += 2;

                // update the test point
                ii = 3*j + 2;
                xi = x[ii];
                yi = y[ii];
                zi = z[ii];

                // rescan with the new point
                j = 0;
                continue;
            }

            j += 1;
        }
    }

    // transfer the remaining points and relable the cells
    nnPts = nPts - nDel;
    nbytes = nnPts*sizeof(coord_t);
    coord_t *nx = (coord_t*)malloc(nbytes);
    coord_t *ny = (coord_t*)malloc(nbytes);
    coord_t *nz = (coord_t*)malloc(nbytes);

    index_t simpleCellArraySize = 3*nSimpleCells;
    nbytes = simpleCellArraySize*sizeof(index_t);
    index_t *oSimpleCells = (index_t*)malloc(nbytes);
    memcpy(oSimpleCells, simpleCells, nbytes);

    nbytes = complexCellArraySize*sizeof(index_t);
    index_t *oComplexCells = (index_t*)malloc(nbytes);
    memcpy(oComplexCells, complexCells, nbytes);

    for (index_t i = 0, j = 0; i < nPts; ++i)
    {
        if (del[i])
        {
            // every point id greater than this is shifted
            for (index_t q = 0; q < simpleCellArraySize; ++q)
            {
                index_t shift =  oSimpleCells[q] > i ? 1 : 0;
                simpleCells[q] -= shift;
            }

            for (index_t q = 0; q < complexCellArraySize; ++q)
            {
                index_t shift =  oComplexCells[q] > i ? 1 : 0;
                complexCells[q] -= shift;
            }
        }
        else
        {
            // copy this point
            nx[j] = x[i];
            ny[j] = y[i];
            nz[j] = z[i];
            j += 1;
        }
    }
    free(oSimpleCells);
    free(oComplexCells);
    free(visit);
    free(del);

    // update the points arrays
    free(x);
    free(y);
    free(z);

    x = nx;
    y = ny;
    z = nz;

    nPts = nnPts;

#if defined(MEM_CHECK)
    // this dumps the state, and in the process touches all
    // values. if there is a bad cell this should flag it
    print(nSimpleCells, simpleCells, nComplexCells, complexCells,
        complexCellLens, complexCellLocs, complexCellArraySize,
        nPts, x, y, z);
#endif

    return 0;
}

// --------------------------------------------------------------------------
template<typename index_t, typename coord_t>
int initialize(index_t &nSimpleCells, coord_t *p0, coord_t *p5,
     coord_t*p1, index_t &nPts, coord_t *&x, coord_t *&y, coord_t *&z,
    index_t *&simpleCells)
{
    // each simple cell has 3 points. first in p0 second in p5 and last in p1
    nPts = 3*nSimpleCells;

    // arrange x-coords for each simple cell
    x = (coord_t*)malloc(nPts*sizeof(coord_t));
    coord_t *dst = x;
    coord_t *src = p0;
    for (index_t i = 0; i < nSimpleCells; ++i)
        dst[3*i] = src[i];

    dst = x+1;
    src = p5;
    for (index_t i = 0; i < nSimpleCells; ++i)
        dst[3*i] = src[i];

    dst = x+2;
    src = p1;
    for (index_t i = 0; i < nSimpleCells; ++i)
        dst[3*i] = src[i];

    // arrange y-coords for each simple cell
    y = (coord_t*)malloc(nPts*sizeof(coord_t));
    dst = y;
    src = p0 + nSimpleCells;
    for (index_t i = 0; i < nSimpleCells; ++i)
        dst[3*i] = src[i];

    dst = y+1;
    src = p5 + nSimpleCells;
    for (index_t i = 0; i < nSimpleCells; ++i)
        dst[3*i] = src[i];

    dst = y+2;
    src = p1 + nSimpleCells;
    for (index_t i = 0; i < nSimpleCells; ++i)
        dst[3*i] = src[i];

    // arrange y-coords for each simple cell
    z = (coord_t*)malloc(nPts*sizeof(coord_t));
    dst = z;
    src = p0 + 2*nSimpleCells;
    for (index_t i = 0; i < nSimpleCells; ++i)
        dst[3*i] = src[i];

    dst = z+1;
    src = p5 + 2*nSimpleCells;
    for (index_t i = 0; i < nSimpleCells; ++i)
        dst[3*i] = src[i];

    dst = z+2;
    src = p1 + 2*nSimpleCells;
    for (index_t i = 0; i < nSimpleCells; ++i)
        dst[3*i] = src[i];

    // with the above ordering initialize the simple cells
    simpleCells = (index_t*)malloc(3*nSimpleCells*sizeof(index_t));
    for (index_t i = 0; i < nPts; ++i)
        simpleCells[i] = i;

    return 0;
}

// --------------------------------------------------------------------------
template<typename index_t, typename coord_t>
int packageVtk(int id, index_t nPts, coord_t *x, coord_t *y, coord_t *z,
    index_t nCells, index_t *cells, index_t *cellLens, index_t *cellLocs,
    index_t cellArraySize, coord_t *dist, coord_t *thick,
    vtkPolyData *&dataset)
{
    // convert cells into VTK's format
    vtkIdTypeArray *cellIds = vtkIdTypeArray::New();
    cellIds->SetName("cellIds");
    cellIds->SetNumberOfTuples(cellArraySize + nCells);
    vtkIdType *dst = cellIds->GetPointer(0);
    for (index_t i = 0; i < nCells; ++i)
    {
        index_t *src = cells + cellLocs[i];
        index_t cellLen = cellLens[i];
        dst[0] = cellLen;
        dst += 1;
        for (index_t j = 0; j < cellLen; ++j)
            dst[j] = src[j];
        dst += cellLen;
    }

    vtkCellArray *cellArray = vtkCellArray::New();
    cellArray->SetCells(nCells, cellIds);
    cellIds->Delete();

#if defined(ZERO_COPY_VTK)
    // zero copy points
    vtkSOADataArrayTemplate<coord_t> *coords =
        vtkSOADataArrayTemplate<coord_t>::New();

    coords->SetNumberOfComponents(3);
    coords->SetArray(0, x, nPts, true, true, 0);
    coords->SetArray(1, y, nPts, false, true, 0);
    coords->SetArray(2, z, nPts, false, true, 0);
    coords->SetName("coords");
#else
    // deep copy points
    vtkAOSDataArrayTemplate<coord_t> *coords =
        vtkAOSDataArrayTemplate<coord_t>::New();

    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(nPts);
    coord_t *pcoords = coords->GetPointer(0);
    for (index_t i = 0; i < nPts; ++i)
    {
        index_t ii = 3*i;
        pcoords[ii  ] = x[i];
        pcoords[ii+1] = y[i];
        pcoords[ii+2] = z[i];
    }
#endif

    vtkPoints *points = vtkPoints::New();
    points->SetData(coords);
    coords->Delete();

    dataset = vtkPolyData::New();
    dataset->SetPoints(points);
    dataset->SetLines(cellArray);
    points->Delete();
    cellArray->Delete();

    vtkIntArray *nid = vtkIntArray::New();
    nid->SetNumberOfTuples(nCells);
    nid->FillValue(id);
    nid->SetName("neuron");
    dataset->GetCellData()->AddArray(nid);
    nid->Delete();

    nid = vtkIntArray::New();
    nid->SetNumberOfTuples(nPts);
    nid->FillValue(id);
    nid->SetName("neuron");
    dataset->GetPointData()->AddArray(nid);
    nid->Delete();

    vtkAOSDataArrayTemplate<coord_t> *distance =
        vtkAOSDataArrayTemplate<coord_t>::New();
    distance->SetName("distance");
#if defined(ZERO_COPY_VTK)
    distance->SetArray(dist, nPts, true);
#else
    distance->SetNumberOfTuples(nPts);
    memcpy(distance->GetPointer(0), dist, nPts*sizeof(coord_t));
#endif
    dataset->GetPointData()->AddArray(distance);
    distance->Delete();

    vtkAOSDataArrayTemplate<coord_t> *thickness =
        vtkAOSDataArrayTemplate<coord_t>::New();
    thickness->SetName("thickness");
#if defined(ZERO_COPY_VTK)
    thickness->SetArray(dist, nPts, true);
#else
    thickness->SetNumberOfTuples(nPts);
    memcpy(thickness->GetPointer(0), thick, nPts*sizeof(coord_t));
#endif
    dataset->GetPointData()->AddArray(thickness);
    thickness->Delete();

    return 0;
}

template<typename index_t, typename coord_t>
struct ImportNeuron
{
    ImportNeuron() :
        blockId(-1), neuronId(-1), baseDir(nullptr), nSimpleCells(0),
        p0(nullptr), p5(nullptr), p1(nullptr), nPts(0), x(nullptr),
        y(nullptr), z(nullptr), simpleCells(nullptr), nComplexCells(0),
        complexCells(nullptr), complexCellLens(nullptr),
        complexCellLocs(nullptr), complexCellArraySize(0),
        dist(nullptr), thick(nullptr)
    {}

    ImportNeuron(int bid, int nid, const char *dir) :
        blockId(bid), neuronId(nid), baseDir(dir), nSimpleCells(0),
        p0(nullptr), p5(nullptr), p1(nullptr), nPts(0), x(nullptr),
        y(nullptr), z(nullptr), simpleCells(nullptr), nComplexCells(0),
        complexCells(nullptr), complexCellLens(nullptr),
        complexCellLocs(nullptr), complexCellArraySize(0),
        dist(nullptr), thick(nullptr)
    {}

    ~ImportNeuron() {}

    ImportNeuron(const ImportNeuron &) = delete;
    ImportNeuron(ImportNeuron &&) = default;

    ImportNeuron &operator=(const ImportNeuron &) = delete;
    ImportNeuron &operator=(ImportNeuron &&) = default;

    void operator()()
    {
        // read the data for this instance
        if (read(baseDir, neuronId, p0, p5, p1, nSimpleCells))
        {
            ERROR("Failed to load neuron " << neuronId)
            return;
        }

        // create the points and simple cells
        initialize(nSimpleCells, p0, p5, p1, nPts, x, y, z, simpleCells);

        // remove duplicate points and merge segements
        clean(nPts, x, y, z, nSimpleCells, simpleCells, nComplexCells,
            complexCells, complexCellLens, complexCellLocs,
            complexCellArraySize);

        // compute distance field
        distance(nComplexCells, complexCells, complexCellLens,
            complexCellLocs,  nPts, x, y, z, dist);

        thickness(nComplexCells, complexCells, complexCellLens,
            complexCellLocs,  nPts, x, y, z, dist, thick);
    }

    void freeMem()
    {
        free(p0);
        free(p5);
        free(p1);
        free(x);
        free(y);
        free(z);
        free(simpleCells);
        free(complexCells);
        free(complexCellLens);
        free(complexCellLocs);
        free(dist);
        free(thick);
        nSimpleCells = 0;
        p0 = nullptr;
        p5 = nullptr;
        p1 = nullptr;
        nPts = 0;
        x = nullptr;
        y = nullptr;
        z = nullptr;
        simpleCells = nullptr;
        nComplexCells = 0;
        complexCells = nullptr;
        complexCellLens = nullptr;
        complexCellLocs = nullptr;
        complexCellArraySize = 0;
        dist = nullptr;
        thick = nullptr;
    }

    void package(vtkPolyData *&block)
    {
        // package in VTK format
        packageVtk(neuronId, nPts, x, y, z, nComplexCells, complexCells,
            complexCellLens, complexCellLocs, complexCellArraySize,
            dist, thick, block);
    }

    int blockId;
    int neuronId;
    const char *baseDir;
    index_t nSimpleCells;
    coord_t *p0;
    coord_t *p5;
    coord_t *p1;
    index_t nPts;
    coord_t *x;
    coord_t *y;
    coord_t *z;
    index_t *simpleCells;
    index_t nComplexCells;
    index_t *complexCells;
    index_t *complexCellLens;
    index_t *complexCellLocs;
    index_t complexCellArraySize;
    coord_t *dist;
    coord_t *thick;
};



template<typename index_t, typename coord_t>
struct Importer
{
    Importer() = delete;

    Importer(int n0, const char *dir,
        std::vector<ImportNeuron<index_t,coord_t>> *n) :
            neuron0(n0), baseDir(dir), neurons(n)
    {}

    void operator()(const tbb::blocked_range<int> &r) const
    {
        for (int i = r.begin(); i != r.end(); ++i)
        {
            ImportNeuron<index_t,coord_t> import(i, neuron0 + i, baseDir);
            import();
            (*neurons)[i] = std::move(import);
#if defined(MEM_CHECK)
            std::cerr << ".";
#endif
        }
    }

    int neuron0;
    const char *baseDir;
    std::vector<ImportNeuron<index_t,coord_t>> *neurons;
};


using integer_t = vtkIdType;
using float_t = float;


int main(int argc, char **argv)
{
    if (argc < 5)
    {
        std::cerr << "Usage: cvt [base dir] [first neuron] "
            << "[last neuron] [out file] [n threads]" << std::endl;
        return -1;
    }

    const char *baseDir = argv[1];
    int neuron0 = atoi(argv[2]);
    int neuron1 = atoi(argv[3]);
    int nNeuron = neuron1 - neuron0 + 1;
    const char *ofBase = argv[4];
    int nThreads = argc > 5 ? atoi(argv[5]) : -1;

    tbb::task_scheduler_init init(nThreads);

    std::cerr << "importing " << nNeuron
        << " neurons from " << baseDir << "...";

    // import the data
    std::vector<ImportNeuron<integer_t,float_t>> neurons(nNeuron);

    tbb::parallel_for(tbb::blocked_range<int>(0,nNeuron),
        Importer<integer_t,float_t>(neuron0, baseDir, &neurons));

    std::cerr << "done!" << std::endl
        << "packaging as VTK...";

    // package as VTK
    vtkMultiBlockDataSet *mbds = vtkMultiBlockDataSet::New();
    mbds->SetNumberOfBlocks(nNeuron);

    for (int i = 0; i < nNeuron; ++i)
    {
        vtkPolyData *block = nullptr;
        neurons[i].package(block);
        mbds->SetBlock(i, block);
        block->Delete();
        neurons[i].freeMem();
    }

    // write VTK
    vtkXMLMultiBlockDataWriter *writer = vtkXMLMultiBlockDataWriter::New();
    writer->SetInputData(mbds);

    char outFileName[256];
    snprintf(outFileName, 256, "%s_%d_%d.%s", ofBase,
        neuron0, neuron1, writer->GetDefaultFileExtension());

    std::cerr << "done!" << std::endl
        << "writing to VTK at " << outFileName << "...";

    writer->SetFileName(outFileName);
    writer->Write();

    mbds->Delete();
    writer->Delete();

    std::cerr << "done!" << std::endl;

    return 0;
}
