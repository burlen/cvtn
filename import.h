#ifndef import_h
#define import_h

#include <cstring>
#include <iostream>
#include <limits>

#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

#include <hdf5.h>
#include <hdf5_hl.h>

//#include <tbb/tbb.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/queuing_mutex.h>

#include <vtkType.h>
#include <vtkSOADataArrayTemplate.h>
#include <vtkAOSDataArrayTemplate.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkIdTypeArray.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkPolyDataWriter.h>

/*#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkPolyDataWriter.h>*/

// make memory aligned to cahce lines/avx register size
#define ALIGN_TO __BIGGEST_ALIGNMENT__
#define MALLOC(_nbytes) \
    aligned_alloc(ALIGN_TO, _nbytes)


namespace import
{

// defining MEM_CHECK will result in a walk of the dataset
// every point of every cell is sent to stderr stream used
// in conjunction with valgrind this validates the cells
//#define MEM_CHECK

// this is meant to prevent multiple concurrent calls into HDF5
// and it is also used to serialize debug output when MEM_CHECK
// is defined
extern tbb::queuing_mutex ioMutex;

// --------------------------------------------------------------------------
#define ERROR(_msg)                                                 \
{                                                                   \
    std::cerr << "Error: ["                                         \
        << __FILE__ << ":" << __LINE__ << "] " _msg << std::endl;   \
}

// --------------------------------------------------------------------------
#define EYESNO(_er) (_er == 0 ? "(yes)" : "(no)")

// --------------------------------------------------------------------------
int readDimensions(hid_t fh, const char *dsname, hsize_t *dims,
    int nDimsReq=1, hsize_t *dimsReq=nullptr)
{
    // get size of buffers for time series
    int nDims = 0;
    if (H5LTget_dataset_ndims(fh, dsname, &nDims) < 0)
    {
        ERROR("Failed to get number of \"" << dsname << "\" dimensions")
        return -1;
    }

    // verify that this dataset has the expected dimensionality
    if (nDims != nDimsReq)
    {
        ERROR("\"" << dsname << "\" has " << nDims << " but we require " << nDimsReq)
        return -1;
    }

    H5T_class_t elemClass;
    size_t elemSize = 0;
    if (H5LTget_dataset_info(fh, dsname, dims, &elemClass, &elemSize) < 0)
    {
        ERROR("Failed to get dataset info for \"" << dsname << "\"")
        return -1;
    }

    // user passed expected dimensions verify the dimensions
    if (dimsReq)
    {
        for (int i =0; i < nDimsReq; ++i)
        {
            if (dims[i] != dimsReq[i])
            {
                ERROR("In \"" << dsname << "\" expected "
                    << dimsReq[i] << " but found " << dims[i]
                    << " at dimension " << i)
                return -1;
            }
        }
    }

    return 0;
}

// --------------------------------------------------------------------------
int readDimensions(hid_t fh, const char *dsname, hsize_t *dims,
    int nDimsReq, hsize_t dimsReq)
{
    return readDimensions(fh, dsname, dims, nDimsReq, &dimsReq);
}

// --------------------------------------------------------------------------
template<typename cppT>
struct cppH5Tt {};

#define declareCppH5Tt(_cppT, _h5Name, _h5Type)                         \
template<> struct cppH5Tt<_cppT>                                        \
{                                                                       \
    static                                                              \
    hid_t typeCode() { return _h5Type; }                                \
                                                                        \
    static                                                              \
    int read(hid_t fh, const char *dsn, _cppT *buf)                     \
    {                                                                   \
        herr_t ierr = H5LTread_dataset_ ## _h5Name (fh, dsn, buf);      \
        if (ierr < 0)                                                   \
            return -1;                                                  \
        return 0;                                                       \
    }                                                                   \
};

declareCppH5Tt(char, char, H5T_NATIVE_CHAR)
declareCppH5Tt(short, short, H5T_NATIVE_SHORT)
declareCppH5Tt(int, int, H5T_NATIVE_INT)
declareCppH5Tt(long, long, H5T_NATIVE_LONG)
declareCppH5Tt(float, float, H5T_NATIVE_FLOAT)
declareCppH5Tt(double, double, H5T_NATIVE_DOUBLE)

// specialization for reading strings
template<> struct cppH5Tt<char **>
{
    static
    int read(hid_t fh, const char *dsn, char *&str, size_t maxlen=64)
    {
        hid_t dset = H5Dopen(fh, dsn, H5P_DEFAULT);
        if (dset < 0)
        {
            ERROR("Failed to open \"" << dsn << "\"")
            return -1;
        }

        hid_t ft = H5Dget_type(dset);
        hid_t dsp = H5Dget_space(dset);
        int ndims = H5Sget_simple_extent_ndims(dsp);
        if (ndims != 0)
        {
            ERROR("\"" << dsn << "\" not a scalar")
            return -1;
        }
        hid_t mt = H5Tcopy(ft);

        char *pstr = nullptr;
        hid_t ierr = H5Dread(dset, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pstr);
        if (ierr < 0)
        {
            ERROR("Failed to reqd \"" << dsn << "\"")
            return -1;
        }

        str = strndup(pstr, maxlen);

        H5Dvlen_reclaim(mt, dsp, H5P_DEFAULT, &pstr);
        H5Dclose(dset);
        H5Sclose(dsp);
        H5Tclose(ft);
        H5Tclose(mt);

        return 0;
    }
};

// --------------------------------------------------------------------------
template <typename index_t, typename data_t>
int readTimeStep(hid_t fh, index_t nSteps, index_t stepSize, index_t step,
    data_t *&data)
{
    if (fh < 0)
    {
        ERROR("Invalid file handle")
        return -1;
    }

    if (step >= nSteps)
    {
        H5Fclose(fh);
        fh = -1;
        ERROR("Invalid step " << step << " file has " << nSteps << " steps")
        return -1;
    }

    hid_t dset = H5Dopen(fh, "/data", H5P_DEFAULT);

    // select the time step
    hsize_t foffs[2] = {hsize_t(step), 0};
    hsize_t fsize[2] = {1, hsize_t(stepSize)};
    hid_t fdsp = H5Dget_space(dset);
    if (H5Sselect_hyperslab(fdsp, H5S_SELECT_SET, foffs, nullptr, fsize, nullptr) < 0)
    {
        H5Fclose(fh);
        fh = -1;
        ERROR("Failed to select time step" << step)
        return -1;
    }

    // allocate a buffer of the right side if needed
    if (data == nullptr)
        data = (data_t*)MALLOC(stepSize*sizeof(data_t));

    // read it in
    hid_t mdsp = H5Screate_simple(1, &fsize[1], nullptr);
    if (H5Dread(dset, cppH5Tt<data_t>::typeCode(), mdsp, fdsp, H5P_DEFAULT, data) < 0)
    {
        free(data);
        data = nullptr;
        H5Sclose(mdsp);
        H5Sclose(fdsp);
        H5Dclose(dset);
        H5Fclose(fh);
        fh = -1;
        ERROR("Failed to read time step  " << step)
        return -1;
    }

    H5Sclose(mdsp);
    H5Sclose(fdsp);
    H5Dclose(dset);
    //H5Fclose(fh);

    return 0;
}

// --------------------------------------------------------------------------
template <typename index_t, typename data_t>
int readTimeSeriesMetadata(const std::string &baseDir,
    hid_t &fh,
    index_t &nNeurons, index_t &nSteps, index_t &stepSize,
    data_t &t0, data_t &t1, data_t &dt,
    index_t *&neuronIds, index_t *&neuronOffs, index_t *&neuronSize)
{
    tbb::queuing_mutex::scoped_lock lock(ioMutex);

    if (fh < 0)
    {
        char fn[256];
        snprintf(fn, 256, "%s/im.h5", baseDir.c_str());

        fh = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (fh < 0)
        {
            ERROR("Failed to open " << fn)
            return -1;
        }
    }

    // read time metadata
    hsize_t dims[2] = {0};
    data_t tmd[3] = {data_t(0)};
    if (readDimensions(fh, "/mapping/time", dims, 1, 3) ||
        cppH5Tt<double>::read(fh, "/mapping/time", tmd))
    {
        H5Fclose(fh);
        fh = -1;
        return -1;
    }
    t0 = tmd[0];
    t1 = tmd[1];
    dt = tmd[2];

    // get size of buffers for time series
    if (readDimensions(fh, "/data", dims, 2))
    {
        H5Fclose(fh);
        fh = -1;
        return -1;
    }
    nSteps = dims[0];
    stepSize = dims[1];

    // get size of buffers for maps
    if (readDimensions(fh, "/mapping/gids", dims, 1))
    {
        H5Fclose(fh);
        fh = -1;
        return -1;
    }
    nNeurons = dims[0];

    // allocate buffers and read maps
    int egi=0, eei=0;
    neuronIds = (index_t*)MALLOC(nNeurons*sizeof(index_t));
    index_t *eids = (index_t*)MALLOC(stepSize*sizeof(index_t));
    if ((egi = cppH5Tt<index_t>::read(fh, "/mapping/gids", neuronIds)) ||
        (eei = cppH5Tt<index_t>::read(fh, "/mapping/element_id", eids)))
    {
        free(neuronIds);
        H5Fclose(fh);
        fh = -1;
        ERROR("Failed to read. neuronIds" << EYESNO(egi) << ", eids" << EYESNO(eei))
        return -1;
    }

    // calculate the offset and size of the neurons
    neuronOffs = (index_t*)MALLOC(nNeurons*sizeof(index_t));
    neuronSize = (index_t*)MALLOC(nNeurons*sizeof(index_t));
    for (index_t i = 0, q = 0; i < nNeurons; ++i)
    {
        neuronOffs[i] = q;
        while (q < stepSize)
        {
             ++q;
             if ((q == stepSize) || (eids[q] == 0))
                 break;
        }
        neuronSize[i] = eids[q-1] + 1;
    }

    free(eids);
    //H5Fclose(fh);
    return 0;
}


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
int readNeuron(const std::string &baseDir, int nId,
    coord_t *& __restrict__ p0, coord_t *& __restrict__ p5,
    coord_t *& __restrict__ p1, coord_t *& __restrict__ d0,
    coord_t *& __restrict__ d5, coord_t *& __restrict__ d1,
    coord_t spos[3], char *& e_type, char *& m_type, char *&ei,
    int &layer, index_t & nSimpleCells)
{
    tbb::queuing_mutex::scoped_lock lock(ioMutex);

    char fn[256];
    snprintf(fn, 256, "%s/seg_coords/%d.h5", baseDir.c_str(), nId);

    hid_t fh = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fh < 0)
    {
        ERROR("Failed to open " << fn)
        return -1;
    }

    hsize_t dims[2] = {0};
    if (readDimensions(fh, "/p0", dims, 2))
    {
        H5Fclose(fh);
        ERROR("Failed to get coordinate dimensions")
        return -1;
    }

    index_t ptSize = dims[0];
    nSimpleCells = dims[1];
    index_t bufSize = ptSize*nSimpleCells;

    int ep0=0,ep5=0,ep1=0,ed0=0,ed1=0,eet=0,emt=0,eei=0,esp=0,el=0;
    p0 = (coord_t*)MALLOC(bufSize*sizeof(coord_t));
    p5 = (coord_t*)MALLOC(bufSize*sizeof(coord_t));
    p1 = (coord_t*)MALLOC(bufSize*sizeof(coord_t));
    d0 = (coord_t*)MALLOC(bufSize*sizeof(coord_t));
    d5 = (coord_t*)MALLOC(bufSize*sizeof(coord_t));
    d1 = (coord_t*)MALLOC(bufSize*sizeof(coord_t));
    if ((ep0 = cppH5Tt<coord_t>::read(fh, "/p0", p0)) ||
        (ep5 = cppH5Tt<coord_t>::read(fh, "/p05", p5)) ||
        (ep1 = cppH5Tt<coord_t>::read(fh, "/p1", p1)) ||
        (ed0 = cppH5Tt<coord_t>::read(fh, "/d0", d0)) ||
        (ed1 = cppH5Tt<coord_t>::read(fh, "/d1", d1)) ||
        (esp = cppH5Tt<coord_t>::read(fh, "/soma_pos", spos)) ||
        (eet = cppH5Tt<char**>::read(fh, "/e_type", e_type)) ||
        (emt = cppH5Tt<char**>::read(fh, "/m_type", m_type)) ||
        (eei = cppH5Tt<char**>::read(fh, "/ei", ei)) ||
        (el = cppH5Tt<int>::read(fh, "/layer", &layer)))
    {
        free(p0);
        free(p5);
        free(p1);
        free(d0);
        free(d5);
        free(d1);
        H5Fclose(fh);
        ERROR("Failed to read. p0" << EYESNO(ep0) << ", p5" << EYESNO(ep5)
            << ", p1" << EYESNO(ep1) << ", d0" << EYESNO(ed0) << ", d1" <<EYESNO(ed1)
            << ", spos" << EYESNO(esp) << ", e_tyoe" << EYESNO(eet)
            << ", m_type" << EYESNO(emt) << ", ei" << EYESNO(eei)
            << ", layer" << EYESNO(el))
        return -1;
    }

    for (index_t i = 0; i < nSimpleCells; ++i)
        d5[i] = (d0[i] + d1[i])*coord_t(0.5);
/*
    // diameter to radius
    for (index_t i = 0; i < nSimpleCells; ++i)
        d0[i] *= 0.5;

    for (index_t i = 0; i < nSimpleCells; ++i)
        d5[i] *= 0.5;

    for (index_t i = 0; i < nSimpleCells; ++i)
        d1[i] *= 0.5;
*/
    H5Fclose(fh);
    return 0;
}

// --------------------------------------------------------------------------
template<typename index_t, typename coord_t>
int thickness(index_t nComplexCells, index_t * __restrict__ complexCells,
    index_t * __restrict__ complexCellLens,
    index_t * __restrict__ complexCellLocs,
    index_t nPts, coord_t *__restrict__ dist,
    coord_t *& __restrict__ thick)
{
    thick = (coord_t*)MALLOC(nPts*sizeof(coord_t));

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
int distance(index_t * __restrict__ complexCells,
    index_t nPts, coord_t * __restrict__ x, coord_t * __restrict__ y,
    coord_t * __restrict__ z, coord_t *& __restrict__ dist)
{
    dist = (coord_t*)MALLOC(nPts*sizeof(coord_t));

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

        dist[i] = r;
    }

    coord_t max_dist = 0.0;
    for (index_t i = 0; i < nPts; ++i)
        max_dist = dist[i] > max_dist ? dist[i] : max_dist;

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
    coord_t *x, coord_t *y, coord_t *z, coord_t *d)
{
    tbb::queuing_mutex::scoped_lock lock(ioMutex);

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
    for (index_t i = 0; i < nPts; ++i)
        std::cerr << "(" << x[i] << ", " << y[i] << ", "
            << z[i] << ", " << d[i] << "), ";
    std::cerr << std::endl;

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
int clean(index_t &nPts, coord_t *& __restrict__ x, coord_t *& __restrict__ y,
    coord_t *& __restrict__ z, coord_t *& __restrict__ d,
    index_t &nSimpleCells, index_t *& __restrict__ simpleCells,
    index_t &nComplexCells, index_t *& __restrict__ complexCells,
    index_t *& __restrict__ complexCellLens,
    index_t *& __restrict__ complexCellLocs,
    index_t &complexCellArraySize, coord_t tol=1.e-3)
{
    // allocate space for the complex cells.
    nComplexCells = 0;
    complexCellArraySize = 3*3*nSimpleCells;
    index_t nbytes = complexCellArraySize*sizeof(index_t);
    complexCells = (index_t*)MALLOC(nbytes);
    memset(complexCells, 0, nbytes);
    complexCellArraySize = 0;

    // random access structures
    nbytes = nSimpleCells*sizeof(index_t);
    complexCellLens = (index_t*)MALLOC(nbytes);
    memset(complexCellLens, 0, nbytes);

    // random access structure
    complexCellLocs = (index_t*)MALLOC(nbytes);
    memset(complexCellLocs, 0, nbytes);

    // this mask indicates if we deleted a given point already
    index_t nnPts = 0;
    index_t nDel = 0;
    nbytes = nPts*sizeof(int);
    int *del = (int*)MALLOC(nbytes);
    memset(del, 0, nbytes);

    // this mask indicates if the cell has been checked
    nbytes = nSimpleCells*sizeof(int);
    int *visit = (int*)MALLOC(nbytes);
    memset(visit, 0, nbytes);

    // check each point for duplicates
    index_t * __restrict__ cc = complexCells;
    for (index_t i = 0; i < nSimpleCells; ++i)
    {
        // skip cells we already checked
        if (visit[i])
            continue;

        // start the new complex cell
        complexCellLocs[nComplexCells] = complexCellArraySize;

        index_t * __restrict__ ccl = complexCellLens + nComplexCells;
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
    coord_t * __restrict__ nx = (coord_t*)MALLOC(nbytes);
    coord_t * __restrict__ ny = (coord_t*)MALLOC(nbytes);
    coord_t * __restrict__ nz = (coord_t*)MALLOC(nbytes);
    coord_t * __restrict__ nd = (coord_t*)MALLOC(nbytes);

    index_t simpleCellArraySize = 3*nSimpleCells;
    nbytes = simpleCellArraySize*sizeof(index_t);
    index_t * __restrict__ oSimpleCells = (index_t*)MALLOC(nbytes);
    memcpy(oSimpleCells, simpleCells, nbytes);

    nbytes = complexCellArraySize*sizeof(index_t);
    index_t * __restrict__ oComplexCells = (index_t*)MALLOC(nbytes);
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
            nd[j] = d[i];
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
    free(d);

    x = nx;
    y = ny;
    z = nz;
    d = nd;

    nPts = nnPts;

#if defined(MEM_CHECK)
    // this dumps the state, and in the process touches all
    // values. if there is a bad cell this should flag it
    print(nSimpleCells, simpleCells, nComplexCells, complexCells,
        complexCellLens, complexCellLocs, complexCellArraySize,
        nPts, x, y, z, d);
#endif

    return 0;
}

// --------------------------------------------------------------------------
template<typename index_t, typename coord_t>
int initialize(index_t & __restrict__ nSimpleCells,
    coord_t * __restrict__ p0, coord_t * __restrict__ p5,
    coord_t * __restrict__ p1, coord_t * __restrict__ d0,
    coord_t * __restrict__ d5, coord_t * __restrict__ d1,
    index_t & nPts,
    coord_t *& __restrict__ x, coord_t *& __restrict__ y,
    coord_t *& __restrict__ z, coord_t *& __restrict__ d,
    index_t *&simpleCells)
{
    // each simple cell has 3 points. first in p0 second in p5 and last in p1
    nPts = 3*nSimpleCells;

    // arrange x-coords for each simple cell
    x = (coord_t*)MALLOC(nPts*sizeof(coord_t));
    coord_t * __restrict__ dst = x;
    coord_t * __restrict__ src = p0;
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
    y = (coord_t*)MALLOC(nPts*sizeof(coord_t));
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
    z = (coord_t*)MALLOC(nPts*sizeof(coord_t));
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

    // arrange diameter for each simple cell
    d = (coord_t*)MALLOC(nPts*sizeof(coord_t));
    dst = d;
    src = d0;
    for (index_t i = 0; i < nSimpleCells; ++i)
        dst[3*i] = src[i];

    dst = d+1;
    src = d5;
    for (index_t i = 0; i < nSimpleCells; ++i)
        dst[3*i] = src[i];

    dst = d+2;
    src = d1;
    for (index_t i = 0; i < nSimpleCells; ++i)
        dst[3*i] = src[i];

    // with the above ordering initialize the simple cells
    simpleCells = (index_t*)MALLOC(3*nSimpleCells*sizeof(index_t));
    for (index_t i = 0; i < nPts; ++i)
        simpleCells[i] = i;

    return 0;
}

// --------------------------------------------------------------------------
template<typename index_t>
int pointCount(index_t * __restrict__ simpleCells, index_t nSimpleCells,
    index_t nPts, index_t *& __restrict__ ptCount)
{
    // allocate buffer
    if (!ptCount)
        ptCount = (index_t*)MALLOC(nPts*sizeof(index_t));

    // zero it out
    memset(ptCount, 0, nPts*sizeof(index_t));

    // because the simples cells have duplicated points that have been removed
    // count keeps track of how many values contributre to each point
    for (index_t i = 0; i < nSimpleCells; ++i)
    {
        index_t ii = 3*i;
        for (index_t j = 0; j < 3; ++j)
        {
            index_t idx = simpleCells[ii+j];
            ptCount[idx] += 1;
        }
    }

    return 0;
}

// --------------------------------------------------------------------------
template<typename index_t, typename data_t>
int cellToPoint(index_t * __restrict__ simpleCells, index_t nSimpleCells,
    index_t nPts, index_t *& __restrict__ ptCount, data_t * __restrict__ cellScalar,
    data_t *& __restrict__ ptScalar)
{
    // allocate buffers
    if (!ptScalar)
        ptScalar = (data_t*)MALLOC(nPts*sizeof(data_t));

    // zero it out
    memset(ptScalar, 0, nPts*sizeof(data_t));

    // transfer values from cell to points.
    for (index_t i = 0; i < nSimpleCells; ++i)
    {
        index_t ii = 3*i;
        data_t val = cellScalar[i];
        for (index_t j = 0; j < 3; ++j)
        {
            index_t idx = simpleCells[ii+j];
            ptScalar[idx] += val;
        }
    }

    // finish the calculation by scalaing by the count.
    for (index_t i = 0; i < nPts; ++i)
    {
        ptScalar[i] /= ptCount[i];
    }

    return 0;
}

// --------------------------------------------------------------------------
template<typename index_t, typename coord_t, typename data_t>
int packageVtk(index_t *nc, coord_t *x0, coord_t *dx,
    int nCellScalars, data_t **cellScalar, const char **cellScalarName,
    vtkImageData *&imd)
{
    index_t nCells = nc[0]*nc[1]*nc[2];
    index_t ptDims[3] = {nc[0] + 1, nc[1] + 1, nc[2] + 1};

    imd = vtkImageData::New();
    imd->SetDimensions(ptDims[0], ptDims[1], ptDims[2]);
    imd->SetOrigin(x0[0], x0[1], x0[2]);
    imd->SetSpacing(dx[0], dx[1], dx[2]);

    for (int i = 0; i < nCellScalars; ++i)
    {
        vtkAOSDataArrayTemplate<data_t> *scalar =
            vtkAOSDataArrayTemplate<data_t>::New();
        scalar->SetName(cellScalarName[i]);
#if defined(ZERO_COPY_VTK)
        scalar->SetArray(dist, nCells, true);
#else
        scalar->SetNumberOfTuples(nCells);
        memcpy(scalar->GetPointer(0), cellScalar[i], nCells*sizeof(data_t));
#endif
        imd->GetCellData()->AddArray(scalar);
        scalar->Delete();
    }

    return 0;
}

// --------------------------------------------------------------------------
template<typename index_t, typename coord_t, typename data_t>
int packageVtk(int id, index_t nPts, coord_t * __restrict__ x,
    coord_t * __restrict__ y, coord_t * __restrict__ z,
    coord_t * __restrict__ d, index_t nCells, index_t * __restrict__ cells,
    index_t * __restrict__ cellLens, index_t * __restrict__ cellLocs,
    index_t cellArraySize, coord_t * __restrict__ dist,
    coord_t * __restrict__ thick, data_t * __restrict__ scalar, coord_t *pos,
    const char *e_type, const char *m_type, const char *ei, int layer,
    vtkPolyData *&dataset)
{
    // TODO -- how should these be packaged?
    (void)e_type;
    (void)m_type;
    (void)pos;

    // convert cells into VTK's format
    vtkIdTypeArray *cellIds = vtkIdTypeArray::New();
    cellIds->SetName("cellIds");
    cellIds->SetNumberOfTuples(cellArraySize + nCells);
    vtkIdType * __restrict__ dst = cellIds->GetPointer(0);
    for (index_t i = 0; i < nCells; ++i)
    {
        index_t * __restrict__ src = cells + cellLocs[i];
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
    coord_t * __restrict__ pcoords = coords->GetPointer(0);
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

    // gid
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

    // layer
    vtkIntArray *nla = vtkIntArray::New();
    nla->SetNumberOfTuples(nCells);
    nla->FillValue(layer);
    nla->SetName("layer");
    dataset->GetCellData()->AddArray(nla);
    nla->Delete();

    nla = vtkIntArray::New();
    nla->SetNumberOfTuples(nPts);
    nla->FillValue(layer);
    nla->SetName("layer");
    dataset->GetPointData()->AddArray(nla);
    nla->Delete();

    // e/i
    vtkIntArray *nei = vtkIntArray::New();
    nei->SetNumberOfTuples(nCells);
    nei->FillValue((ei[0] == 'e' ? 1 : 0));
    nei->SetName("excite");
    dataset->GetCellData()->AddArray(nei);
    nei->Delete();

    nei = vtkIntArray::New();
    nei->SetNumberOfTuples(nPts);
    nei->FillValue((ei[0] == 'e' ? 1 : 0));
    nei->SetName("excite");
    dataset->GetPointData()->AddArray(nei);
    nei->Delete();

    // distance from first point
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

    // thickness
    vtkAOSDataArrayTemplate<coord_t> *thickness =
        vtkAOSDataArrayTemplate<coord_t>::New();
    thickness->SetName("thickness");
#if defined(ZERO_COPY_VTK)
    thickness->SetArray(thick, nPts, true);
#else
    thickness->SetNumberOfTuples(nPts);
    memcpy(thickness->GetPointer(0), thick, nPts*sizeof(coord_t));
#endif
    dataset->GetPointData()->AddArray(thickness);
    thickness->Delete();

    // diameter
    vtkAOSDataArrayTemplate<coord_t> *diameter =
        vtkAOSDataArrayTemplate<coord_t>::New();
    diameter->SetName("diameter");
#if defined(ZERO_COPY_VTK)
    diameter->SetArray(d, nPts, true);
#else
    diameter->SetNumberOfTuples(nPts);
    memcpy(diameter->GetPointer(0), d, nPts*sizeof(coord_t));
#endif
    dataset->GetPointData()->AddArray(diameter);
    diameter->Delete();

    // scalar
    if (scalar)
    {
        vtkAOSDataArrayTemplate<data_t> *sa =
            vtkAOSDataArrayTemplate<data_t>::New();
        sa->SetName("scalar");
#if defined(ZERO_COPY_VTK)
        sa->SetArray(scalar, nPts, true);
#else
        sa->SetNumberOfTuples(nPts);
        memcpy(sa->GetPointer(0), scalar, nPts*sizeof(data_t));
#endif
        dataset->GetPointData()->AddArray(sa);
        sa->Delete();
    }

    return 0;
}

//--------------------------------------------------------------------------
template<typename index_t, typename coord_t>
int getBounds(index_t nPts, const coord_t * __restrict__ x,
    const coord_t * __restrict__ y, const coord_t * __restrict__ z,
    coord_t * __restrict__ bds)
{
    // find the spatial bounds of the data
    coord_t mnx = std::numeric_limits<coord_t>::max();
    for (index_t i = 0; i < nPts; ++i)
        mnx = mnx > x[i] ? x[i] : mnx;

    coord_t mxx = std::numeric_limits<coord_t>::lowest();
    for (index_t i = 0; i < nPts; ++i)
        mxx = mxx < x[i] ? x[i] : mxx;

    coord_t mny = std::numeric_limits<coord_t>::max();
    for (index_t i = 0; i < nPts; ++i)
        mny = mny > y[i] ? y[i] : mny;

    coord_t mxy = std::numeric_limits<coord_t>::lowest();
    for (index_t i = 0; i < nPts; ++i)
        mxy = mxy < y[i] ? y[i] : mxy;

    coord_t mnz = std::numeric_limits<coord_t>::max();
    for (index_t i = 0; i < nPts; ++i)
        mnz = mnz > z[i] ? z[i] : mnz;

    coord_t mxz = std::numeric_limits<coord_t>::lowest();
    for (index_t i = 0; i < nPts; ++i)
        mxz = mxz < z[i] ? z[i] : mxz;

    bds[0] = mnx;
    bds[1] = mxx;
    bds[2] = mny;
    bds[3] = mxy;
    bds[4] = mnz;
    bds[5] = mxz;

    return 0;
}

//--------------------------------------------------------------------------
template<typename coord_t, typename index_t>
int meshParams(coord_t * __restrict__ bounds, index_t nCellsMajorAxis,
    index_t * __restrict__ nc, coord_t * __restrict__ x0,
    coord_t * __restrict__ dx)
{
    // mesh resolution
    coord_t dxt[3] = {bounds[1] - bounds[0],
        bounds[3] - bounds[2], bounds[5] - bounds[4]};

    coord_t dxm = 0.0;

    if ((dxt[0] >= dxt[1]) && (dxt[0] >= dxt[2]))
    {
        // x axis major
        dxm = dxt[0]/(nCellsMajorAxis - 1);
        nc[0] = nCellsMajorAxis;
        nc[1] = dxt[1]/dxm + 1;
        nc[2] = dxt[2]/dxm + 1;
    }
    else if ((dxt[1] >= dxt[0]) && (dxt[1] >= dxt[2]))
    {
        // y axis major
        dxm = dxt[1]/(nCellsMajorAxis - 1);
        nc[0] = dxt[0]/dxm + 1;
        nc[1] = nCellsMajorAxis;
        nc[2] = dxt[2]/dxm + 1;
    }
    else
    {
        // z axis major
        dxm = dxt[2]/(nCellsMajorAxis - 1);
        nc[0] = dxt[0]/dxm + 1;
        nc[1] = dxt[1]/dxm + 1;
        nc[2] = nCellsMajorAxis;
    }

    // spacing
    dx[0] = dxm;
    dx[1] = dxm;
    dx[2] = dxm;

    // set the origin
    x0[0] = bounds[0]; //- dx[0]/2.0;
    x0[1] = bounds[2]; //- dx[1]/2.0;
    x0[2] = bounds[4]; //- dx[2]/2.0;

    return 0;
}

//--------------------------------------------------------------------------
template<typename index_t, typename coord_t>
void meshIds(
    // point set
    const coord_t * __restrict__ x, const coord_t * __restrict__ z, index_t nPts,
    // origin, spacing and local extent
    coord_t x0, coord_t z0, coord_t dx, coord_t dz, index_t nxc, index_t nzc,
    // mesh ids that each point falls in
    index_t * __restrict__ ids)
{
    (void)nzc;
    coord_t dxi = coord_t(1)/dx;
    coord_t dzi = coord_t(1)/dz;

    for (index_t q = 0; q < nPts; ++q)
    {
        ids[q] = static_cast<index_t>((z[q] - z0)*dzi)*nxc
            + static_cast<index_t>((x[q] - x0)*dxi);
    }
}

#include <cassert>
//--------------------------------------------------------------------------
template<typename index_t, typename coord_t>
void meshIds(
    // point set
    index_t nPts,
    const coord_t * __restrict__ x,
    const coord_t * __restrict__ y,
    const coord_t * __restrict__ z,
    // origin, spacing and local extent
    index_t nxc, index_t nyc, index_t nzc,
    coord_t x0, coord_t y0, coord_t z0,
    coord_t dx, coord_t dy, coord_t dz,
    // mesh ids that each point falls in
    index_t *& __restrict__ ids)
{
    (void)nzc;
    ids = (index_t*)MALLOC(nPts*sizeof(index_t));

    coord_t dxi = coord_t(1)/dx;
    coord_t dyi = coord_t(1)/dy;
    coord_t dzi = coord_t(1)/dz;
    index_t nxyc = nxc*nyc;

    for (index_t q = 0; q < nPts; ++q)
    {
        ids[q] = static_cast<index_t>((z[q] - z0)*dzi)*nxyc
            + static_cast<index_t>((y[q] - y0)*dyi)*nxc
            + static_cast<index_t>((x[q] - x0)*dxi);

        assert((ids[q] >= 0) && (ids[q] < nxc*nyc*nzc));
    }
}

template<typename index_t, typename data_t>
struct ReduceDensity
{
    ReduceDensity() : nc(0), result(nullptr) {}

    void initialize(index_t n)
    {
        if (nc < n)
        {
            nc = n;
            result = (data_t*)MALLOC(nc*sizeof(data_t));
        }
        memset(result, 0, nc*sizeof(data_t));
    }

    void reduce(index_t nPts, const index_t * __restrict__ ids,
        const data_t * __restrict__ scalar)
    {
        for (index_t q = 0; q < nPts; ++q)
            result[ids[q]] += data_t(1);
    }

    void finalize()
    {
    }

    data_t *getResult()
    {
        return result;
    }

    const char *getName()
    {
        return "density";
    }

    void freeMem()
    {
        free(result);
    }

    index_t nc;
    data_t *result;
};

template<typename index_t, typename data_t>
struct ReduceAverage
{
    ReduceAverage() : nc(0), result(nullptr), count(nullptr) {}

    void initialize(index_t n)
    {
        if (nc < n)
        {
            nc = n;
            result = (data_t*)MALLOC(nc*sizeof(data_t));
            count = (index_t*)MALLOC(nc*sizeof(index_t));
        }
        memset(result, 0, nc*sizeof(data_t));
        memset(count, 0, nc*sizeof(index_t));
    }

    void reduce(index_t nPts, const index_t * __restrict__ meshIds,
        const data_t * __restrict__ scalar)
    {
        for (index_t q = 0; q < nPts; ++q)
        {
            index_t qq = meshIds[q];
            result[qq] += scalar[q];
            count[qq] += index_t(1);
        }
    }

    void finalize()
    {
        for (index_t i = 0; i < nc; ++i)
            result[i] /= (count[i] ? data_t(count[i]) : data_t(1));
    }

    data_t *getResult()
    {
        return result;
    }

    const char *getName()
    {
        return "average";
    }

    void freeMem()
    {
        nc = 0;
        free(result);
        free(count);
    }

    index_t nc;
    data_t *result;
    index_t *count;
};


template<typename index_t, typename data_t>
struct ReduceMaxAbs
{
    ReduceMaxAbs() : nc(0), result(nullptr) {}

    void initialize(index_t n)
    {
        if (nc < n)
        {
            nc = n;
            result = (data_t*)MALLOC(nc*sizeof(data_t));
        }
        data_t ival = data_t(0);
        for (index_t q = 0; q < nc; ++q)
            result[q] = ival;
    }

    void reduce(index_t nPts, const index_t * __restrict__ meshIds,
        const data_t * __restrict__ scalar)
    {
        for (index_t q = 0; q < nPts; ++q)
        {
            index_t qq = meshIds[q];
            result[qq] = fabs(scalar[q]) > fabs(result[qq]) ?
                scalar[q] : result[qq];
        }
    }

    void finalize()
    {
    }

    data_t *getResult()
    {
        return result;
    }

    const char *getName()
    {
        return "maxabs";
    }

    void freeMem()
    {
        nc = 0;
        free(result);
    }

    index_t nc;
    data_t *result;
};


template<typename index_t, typename data_t>
struct ReduceMax
{
    ReduceMax() : nc(0), result(nullptr) {}

    void initialize(index_t n)
    {
        if (nc < n)
        {
            nc = n;
            result = (data_t*)MALLOC(nc*sizeof(data_t));
        }
        data_t ival = std::numeric_limits<data_t>::lowest();
        for (index_t q = 0; q < nc; ++q)
            result[q] = ival;
    }

    void reduce(index_t nPts, const index_t * __restrict__ meshIds,
        const data_t * __restrict__ scalar)
    {
        for (index_t q = 0; q < nPts; ++q)
        {
            index_t qq = meshIds[q];
            result[qq] = scalar[q] > result[qq] ?  scalar[q] : result[qq];
        }
    }

    void finalize()
    {
        data_t ival = 10.0*std::numeric_limits<data_t>::lowest();
        for (index_t q = 0; q < nc; ++q)
            result[q] = result[q] <= ival ? 0.0 : result[q];
    }

    data_t *getResult()
    {
        return result;
    }

    const char *getName()
    {
        return "max";
    }


    void freeMem()
    {
        nc = 0;
        free(result);
    }

    index_t nc;
    data_t *result;
};




// helper class that does the work of importing a single neuron
// from an hdf5 file named [N].h5, where N is the global id of
// the neuron
template<typename index_t, typename coord_t, typename data_t>
struct Neuron
{
    Neuron() :
        blockId(-1), neuronId(-1), baseDir(nullptr), nSimpleCells(0),
        p0(nullptr), p5(nullptr), p1(nullptr), d0(nullptr), d5(nullptr),
        d1(nullptr), spos{-1,-1,-1}, e_type(nullptr), m_type(nullptr),
        ei(nullptr), layer(-1), nPts(0), x(nullptr), y(nullptr), z(nullptr),
        d(nullptr), simpleCells(nullptr), nComplexCells(0),
        complexCells(nullptr), complexCellLens(nullptr),
        complexCellLocs(nullptr), complexCellArraySize(0), dist(nullptr),
        thick(nullptr), scalar(nullptr), count(nullptr), meshIds(nullptr)
    {}

    Neuron(int bid, int nid, const char *dir) :
        blockId(bid), neuronId(nid), baseDir(dir), nSimpleCells(0),
        p0(nullptr), p5(nullptr), p1(nullptr), d0(nullptr), d5(nullptr),
        d1(nullptr), spos{-1,-1,-1}, e_type(nullptr), m_type(nullptr),
        ei(nullptr), layer(-1), nPts(0), x(nullptr), y(nullptr), z(nullptr),
        d(nullptr), simpleCells(nullptr), nComplexCells(0),
        complexCells(nullptr), complexCellLens(nullptr),
        complexCellLocs(nullptr), complexCellArraySize(0), dist(nullptr),
        thick(nullptr), scalar(nullptr), count(nullptr), meshIds(nullptr)
    {}

    ~Neuron() {}

    Neuron(const Neuron &) = delete;
    Neuron(Neuron &&) = default;

    Neuron &operator=(const Neuron &) = delete;
    Neuron &operator=(Neuron &&) = default;

    int import()
    {
        // read the data for this instance
        if (readNeuron(baseDir, neuronId, p0, p5, p1,
            d0, d5, d1, spos, e_type, m_type, ei, layer,
            nSimpleCells))
        {
            ERROR("Failed to load neuron " << neuronId)
            return -1;
        }

        // create the points and simple cells
        initialize(nSimpleCells, p0, p5, p1, d0, d5, d1,
            nPts, x, y, z, d, simpleCells);

        // remove duplicate points and merge segements
        clean(nPts, x, y, z, d, nSimpleCells, simpleCells,
            nComplexCells, complexCells, complexCellLens,
            complexCellLocs, complexCellArraySize);

        // get the spatial bounds
        getBounds(nPts, x, y, z, bounds);

        // compute distance field
        distance(complexCells, nPts, x, y, z, dist);

        // thickenss as a function of distance
        thickness(nComplexCells, complexCells, complexCellLens,
            complexCellLocs, nPts, dist, thick);

        // generate counts
        pointCount(simpleCells, nSimpleCells, nPts, count);

        return 0;
    }

    int setMeshParams(index_t *nc, coord_t *x0, coord_t *dx)
    {
        import::meshIds(nPts, x, y, z, nc[0], nc[1], nc[2],
            x0[0], x0[1], x0[2], dx[0], dx[1], dx[2], meshIds);
        return 0;
    }

    index_t *getMeshIds()
    {
        return meshIds;
    }

    int setScalar(data_t * __restrict__ cellScalar)
    {
        // interpolate from simple cells to reduced point set
        return cellToPoint(simpleCells, nSimpleCells, nPts,
            count, cellScalar, scalar);
    }

    void freeMem()
    {
        free(p0);
        free(p5);
        free(p1);
        free(d0);
        free(d5);
        free(d1);
        free(e_type);
        free(m_type);
        free(ei);
        free(x);
        free(y);
        free(z);
        free(d);
        free(simpleCells);
        free(complexCells);
        free(complexCellLens);
        free(complexCellLocs);
        free(dist);
        free(thick);
        free(scalar);
        free(count);
        free(meshIds);
        nSimpleCells = 0;
        p0 = nullptr;
        p5 = nullptr;
        p1 = nullptr;
        d0 = nullptr;
        d5 = nullptr;
        d1 = nullptr;
        memset(spos, -1, sizeof(spos));
        e_type = nullptr;
        m_type = nullptr;
        ei = nullptr;
        layer = -1;
        nPts = 0;
        x = nullptr;
        y = nullptr;
        z = nullptr;
        d = nullptr;
        simpleCells = nullptr;
        nComplexCells = 0;
        complexCells = nullptr;
        complexCellLens = nullptr;
        complexCellLocs = nullptr;
        complexCellArraySize = 0;
        dist = nullptr;
        thick = nullptr;
        scalar = nullptr;
        count = nullptr;
        meshIds = nullptr;
    }

    void package(vtkPolyData *&block)
    {
        // package in VTK format
        packageVtk(neuronId, nPts, x, y, z, d, nComplexCells, complexCells,
            complexCellLens, complexCellLocs, complexCellArraySize,
            dist, thick, scalar, spos, e_type, m_type, ei, layer, block);
    }

    int blockId;
    int neuronId;
    const char *baseDir;
    index_t nSimpleCells;
    coord_t *p0;
    coord_t *p5;
    coord_t *p1;
    coord_t *d0;
    coord_t *d5;
    coord_t *d1;
    coord_t spos[3];
    char *e_type;
    char *m_type;
    char *ei;
    int layer;
    index_t nPts;
    coord_t *x;
    coord_t *y;
    coord_t *z;
    coord_t *d;
    coord_t bounds[6];
    index_t *simpleCells;
    index_t nComplexCells;
    index_t *complexCells;
    index_t *complexCellLens;
    index_t *complexCellLocs;
    index_t complexCellArraySize;
    coord_t *dist;
    coord_t *thick;
    data_t *scalar;
    index_t *count;
    index_t *meshIds;
};

template<typename index_t, typename coord_t, typename data_t>
struct Importer
{
    Importer() = delete;

    Importer(int n0, const char *dir,
        std::vector<Neuron<index_t,coord_t,data_t>> *n) :
            neuron0(n0), baseDir(dir), neurons(n)
    {}

    void operator()(const tbb::blocked_range<int> &r) const
    {
        for (int i = r.begin(); i != r.end(); ++i)
        {
            Neuron<index_t,coord_t,data_t> import(i, neuron0 + i, baseDir);
            import.import();
            (*neurons)[i] = std::move(import);
#if defined(MEM_CHECK)
            std::cerr << ".";
#endif
        }
    }

    int neuron0;
    const char *baseDir;
    std::vector<Neuron<index_t,coord_t,data_t>> *neurons;
};


template<typename index_t>
struct lessThan
{
    lessThan(index_t *d) : data(d) {}
    bool operator()(index_t i, index_t j)
    {
        return data[i] < data[j];
    }
    index_t *data;
};

template <typename index_t, typename coord_t, typename data_t>
struct TimeSeries
{
    TimeSeries() : baseDir(nullptr), fh(-1), t0(-1.), t1(-1.), dt(-1.),
        nNeurons(0), nSteps(0), stepSize(0), neuronIds(nullptr),
        neuronOffs(nullptr), neuronSize(nullptr), raMap(nullptr),
        scalar(nullptr), neurons(nullptr)
    {}

    TimeSeries(const char *baseDir,
        std::vector<Neuron<index_t,coord_t,data_t>> *nrns) :
        baseDir(baseDir), fh(-1), t0(-1.), t1(-1.), dt(-1.),
        nNeurons(0), nSteps(0), stepSize(0), neuronIds(nullptr),
        neuronOffs(nullptr), neuronSize(nullptr), raMap(nullptr),
        scalar(nullptr), neurons(nrns)
    {}

    ~TimeSeries() {}

    TimeSeries(const TimeSeries &) = delete;
    TimeSeries(TimeSeries &&) = default;

    TimeSeries &operator=(const TimeSeries &) = delete;
    TimeSeries &operator=(TimeSeries &&) = default;

    int initialize()
    {
        if (readTimeSeriesMetadata(baseDir, fh, nNeurons, nSteps,
            stepSize, t0, t1, dt, neuronIds, neuronOffs, neuronSize))
        {
            ERROR("Failed to import time series metadata")
            return -1;
        }

        // make random access structure. neuronIds are not in
        // ascending order, so given a neuronId the map tells you
        // where in the offset and size array its data lives
        raMap = (index_t*)MALLOC(nNeurons*sizeof(index_t));
        for (index_t i = 0; i < nNeurons; ++i)
            raMap[i] = i;

        std::sort(raMap, raMap+nNeurons, lessThan(neuronIds));

        return 0;
    }

    int importTimeStep(index_t step)
    {
        // read this time step
        if (readTimeStep(fh, nSteps, stepSize, step, scalar))
        {
            ERROR("Failed to read time step " << step)
            return -1;
        }

        // update the neurons
        index_t nnrn = neurons->size();
        for (index_t i = 0; i < nnrn; ++i)
        {
            import::Neuron<index_t, coord_t, data_t> &ni = (*neurons)[i];
            ni.setScalar(getScalar(ni.neuronId));
        }

        return 0;
    }

    data_t *getScalar(index_t nid)
    {
        // return a pointer to the start of this neurons data
        return scalar + neuronOffs[raMap[nid]];
    }

    void freeMem()
    {
        free(neuronIds);
        free(neuronOffs);
        free(neuronSize);
        free(raMap);
        free(scalar);
        H5Fclose(fh);
        fh = -1;
    }

    void print()
    {
        std::cerr << "nNeurons=" << nNeurons << ", nSteps=" << nSteps
            << ", stepSize=" << stepSize << " t0=" << t0
            << ", t1=" << t1 << ", dt=" <<  dt << std::endl;
#if defined(MEM_CHECK)
        std::cerr << "access structs (" << nNeurons << ") = ";
        for (index_t i = 0; i < nNeurons; ++i)
        {
            index_t ii = raMap[i];
            std::cerr <<"(" <<  neuronIds[ii] << ", "
                << neuronOffs[ii] << ", " << neuronSize[ii] << ") ";
        }
        std::cerr << std::endl;
#endif
    }

    const char *baseDir;
    hid_t fh;
    double t0;
    double t1;
    double dt;
    index_t nNeurons;
    index_t nSteps;
    index_t stepSize;
    index_t *neuronIds;
    index_t *neuronOffs;
    index_t *neuronSize;
    index_t *raMap;
    data_t *scalar;
    std::vector<Neuron<index_t,coord_t,data_t>> *neurons;
};


template<typename index_t, typename coord_t, typename data_t>
struct Mesher
{
    Mesher() : nCells(0), bounds{0.}, x0{0.},
        dx{0.}, nc{0}, neurons(nullptr)
    {}

    Mesher(std::vector<Neuron<index_t,coord_t,data_t>> *nrns) :
            nCells(0), bounds{0.}, x0{0.}, dx{0.}, nc{0}, neurons(nrns)
    {}

    ~Mesher() {}

    Mesher(const Mesher &) = delete;
    Mesher(Mesher &&) = default;

    Mesher &operator=(const Mesher &) = delete;
    Mesher &operator=(Mesher &&) = default;

    void freeMem()
    {
    }

    void print()
    {
        std::cerr << "bounds=" << bounds[0] << ", " << bounds[1]
            << ", " << bounds[2] << ", " << bounds[3] << ", "
            << bounds[4] << ", " << bounds[5] << std::endl
            << "x0=" << x0[0] << ", " << x0[1] << ", " << x0[2] << std::endl
            << "x1=" << x0[0] + dx[0]*(nc[0]+1) << ", " << x0[1] + dx[1]*(nc[1]+1)
            << ", " << x0[2] + dx[2]*(nc[2]+1) << std::endl
            << "dx=" << dx[0] << ", " << dx[1] << ", " << dx[2] << std::endl
            << "nc=" << nc[0] << ", " << nc[1] << ", " << nc[2] << std::endl;
    }

    void getMeshParams(index_t *anc, coord_t *ax0, coord_t *adx, index_t &anct)
    {
        anct = 1;
        for (int i = 0; i < 3; ++i)
        {
            anc[i] = nc[i];
            ax0[i] = x0[i];
            adx[i] = dx[i];
            anct *= nc[i];
        }
    }

    int initialize(index_t longSideCells)
    {
        // calculate bounding box
        for (int i = 0; i < 3; ++i)
        {
            bounds[2*i  ] = std::numeric_limits<coord_t>::max();
            bounds[2*i+1] = std::numeric_limits<coord_t>::lowest();
        }

        index_t nNeuron = neurons->size();
        for (int i = 0; i < nNeuron; ++i)
        {
            const import::Neuron<index_t, coord_t, data_t> &ni = (*neurons)[i];
            bounds[0] = std::min(bounds[0], ni.bounds[0]);
            bounds[1] = std::max(bounds[1], ni.bounds[1]);
            bounds[2] = std::min(bounds[2], ni.bounds[2]);
            bounds[3] = std::max(bounds[3], ni.bounds[3]);
            bounds[4] = std::min(bounds[4], ni.bounds[4]);
            bounds[5] = std::max(bounds[5], ni.bounds[5]);
        }

        // derive mesh parameters
        import::meshParams(bounds, longSideCells, nc, x0, dx);

        nCells = nc[0]*nc[1]*nc[2];

        // initialize neurons for meshing, this generates mesh ids
        // needed for reductions
        //for (int i = 0; i < nNeuron; ++i)
        tbb::parallel_for(size_t(0), size_t(nNeuron), size_t(1),
            [&](size_t i)
            {
                (*neurons)[i].setMeshParams(nc, x0, dx);
            });

        return 0;
    }

    template<template<typename,typename> typename op_t>
    void mesh(op_t<index_t, coord_t> &op)
    {
        op.initialize(nCells);

        index_t nNeuron = neurons->size();
        for (int i = 0; i < nNeuron; ++i)
        {
            const Neuron<index_t,coord_t,data_t> &ni = (*neurons)[i];
            op.reduce(ni.nPts, ni.meshIds, ni.scalar);
        }

       op.finalize();
    }

    index_t nCells;
    coord_t bounds[6];
    coord_t x0[3];
    coord_t dx[3];
    index_t nc[3];
    std::vector<Neuron<index_t,coord_t,data_t>> *neurons;
};

template <typename index_t, typename coord_t, typename data_t>
struct Exporter
{
    Exporter() : baseDir("./"),
        dsetName("neurons"), nct(0), nc{0}, x0{0.}, dx{0.}
    {}

    template <typename mesher_t>
    int initialize(const char *bdir, const char *dset, mesher_t &mesher)
    {
        baseDir = bdir;
        dsetName = dset;
        mesher.getMeshParams(nc, x0, dx, nct);

        // make the output directory
        DIR* de = opendir(baseDir);
        if (de)
        {
            // all set
            closedir(de);
        }
        else if (errno == ENOENT)
        {
            // doesn't yet exist
            if (mkdir(baseDir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH))
            {
                ERROR("failed to create the output directory \"" << baseDir << "\"")
                return -1;
            }
        }
        else
        {
            // other error
            const char *estr = strerror(errno);
            ERROR("failed to create output directory \"" << baseDir << "\" " << estr)
            return -1;
        }

        return 0;
    }

    int write(index_t step, index_t nArrays, data_t **arrays, const char **names)
    {
        // package as VTK
        vtkImageData *mesh = nullptr;
        packageVtk(nc, x0, dx, nArrays, arrays, names, mesh);

        vtkXMLImageDataWriter *writer = vtkXMLImageDataWriter::New();

        // write VTK
        char outFileName[1024];
        snprintf(outFileName, 1023, "%s/%s_mesh_%06d.%s", baseDir,
            dsetName, int(step), writer->GetDefaultFileExtension());

        writer->SetInputData(mesh);
        writer->SetFileName(outFileName);

        int ierr = 0;
        writer->SetFileName(outFileName);
        if (writer->Write() == 0)
        {
            ERROR("vtk writer failed to write mesh \"" << outFileName << "\"")
            ierr = -1;
        }

        writer->Delete();
        mesh->Delete();

        return ierr;
    }

    int write(index_t step,
        std::vector<Neuron<index_t,coord_t,data_t>> *neurons)
    {
        index_t nNeuron = neurons->size();

        // package as VTK
        vtkMultiBlockDataSet *mbds = vtkMultiBlockDataSet::New();
        mbds->SetNumberOfBlocks(nNeuron);

        //for (int i = 0; i < nNeuron; ++i)
        tbb::parallel_for(size_t(0), size_t(nNeuron), size_t(1),
            [&](size_t i)
            {
                vtkPolyData *block = nullptr;
                (*neurons)[i].package(block);
                mbds->SetBlock(i, block);
                block->Delete();
            });

        // write VTK
        vtkXMLMultiBlockDataWriter *writer = vtkXMLMultiBlockDataWriter::New();
        writer->SetInputData(mbds);

        char outFileName[1024];
        snprintf(outFileName, 1023, "%s/%s_neuron_%06d.%s", baseDir,
            dsetName, int(step), writer->GetDefaultFileExtension());

        int ierr = 0;
        writer->SetFileName(outFileName);
        if (writer->Write() == 0)
        {
            ERROR("vtk writer failed to write neurons \"" << outFileName << "\"")
            ierr = -1;
        }

        mbds->Delete();
        writer->Delete();

        return ierr;
    }

    void freeMem()
    {
    }

    const char *baseDir;
    const char *dsetName;
    index_t nct;
    index_t nc[3];
    coord_t x0[3];
    coord_t dx[3];
};

};

#endif
