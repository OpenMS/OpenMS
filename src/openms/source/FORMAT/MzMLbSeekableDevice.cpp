//
// $Id: MzMLbSeekableDevice.cpp
//
//
// Original authors: Andrew Dowsey <andrew.dowsey@bristol.ac.uk>
//
// Copyright 2017 biospi Laboratory,
//                University of Bristol, UK
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#include <OpenMS/FORMAT/MzMLbSeekableDevice.h>
#include <stdexcept>
#include <iostream>
#include <sstream>


//using namespace std;
using namespace boost::iostreams;

#define CURRENT_VERSION "mzMLb 1.0"

namespace OpenMS {
/*
namespace pwiz {
namespace msdata {
namespace mzmlb {
*/
MzMLbSeekableDevice::MzMLbSeekableDevice(const std::string& id, bool identifyOnly)
{
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    // open HDF5 file for reading
    file_ = H5Fopen(id.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_ < 0)
        throw std::runtime_error("[MzMLbSeekableDevice::open()] Could not open mzMLb file for reading.");

    // open dataset with stored mzML XML to find chunk size
    hsize_t chunk_size;
    hid_t dataset = H5Dopen(file_, "mzML", H5P_DEFAULT);
    {
        if (dataset < 0)
        {
            H5Fclose(file_);
            throw std::runtime_error("[MzMLbSeekableDevice::open()] Could not open mzML dataset for reading.");
        }
        
        hid_t dcpl = H5Dget_create_plist(dataset);
        {
            H5Pget_chunk(dcpl, 1, &chunk_size);
        }
        H5Pclose(dcpl);
    }
    H5Dclose(dataset);
    
    // open again with appropraite dataset cache
    hid_t dapl = H5Pcreate(H5P_DATASET_ACCESS);
    {
        size_t nslots; // Number of slots in the hash table
        size_t nbytes; // Size of chunk cache in bytes
        double w0; // Chunk preemption policy
        H5Pget_chunk_cache(dapl, &nslots, &nbytes, &w0);
        nbytes = (nbytes > chunk_size) ? nbytes : chunk_size;
        w0 = 1.0; // since pwiz only writes a spectrum once
        H5Pset_chunk_cache(dapl, nslots, nbytes, w0);
        mzML_.dataset = H5Dopen(file_, "mzML", dapl); // open XML part
        mzML_.space = H5Dget_space(mzML_.dataset); // get handle for dataset retrieval
        hsize_t size, maxdims;
        H5Sget_simple_extent_dims(mzML_.space, &size, &maxdims); // retrieve size and number of dimensions (TODO: check if more than one dimension)
        mzML_.size = size; // size of XML part

        // get version number
        hid_t aid = H5Aopen(mzML_.dataset, "version", H5P_DEFAULT);
        {
            if (aid < 0)
            {
                H5Aclose(aid);
                close();
                throw std::runtime_error("[MzMLbSeekableDevice::open()] This does not look like an mzMLb file.");
            }
            
            hid_t atype = H5Aget_type(aid);
            // at this point it's definitely an mzMLb file and if version is wrong it should not throw an exception that would cause Reader::identify to not identify the file
            if (!identifyOnly)
            {
                hid_t atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);
                {
                    char* ver = new char[H5Tget_size(atype)];
                    H5Aread(aid, atype_mem, ver);
                    std::string version = ver;
                    std::cout << "mzMLb version: " << version << std::endl;
                    if (version != CURRENT_VERSION)
                    {
                        H5Aclose(aid);
                        H5Aclose(atype_mem);
                        close();
                        throw std::runtime_error("[MzMLbSeekableDevice::open()] Cannot read this version of mzMLb: \"" + version + "\" (or version is not fixed-length string); only " CURRENT_VERSION " is supported");
                    }                   
                }
                H5Aclose(atype_mem);
            }
            H5Aclose(atype);
        }
        H5Aclose(aid);        
    }
    H5Pclose(dapl);

    if (identifyOnly)
    {
        close();
        return;
    }

    opaque_id_ = H5Tcreate(H5T_OPAQUE, 1);
 }


MzMLbSeekableDevice::MzMLbSeekableDevice(const std::string& id, int chunk_size, int compression_level) :
    chunk_size_(chunk_size),
    compression_level_(compression_level)
{    
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    // create/truncate HDF5 file for writing
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    {
        int nelemts; // Dummy parameter in API, no longer used
        size_t nslots; // Number of slots in the hash table
        size_t nbytes; // Size of chunk cache in bytes
        double w0; // Chunk preemption policy
        H5Pget_cache(fapl, &nelemts, &nslots, &nbytes, &w0);
        nbytes = (nbytes > chunk_size) ? nbytes : chunk_size; // Set per dataset cache to twice the chunk size please
        w0 = 1.0; // since pwiz only writes a spectrum once
        H5Pset_cache(fapl, nelemts, nslots, nbytes, w0);
        file_ = H5Fcreate(id.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
        if (file_ < 0)
        {
            H5Pclose(fapl);
            throw std::runtime_error("[MzMLbSeekableDevice::MzMLbSeekableDevice()] Could not open or create mzMLb file for writing.");
        }

        // create dataset to store mzML XML
        hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        {        
            hsize_t cdims = chunk_size;
            H5Pset_chunk(dcpl, 1, &cdims);
            if (compression_level > 0)
            {
                hsize_t level = compression_level;
                H5Pset_deflate(dcpl, level);
            }
            H5Pset_fletcher32(dcpl);
            hsize_t maxdims = H5S_UNLIMITED;
            mzML_.space = H5Screate_simple(1, &mzML_.size, &maxdims);
            mzML_.dataset = H5Dcreate(file_, "mzML", H5T_NATIVE_CHAR, mzML_.space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
        }
        H5Pclose(dcpl);

        // write version string
        hid_t aid = H5Screate(H5S_SCALAR);
        {
            hid_t atype = H5Tcopy(H5T_C_S1);
            {
                H5Tset_size(atype, 10);
                H5Tset_strpad(atype, H5T_STR_NULLTERM);
                hid_t attrid = H5Acreate(mzML_.dataset, "version", atype, aid, H5P_DEFAULT, H5P_DEFAULT);
                {
                    std::string version_out = CURRENT_VERSION;
                    H5Awrite(attrid, atype, version_out.c_str());
                }
                H5Aclose(attrid);
            }
            H5Tclose(atype);
        }
        H5Sclose(aid);
    }
    H5Pclose(fapl);
    
    opaque_id_ = H5Tcreate(H5T_OPAQUE, 1);
}


void MzMLbSeekableDevice::close()
{
    H5Tclose(opaque_id_);

    H5Dclose(mzML_.dataset);
    H5Sclose(mzML_.space);

    for (std::map<std::string, Stream>::iterator it = binary_.begin(); it != binary_.end(); ++it)
    {
        H5Dclose(it->second.dataset);
        H5Sclose(it->second.space);
    }

    H5Fclose(file_);
}

// read mzMLb "mzML" dataset
std::streamsize MzMLbSeekableDevice::read(char* s, std::streamsize n)
{
    // don't read past end of dataset
    if (mzML_.pos + n > mzML_.size)
    {
        n = mzML_.size - mzML_.pos;
    }

    if (n > 0)
    {
        // read
        hsize_t count = n;
        H5Sselect_hyperslab(mzML_.space, H5S_SELECT_SET, &mzML_.pos, NULL, &count, NULL);
        hid_t mspace = H5Screate_simple(1, &count, &count);
        {
            H5Dread(mzML_.dataset, H5T_NATIVE_CHAR, mspace, mzML_.space, H5P_DEFAULT, s);
        }
        H5Sclose(mspace);

        mzML_.pos += n;

        return n;
    }
    else
    {
        return -1;
    }
}


// write mzMLb "mzML" dataset
std::streamsize MzMLbSeekableDevice::write(const char* s, std::streamsize n)
{
    // extend dataset size if needed
    if (mzML_.pos + n > mzML_.size)
    {
        mzML_.size = mzML_.pos + n;
        H5Dset_extent(mzML_.dataset, &mzML_.size);
        H5Sclose(mzML_.space);
        mzML_.space = H5Dget_space(mzML_.dataset);
    }

    // write
    hsize_t count = n;
    H5Sselect_hyperslab(mzML_.space, H5S_SELECT_SET, &mzML_.pos, NULL, &count, NULL);
    hid_t mspace = H5Screate_simple(1, &count, &count);
    {
        H5Dwrite(mzML_.dataset, H5T_NATIVE_CHAR, mspace, mzML_.space, H5P_DEFAULT, s);
    }
    H5Sclose(mspace);
    
    mzML_.pos += n;

    return n;
}


// seek mzMLb "mzML" dataset
stream_offset MzMLbSeekableDevice::seek(stream_offset off, std::ios_base::seekdir way)
{   
    switch (way)
    {
    case std::ios_base::beg:
        mzML_.pos = off;
        break;
    case std::ios_base::cur:
        mzML_.pos += off;
        break;
    case std::ios_base::end:
    default:
        mzML_.pos = mzML_.size - off;
        break;
    }
    
 
    return mzML_.pos;
}


// read mzMLb mzML index
/*void MzMLbSeekableDevice::readIndex(const std::string& id, std::vector<stream_offset>& positions)
{
    // open dataset to find chunk size
    hsize_t chunk_size;
    hid_t dataset = H5Dopen(file_, id.c_str(), H5P_DEFAULT);
    {
        if (dataset < 0)
        {
            throw std::runtime_error("[MzMLbSeekableDevice::read()] Could not open dataset " + id + " for reading.");
        }
        hid_t dcpl = H5Dget_create_plist(dataset);
        {
            H5Pget_chunk(dcpl, 1, &chunk_size);
        }
        H5Pclose(dcpl);
    }
    H5Dclose(dataset);
   
    // open dataset
    hid_t dapl = H5Pcreate(H5P_DATASET_ACCESS);
    {
        size_t nslots; // Number of slots in the hash table
        size_t nbytes; // Size of chunk cache in bytes
        double w0; // Chunk preemption policy
        H5Pget_chunk_cache(dapl, &nslots, &nbytes, &w0);
        nbytes = (nbytes > chunk_size) ? nbytes : chunk_size; // Set per dataset cache to twice the chunk size please
        w0 = 1.0; // since pwiz only writes a spectrum once
        H5Pset_chunk_cache(dapl, nslots, nbytes, w0);
        hid_t dataset = H5Dopen(file_, id.c_str(), dapl);   
        {
            hid_t space = H5Dget_space(dataset);
            {
                hsize_t size, maxdims;
                H5Sget_simple_extent_dims(space, &size, &maxdims);
        
                // read
                hid_t mspace = H5Screate_simple(1, &size, &size);
                {
                    std::vector<long long> positions_(size);
                    H5Dread(dataset, H5T_NATIVE_LLONG, mspace, space, H5P_DEFAULT, &positions_[0]);
                    positions.assign(positions_.begin(), positions_.end());
                }
                H5Sclose(mspace);
            }
            H5Sclose(space);
        }
        H5Dclose(dataset);
    }
    H5Pclose(dapl);
}


// write mzMLb mzML index
void MzMLbSeekableDevice::writeIndex(const std::string& id, const std::vector<stream_offset>& positions)
{
    // create dataset
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    {
        hsize_t cdims = chunk_size_;
        H5Pset_chunk(dcpl, 1, &cdims);
        if (compression_level_ > 0)
        {
            hsize_t level = compression_level_;
            H5Pset_shuffle(dcpl);
            H5Pset_deflate(dcpl, level);
        }
        H5Pset_fletcher32(dcpl);
        hsize_t maxdims = H5S_UNLIMITED;
        hsize_t count = positions.size();
        hid_t space = H5Screate_simple(1, &count, &maxdims);
        {
            hid_t dataset = H5Dcreate(file_, id.c_str(), H5T_NATIVE_LLONG, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
            { 
                // write
                std::vector<long long> positions_(positions.begin(), positions.end());
                hid_t mspace = H5Screate_simple(1, &count, &count);
                {
                    H5Dwrite(dataset, H5T_NATIVE_LLONG, mspace, space, H5P_DEFAULT, &positions_[0]);
                }
                H5Sclose(mspace);
            }
            H5Dclose(dataset);
        }
        H5Sclose(space);
    }
    H5Pclose(dcpl);
 }*/


bool MzMLbSeekableDevice::exists(const std::string& id)
{
    return H5Lexists(file_, id.c_str(), H5P_DEFAULT) > 0;
}


std::streamsize MzMLbSeekableDevice::size(const std::string& id)
{
    hid_t dataset = H5Dopen(file_, id.c_str(), H5P_DEFAULT);

    if (dataset < 0)
        throw std::runtime_error("[MzMLbSeekableDevice::read()] Could not open dataset " + id + ".");
    
    hid_t space = H5Dget_space(dataset);
    
    hsize_t size;
    H5Sget_simple_extent_dims(space, &size, 0);
    
    H5Sclose(space);
    H5Dclose(dataset);

    return size;
}


std::streamsize MzMLbSeekableDevice::read_opaque(const std::string& id, void* buf, std::streamsize n)
{
    return read(id, buf, n, opaque_id_);
}


std::streamsize MzMLbSeekableDevice::read(const std::string& id, char* buf, std::streamsize n)
{
    return read(id, buf, n, H5T_NATIVE_CHAR);
}


std::streamsize MzMLbSeekableDevice::read(const std::string& id, double* buf, std::streamsize n)
{
    return read(id, buf, n, H5T_NATIVE_DOUBLE);
}


std::streamsize MzMLbSeekableDevice::read(const std::string& id,long* buf, std::streamsize n)
{
    return read(id, buf, n, H5T_NATIVE_LONG);
}


std::streamsize MzMLbSeekableDevice::read(const std::string& id, long long* buf, std::streamsize n)
{
    return read(id, buf, n, H5T_NATIVE_LLONG);
}


std::streamsize MzMLbSeekableDevice::read(const std::string& id, void* buf, std::streamsize n, hid_t native_format)
{
    Stream& s_ = binary_[id];
    if (!s_.dataset)
    {
        // open dataset to find chunk size
        hsize_t chunk_size;
        hid_t dataset = H5Dopen(file_, id.c_str(), H5P_DEFAULT);
        {
            if (dataset < 0)
            {
                throw std::runtime_error("[MzMLbSeekableDevice::read()] Could not open dataset " + id + " for reading.");
            }
            hid_t dcpl = H5Dget_create_plist(dataset);
            {
                H5Pget_chunk(dcpl, 1, &chunk_size);
            }
            H5Pclose(dcpl);
        }
        H5Dclose(dataset);
       
        // open dataset
        hid_t dapl = H5Pcreate(H5P_DATASET_ACCESS);
        {
            size_t nslots; // Number of slots in the hash table
            size_t nbytes; // Size of chunk cache in bytes
            double w0; // Chunk preemption policy
            H5Pget_chunk_cache(dapl, &nslots, &nbytes, &w0);
             nbytes = (nbytes > chunk_size) ? nbytes : chunk_size; // Set per dataset cache to twice the chunk size please
            w0 = 1.0; // since pwiz only writes a spectrum once
            H5Pset_chunk_cache(dapl, nslots, nbytes, w0);
            s_.dataset = H5Dopen(file_, id.c_str(), dapl);        
            s_.space = H5Dget_space(s_.dataset);
            hsize_t size;
            H5Sget_simple_extent_dims(s_.space, &size, 0);
            s_.size = size;
        }
        H5Pclose(dapl);
    }

    // don't read past end of dataset
    if (s_.pos + n > s_.size)
    {
        n = s_.size - s_.pos;
    }

    if (n > 0)
    {
        // read
        hsize_t count = n;
        H5Sselect_hyperslab(s_.space, H5S_SELECT_SET, &s_.pos, NULL, &count, NULL);
        hid_t mspace = H5Screate_simple(1, &count, &count);
        {
            H5Dread(s_.dataset, native_format, mspace, s_.space, H5P_DEFAULT, buf);
        }
        H5Sclose(mspace);
        
        s_.pos += n;

        return n;
    }
    else
    {
       return -1;
    }
}


std::streamsize MzMLbSeekableDevice::write_opaque(const std::string& id, const void* buf, std::streamsize n)
{
    return write(id, buf, n, opaque_id_, opaque_id_, 1);
}


std::streamsize MzMLbSeekableDevice::write(const std::string& id, const char* buf, std::streamsize n)
{
    return write(id, buf, n, H5T_NATIVE_CHAR, H5T_NATIVE_CHAR, sizeof(char));
}


std::streamsize MzMLbSeekableDevice::write(const std::string& id, const float* buf, std::streamsize n)
{
    return write(id, buf, n, H5T_NATIVE_FLOAT, H5T_NATIVE_FLOAT, sizeof(float));
}


std::streamsize MzMLbSeekableDevice::write(const std::string& id, const double* buf, std::streamsize n)
{
    return write(id, buf, n, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, sizeof(double));
}


std::streamsize MzMLbSeekableDevice::write(const std::string& id, const long* buf, std::streamsize n)
{
    return write(id, buf, n, H5T_NATIVE_LONG, H5T_NATIVE_LONG, sizeof(long));
}


std::streamsize MzMLbSeekableDevice::write(const std::string& id, const long long* buf, std::streamsize n)
{
    return write(id, buf, n, H5T_NATIVE_LLONG, H5T_NATIVE_LLONG, sizeof(long long));
}


std::streamsize MzMLbSeekableDevice::write(const std::string& id, const void* buf, std::streamsize n, hid_t native_format, hid_t format, size_t bytes)
{
    Stream& stream = binary_[id];
    if (!stream.dataset)
    {
        // create dataset
        hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        {
            hsize_t cdims = chunk_size_ / bytes;
            H5Pset_chunk(dcpl, 1, &cdims);

            if (compression_level_ > 0)
            {
                hsize_t level = compression_level_;
                if (bytes > 1) H5Pset_shuffle(dcpl);
                H5Pset_deflate(dcpl, level);
            }
            H5Pset_fletcher32(dcpl);
            hsize_t maxdims = H5S_UNLIMITED;
            stream.size = n;
            stream.space = H5Screate_simple(1, &stream.size, &maxdims);
            stream.dataset = H5Dcreate(file_, id.c_str(), format, stream.space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
            stream.format = format;
        }
        H5Pclose(dcpl);
    }
    else
    {
        // extend dataset size if needed
        if (stream.pos + n > stream.size)
        {
            stream.size = stream.pos + n;
            H5Dset_extent(stream.dataset, &stream.size);
            H5Sclose(stream.space);
            stream.space = H5Dget_space(stream.dataset);
       }
    }
    
 
    // write
    hsize_t count = n;
    H5Sselect_hyperslab(stream.space, H5S_SELECT_SET, &stream.pos, NULL, &count, NULL);
    hid_t mspace = H5Screate_simple(1, &count, &count);
    {
        H5Dwrite(stream.dataset, native_format, mspace, stream.space, H5P_DEFAULT, buf);
    }
    H5Sclose(mspace);

    stream.pos += n;

    return n;
}


stream_offset MzMLbSeekableDevice::seek(const std::string& id, stream_offset off, std::ios_base::seekdir way)
{
    Stream& stream = binary_[id];

    switch (way)
    {
    case std::ios_base::beg:
        stream.pos = off;
        break;
    case std::ios_base::cur:
        stream.pos += off;
        break;
    case std::ios_base::end:
    default:
        stream.pos = stream.size - off;
        break;
    }
    
    return stream.pos;
}

}
/*
} // mzmlb
} // msdata
} // pwiz
*/
