//
// $Id: Connection_mzMLb.hpp
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

#ifndef CONNECTION_MZMLB_HPP_
#define CONNECTION_MZMLB_HPP_


#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/concepts.hpp>  // seekable_device
#include <fstream>
#include <vector>
#include <map>
#include <hdf5.h>


namespace pwiz {
namespace msdata {
namespace mzmlb {

using namespace boost::iostreams;

class Connection_mzMLb : public device<seekable> {
public:   
    Connection_mzMLb(const std::string& filename, int chunk_size, int compression_level); // open for writing
    Connection_mzMLb(const std::string& filename, bool identifyOnly = false); // open for reading or identify()
    void close(); // close (called by boost stream on final destruction)

    // boost device<seekable> methods for reading/writing text to "mzML" dataset
    std::streamsize read(char* s, std::streamsize n);
    std::streamsize write(const char* s, std::streamsize n);
    stream_offset seek(stream_offset off, std::ios_base::seekdir way);    

    // methods for reading/writing other binary datasets
    bool exists(const std::string& id);
    std::streamsize size(const std::string& id);
    
    std::streamsize read_opaque(const std::string& id, void* buf, std::streamsize n);
    std::streamsize read(const std::string& id, char* buf, std::streamsize n);
    std::streamsize read(const std::string& id, double* buf, std::streamsize n);
    std::streamsize read(const std::string& id, long* buf, std::streamsize n);
    std::streamsize read(const std::string& id, long long* buf, std::streamsize n);

    std::streamsize write_opaque(const std::string& id, const void* buf, std::streamsize n);
    std::streamsize write(const std::string& id, const char* buf, std::streamsize n);
    std::streamsize write(const std::string& id, const float* buf, std::streamsize n);
    std::streamsize write(const std::string& id, const double* buf, std::streamsize n);
    std::streamsize write(const std::string& id, const long* buf, std::streamsize n);
    std::streamsize write(const std::string& id, const long long* buf, std::streamsize n);

    stream_offset seek(const std::string& id, stream_offset off, std::ios_base::seekdir way);

    
private:
    std::streamsize read(const std::string& id, void* buf, std::streamsize n, hid_t native_format);
    std::streamsize write(const std::string& path, const void* buf, std::streamsize n, hid_t native_format, hid_t format, size_t bytes);

    hid_t opaque_id_;
    hid_t file_;
    unsigned long chunk_size_;
    unsigned long compression_level_;

    struct Stream {
        hid_t dataset;
        hid_t space;
        hsize_t pos;
        hsize_t size;
        hid_t format;

        Stream() : dataset(0), space(0), pos(0), size(0) {}
    };
    
    Stream& open(const std::string& id, hid_t format);
    Stream& create(const std::string& id, hid_t format, size_t bytes);
    
    Stream mzML_; // stream parameters for "mzML" text dataset
    std::map<std::string, Stream> binary_; // stream parameters for binary datasets 
};

} // mzmlb
} // msdata
} // pwiz

#endif /* CONNECTION_MZMLB_HPP_ */
