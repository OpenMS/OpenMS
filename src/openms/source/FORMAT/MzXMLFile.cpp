// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzXMLFile.h>

#include <OpenMS/FORMAT/HANDLERS/MzXMLHandler.h>

using namespace std;

namespace OpenMS
{
  MzXMLFile::MzXMLFile() :
    XMLFile("/SCHEMAS/mzXML_idx_3.1.xsd", "3.1")
  {
  }

  MzXMLFile::~MzXMLFile()
  {
  }

  PeakFileOptions & MzXMLFile::getOptions()
  {
    return options_;
  }

  const PeakFileOptions & MzXMLFile::getOptions() const
  {
    return options_;
  }

  void MzXMLFile::setOptions(const PeakFileOptions & options)
  {
      options_ = options;
  }

  void MzXMLFile::load(const String & filename, MapType & map)
  {
    map.reset();

    //set DocumentIdentifier
    map.setLoadedFileType(filename);
    map.setLoadedFilePath(filename);

    Internal::MzXMLHandler handler(map, filename, schema_version_, *this);
    handler.setOptions(options_);
    parse_(filename, &handler);
  }

  void MzXMLFile::store(const String & filename, const MapType & map) const
  {
    Internal::MzXMLHandler handler(map, filename, schema_version_, *this);
    handler.setOptions(options_);
    save_(filename, &handler);
  }

  void MzXMLFile::transform(const String& filename_in, Interfaces::IMSDataConsumer * consumer, bool skip_full_count)
  {
    // First pass through the file -> get the meta-data and hand it to the consumer
    transformFirstPass_(filename_in, consumer, skip_full_count);
    
    // Second pass through the data, now read the spectra!
    {
      MapType dummy;
      Internal::MzXMLHandler handler(dummy, filename_in, getVersion(), *this);
      handler.setOptions(options_);
      handler.setMSDataConsumer(consumer);
      parse_(filename_in, &handler);
    }
  }

  void MzXMLFile::transform(const String& filename_in, Interfaces::IMSDataConsumer * consumer, MapType& map, bool skip_full_count)
  {
    // First pass through the file -> get the meta-data and hand it to the consumer
    transformFirstPass_(filename_in, consumer, skip_full_count);

    // Second pass through the data, now read the spectra!
    {
      PeakFileOptions tmp_options(options_);
      Internal::MzXMLHandler handler(map, filename_in, getVersion(), *this);
      tmp_options.setAlwaysAppendData(true);
      handler.setOptions(tmp_options);
      handler.setMSDataConsumer(consumer);

      parse_(filename_in, &handler);
    }
  }

  void MzXMLFile::transformFirstPass_(const String& filename_in, Interfaces::IMSDataConsumer * consumer, bool skip_full_count)
  {
    // Create temporary objects and counters
    PeakFileOptions tmp_options(options_);
    Size scount = 0, ccount = 0;
    MapType experimental_settings;
    Internal::MzXMLHandler handler(experimental_settings, filename_in, getVersion(), *this);

    // set temporary options for handler
    tmp_options.setMetadataOnly( skip_full_count );
    handler.setOptions(tmp_options);
    handler.setLoadDetail(Internal::XMLHandler::LD_RAWCOUNTS);
    parse_(filename_in, &handler);

    // After parsing, collect information
    scount = handler.getScanCount();
    consumer->setExpectedSize(scount, ccount);
    consumer->setExperimentalSettings(experimental_settings);
  }

} // namespace OpenMS

