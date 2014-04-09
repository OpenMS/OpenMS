// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MZXMLFILE_H
#define OPENMS_FORMAT_MZXMLFILE_H

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/FORMAT/HANDLERS/MzXMLHandler.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>

namespace OpenMS
{
  class String;

  /**
      @brief File adapter for MzXML 2.1 files

      @ingroup FileIO
  */
  class OPENMS_DLLAPI MzXMLFile :
    public Internal::XMLFile,
    public ProgressLogger
  {
public:
    ///Default constructor
    MzXMLFile();
    ///Destructor
    ~MzXMLFile();

    /// Mutable access to the options for loading/storing
    PeakFileOptions & getOptions();

    /// Non-mutable access to the options for loading/storing
    const PeakFileOptions & getOptions() const;

    /// set options for loading/storing
    void setOptions(const PeakFileOptions &);

    /**
        @brief Loads a map from a MzXML file.

        @p map has to be a MSExperiment or have the same interface.

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    template <typename MapType>
    void load(const String & filename, MapType & map)
    {
      map.reset();

      //set DocumentIdentifier
      map.setLoadedFileType(filename);
      map.setLoadedFilePath(filename);

      Internal::MzXMLHandler<MapType> handler(map, filename, schema_version_, *this);
      handler.setOptions(options_);
      parse_(filename, &handler);
    }

    /**
        @brief Stores a map in a MzXML file.

        @p map has to be a MSExperiment or have the same interface.

        @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    template <typename MapType>
    void store(const String & filename, const MapType & map) const
    {
      Internal::MzXMLHandler<MapType> handler(map, filename, schema_version_, *this);
      save_(filename, &handler);
    }


    /**
      @brief Transforms a map while loading using the supplied MSDataConsumer.

      The result will not be stored by this class.
    */
    template <typename MapType>
    void transform(const String& filename_in, /* const String& filename_out, */ Interfaces::IMSDataConsumer<MapType> * consumer/* , const MapType& map */)
    {
      
      // First pass through the file -> get the meta-data and hand it to the consumer
      {
        Size scount = 0, ccount = 0;
        MapType experimental_settings;
        bool size_only_before_ = options_.getSizeOnly();
        options_.setSizeOnly(true);
        Internal::MzXMLHandler<MapType> handler(experimental_settings, filename_in, getVersion(), *this);
        handler.setOptions(options_);
        parse_(filename_in, &handler);
        scount = handler.getScanCount();
        options_.setSizeOnly(size_only_before_);
        consumer->setExpectedSize(scount, ccount);
        consumer->setExperimentalSettings(experimental_settings);
      }
      
      // Second pass through the data, now read the spectra!
      {
        MapType dummy;
        Internal::MzXMLHandler<MapType> handler(dummy, filename_in, getVersion(), *this);
        handler.setOptions(options_);
        handler.setMSDataConsumer(consumer);
        // TODO catch errors as above ?
        parse_(filename_in, &handler);
      }
    }

    /**
      @brief Transforms a map while loading using the supplied MSDataConsumer

      The result will be stored in the supplied map.
    */
    template <typename MapType>
    void transform(const String& filename_in, /* const String& filename_out, */ Interfaces::IMSDataConsumer<MapType> * consumer, MapType& map)
    {
      // First pass through the file -> get the meta-data and hand it to the consumer
      {
        Size scount = 0, ccount = 0;
        MapType experimental_settings;
        bool size_only_before_ = options_.getSizeOnly();
        options_.setSizeOnly(true);
        Internal::MzXMLHandler<MapType> handler(experimental_settings, filename_in, getVersion(), *this);
        handler.setOptions(options_);
        parse_(filename_in, &handler);
        scount = handler.getScanCount();
        options_.setSizeOnly(size_only_before_);
        consumer->setExpectedSize(scount, ccount);
        consumer->setExperimentalSettings(experimental_settings);
      }

      // Second pass through the data, now read the spectra!
      {
        Internal::MzXMLHandler<MapType> handler(map, filename_in, getVersion(), *this);
        bool always_append_data = options_.getAlwaysAppendData();
        options_.setAlwaysAppendData(true);
        handler.setOptions(options_);
        handler.setMSDataConsumer(consumer);

        // TODO catch errors as above ?
        parse_(filename_in, &handler);
        options_.setAlwaysAppendData(always_append_data);
      }
    }

private:
    PeakFileOptions options_;
  };
} // namespace OpenMS

#endif // OPENMS_FOMAT_MZXMLFILE_H
