// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#ifndef OPENMS_FORMAT_MZXMLFILE_H
#define OPENMS_FORMAT_MZXMLFILE_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>

namespace OpenMS
{
  class String;

  /**
      @brief File adapter for MzXML 3.1 files

      @ingroup FileIO
  */
  class OPENMS_DLLAPI MzXMLFile :
    public Internal::XMLFile,
    public ProgressLogger
  {
	typedef PeakMap MapType;

public:
    ///Default constructor
    MzXMLFile();
    ///Destructor
    ~MzXMLFile() override;

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
    void load(const String & filename, MapType & map);

    /**
        @brief Stores a map in a MzXML file.

        @p map has to be a MSExperiment or have the same interface.

        @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String & filename, const MapType & map) const;

    /**
      @brief Transforms a map while loading using the supplied MSDataConsumer.

      The result will not be stored directly but is available through the
      events triggered by the parser and caught by the provided IMSDataConsumer
      object.

      This function should be used if processing and storage of the result can
      be performed directly in the provided IMSDataConsumer object.

      @note Transformation can be speed up by setting skip_full_count which
      does not require a full first pass through the file to compute the
      correct number of spectra and chromatograms in the input file.

      @param filename_in Filename of input mzML file to transform 
      @param consumer Consumer class to operate on the input filename (implementing a transformation)
      @param skip_full_count Whether to skip computing the correct number of spectra and chromatograms in the input file 
    */
    void transform(const String& filename_in, Interfaces::IMSDataConsumer * consumer, bool skip_full_count = false);

    /**
      @brief Transforms a map while loading using the supplied MSDataConsumer

      The result will be stored in the provided map.

      This function should be used if a specific pre-processing should be
      applied to the data before storing them in a map (e.g. if data-reduction
      should be applied to the data before loading all data into memory).

      @param filename_in Filename of input mzXML file to transform 
      @param consumer Consumer class to operate on the input filename (implementing a transformation)
      @param map Map to store the resulting spectra and chromatograms
      @param skip_full_count Whether to skip computing the correct number of spectra and chromatograms in the input file 
    */
    void transform(const String& filename_in, Interfaces::IMSDataConsumer * consumer, MapType& map, bool skip_full_count = false);

protected:

    /// Perform first pass through the file and retrieve the meta-data to initialize the consumer
    void transformFirstPass_(const String& filename_in, Interfaces::IMSDataConsumer * consumer, bool skip_full_count);

private:

    PeakFileOptions options_;
  };
} // namespace OpenMS

#endif // OPENMS_FORMAT_MZXMLFILE_H
