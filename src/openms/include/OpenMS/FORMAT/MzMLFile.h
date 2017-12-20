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

#ifndef OPENMS_FORMAT_MZMLFILE_H
#define OPENMS_FORMAT_MZMLFILE_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h> // StringList
#include <OpenMS/INTERFACES/IMSDataConsumer.h>

namespace OpenMS
{
  /**
    @brief File adapter for MzML files

    This implementation does currently not support the whole functionality of MzML.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI MzMLFile :
    public Internal::XMLFile,
    public ProgressLogger
  {
public:
    ///Default constructor
    MzMLFile();
    ///Destructor
    ~MzMLFile() override;

    /// Mutable access to the options for loading/storing
    PeakFileOptions& getOptions();

    /// Non-mutable access to the options for loading/storing
    const PeakFileOptions& getOptions() const;

    /// set options for loading/storing
    void setOptions(const PeakFileOptions &);

    /**
      @brief Loads a map from a MzML file. Spectra and chromatograms are sorted by default (this can be disabled using PeakFileOptions).

      @p filename The filename with the data
      @p map Is an MSExperiment

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename, PeakMap& map);

    /**
      @brief Loads a map from a MzML file stored in a buffer (in memory).

      @p filename The buffer with the data
      @p map Is an MSExperiment

      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void loadBuffer(const std::string& buffer, PeakMap& map);

    /**
      @brief Only count the number of spectra and chromatograms from a file
    */
    void loadSize(const String & filename, Size& scount, Size& ccount);

    /**
      @brief Stores a map in an MzML file.

      @p map has to be an MSExperiment or have the same interface.

      @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String& filename, const PeakMap& map) const;

    /**
      @brief Stores a map in an output string.

      @p output An empty string to store the result
      @p map has to be an MSExperiment
    */
    void storeBuffer(std::string & output, const PeakMap& map) const;

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
    void transform(const String& filename_in, Interfaces::IMSDataConsumer * consumer, bool skip_full_count = false, bool skip_first_pass = false);

    /**
      @brief Transforms a map while loading using the supplied MSDataConsumer

      The result will be stored in the provided map.

      This function should be used if a specific pre-processing should be
      applied to the data before storing them in a map (e.g. if data-reduction
      should be applied to the data before loading all data into memory).

      @param filename_in Filename of input mzML file to transform
      @param consumer Consumer class to operate on the input filename (implementing a transformation)
      @param map Map to store the resulting spectra and chromatograms
      @param skip_full_count Whether to skip computing the correct number of spectra and chromatograms in the input file
    */
    void transform(const String& filename_in, Interfaces::IMSDataConsumer * consumer, PeakMap& map, bool skip_full_count = false, bool skip_first_pass = false);

    /**
      @brief Checks if a file validates against the XML schema.

      @exception Exception::FileNotFound is thrown if the file cannot be found.
    */
    bool isValid(const String& filename, std::ostream& os = std::cerr);

    /**
      @brief Checks if a file is valid with respect to the mapping file and the controlled vocabulary.

      @param filename File name of the file to be checked.
      @param errors Errors during the validation are returned in this output parameter.
      @param warnings Warnings during the validation are returned in this output parameter.

      @exception Exception::FileNotFound is thrown if the file could not be opened
    */
    bool isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings);

protected:

    /// Perform first pass through the file and retrieve the meta-data to initialize the consumer
    void transformFirstPass_(const String& filename_in, Interfaces::IMSDataConsumer * consumer, bool skip_full_count);

    /// Safe parse that catches exceptions and handles them accordingly
    void safeParse_(const String & filename, Internal::XMLHandler * handler);

private:

    /// Options for loading / storing
    PeakFileOptions options_;

    /// Location of indexed mzML schema
    String indexed_schema_location_;
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_MZMLFILE_H

