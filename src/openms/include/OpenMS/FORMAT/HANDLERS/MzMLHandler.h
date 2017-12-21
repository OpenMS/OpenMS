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

#ifndef OPENMS_FORMAT_HANDLERS_MZMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZMLHANDLER_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/VersionInfo.h>

#include <OpenMS/DATASTRUCTURES/CVMappings.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLHandlerHelper.h>
#include <OpenMS/FORMAT/VALIDATORS/MzMLValidator.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/CONCEPT/Helpers.h>

#include <OpenMS/SYSTEM/File.h>

#include <sstream>
#include <boost/shared_ptr.hpp>
#include <iostream>

#include <QRegExp>

//MISSING:
// - more than one selected ion per precursor (warning if more than one)
// - scanWindowList for each acquisition separately (currently for the whole spectrum only)
// - instrumentConfigurationRef attribute for scan (why should the instrument change between scans? - warning if used)
// - scanSettingsRef attribute for instrumentConfiguration tag (currently no information there because of missing mapping file entry - warning if used)

// xs:id/xs:idref prefix list
// - sf_ru : sourceFile (run)
// - sf_sp : sourceFile (spectrum)
// - sf_pr : sourceFile (precursor)
// - sf_ac : sourceFile (acquisition)
// - sa    : sample
// - ic    : instrumentConfiguration
// - so_dp : software (data processing)
// - so_in : software (instrument)
// - dp_sp : dataProcessing (spectrum)
// - dp_bi : dataProcessing (binary data array)
// - dp_ch : dataProcessing (chromatogram)

namespace OpenMS
{
  class ControlledVocabulary;
  namespace Internal
  {

    /**
        @brief XML handler for MzMLFile

        MapType has to be an MSExperiment or have the same interface. In
        read-mode, this class will parse an MzML XML file and append the input
        spectra to the provided MapType object or (if provided separately
        through setMSDataConsumer) to the provided IMSDataConsumer Interface.

        @note Do not use this class. It is only needed in MzMLFile.

        @note Only upon destruction of this class it can be guaranteed that all
        data has been appended to the appropriate consumer of the data. Do not
        try to access the data before that.

        @todo replace hardcoded cv stuff with more flexible handling via obo r/w.
    */

	typedef PeakMap MapType;
	typedef MSSpectrum SpectrumType;
	typedef MSChromatogram ChromatogramType;

    class OPENMS_DLLAPI MzMLHandler :
      public XMLHandler
    {
public:
      /**@name Constructors and destructor */
      //@{

      /// Constructor for a read-only  handler
      MzMLHandler(MapType& exp, const String& filename, const String& version, ProgressLogger& logger);

      /// Constructor for a write-only handler
      MzMLHandler(const MapType& exp, const String& filename, const String& version, const ProgressLogger& logger);

      /// Destructor
      ~MzMLHandler() override;
      //@}

      /**@name XML Handling functions and output writing */
      //@{

      // Docu in base class
      void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname) override;

      // Docu in base class
      void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

      // Docu in base class
      void characters(const XMLCh* const chars, const XMLSize_t length) override;

      //Docu in base class
      void writeTo(std::ostream& os) override;

      //@}

      /**@name PeakFileOptions setters/getters

        The PeakFileOptions object determine the reading and writing of the
        MzML file. In read-mode the lazy-loading options determine whether
        meta-data only or the full raw data is read into memory and how this
        data should be handled. The MS-level, m/z, RT and Intensity range
        options determine which part of the MzML file is read into memory.

       */
      //@{

      /// Set the peak file options
      void setOptions(const PeakFileOptions& opt)
      {
        options_ = opt;
        spectrum_data_.reserve(options_.getMaxDataPoolSize());
        chromatogram_data_.reserve(options_.getMaxDataPoolSize());
      }

      /// Get the peak file options
      PeakFileOptions& getOptions()
      {
        return options_;
      }

      //@}

      /// Get the spectra and chromatogram counts of a file
      void getCounts(Size& spectra_counts, Size& chromatogram_counts)
      {
        spectra_counts = scan_count;
        chromatogram_counts = chromatogram_count;
      }

      /// Set the IMSDataConsumer consumer which will consume the read data
      void setMSDataConsumer(Interfaces::IMSDataConsumer* consumer)
      {
        consumer_ = consumer;
      }

protected:

      /// Peak type
      typedef MapType::PeakType PeakType;
      /// Chromatogram peak type
      typedef MapType::ChromatogramPeakType ChromatogramPeakType;
      /// Spectrum type
      typedef MSSpectrum SpectrumType;
      /// Spectrum type
      typedef MSChromatogram ChromatogramType;

      typedef MzMLHandlerHelper::BinaryData BinaryData;

      void writeSpectrum_(std::ostream& os, const SpectrumType& spec, Size s,
                          Internal::MzMLValidator& validator, bool renew_native_ids,
                          std::vector<std::vector< ConstDataProcessingPtr > >& dps);

      void writeChromatogram_(std::ostream& os, const ChromatogramType& chromatogram, Size c, Internal::MzMLValidator& validator);

      template <typename ContainerT>
      void writeContainerData(std::ostream& os, const PeakFileOptions& pf_options_, const ContainerT& container, String array_type);

      /**
          @brief Populate all spectra on the stack with data from input

          Will populate all spectra on the current work stack with data (using
          multiple threads if available) and append them to the result.
      */
      void populateSpectraWithData();

      /**
          @brief Populate all chromatograms on the stack with data from input

          Will populate all chromatograms on the current work stack with data (using
          multiple threads if available) and append them to the result.
      */
      void populateChromatogramsWithData();

      void addSpectrumMetaData_(const std::vector<MzMLHandlerHelper::BinaryData>& input_data, 
                                const Size n, SpectrumType& spectrum) const;

      /**
          @brief Fill a single spectrum with data from input

          @note Do not modify any internal state variables of the class since
          this function will be executed in parallel.

          Speed: this function takes about 50 % of total load time with a
          single thread and parallelizes linearly up to at least 10 threads.

      */
      void populateSpectraWithData_(std::vector<MzMLHandlerHelper::BinaryData>& input_data,
                                    Size& default_arr_length, const PeakFileOptions& peak_file_options,
                                    SpectrumType& spectrum);

      /**
          @brief Fill a single chromatogram with data from input

          @note Do not modify any internal state variables of the class since
          this function will be executed in parallel.

      */
      void populateChromatogramsWithData_(std::vector<MzMLHandlerHelper::BinaryData>& input_data,
                                          Size& default_arr_length, const PeakFileOptions& peak_file_options,
                                          ChromatogramType& inp_chromatogram);

      template <typename DataType>
      void writeBinaryDataArray(std::ostream& os, const PeakFileOptions& pf_options_, std::vector<DataType> data_to_encode, bool is32bit, String array_type);

      void writeHeader_(std::ostream& os, const MapType& exp, std::vector<std::vector< ConstDataProcessingPtr > >& dps, Internal::MzMLValidator& validator);

      /// map pointer for reading
      MapType* exp_;
      /// map pointer for writing
      const MapType* cexp_;

      /// Options that can be set for loading/storing
      PeakFileOptions options_;

      /**@name temporary data structures to hold parsed data */
      //@{
      /// The current spectrum
      SpectrumType spec_;
      /// The current chromatogram
      ChromatogramType chromatogram_;
      /// The spectrum data (or chromatogram data)
      std::vector<BinaryData> data_;
      /// The default number of peaks in the current spectrum
      Size default_array_length_;
      /// Flag that indicates that we're inside a spectrum (in contrast to a chromatogram)
      bool in_spectrum_list_;
      /// Id of the current list. Used for referencing param group, source file, sample, software, ...
      String current_id_;
      /// The referencing param groups: id => array (accession, value)
      Map<String, std::vector<SemanticValidator::CVTerm> > ref_param_;
      /// The source files: id => SourceFile
      Map<String, SourceFile> source_files_;
      /// The sample list: id => Sample
      Map<String, Sample> samples_;
      /// The software list: id => Software
      Map<String, Software> software_;
      /// The data processing list: id => Instrument
      Map<String, Instrument> instruments_;
      /// The data processing list: id => Instrument
      Map<String, std::vector< DataProcessingPtr > > processing_;
      /// id of the default data processing (used when no processing is defined)
      String default_processing_;

      /**
          @brief Data necessary to generate a single spectrum

          Small struct holds all data necessary to populate a spectrum at a
          later timepoint (since reading of the base64 data and generation of
          spectra can be done at distinct timepoints).
      */
      struct SpectrumData
      {
        std::vector<BinaryData> data;
        Size default_array_length;
        SpectrumType spectrum;
        bool skip_data;
      };

      /// Vector of spectrum data stored for later parallel processing
      std::vector<SpectrumData> spectrum_data_;

      /**
          @brief Data necessary to generate a single chromatogram

          Small struct holds all data necessary to populate a chromatogram at a
          later timepoint (since reading of the base64 data and generation of
          chromatogram can be done at distinct timepoints).
      */
      struct ChromatogramData
      {
        std::vector<BinaryData> data;
        Size default_array_length;
        ChromatogramType chromatogram;
      };

      /// Vector of chromatogram data stored for later parallel processing
      std::vector<ChromatogramData> chromatogram_data_;

      //@}
      /**@name temporary data structures to hold written data */
      //@{
      std::vector<std::pair<std::string, long> > spectra_offsets;
      std::vector<std::pair<std::string, long> > chromatograms_offsets;
      //@}

      /// Decoder/Encoder for Base64-data in MzML
      Base64 decoder_;

      /// Progress logger
      const ProgressLogger& logger_;

      /// Consumer class to work on spectra
      Interfaces::IMSDataConsumer* consumer_;

      /// Counting spectra and chromatograms
      UInt scan_count;
      UInt chromatogram_count;

      /// Flag that indicates whether this spectrum should be skipped (due to options)
      bool skip_chromatogram_;
      bool skip_spectrum_;

      // Remember whether the RT of the spectrum was set or not
      bool rt_set_;

      ///Controlled vocabulary (psi-ms from OpenMS/share/OpenMS/CV/psi-ms.obo)
      ControlledVocabulary cv_;
      CVMappings mapping_;
      //~ Internal::MzMLValidator validator_;

      ///Count of selected ions
      UInt selected_ion_count_;

      /*
      /// Fills the current spectrum with peaks and meta data
      void fillData_();
      */

      /// Fills the current chromatogram with data points and meta data
      void fillChromatogramData_();

      /// Handles CV terms
      void handleCVParam_(const String& parent_parent_tag, const String& parent_tag, /*  const String & cvref, */ const String& accession, const String& name, const String& value, const String& unit_accession = "");

      /// Handles user terms
      void handleUserParam_(const String& parent_parent_tag, const String& parent_tag, const String& name, const String& type, const String& value);

      /// Writes user terms
      void writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, UInt indent, String path, Internal::MzMLValidator& validator) const;

      /// Looks up a child CV term of @p parent_accession with the name @p name. If no such term is found, an empty term is returned.
      ControlledVocabulary::CVTerm getChildWithName_(const String& parent_accession, const String& name) const;

      /// Helper method that writes a software
      void writeSoftware_(std::ostream& os, const String& id, const Software& software, Internal::MzMLValidator& validator);

      /// Helper method that writes a source file
      void writeSourceFile_(std::ostream& os, const String& id, const SourceFile& software, Internal::MzMLValidator& validator);

      /// Helper method that writes a data processing list
      void writeDataProcessing_(std::ostream& os, const String& id, const std::vector< ConstDataProcessingPtr >& dps, Internal::MzMLValidator& validator);

      /// Helper method that write precursor information from spectra and chromatograms
      void writePrecursor_(std::ostream& os, const Precursor& precursor, Internal::MzMLValidator& validator);

      /// Helper method that write precursor information from spectra and chromatograms
      void writeProduct_(std::ostream& os, const Product& product, Internal::MzMLValidator& validator);

      /// Helper method to write an CV based on a meta value
      String writeCV_(const ControlledVocabulary::CVTerm& c, const DataValue& metaValue) const;

      /// Helper method to validate if the given CV is allowed in the current location (path)
      bool validateCV_(const ControlledVocabulary::CVTerm& c, const String& path, const Internal::MzMLValidator& validator) const;
    };

    //--------------------------------------------------------------------------------

  } // namespace Internal
} // namespace OpenMS

#endif
