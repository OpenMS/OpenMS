// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Chris Bielow, Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Helpers.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/DATASTRUCTURES/CVMappings.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLHandlerHelper.h>

#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>

#include <map>


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
  namespace Interfaces
  {
    class IMSDataConsumer;
  }

  namespace Internal
  {
    class MzMLValidator;

	  typedef PeakMap MapType;
	  typedef MSSpectrum SpectrumType;
	  typedef MSChromatogram ChromatogramType;

    /**@brief Handler for mzML file format
     *
     * This class handles parsing and writing of the mzML file format. It
     * supports reading data directly into memory or parsing on-the-fly using a
     * consumer (see @ref consumer_a "Setting a consumer").  In read-mode, this
     * class will parse an MzML XML file and append the input spectra to the
     * provided MSExperiment object or to the provided
     * Interfaces::IMSDataConsumer (needs to be provided separately through
     * setMSDataConsumer()).
     *
     * Functions constituting the XML reading/writing interface can be found
     * under @ref xml_handling "XML Handling functions", helper functions
     * specifically used for writing out to XML are organized under @ref
     * helper_write "Writing functions" and helper functions used for reading
     * in XML from disk are organized under @ref helper_read "Reading
     * functions".
     *
     * See the MzMLHandlerHelper for additional helper functions that are
     * independent of state.
     *
     * @note Do not use this class directly. It is only needed in MzMLFile.
     *
     * @note Only upon destruction of this class it can be guaranteed that all
     * data has been appended to the appropriate consumer of the data. Do not
     * try to access the data before that.
     *
     * @todo replace hardcoded cv stuff with more flexible handling via obo r/w.
     *
     **/
    class OPENMS_DLLAPI MzMLHandler :
      public XMLHandler
    {
public:

      /**@name Constructors and destructor */
      //@{

      /// Constructor for a read-only  handler
      MzMLHandler(MapType& exp, const String& filename, const String& version, const ProgressLogger& logger);

      /// Constructor for a write-only handler
      MzMLHandler(const MapType& exp, const String& filename, const String& version, const ProgressLogger& logger);

      /// Destructor
      ~MzMLHandler() override;
      //@}

      /**
       * @anchor xml_handling
       * @name XML Handling and data parsing
       **/
      //@{

      /// Docu in base class XMLHandler::endElement
      void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname) override;

      /// Docu in base class XMLHandler::startElelement
      void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

      /// Docu in base class XMLHandler::characters
      void characters(const XMLCh* const chars, const XMLSize_t length) override;

      /// Docu in base class XMLHandler::writeTo
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
      void setOptions(const PeakFileOptions& opt);

      /// Get the peak file options
      PeakFileOptions& getOptions();

      //@}

      /// Get the spectra and chromatogram counts of a file
      void getCounts(Size& spectra_counts, Size& chromatogram_counts);

      /**@name IMSDataConsumer setter
         @anchor consumer_a

        The IMSDataConsumer object allows the user to specify a callback object
        which can consume spectra and chromatograms on the fly. The consumer
        does not have to wait until data is read fully into memory, but will
        start receiving data as soon as it is available (read from disk).
       */
      //@{

      /// Set the IMSDataConsumer consumer which will consume the read data
      void setMSDataConsumer(Interfaces::IMSDataConsumer* consumer);
      //@}

      /// handler which support partial loading, implement this method
      LOADDETAIL getLoadDetail() const override;

      /// handler which support partial loading, implement this method
      void setLoadDetail(const LOADDETAIL d) override;

protected:

      /// delegated constructor for the two public versions
      MzMLHandler(const String& filename, const String& version, const ProgressLogger& logger);

      /// Peak type
      typedef MapType::PeakType PeakType;
      /// Chromatogram peak type
      typedef MapType::ChromatogramPeakType ChromatogramPeakType;
      /// Spectrum type
      typedef MSSpectrum SpectrumType;
      /// Spectrum type
      typedef MSChromatogram ChromatogramType;

      typedef MzMLHandlerHelper::BinaryData BinaryData;

      /**@name Helper functions for storing data in memory
       * @anchor helper_read
       */
      //@{

      /**
          @brief Populate all spectra on the stack with data from input

          Will populate all spectra on the current work stack with data (using
          multiple threads if available) and append them to the result.
      */
      void populateSpectraWithData_();

      /**
          @brief Populate all chromatograms on the stack with data from input

          Will populate all chromatograms on the current work stack with data (using
          multiple threads if available) and append them to the result.
      */
      void populateChromatogramsWithData_();

      /**
          @brief Add extra data arrays to a spectrum

          Add the float, integer and string data arrays to a spectrum.
      */
      void addSpectrumMetaData_(const std::vector<MzMLHandlerHelper::BinaryData>& input_data,
                                const Size n,
                                SpectrumType& spectrum) const;

      /**
          @brief Fill a single spectrum with data from input

          @note Do not modify any internal state variables of the class since
          this function will be executed in parallel.

          @note This function takes about 50 % of total load time with a
          single thread and parallelizes linearly up to at least 10 threads.

          @param input_data The input data with which to fill the spectra
          @param length The input data length (number of data points)
          @param peak_file_options Will be used if only part of the data should be copied (RT, mz or intensity range)
          @param spectrum The output spectrum

      */
      void populateSpectraWithData_(std::vector<MzMLHandlerHelper::BinaryData>& input_data,
                                    Size& length,
                                    const PeakFileOptions& peak_file_options,
                                    SpectrumType& spectrum);

      /**
          @brief Fill a single chromatogram with data from input

          @note Do not modify any internal state variables of the class since
          this function will be executed in parallel.

          @param input_data The input data with which to fill the spectra
          @param length The input data length (number of data points)
          @param peak_file_options Will be used if only part of the data should be copied (RT, mz or intensity range)
          @param chromatogram The output chromatogram

      */
      void populateChromatogramsWithData_(std::vector<MzMLHandlerHelper::BinaryData>& input_data,
                                          Size& length,
                                          const PeakFileOptions& peak_file_options,
                                          ChromatogramType& chromatogram);

      /// Fills the current chromatogram with data points and meta data
      void fillChromatogramData_();

      /// Handles CV terms
      void handleCVParam_(const String& parent_parent_tag,
                          const String& parent_tag,
                          const String& accession,
                          const String& name,
                          const String& value,
                          const String& unit_accession = "");

      /// Handles user terms
      void handleUserParam_(const String& parent_parent_tag,
                            const String& parent_tag,
                            const String& name,
                            const String& type,
                            const String& value,
                            const String& unit_accession = "");
      //@}

      /**
       * @anchor helper_write
       * @name Helper functions for writing data
       */
      //@{

      /// Write out XML header including (everything up to spectrumList / chromatogramList
      void writeHeader_(std::ostream& os,
                        const MapType& exp,
                        std::vector<std::vector< ConstDataProcessingPtr > >& dps,
                        const Internal::MzMLValidator& validator);


      /// Write out a single spectrum
      void writeSpectrum_(std::ostream& os,
                          const SpectrumType& spec,
                          Size spec_idx,
                          const Internal::MzMLValidator& validator,
                          bool renew_native_ids,
                          std::vector<std::vector< ConstDataProcessingPtr > >& dps);

      /// Write out a single chromatogram
      void writeChromatogram_(std::ostream& os,
                              const ChromatogramType& chromatogram,
                              Size chrom_idx,
                              const Internal::MzMLValidator& validator);

      template <typename ContainerT>
      void writeContainerData_(std::ostream& os, const PeakFileOptions& pf_options_, const ContainerT& container, const String& array_type);

      /**
          @brief Write a single \<binaryDataArray\> element to the output

          @param os The stream into which to write
          @param options The PeakFileOptions which determines the compression type to use
          @param data The data to write (32bit float or 64 bit double)
          @param is32bit Whether data is 32bit
          @param array_type Which type of data array is written (mz, time, intensity or float_data)

          @note The data argument may be modified by the function (see Base64 for reasons why)

      */
      template <typename DataType>
      void writeBinaryDataArray_(std::ostream& os,
                                 const PeakFileOptions& options,
                                 std::vector<DataType>& data,
                                 bool is32bit,
                                 String array_type);

      /**
          @brief Write a single \<binaryDataArray\> element for a float data array to the output

          This is only for non-standard data arrays which are treated slightly
          differently by the standard.

          @param os The stream into which to write
          @param options The PeakFileOptions which determines the compression type to use
          @param array The data to write
          @param spec_chrom_idx The index of the current spectrum or chromatogram
          @param array_idx The index of the current float data array
          @param is_spectrum Whether data is associated with a spectrum (if false, a chromatogram is assumed)
          @param validator Validator object
      */
      void writeBinaryFloatDataArray_(std::ostream& os,
                                      const PeakFileOptions& options,
                                      const OpenMS::DataArrays::FloatDataArray& array,
                                      const Size spec_chrom_idx,
                                      const Size array_idx,
                                      bool is_spectrum,
                                      const Internal::MzMLValidator& validator);

      /// Writes user terms
      void writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, UInt indent, const String& path, const Internal::MzMLValidator& validator, const std::set<String>& exclude = {}) const;

      /// Helper method that writes a software
      void writeSoftware_(std::ostream& os, const String& id, const Software& software, const Internal::MzMLValidator& validator);

      /// Helper method that writes a source file
      void writeSourceFile_(std::ostream& os, const String& id, const SourceFile& software, const Internal::MzMLValidator& validator);

      /// Helper method that writes a data processing list
      void writeDataProcessing_(std::ostream& os, const String& id, const std::vector< ConstDataProcessingPtr >& dps, const Internal::MzMLValidator& validator);

      /// Helper method that write precursor information from spectra and chromatograms
      void writePrecursor_(std::ostream& os, const Precursor& precursor, const Internal::MzMLValidator& validator);

      /// Helper method that write precursor information from spectra and chromatograms
      void writeProduct_(std::ostream& os, const Product& product, const Internal::MzMLValidator& validator);

      /// Helper method to write an CV based on a meta value
      String writeCV_(const ControlledVocabulary::CVTerm& c, const DataValue& metaValue) const;

      /// Helper method to validate if the given CV is allowed in the current location (path)
      bool validateCV_(const ControlledVocabulary::CVTerm& c, const String& path, const Internal::MzMLValidator& validator) const;

      /// Helper method to look up a child CV term of @p parent_accession with the name @p name. If no such term is found, an empty term is returned.
      ControlledVocabulary::CVTerm getChildWithName_(const String& parent_accession, const String& name) const;

      //@}

      // MEMBERS

      /// map pointer for reading
      MapType* exp_{ nullptr };

      /// map pointer for writing
      const MapType* cexp_{ nullptr };

      /// Options that can be set for loading/storing
      PeakFileOptions options_;

      /**@name temporary data structures to hold parsed data */
      //@{
      /// The current spectrum
      SpectrumType spec_;
      /// The current chromatogram
      ChromatogramType chromatogram_;
      /// The spectrum data (or chromatogram data)
      std::vector<BinaryData> bin_data_;
      /// The default number of peaks in the current spectrum
      Size default_array_length_;
      /// Flag that indicates that we're inside a spectrum (in contrast to a chromatogram)
      bool in_spectrum_list_{ false };
      /// Flag that indicates whether this spectrum should be skipped (e.g. due to options)
      bool skip_spectrum_{ false };
      /// Flag that indicates whether this chromatogram should be skipped (e.g. due to options)
      bool skip_chromatogram_{ false };
      /// Remember whether the RT of the spectrum was set or not
      bool rt_set_{ false };
      /// Id of the current list. Used for referencing param group, source file, sample, software, ...
      String current_id_;
      /// The referencing param groups: id => array (accession, value)
      std::map<String, std::vector<SemanticValidator::CVTerm> > ref_param_;
      /// The source files: id => SourceFile
      std::map<String, SourceFile> source_files_;
      /// The sample list: id => Sample
      std::map<String, Sample> samples_;
      /// The software list: id => Software
      std::map<String, Software> software_;
      /// The data processing list: id => Instrument
      std::map<String, Instrument> instruments_;
      /// CV terms-path-combinations that have been checked in validateCV_()
      mutable std::map<std::pair<String, String>, bool> cached_terms_;
      /// The data processing list: id => Instrument
      std::map<String, std::vector< DataProcessingPtr > > processing_;
      /// id of the default data processing (used when no processing is defined)
      String default_processing_;
      /// Count of selected ions
      UInt selected_ion_count_{ 0 };

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
      
      /**@name temporary data structures to hold written data
       *
       * These data structures are used to store binary offsets required by the
       * indexedMzML format, specifically the start of each \<spectrum\> and
       * \<chromatogram\> tag is stored and will then be stored at the end of the file.
       **/
      //@{
      std::vector<std::pair<std::string, Int64> > spectra_offsets_; ///< Stores binary offsets for each \<spectrum\> tag
      std::vector<std::pair<std::string, Int64> > chromatograms_offsets_; ///< Stores binary offsets for each \<chromatogram\> tag
      //@}

      /// Progress logger
      const ProgressLogger& logger_;

      /// Consumer class to work on spectra
      Interfaces::IMSDataConsumer* consumer_{ nullptr };

      /**@name temporary data structures for counting spectra and chromatograms */
      UInt scan_count_{ 0 };  ///< number of scans which pass the options-filter
      UInt chromatogram_count_{ 0 }; ///< number of chromatograms which pass the options-filter
      Int scan_count_total_{ -1 }; ///< total number of scans in mzML file (according to 'count' attribute)
      Int chrom_count_total_{ -1 }; ///< total number of chromatograms in mzML file (according to 'count' attribute)
      //@}

      ///Controlled vocabulary (psi-ms from OpenMS/share/OpenMS/CV/psi-ms.obo)
      const ControlledVocabulary& cv_;
      CVMappings mapping_;

    };

    //--------------------------------------------------------------------------------

  } // namespace Internal
} // namespace OpenMS

