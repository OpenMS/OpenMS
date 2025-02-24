// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{
  class MetaInfoInterface;
  namespace Interfaces
  {
    class IMSDataConsumer;
  }

  namespace Internal
  {
    

    /**
      @brief XML handlers for MzXMLFile

          MapType has to be a MSExperiment or have the same interface.
          Do not use this class. It is only needed in MzXMLFile.
    */

    typedef PeakMap MapType;
    typedef MSSpectrum SpectrumType;

    class OPENMS_DLLAPI MzXMLHandler :
      public XMLHandler
    {
public:
      /**@name Constructors and destructor */
      //@{
      /// Constructor for a read-only handler
      MzXMLHandler(MapType& exp, const String& filename, const String& version, ProgressLogger& logger);

      /// Constructor for a write-only handler
      MzXMLHandler(const MapType& exp, const String& filename, const String& version, const ProgressLogger& logger);

      /// Destructor
      ~MzXMLHandler() override {}
      //@}

      /// handler which support partial loading, implement this method
      LOADDETAIL getLoadDetail() const override;

      /// handler which support partial loading, implement this method
      void setLoadDetail(const LOADDETAIL d) override;

      // Docu in base class
      void endElement(const XMLCh* const uri, const XMLCh* const local_name, const XMLCh* const qname) override;

      // Docu in base class
      void startElement(const XMLCh* const uri, const XMLCh* const local_name, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

      // Docu in base class
      void characters(const XMLCh* const chars, const XMLSize_t length) override;

      /// Write the contents to a stream
      void writeTo(std::ostream& os) override;

      /// Sets the options
      void setOptions(const PeakFileOptions& options)
      {
        options_ = options;
      }

      ///Gets the scan count
      UInt getScanCount() const
      {
        return scan_count_;
      }

      /// Set the IMSDataConsumer consumer which will consume the read data
      void setMSDataConsumer(Interfaces::IMSDataConsumer * consumer)
      {
        consumer_ = consumer;
      }

protected:

      /// Peak type
      typedef MapType::PeakType PeakType;
      /// Spectrum type
      typedef MSSpectrum SpectrumType;

      /// map pointer for reading
      MapType* exp_;
      /// map pointer for writing
      const MapType* cexp_;

      /// Options for loading and storing
      PeakFileOptions options_;

      /**@name temporary data structures to hold parsed data */
      //@{
      Int nesting_level_;

      /**
          @brief Data necessary to generate a single spectrum

          Small struct holds all data necessary to populate a spectrum at a
          later timepoint (since reading of the base64 data and generation of
          spectra can be done at distinct timepoints).
      */
      struct SpectrumData
      {
        UInt peak_count_;
        String precision_;
        String compressionType_;
        String char_rest_;
        SpectrumType spectrum;
        bool skip_data;
      };

      /// Vector of spectrum data stored for later parallel processing
      std::vector< SpectrumData > spectrum_data_;
      //@}

      /// Flag that indicates whether this spectrum should be skipped (due to options)
      bool skip_spectrum_;

      /// spectrum counter (spectra without peaks are not written)
      UInt spec_write_counter_;

      /// Consumer class to work on spectra
      Interfaces::IMSDataConsumer* consumer_;

      /// Consumer class to work on spectra
      UInt scan_count_;

      /// Progress logging class
      const ProgressLogger& logger_;


      /// write metaInfo to xml (usually in nameValue-tag)
      /// returns true if metavalue existed and data was written to the stream
      inline bool writeAttributeIfExists_(std::ostream& os, const MetaInfoInterface& meta, const String& metakey, const String& attname);

      /// write metaInfo to xml (usually in nameValue-tag)
      inline void writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, int indent = 4, const String& tag = "nameValue");

      /**
      @brief Fill a single spectrum with data from input

      @note Do not modify any internal state variables of the class since
      this function will be executed in parallel.
      */
      void doPopulateSpectraWithData_(SpectrumData & spectrum_data);

      /**
      @brief Populate all spectra on the stack with data from input

      Will populate all spectra on the current work stack with data (using
      multiple threads if available) and append them to the result.
      */
      void populateSpectraWithData_();

      /// data processing auxiliary variable
      std::vector< boost::shared_ptr< DataProcessing> > data_processing_;


private:
      /// Not implemented
      MzXMLHandler();
      
      /// initialize members (call from C'tor)
      void init_();
    };

  } // namespace Internal

} // namespace OpenMS

