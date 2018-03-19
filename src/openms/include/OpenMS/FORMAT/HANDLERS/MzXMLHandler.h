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
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_MZXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZXMLHANDLER_H

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>

#include <stack>

namespace OpenMS
{
  class MetaInfoInterface;

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
      UInt getScanCount()
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
      Base64 decoder_;
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
      inline std::ostream& writeAttributeIfExists_(std::ostream& os, const MetaInfoInterface& meta, const String& metakey, const String& attname);

      /// write metaInfo to xml (usually in nameValue-tag)
      inline void writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, int indent = 4, String tag = "nameValue");

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

      // init all the static members, which is necessary because otherwise the undefined order will cause problems
      void initStaticMembers_();

      static const XMLCh* s_value_;
      static const XMLCh* s_count_;
      static const XMLCh* s_type_;
      static const XMLCh* s_name_;
      static const XMLCh* s_version_;
      static const XMLCh* s_filename_;
      static const XMLCh* s_filetype_;
      static const XMLCh* s_filesha1_;
      static const XMLCh* s_completiontime_;
      static const XMLCh* s_precision_;
      static const XMLCh* s_byteorder_;
      static const XMLCh* s_contentType_;
      static const XMLCh* s_compressionType_;
      static const XMLCh* s_precursorintensity_;
      static const XMLCh* s_precursorcharge_;
      static const XMLCh* s_windowwideness_;
      static const XMLCh* s_mslevel_;
      static const XMLCh* s_peakscount_;
      static const XMLCh* s_polarity_;
      static const XMLCh* s_scantype_;
      static const XMLCh* s_filterline_;
      static const XMLCh* s_retentiontime_;
      static const XMLCh* s_startmz_;
      static const XMLCh* s_endmz_;
      static const XMLCh* s_first_;
      static const XMLCh* s_last_;
      static const XMLCh* s_phone_;
      static const XMLCh* s_email_;
      static const XMLCh* s_uri_;
      static const XMLCh* s_num_;
      static const XMLCh* s_intensitycutoff_;
      static const XMLCh* s_centroided_;
      static const XMLCh* s_deisotoped_;
      static const XMLCh* s_chargedeconvoluted_;
    };

  } // namespace Internal

} // namespace OpenMS

#endif
