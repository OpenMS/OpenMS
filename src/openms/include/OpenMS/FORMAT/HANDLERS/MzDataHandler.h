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

#ifndef OPENMS_FORMAT_HANDLERS_MZDATAHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZDATAHANDLER_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <sstream>

namespace OpenMS
{
  namespace Internal
  {
    /**
        @brief XML handler for MzDataFile

        MapType has to be a MSExperiment or have the same interface.
        Do not use this class. It is only needed in MzDataFile.

        @improvement Add implementation and tests of 'supDataArray' to store IntegerDataArray and StringDataArray of MSSpectrum (Hiwi)
    */

	typedef PeakMap MapType;
	typedef MSSpectrum SpectrumType;
	typedef MSChromatogram ChromatogramType;

    class OPENMS_DLLAPI MzDataHandler :
      public XMLHandler
    {
public:
      /**@name Constructors and destructor */
      //@{
      /// Constructor for a write-only handler
      MzDataHandler(MapType & exp, const String & filename, const String & version, ProgressLogger & logger);

      /// Constructor for a read-only handler
      MzDataHandler(const MapType & exp, const String & filename, const String & version, const ProgressLogger & logger);

      /// Destructor
      ~MzDataHandler() override
      {
      }

      //@}


      // Docu in base class
      void endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname) override;

      // Docu in base class
      void startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes) override;

      // Docu in base class
      void characters(const XMLCh * const chars, const XMLSize_t length) override;

      /// Writes the contents to a stream
      void writeTo(std::ostream & os) override;

      ///Sets the options
      void setOptions(const PeakFileOptions & options)
      {
        options_ = options;
      }

private:
      void init_();

protected:

      /// Peak type
      typedef MapType::PeakType PeakType;
      /// Spectrum type
      typedef MSSpectrum SpectrumType;

      /// map pointer for reading
      MapType * exp_;
      /// map pointer for writing
      const MapType * cexp_;

      ///Options that can be set for loading/storing
      PeakFileOptions options_;

      /**@name temporary datastructures to hold parsed data */
      //@{
      /// The number of peaks in the current spectrum
      UInt peak_count_;
      /// The current spectrum
      SpectrumType spec_;
      /// An array of pairs MetaInfodescriptions and their ids
      std::vector<std::pair<String, MetaInfoDescription> > meta_id_descs_;
      /// encoded data which is read and has to be decoded
      std::vector<String> data_to_decode_;
      /// floating point numbers which have to be encoded and written
      std::vector<float> data_to_encode_;
      std::vector<std::vector<float> > decoded_list_;
      std::vector<std::vector<double> > decoded_double_list_;
      std::vector<String> precisions_;
      std::vector<String> endians_;
      //@}

      /// Decoder/Encoder for Base64-data in MzData
      Base64 decoder_;

      /// Flag that indicates whether this spectrum should be skipped (due to options)
      bool skip_spectrum_;

      /// Progress logger
      const ProgressLogger & logger_;

      /// fills the current spectrum with peaks and meta data
      void fillData_();

      ///@name cvParam and userParam handling methods (for mzData and featureXML)
      //@{
      /**
          @brief write cvParam containing strings to stream

          @p value string value
          @p acc accession number defined by ontology
          @p name term defined by ontology
          @p indent number of tabs used in front of tag

          Example:
          &lt;cvParam cvLabel="psi" accession="PSI:@p acc" name="@p name" value="@p value"/&gt;
      */
      inline void writeCVS_(std::ostream & os, double value, const String & acc, const String & name, UInt indent = 4) const;

      /**
          @brief write cvParam containing strings to stream

          @p value string value
          @p acc accession number defined by ontology
          @p name term defined by ontology
          @p indent number of tabs used in front of tag

          Example:
          &lt;cvParam cvLabel="psi" accession="PSI:@p acc" name="@p name" value="@p value"/&gt;
      */
      inline void writeCVS_(std::ostream & os, const String & value, const String & acc, const String & name, UInt indent = 4) const;

      /**
          @brief write cvParam element to stream

          @p os Output stream
          @p value enumeration value
          @p map index if the terms in cv_terms_
          @p acc accession number defined by ontology
          @p name term defined by ontology
          @p indent number of tabs used in front of tag

          Example:
          &lt;cvParam cvLabel="psi" accession="PSI:@p acc" name="@p name" value=""/&gt;
      */
      inline void writeCVS_(std::ostream & os, UInt value, UInt map, const String & acc, const String & name, UInt indent = 4);

      ///Writing the MetaInfo as UserParam to the file
      inline void writeUserParam_(std::ostream & os, const MetaInfoInterface & meta, UInt indent = 4);

      /**
          @brief read attributes of MzData's cvParamType

          Example:
          &lt;cvParam cvLabel="psi" accession="PSI:1000001" name="@p name" value="@p value"/&gt;
          @p name and sometimes @p value are defined in the MzData ontology.
      */
      void cvParam_(const String & name, const String & value);
      //@}

      /**
          @brief write binary data to stream (first one)

          The @p name and @p id are only used if the @p tag is @em supDataArrayBinary or @em supDataArray.
      */
      inline void writeBinary_(std::ostream & os, Size size, const String & tag, const String & name = "", SignedSize id = -1);

      //Data processing auxiliary variable
      boost::shared_ptr< DataProcessing > data_processing_;

    };

    //--------------------------------------------------------------------------------
  }   // namespace Internal

} // namespace OpenMS

#endif
