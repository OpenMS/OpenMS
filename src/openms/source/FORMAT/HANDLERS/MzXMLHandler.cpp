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

#include <OpenMS/FORMAT/HANDLERS/MzXMLHandler.h>

namespace OpenMS
{
  namespace Internal
  {

    // class holding the byte offsets to '<scan>' tags in the mzXML file; req. to create the index at the end
    struct IndexPos
    {
      Size id_;
      std::ostream::streampos pos_;
      IndexPos(const Size id, const std::ostream::streampos pos)
        : id_(id),
        pos_(pos) {}
    };


    //--------------------------------------------------------------------------------

    // this cannot be moved into a function as VS2008 does not allow more than 31 static members in a function .. don't ask...
    const XMLCh * MzXMLHandler::s_value_ = 0;
    const XMLCh * MzXMLHandler::s_count_ = 0;
    const XMLCh * MzXMLHandler::s_type_ = 0;
    const XMLCh * MzXMLHandler::s_name_ = 0;
    const XMLCh * MzXMLHandler::s_version_ = 0;
    const XMLCh * MzXMLHandler::s_filename_ = 0;
    const XMLCh * MzXMLHandler::s_filetype_ = 0;
    const XMLCh * MzXMLHandler::s_filesha1_ = 0;
    const XMLCh * MzXMLHandler::s_completiontime_ = 0;
    const XMLCh * MzXMLHandler::s_precision_ = 0;
    const XMLCh * MzXMLHandler::s_byteorder_ = 0;
    const XMLCh * MzXMLHandler::s_contentType_ = 0;
    const XMLCh * MzXMLHandler::s_compressionType_ = 0;
    const XMLCh * MzXMLHandler::s_precursorintensity_ = 0;
    const XMLCh * MzXMLHandler::s_precursorcharge_ = 0;
    const XMLCh * MzXMLHandler::s_windowwideness_ = 0;
    const XMLCh * MzXMLHandler::s_mslevel_ = 0;
    const XMLCh * MzXMLHandler::s_peakscount_ = 0;
    const XMLCh * MzXMLHandler::s_polarity_ = 0;
    const XMLCh * MzXMLHandler::s_scantype_ = 0;
    const XMLCh * MzXMLHandler::s_filterline_ = 0;
    const XMLCh * MzXMLHandler::s_retentiontime_ = 0;
    const XMLCh * MzXMLHandler::s_startmz_ = 0;
    const XMLCh * MzXMLHandler::s_endmz_ = 0;
    const XMLCh * MzXMLHandler::s_first_ = 0;
    const XMLCh * MzXMLHandler::s_last_ = 0;
    const XMLCh * MzXMLHandler::s_phone_ = 0;
    const XMLCh * MzXMLHandler::s_email_ = 0;
    const XMLCh * MzXMLHandler::s_uri_ = 0;
    const XMLCh * MzXMLHandler::s_num_ = 0;
    const XMLCh * MzXMLHandler::s_intensitycutoff_ = 0;
    const XMLCh * MzXMLHandler::s_centroided_ = 0;
    const XMLCh * MzXMLHandler::s_deisotoped_ = 0;
    const XMLCh * MzXMLHandler::s_chargedeconvoluted_ = 0;

    /// Constructor for a read-only handler
    MzXMLHandler::MzXMLHandler(MapType& exp, const String& filename, const String& version, ProgressLogger& logger) :
      XMLHandler(filename, version),
      exp_(&exp),
      cexp_(0),
      decoder_(),
      nesting_level_(0),
      skip_spectrum_(false),
      spec_write_counter_(1),
      consumer_(NULL),
      scan_count_(0),
      logger_(logger)
    {
      init_();
    }

    /// Constructor for a write-only handler
    MzXMLHandler::MzXMLHandler(const MapType& exp, const String& filename, const String& version, const ProgressLogger& logger) :
      XMLHandler(filename, version),
      exp_(0),
      cexp_(&exp),
      decoder_(),
      nesting_level_(0),
      skip_spectrum_(false),
      spec_write_counter_(1),
      consumer_(NULL),
      scan_count_(0),
      logger_(logger)
    {
      init_();
    }

    void MzXMLHandler::startElement(const XMLCh* const /*uri*/,
      const XMLCh* const /*local_name*/, const XMLCh* const qname,
      const xercesc::Attributes& attributes)
    {
      OPENMS_PRECONDITION(nesting_level_ >= 0, "Nesting level needs to be zero or more")

        static bool init_static_members(false);
      if (!init_static_members)
      {
        initStaticMembers_();
      }

      String tag = sm_.convert(qname);
      open_tags_.push_back(tag);
      //std::cout << " -- Start -- "<< tag << " -- " << "\n";

      //Skip all tags until the the next scan
      if (skip_spectrum_ && tag != "scan")
        return;

      if (tag == "msRun")
      {
        Int count = 0;
        optionalAttributeAsInt_(count, attributes, s_count_);
        exp_->reserve(count);
        logger_.startProgress(0, count, "loading mzXML file");
        scan_count_ = 0;
        data_processing_.clear();
        //start and end time are xs:duration. This makes no sense => ignore them
      }
      else if (tag == "parentFile")
      {
        SourceFile sf;
        sf.setNameOfFile(attributeAsString_(attributes, s_filename_));
        sf.setFileType(attributeAsString_(attributes, s_filetype_));
        sf.setChecksum(attributeAsString_(attributes, s_filesha1_), SourceFile::SHA1);
        exp_->getSourceFiles().push_back(sf);
      }
      else if (tag == "software")
      {
        String& parent_tag = *(open_tags_.end() - 2);
        if (parent_tag == "dataProcessing")
        {
          data_processing_.back()->getSoftware().setVersion(attributeAsString_(attributes, s_version_));
          data_processing_.back()->getSoftware().setName(attributeAsString_(attributes, s_name_));
          data_processing_.back()->setMetaValue("#type", String(attributeAsString_(attributes, s_type_)));

          String time;
          optionalAttributeAsString_(time, attributes, s_completiontime_);
          data_processing_.back()->setCompletionTime(asDateTime_(time));
        }
        else if (parent_tag == "msInstrument")
        {
          exp_->getInstrument().getSoftware().setVersion(attributeAsString_(attributes, s_version_));
          exp_->getInstrument().getSoftware().setName(attributeAsString_(attributes, s_name_));
        }
      }
      else if (tag == "peaks")
      {
        //precision
        spectrum_data_.back().precision_ = "32";
        optionalAttributeAsString_(spectrum_data_.back().precision_, attributes, s_precision_);
        if (spectrum_data_.back().precision_ != "32" && spectrum_data_.back().precision_ != "64")
        {
          error(LOAD, String("Invalid precision '") + spectrum_data_.back().precision_ + "' in element 'peaks'");
        }
        //byte order
        String byte_order = "network";
        optionalAttributeAsString_(byte_order, attributes, s_byteorder_);
        if (byte_order != "network")
        {
          error(LOAD, String("Invalid or missing byte order '") + byte_order + "' in element 'peaks'. Must be 'network'!");
        }
        //pair order
        String pair_order = "m/z-int";
        optionalAttributeAsString_(pair_order, attributes, s_contentType_);
        if (pair_order != "m/z-int")
        {
          error(LOAD, String("Invalid or missing pair order '") + pair_order + "' in element 'peaks'. Must be 'm/z-int'!");
        }
        //compressionType
        spectrum_data_.back().compressionType_ = "none";
        optionalAttributeAsString_(spectrum_data_.back().compressionType_, attributes, s_compressionType_);
        if (spectrum_data_.back().compressionType_ != "none" && spectrum_data_.back().compressionType_ != "zlib")
        {
          error(LOAD, String("Invalid compression type ") + spectrum_data_.back().compressionType_ + "in elements 'peaks'. Must be 'none' or 'zlib'! ");
        }
      }
      else if (tag == "precursorMz")
      {
        //add new precursor
        spectrum_data_.back().spectrum.getPrecursors().push_back(Precursor());
        //intensity
        try
        {
          spectrum_data_.back().spectrum.getPrecursors().back().setIntensity(attributeAsDouble_(attributes, s_precursorintensity_));
        }
        catch (Exception::ParseError& /*e*/)
        {
          error(LOAD, "Mandatory attribute 'precursorIntensity' of tag 'precursorMz' not found! Setting precursor intensity to zero!");
        }
        //charge
        Int charge = 0;
        if (optionalAttributeAsInt_(charge, attributes, s_precursorcharge_))
        {
          spectrum_data_.back().spectrum.getPrecursors().back().setCharge(charge);
        }
        //window bounds (here only the width is stored in both fields - this is corrected when we parse the m/z position)
        double window = 0.0;
        if (optionalAttributeAsDouble_(window, attributes, s_windowwideness_))
        {
          spectrum_data_.back().spectrum.getPrecursors().back().setIsolationWindowLowerOffset(window);
        }
      }
      else if (tag == "scan")
      {
        skip_spectrum_ = false;
        nesting_level_++;

        if (options_.getMetadataOnly())
          throw EndParsingSoftly(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);

        // check if the scan is in the desired MS / RT range
        UInt ms_level = attributeAsInt_(attributes, s_mslevel_);
        if (ms_level == 0)
        {
          warning(LOAD, String("Invalid 'msLevel' attribute with value '0' in 'scan' element found. Assuming ms level 1!"));
          ms_level = 1;
        }

        //parse retention time and convert it from xs:duration to seconds
        double retention_time = 0.0;
        String time_string = "";
        if (optionalAttributeAsString_(time_string, attributes, s_retentiontime_))
        {
          time_string = time_string.suffix('T');
          //std::cout << "Initial trim: " << time_string << "\n";
          if (time_string.has('H'))
          {
            retention_time += 3600 * asDouble_(time_string.prefix('H'));
            time_string = time_string.suffix('H');
            //std::cout << "After H: " << time_string << "\n";
          }
          if (time_string.has('M'))
          {
            retention_time += 60 * asDouble_(time_string.prefix('M'));
            time_string = time_string.suffix('M');
            //std::cout << "After M: " << time_string << "\n";
          }
          if (time_string.has('S'))
          {
            retention_time += asDouble_(time_string.prefix('S'));
            time_string = time_string.suffix('S');
            //std::cout << "After S: " << time_string << "\n";
          }
        }

        logger_.setProgress(scan_count_);

        if ((options_.hasRTRange() && !options_.getRTRange().encloses(DPosition<1>(retention_time)))
          || (options_.hasMSLevels() && !options_.containsMSLevel(ms_level))
          || options_.getSizeOnly())
        {
          // skip this tag
          skip_spectrum_ = true;
          ++scan_count_;
          return;
        }

        // Add a new spectrum, initialize and set MS level and RT
        spectrum_data_.resize(spectrum_data_.size() + 1);
        spectrum_data_.back().peak_count_ = 0;

        spectrum_data_.back().spectrum.setMSLevel(ms_level);
        spectrum_data_.back().spectrum.setRT(retention_time);
        spectrum_data_.back().spectrum.setNativeID(String("scan=") + attributeAsString_(attributes, s_num_));
        //peak count == twice the scan size
        spectrum_data_.back().peak_count_ = attributeAsInt_(attributes, s_peakscount_);
        spectrum_data_.back().spectrum.reserve(spectrum_data_.back().peak_count_ / 2 + 1);
        spectrum_data_.back().spectrum.setDataProcessing(data_processing_);

        //centroided, chargeDeconvoluted, deisotoped, collisionEnergy are ignored

        //other optional attributes
        ScanWindow window;
        optionalAttributeAsDouble_(window.begin, attributes, s_startmz_);
        optionalAttributeAsDouble_(window.end, attributes, s_endmz_);
        if (window.begin != 0.0 || window.end != 0.0)
        {
          spectrum_data_.back().spectrum.getInstrumentSettings().getScanWindows().push_back(window);
        }

        String polarity = "any";
        optionalAttributeAsString_(polarity, attributes, s_polarity_);
        spectrum_data_.back().spectrum.getInstrumentSettings().setPolarity((IonSource::Polarity) cvStringToEnum_(0, polarity, "polarity"));

        // Filter string (see CV term MS:1000512 in mzML)
        String filterLine = "";
        optionalAttributeAsString_(filterLine, attributes, s_filterline_);
        if (!filterLine.empty())
        {
          spectrum_data_.back().spectrum.setMetaValue("filter string", filterLine);
        }

        String type = "";
        optionalAttributeAsString_(type, attributes, s_scantype_);
        if (type == "")
        {
          //unknown/unset => do nothing here => no warning in the end
        }
        else if (type == "zoom")
        {
          spectrum_data_.back().spectrum.getInstrumentSettings().setZoomScan(true);
          spectrum_data_.back().spectrum.getInstrumentSettings().setScanMode(InstrumentSettings::MASSSPECTRUM);
        }
        else if (type == "Full")
        {
          if (ms_level > 1)
            spectrum_data_.back().spectrum.getInstrumentSettings().setScanMode(InstrumentSettings::MSNSPECTRUM);
          else
            spectrum_data_.back().spectrum.getInstrumentSettings().setScanMode(InstrumentSettings::MASSSPECTRUM);
        }
        else if (type == "SIM")
        {
          spectrum_data_.back().spectrum.getInstrumentSettings().setScanMode(InstrumentSettings::SIM);
        }
        else if (type == "SRM" || type == "MRM")
        {
          spectrum_data_.back().spectrum.getInstrumentSettings().setScanMode(InstrumentSettings::SRM);
        }
        else if (type == "CRM")
        {
          spectrum_data_.back().spectrum.getInstrumentSettings().setScanMode(InstrumentSettings::CRM);
        }
        else if (type == "Q1")
        {
          spectrum_data_.back().spectrum.getInstrumentSettings().setScanMode(InstrumentSettings::MASSSPECTRUM);
        }
        else if (type == "Q3")
        {
          spectrum_data_.back().spectrum.getInstrumentSettings().setScanMode(InstrumentSettings::MASSSPECTRUM);
        }
        else if (type == "EMS") //Non-standard type: Enhanced MS (ABI - Sashimi converter)
        {
          spectrum_data_.back().spectrum.getInstrumentSettings().setScanMode(InstrumentSettings::MASSSPECTRUM);
        }
        else if (type == "EPI") //Non-standard type: Enhanced Product Ion (ABI - Sashimi converter)
        {
          spectrum_data_.back().spectrum.getInstrumentSettings().setScanMode(InstrumentSettings::MASSSPECTRUM);
          spectrum_data_.back().spectrum.setMSLevel(2);
        }
        else if (type == "ER") // Non-standard type: Enhanced Resolution (ABI - Sashimi converter)
        {
          spectrum_data_.back().spectrum.getInstrumentSettings().setZoomScan(true);
          spectrum_data_.back().spectrum.getInstrumentSettings().setScanMode(InstrumentSettings::MASSSPECTRUM);
        }
        else
        {
          spectrum_data_.back().spectrum.getInstrumentSettings().setScanMode(InstrumentSettings::MASSSPECTRUM);
          warning(LOAD, String("Unknown scan mode '") + type + "'. Assuming full scan");
        }

        ++scan_count_;
      }
      else if (tag == "operator")
      {
        exp_->getContacts().resize(1);
        exp_->getContacts().back().setFirstName(attributeAsString_(attributes, s_first_));
        exp_->getContacts().back().setLastName(attributeAsString_(attributes, s_last_));

        String tmp = "";
        optionalAttributeAsString_(tmp, attributes, s_email_);
        exp_->getContacts().back().setEmail(tmp);

        tmp = "";
        optionalAttributeAsString_(tmp, attributes, s_phone_);
        if (tmp != "")
        {
          exp_->getContacts().back().setMetaValue("#phone", tmp);
        }

        tmp = "";
        optionalAttributeAsString_(tmp, attributes, s_uri_);
        exp_->getContacts().back().setURL(tmp);
      }
      else if (tag == "msManufacturer")
      {
        exp_->getInstrument().setVendor(attributeAsString_(attributes, s_value_));
      }
      else if (tag == "msModel")
      {
        exp_->getInstrument().setModel(attributeAsString_(attributes, s_value_));
      }
      else if (tag == "msIonisation")
      {
        exp_->getInstrument().getIonSources().resize(1);
        exp_->getInstrument().getIonSources()[0].setIonizationMethod((IonSource::IonizationMethod) cvStringToEnum_(2, attributeAsString_(attributes, s_value_), "msIonization"));
      }
      else if (tag == "msMassAnalyzer")
      {
        exp_->getInstrument().getMassAnalyzers().resize(1);
        exp_->getInstrument().getMassAnalyzers()[0].setType((MassAnalyzer::AnalyzerType) cvStringToEnum_(3, attributeAsString_(attributes, s_value_), "msMassAnalyzer"));
      }
      else if (tag == "msDetector")
      {
        exp_->getInstrument().getIonDetectors().resize(1);
        exp_->getInstrument().getIonDetectors()[0].setType((IonDetector::Type) cvStringToEnum_(4, attributeAsString_(attributes, s_value_), "msDetector"));
      }
      else if (tag == "msResolution")
      {
        exp_->getInstrument().getMassAnalyzers()[0].setResolutionMethod((MassAnalyzer::ResolutionMethod) cvStringToEnum_(5, attributeAsString_(attributes, s_value_), "msResolution"));
      }
      else if (tag == "dataProcessing")
      {
        data_processing_.push_back(DataProcessingPtr(new DataProcessing));

        String boolean = "";
        optionalAttributeAsString_(boolean, attributes, s_deisotoped_);
        if (boolean == "true" || boolean == "1")
        {
          data_processing_.back()->getProcessingActions().insert(DataProcessing::DEISOTOPING);
        }

        boolean = "";
        optionalAttributeAsString_(boolean, attributes, s_chargedeconvoluted_);
        if (boolean == "true" || boolean == "1")
        {
          data_processing_.back()->getProcessingActions().insert(DataProcessing::CHARGE_DECONVOLUTION);
        }

        double cutoff = 0.0;
        optionalAttributeAsDouble_(cutoff, attributes, s_intensitycutoff_);
        if (cutoff != 0.0)
        {
          data_processing_.back()->setMetaValue("#intensity_cutoff", cutoff);
        }

        boolean = "";
        optionalAttributeAsString_(boolean, attributes, s_centroided_);
        if (boolean == "true" || boolean == "1")
        {
          data_processing_.back()->getProcessingActions().insert(DataProcessing::PEAK_PICKING);
        }
      }
      else if (tag == "nameValue")
      {
        String name = "";
        optionalAttributeAsString_(name, attributes, s_name_);
        if (name == "")
          return;

        String value = "";
        optionalAttributeAsString_(value, attributes, s_value_);

        String& parent_tag = *(open_tags_.end() - 2);

        if (parent_tag == "msInstrument")
        {
          exp_->getInstrument().setMetaValue(name, value);
        }
        else if (parent_tag == "scan")
        {
          spectrum_data_.back().spectrum.setMetaValue(name, value);
        }
        else
        {
          std::cout << " Warning: Unexpected tag 'nameValue' in tag '" << parent_tag << "'" << "\n";
        }
      }
      else if (tag == "processingOperation")
      {
        String name = "";
        optionalAttributeAsString_(name, attributes, s_name_);
        if (name == "")
          return;

        String value = "";
        optionalAttributeAsString_(value, attributes, s_value_);

        data_processing_.back()->setMetaValue(name, value);
      }

      //std::cout << " -- !Start -- " << "\n";
    }

    void MzXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
      OPENMS_PRECONDITION(nesting_level_ >= 0, "Nesting level needs to be zero or more")

        //std::cout << " -- End -- " << sm_.convert(qname) << " -- " << "\n";

        static const XMLCh* s_mzxml = xercesc::XMLString::transcode("mzXML");
      static const XMLCh* s_scan = xercesc::XMLString::transcode("scan");

      open_tags_.pop_back();

      if (equal_(qname, s_mzxml))
      {
        // Flush the remaining data
        populateSpectraWithData_();

        // End of mzXML
        logger_.endProgress();
      }
      else if (equal_(qname, s_scan))
      {
        // End of scan: go up one nesting level
        // Check whether to populate spectra when on highest nesting level
        nesting_level_--;
        OPENMS_PRECONDITION(nesting_level_ >= 0, "Nesting level needs to be zero or more")

          if (nesting_level_ == 0 && spectrum_data_.size() >= options_.getMaxDataPoolSize())
          {
          populateSpectraWithData_();
          }
      }
      //std::cout << " -- End -- " << "\n";
    }

    void MzXMLHandler::characters(const XMLCh* const chars, const XMLSize_t length)
    {
      //Abort if this spectrum should be skipped
      if (skip_spectrum_)
        return;

      if (open_tags_.back() == "peaks")
      {
        //chars may be split to several chunks => concatenate them
        if (options_.getFillData())
        {
          // Since we convert a Base64 string here, it can only contain plain ASCII
          sm_.appendASCII(chars, length, spectrum_data_.back().char_rest_);
        }
      }
      else if (open_tags_.back() == "offset" || open_tags_.back() == "indexOffset" || open_tags_.back() == "sha1")
      {

      }
      else if (open_tags_.back() == "precursorMz")
      {
        String transcoded_chars = sm_.convert(chars);
        double mz_pos = asDouble_(transcoded_chars);
        //precursor m/z
        spectrum_data_.back().spectrum.getPrecursors().back().setMZ(mz_pos);
        //update window bounds - center them around the m/z pos
        double window_width = spectrum_data_.back().spectrum.getPrecursors().back().getIsolationWindowLowerOffset();
        if (window_width != 0.0)
        {
          spectrum_data_.back().spectrum.getPrecursors().back().setIsolationWindowLowerOffset(0.5 * window_width);
          spectrum_data_.back().spectrum.getPrecursors().back().setIsolationWindowUpperOffset(0.5 * window_width);
        }
      }
      else if (open_tags_.back() == "comment")
      {
        String transcoded_chars = sm_.convert(chars);
        String parent_tag = *(open_tags_.end() - 2);
        //std::cout << "- Comment of parent " << parent_tag << "\n";

        if (parent_tag == "msInstrument")
        {
          exp_->getInstrument().setMetaValue("#comment", transcoded_chars);
        }
        else if (parent_tag == "dataProcessing")
        {
          //this is currently ignored
        }
        else if (parent_tag == "scan")
        {
          spectrum_data_.back().spectrum.setComment(transcoded_chars);
        }
        else if (transcoded_chars.trim() != "")
        {
          warning(LOAD, String("Unhandled comment '") + transcoded_chars + "' in element '" + open_tags_.back() + "'");
        }
      }
      else
      {
        String transcoded_chars = sm_.convert(chars);
        if (transcoded_chars.trim() != "")
        {
          warning(LOAD, String("Unhandled character content '") + transcoded_chars + "' in element '" + open_tags_.back() + "'");
        }
      }
    }

    void MzXMLHandler::writeTo(std::ostream& os)
    {
      //determine how many spectra there are (count only those with peaks)
      UInt count_tmp_ = 0;
      for (Size s = 0; s < cexp_->size(); s++)
      {
        const SpectrumType& spec = (*cexp_)[s];
        if (spec.size() != 0)
          ++count_tmp_;
      }
      if (count_tmp_ == 0) ++count_tmp_;

      logger_.startProgress(0, cexp_->size(), "storing mzXML file");
      double min_rt(0), max_rt(0);
      if (cexp_->size() > 0)
      {
        min_rt = cexp_->begin()->getRT();
        max_rt = (cexp_->end() - 1)->getRT();
      }
      os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
        << "<mzXML xmlns=\"http://sashimi.sourceforge.net/schema_revision/mzXML_3.1\" \n"
        << " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" \n"
        << " xsi:schemaLocation=\"http://sashimi.sourceforge.net/schema_revision/mzXML_3.1"
        << " http://sashimi.sourceforge.net/schema_revision/mzXML_3.1/mzXML_idx_3.1.xsd\">\n"
        << "\t<msRun scanCount=\"" << count_tmp_ << "\" startTime=\"PT" << min_rt << "S\" endTime=\"PT" << max_rt << "S\" >\n";

      //----------------------------------------------------------------------------------------
      // parent files
      //----------------------------------------------------------------------------------------
      if (cexp_->getSourceFiles().empty())
      {
        os << "\t\t<parentFile fileName=\"\" fileType=\"processedData\" fileSha1=\"0000000000000000000000000000000000000000\"/>\n";
      }
      else
      {
        for (Size i = 0; i < cexp_->getSourceFiles().size(); ++i)
        {
          const SourceFile& sf = cexp_->getSourceFiles()[i];
          os << "\t\t<parentFile fileName=\"" << sf.getNameOfFile() << "\" fileType=\"";
          //file type is an enum in mzXML => search for 'raw' string
          if (String(sf.getFileType()).toLower().hasSubstring("raw"))
          {
            os << "RAWData";
          }
          else
          {
            os << "processedData";
          }
          //Sha1 checksum must have 40 characters => create a fake if it is unknown
          os << "\" fileSha1=\"";
          if (sf.getChecksum().size() != 40 || sf.getChecksumType() != SourceFile::SHA1)
          {
            os << "0000000000000000000000000000000000000000";
          }
          else
          {
            os << sf.getChecksum();
          }
          os << "\"/>\n";
        }
      }

      //----------------------------------------------------------------------------------------
      //instrument
      //----------------------------------------------------------------------------------------
      if (cexp_->getInstrument() != Instrument() || cexp_->getContacts().size() != 0)
      {
        const Instrument& inst = cexp_->getInstrument();
        // the Instrument Manufacturer is paramount for some downstream tools
        // Since the .getVendor() is usually empty, we infer this via the Aquisiton Software, which is unique to Thermo
        String manufacturer = inst.getVendor();
        if (options_.getForceMQCompatability() || 
            (manufacturer.empty() && String(inst.getSoftware().getName()).toLower().hasSubstring("xcalibur")))
        { // MaxQuant's internal parameter defaults require either "Thermo Scientific" (MaxQuant 1.2 - 1.5), or "Thermo Finnigan" (MaxQuant 1.3 - 1.5)
          manufacturer = "Thermo Scientific";
          LOG_INFO << "Detected software '" << inst.getSoftware().getName() << "'. Setting <msManufacturer> as '" << manufacturer << "'." << std::endl;
        }
        os << "\t\t<msInstrument>\n"
          << "\t\t\t<msManufacturer category=\"msManufacturer\" value=\"" << manufacturer << "\"/>\n"
          << "\t\t\t<msModel category=\"msModel\" value=\"" << inst.getModel() << "\"/>\n";

        if (inst.getIonSources().empty() || !inst.getIonSources()[0].getIonizationMethod() || cv_terms_[2][inst.getIonSources()[0].getIonizationMethod()].empty())
        { // can be empty for MaxQuant
          os << "\t\t\t<msIonisation category=\"msIonisation\" value=\"\"/>\n";
        }
        else
        {
          os << "\t\t\t<msIonisation category=\"msIonisation\" value=\"" << cv_terms_[2][inst.getIonSources()[0].getIonizationMethod()] << "\"/>\n";
        }
        
        const std::vector<MassAnalyzer>& analyzers = inst.getMassAnalyzers();
        if (analyzers.empty() || cv_terms_[3][analyzers[0].getType()].empty())
        { // can be empty for MaxQuant
          os << "\t\t\t<msMassAnalyzer category=\"msMassAnalyzer\" value=\"\"/>\n";
        }
        else
        {
          os << "\t\t\t<msMassAnalyzer category=\"msMassAnalyzer\" value=\"" << cv_terms_[3][analyzers[0].getType()] << "\"/>\n";
        }

        if (inst.getIonDetectors().empty() || !inst.getIonDetectors()[0].getType() || cv_terms_[4][inst.getIonDetectors()[0].getType()].empty())
        { // can be empty for MaxQuant
          os << "\t\t\t<msDetector category=\"msDetector\" value=\"\"/>\n";
        }
        else
        {
          os << "\t\t\t<msDetector category=\"msDetector\" value=\"" << cv_terms_[4][inst.getIonDetectors()[0].getType()] << "\"/>\n";
        }
        os << "\t\t\t<software type=\"acquisition\" name=\"" << inst.getSoftware().getName() << "\" version=\"" << inst.getSoftware().getVersion() << "\"/>\n";
        if (!(analyzers.empty() || !analyzers[0].getResolutionMethod() || cv_terms_[5][analyzers[0].getResolutionMethod()].empty()))
        { // must not be empty, otherwise MaxQuant crashes upon loading mzXML
          os << "\t\t\t<msResolution category=\"msResolution\" value=\"" << cv_terms_[5][analyzers[0].getResolutionMethod()] << "\"/>\n";
        }

        if (cexp_->getContacts().size() > 0)
        {
          const ContactPerson& cont = cexp_->getContacts()[0];

          os << "\t\t\t<operator first=\"" << cont.getFirstName() << "\" last=\"" << cont.getLastName() << "\"";

          if (cont.getEmail() != "")
          {
            os << " email=\"" << cont.getEmail() << "\"";
          }

          if (cont.getURL() != "")
          {
            os << " URI=\"" << cont.getURL() << "\"";
          }

          if (cont.metaValueExists("#phone"))
          {
            os << " phone=\"" << writeXMLEscape(cont.getMetaValue("#phone").toString()) << "\"";
          }

          os << "/>\n";
        }
        writeUserParam_(os, inst, 3);

        if (inst.metaValueExists("#comment"))
        {
          os << "\t\t\t<comment>" << writeXMLEscape(inst.getMetaValue("#comment")) << "</comment>\n";
        }

        os << "\t\t</msInstrument>\n";
      }

      //----------------------------------------------------------------------------------------
      // data processing (the information of the first spectrum is assigned to the whole file)
      //----------------------------------------------------------------------------------------
      if (cexp_->size() == 0 || (*cexp_)[0].getDataProcessing().empty())
      {
        os << "\t\t<dataProcessing>\n"
          << "\t\t\t<software type=\"processing\" name=\"\" version=\"\"/>\n"
          << "\t\t</dataProcessing>\n";
      }
      else
      {
        for (Size i = 0; i < (*cexp_)[0].getDataProcessing().size(); ++i)
        {
          const DataProcessing& data_processing = *(*cexp_)[0].getDataProcessing()[i].get();
          os << "\t\t<dataProcessing deisotoped=\""
            << data_processing.getProcessingActions().count(DataProcessing::DEISOTOPING)
            << "\" chargeDeconvoluted=\""
            << data_processing.getProcessingActions().count(DataProcessing::CHARGE_DECONVOLUTION)
            << "\" centroided=\""
            << data_processing.getProcessingActions().count(DataProcessing::PEAK_PICKING)
            << "\"";
          if (data_processing.metaValueExists("#intensity_cutoff"))
          {
            os << " intensityCutoff=\"" << writeXMLEscape(data_processing.getMetaValue("#intensity_cutoff").toString()) << "\"";
          }
          os << ">\n"
            << "\t\t\t<software type=\"";
          if (data_processing.metaValueExists("#type"))
          {
            os << writeXMLEscape(data_processing.getMetaValue("#type").toString());
          }
          else
          {
            os << "processing";
          }

          os << "\" name=\"" << data_processing.getSoftware().getName()
            << "\" version=\"" << data_processing.getSoftware().getVersion();

          if (data_processing.getCompletionTime() != DateTime())
          {
            os << "\" completionTime=\"" << data_processing.getCompletionTime().get().substitute(' ', 'T');
          }
          os << "\"/>\n";
          writeUserParam_(os, data_processing, 3, "processingOperation");

          os << "\t\t</dataProcessing>\n";
        }
      }

      //check if the nativeID of all spectra are numbers or numbers prefixed with 'scan='
      //If not we need to renumber all spectra.
      bool all_numbers = true;
      bool all_empty = true;
      bool all_prefixed_numbers = true;
      for (Size s = 0; s < cexp_->size(); s++)
      {
        String native_id = (*cexp_)[s].getNativeID();
        if (!native_id.hasPrefix("scan="))
        {
          all_prefixed_numbers = false;
        }
        else
        {
          native_id = native_id.substr(5);
        }
        try
        {
          native_id.toInt();
        }
        catch (Exception::ConversionError&)
        {
          all_numbers = false;
          all_prefixed_numbers = false;
          if (native_id != "")
          {
            all_empty = false;
          }
        }
      }
      //If we need to renumber and the nativeIDs were not empty, warn the user
      if (!all_numbers && !all_empty)
      {
        warning(STORE, "Not all spectrum native IDs are numbers or correctly prefixed with 'scan='. The spectra are renumbered and the native IDs are lost!");
      }

      std::vector<IndexPos> scan_index_positions;

      // write scans
      std::stack<UInt> open_scans;
      for (Size s = 0; s < cexp_->size(); s++)
      {
        logger_.setProgress(s);
        const SpectrumType& spec = (*cexp_)[s];

        UInt ms_level = spec.getMSLevel();
        open_scans.push(ms_level);

        Size spectrum_id = s + 1;
        if (all_prefixed_numbers)
        {
          spectrum_id = spec.getNativeID().substr(5).toInt();
        }
        else if (all_numbers)
        {
          spectrum_id = spec.getNativeID().toInt();
        }

        os << String(ms_level + 1, '\t');

        scan_index_positions.push_back(IndexPos(spectrum_id, os.tellp())); // remember scan index
        os << "<scan num=\"" << spectrum_id << "\""
          << " msLevel=\"" << ms_level << "\""
          << " peaksCount=\"" << spec.size() << "\""
          << " polarity=\"";
        if (spec.getInstrumentSettings().getPolarity() == IonSource::POSITIVE)
        {
          os << "+";
        }
        else if (spec.getInstrumentSettings().getPolarity() == IonSource::NEGATIVE)
        {
          os << "-";
        }
        else
        {
          os << "any";
        }
        os << "\"";
        //scan type
        switch (spec.getInstrumentSettings().getScanMode())
        {
        case InstrumentSettings::UNKNOWN:
          break;

        case InstrumentSettings::MASSSPECTRUM:
        case InstrumentSettings::MS1SPECTRUM:
        case InstrumentSettings::MSNSPECTRUM:
          if (spec.getInstrumentSettings().getZoomScan())
          {
            os << " scanType=\"zoom\"";
          }
          else
          {
            os << " scanType=\"Full\"";
          }
          break;

        case InstrumentSettings::SIM:
          os << " scanType=\"SIM\"";
          break;

        case InstrumentSettings::SRM:
          os << " scanType=\"SRM\"";
          break;

        case InstrumentSettings::CRM:
          os << " scanType=\"CRM\"";
          break;

        default:
          os << " scanType=\"Full\"";
          warning(STORE, String("Scan type '") + InstrumentSettings::NamesOfScanMode[spec.getInstrumentSettings().getScanMode()] + "' not supported by mzXML. Using 'Full' scan mode!");
        }

        // filter line
        if (spec.metaValueExists("filter string"))
        {
          os << " filterLine=\"";
          os << writeXMLEscape((String)spec.getMetaValue("filter string"));
          os << "\"";
        }

        // retention time
        os << " retentionTime=\"";
        if (spec.getRT() < 0) os << "-";
        os << "PT" << std::fabs(spec.getRT()) << "S\"";
        if (!spec.getInstrumentSettings().getScanWindows().empty())
        {
          os << " startMz=\"" << spec.getInstrumentSettings().getScanWindows()[0].begin << "\" endMz=\"" << spec.getInstrumentSettings().getScanWindows()[0].end << "\"";
        }
        if (spec.getInstrumentSettings().getScanWindows().size() > 1)
        {
          warning(STORE, "The MzXML format can store only one scan window for each scan. Only the first one is stored!");
        }

        // convert meta values to tags
        writeAttributeIfExists_(os, spec, "lowest observed m/z", "lowMz");
        writeAttributeIfExists_(os, spec, "highest observed m/z", "highMz");
        if (spec.metaValueExists("base peak m/z")) writeAttributeIfExists_(os, spec, "base peak m/z", "basePeakMz");
        else
        { // base peak mz (used by some programs like MAVEN), according to xsd: "m/z of the base peak (most intense peak)"
          os << " basePeakMz=\"";
          double basePeakInt(0), basePeakMz(0);
          for (Size j = 0; j < spec.size(); j++)
          {
            if (spec[j].getIntensity() > basePeakInt)
            {
              basePeakInt = spec[j].getIntensity();
              basePeakMz = spec[j].getMZ();
            }
          }
          os << basePeakMz << "\"";
        }

        writeAttributeIfExists_(os, spec, "base peak intensity", "basePeakIntensity");
        writeAttributeIfExists_(os, spec, "total ion current", "totIonCurrent");

        if (ms_level == 2 &&
          !spec.getPrecursors().empty() &&
          spec.getPrecursors().front().metaValueExists("collision energy"))
        {
          os << " collisionEnergy=\"" << spec.getPrecursors().front().getMetaValue("collision energy") << "\" ";
        }
        // end of "scan" attributes
        os << ">\n";


        for (Size i = 0; i < spec.getPrecursors().size(); ++i)
        {
          const Precursor& precursor = spec.getPrecursors()[i];
          // intensity
          os << String(ms_level + 2, '\t') << "<precursorMz precursorIntensity=\"" << (int)precursor.getIntensity() << "\"";
          // charge
          if (precursor.getCharge() != 0)
          {
            os << " precursorCharge=\"" << precursor.getCharge() << "\"";
          }
          // window size
          if (precursor.getIsolationWindowLowerOffset() + precursor.getIsolationWindowUpperOffset() > 0.0)
          {
            os << " windowWideness=\"" << (precursor.getIsolationWindowUpperOffset() + precursor.getIsolationWindowLowerOffset()) << "\"";
          }
          if (!precursor.getActivationMethods().empty() && !Precursor::NamesOfActivationMethodShort[int(*(precursor.getActivationMethods().begin()))].empty())
          { // must not be empty, but technically only ETD, ECD, CID are allowed in mzXML 3.1
            os << " activationMethod=\"" << Precursor::NamesOfActivationMethodShort[int(*(precursor.getActivationMethods().begin()))] << "\" ";
          }

          //m/z
          double mz = precursor.getMZ();
          if (!spec.getAcquisitionInfo().empty() &&
            spec.getAcquisitionInfo().begin()->metaValueExists("[Thermo Trailer Extra]Monoisotopic M/Z:"))
          { // this value is usually more accurate; the old ReAdw-converter uses it as well
            mz = spec.getAcquisitionInfo().begin()->getMetaValue("[Thermo Trailer Extra]Monoisotopic M/Z:");
          }
          os << ">" << mz << "</precursorMz>\n";
        }

        // Note: Some parsers require the following line breaks (MaxQuants
        // mzXML reader will fail otherwise! -- dont ask..) while others cannot
        // deal with them (mostly TPP tools such as SpectraST).
        String s_peaks;
        if (options_.getForceMQCompatability() )
        {
          s_peaks = "<peaks precision=\"32\"\n byteOrder=\"network\"\n contentType=\"m/z-int\"\n compressionType=\"none\"\n compressedLen=\"0\" ";
        }
        else
        {
          s_peaks = "<peaks precision=\"32\" byteOrder=\"network\" contentType=\"m/z-int\" compressionType=\"none\" compressedLen=\"0\" ";
        }
        if (options_.getForceMQCompatability() && !s_peaks.has('\n'))
        { // internal check against inadvertedly removing line breaks above!
          fatalError(STORE, "Internal error: <peaks> tag does not contain newlines as required by MaxQuant. Please report this as a bug.", __LINE__, 0);
        }
        os << String(ms_level + 2, '\t') << s_peaks;
        if (!spec.empty())
        {
          os << ">";
          // for MaxQuant-compatible mzXML, the data type must be 'float', i.e. precision=32 bit. 64bit will crash MaxQuant!
          std::vector<float> tmp;
          for (Size i = 0; i < spec.size(); i++)
          {
            tmp.push_back(spec[i].getMZ());
            tmp.push_back(spec[i].getIntensity());
          }

          String encoded;
          decoder_.encode(tmp, Base64::BYTEORDER_BIGENDIAN, encoded);
          os << encoded << "</peaks>\n";
        }
        else
        {
          os << " xsi:nil=\"true\" />\n";
        }

        writeUserParam_(os, spec, ms_level + 2);
        if (spec.getComment() != "")
        {
          os << String(ms_level + 2, '\t') << "<comment>" << spec.getComment() << "</comment>\n";
        }

        //check MS level of next scan and close scans (scans can be nested)
        UInt next_ms_level = 0;
        if (s < cexp_->size() - 1)
        {
          next_ms_level = ((*cexp_)[s + 1]).getMSLevel();
        }
        //std::cout << "scan: " << s << " this: " << ms_level << " next: " << next_ms_level << "\n";
        if (next_ms_level <= ms_level)
        {
          for (Size i = 0; i <= ms_level - next_ms_level && !open_scans.empty(); ++i)
          {
            os << String(ms_level - i + 1, '\t') << "</scan>\n";
            open_scans.pop();
          }
        }
      }

      os << "\t</msRun>\n";

      if (options_.getWriteIndex() || options_.getForceMQCompatability())
      { // create scan index (does not take a lot of space and is required for MaxQuant)
        if (!options_.getWriteIndex())
        {
          LOG_INFO << "mzXML: index was not requested, but will be written to maintain MaxQuant compatibility." << std::endl;
        }
        std::ostream::streampos index_offset = os.tellp();
        os << "<index name = \"scan\" >\n";
        for (Size i = 0; i < scan_index_positions.size(); ++i)
        {
          os << "<offset id = \"" << scan_index_positions[i].id_ << "\" >" << scan_index_positions[i].pos_ << "</offset>\n";
        }
        os << "</index>\n";
        os << "<indexOffset>" << index_offset << "</indexOffset>\n";
      }

      os << "</mzXML>\n";

      logger_.endProgress();
      spec_write_counter_ = 1;
    }

    void MzXMLHandler::initStaticMembers_()
    {
      static bool init(false);
      if (!init)
      {
        s_value_ = xercesc::XMLString::transcode("value");
        s_count_ = xercesc::XMLString::transcode("scanCount");
        s_type_ = xercesc::XMLString::transcode("type");
        s_name_ = xercesc::XMLString::transcode("name");
        s_version_ = xercesc::XMLString::transcode("version");
        s_filename_ = xercesc::XMLString::transcode("fileName");
        s_filetype_ = xercesc::XMLString::transcode("fileType");
        s_filesha1_ = xercesc::XMLString::transcode("fileSha1");
        s_completiontime_ = xercesc::XMLString::transcode("completionTime");
        s_precision_ = xercesc::XMLString::transcode("precision");
        s_byteorder_ = xercesc::XMLString::transcode("byteOrder");
        s_contentType_ = xercesc::XMLString::transcode("contentType");
        s_compressionType_ = xercesc::XMLString::transcode("compressionType");
        s_precursorintensity_ = xercesc::XMLString::transcode("precursorIntensity");
        s_precursorcharge_ = xercesc::XMLString::transcode("precursorCharge");
        s_windowwideness_ = xercesc::XMLString::transcode("windowWideness");
        s_mslevel_ = xercesc::XMLString::transcode("msLevel");
        s_peakscount_ = xercesc::XMLString::transcode("peaksCount");
        s_polarity_ = xercesc::XMLString::transcode("polarity");
        s_scantype_ = xercesc::XMLString::transcode("scanType");
        s_filterline_ = xercesc::XMLString::transcode("filterLine");
        s_retentiontime_ = xercesc::XMLString::transcode("retentionTime");
        s_startmz_ = xercesc::XMLString::transcode("startMz");
        s_endmz_ = xercesc::XMLString::transcode("endMz");
        s_first_ = xercesc::XMLString::transcode("first");
        s_last_ = xercesc::XMLString::transcode("last");
        s_phone_ = xercesc::XMLString::transcode("phone");
        s_email_ = xercesc::XMLString::transcode("email");
        s_uri_ = xercesc::XMLString::transcode("URI");
        s_num_ = xercesc::XMLString::transcode("num");
        s_intensitycutoff_ = xercesc::XMLString::transcode("intensityCutoff");
        s_centroided_ = xercesc::XMLString::transcode("centroided");
        s_deisotoped_ = xercesc::XMLString::transcode("deisotoped");
        s_chargedeconvoluted_ = xercesc::XMLString::transcode("chargeDeconvoluted");

        init = true;
      }
      return;
    }

    inline std::ostream& MzXMLHandler::writeAttributeIfExists_(std::ostream& os, const MetaInfoInterface& meta, const String& metakey, const String& attname)
    {
      if (meta.metaValueExists(metakey))
      {
        os << " " << attname << "=\"" << meta.getMetaValue(metakey) << "\"";
      }
      return os;
    }

    inline void MzXMLHandler::writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, int indent, String tag)
    {
      std::vector<String> keys; // Vector to hold keys to meta info
      meta.getKeys(keys);

      for (std::vector<String>::const_iterator it = keys.begin(); it != keys.end(); ++it)
      {
        if ((*it)[0] != '#') // internally used meta info start with '#'
        {
          os << String(indent, '\t') << "<" << tag << " name=\"" << *it << "\" value=\"" << writeXMLEscape(meta.getMetaValue(*it)) << "\"/>\n";
        }
      }
    }

    void MzXMLHandler::doPopulateSpectraWithData_(SpectrumData & spectrum_data)
    {
      typedef SpectrumType::PeakType PeakType;

      //std::cout << "reading scan" << "\n";
      if (spectrum_data.char_rest_ == "") // no peaks
      {
        return;
      }

      //remove whitespaces from binary data
      //this should not be necessary, but linebreaks inside the base64 data are unfortunately no exception
      spectrum_data.char_rest_.removeWhitespaces();

      if (spectrum_data.precision_ == "64")
      {
        std::vector<double> data;
        if (spectrum_data.compressionType_ == "zlib")
        {
          decoder_.decode(spectrum_data.char_rest_, Base64::BYTEORDER_BIGENDIAN, data, true);
        }
        else
        {
          decoder_.decode(spectrum_data.char_rest_, Base64::BYTEORDER_BIGENDIAN, data);
        }
        spectrum_data.char_rest_ = "";
        PeakType peak;
        //push_back the peaks into the container
        for (Size n = 0; n < (2 * spectrum_data.peak_count_); n += 2)
        {
          // check if peak in in the specified m/z  and intensity range
          if ((!options_.hasMZRange() || options_.getMZRange().encloses(DPosition<1>(data[n])))
             && (!options_.hasIntensityRange() || options_.getIntensityRange().encloses(DPosition<1>(data[n + 1]))))
          {
            peak.setMZ(data[n]);
            peak.setIntensity(data[n + 1]);
            spectrum_data.spectrum.push_back(peak);
          }
        }
      }
      else //precision 32
      {
        std::vector<float> data;
        if (spectrum_data.compressionType_ == "zlib")
        {
          decoder_.decode(spectrum_data.char_rest_, Base64::BYTEORDER_BIGENDIAN, data, true);
        }
        else
        {
          decoder_.decode(spectrum_data.char_rest_, Base64::BYTEORDER_BIGENDIAN, data);
        }
        spectrum_data.char_rest_ = "";
        PeakType peak;
        //push_back the peaks into the container
        for (Size n = 0; n < (2 * spectrum_data.peak_count_); n += 2)
        {
          if ((!options_.hasMZRange() || options_.getMZRange().encloses(DPosition<1>(data[n])))
             && (!options_.hasIntensityRange() || options_.getIntensityRange().encloses(DPosition<1>(data[n + 1]))))
          {
            peak.setMZ(data[n]);
            peak.setIntensity(data[n + 1]);
            spectrum_data.spectrum.push_back(peak);
          }
        }
      }
    }

    void MzXMLHandler::populateSpectraWithData_()
    {

      // Whether spectrum should be populated with data
      if (options_.getFillData())
      {
        size_t errCount = 0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (SignedSize i = 0; i < (SignedSize)spectrum_data_.size(); i++)
        {
          // parallel exception catching and re-throwing business
          if (!errCount) // no need to parse further if already an error was encountered
          {
            try
            {
              doPopulateSpectraWithData_(spectrum_data_[i]);
              if (options_.getSortSpectraByMZ() && !spectrum_data_[i].spectrum.isSorted())
              {
                spectrum_data_[i].spectrum.sortByPosition();
              }
            }
            catch (...)
            {
              #pragma omp critical(HandleException)
              ++errCount;
            }
          }
        }
        if (errCount != 0)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file_, "Error during parsing of binary data.");
        }
      }

      // Append all spectra
      for (Size i = 0; i < spectrum_data_.size(); i++)
      {
        if (consumer_ != NULL)
        {
          consumer_->consumeSpectrum(spectrum_data_[i].spectrum);
          if (options_.getAlwaysAppendData())
          {
            exp_->addSpectrum(spectrum_data_[i].spectrum);
          }
        }
        else
        {
          exp_->addSpectrum(spectrum_data_[i].spectrum);
        }
      }

      // Delete batch
      spectrum_data_.clear();
    }

    void MzXMLHandler::init_()
    {
      cv_terms_.resize(6);
      //Polarity
      String("any;+;-").split(';', cv_terms_[0]);
      //Scan type
      // is no longer used cv_terms_[1] is empty now

      //Ionization method
      String(";ESI;EI;CI;FAB;;;;;;;;;;;;;APCI;;;NSI;;SELDI;;;MALDI").split(';', cv_terms_[2]);
      cv_terms_[2].resize(IonSource::SIZE_OF_IONIZATIONMETHOD);

      //Mass analyzer
      String(";Quadrupole;Quadrupole Ion Trap;;;TOF;Magnetic Sector;FT-ICR;;;;;;FTMS").split(';', cv_terms_[3]);
      cv_terms_[3].resize(MassAnalyzer::SIZE_OF_ANALYZERTYPE);

      //Detector
      String(";EMT;;;Faraday Cup;;;;;Channeltron;Daly;Microchannel plate").split(';', cv_terms_[4]);
      cv_terms_[4].resize(IonDetector::SIZE_OF_TYPE);

      //Resolution method
      String(";FWHM;TenPercentValley;Baseline").split(';', cv_terms_[5]);
      cv_terms_[5].resize(MassAnalyzer::SIZE_OF_RESOLUTIONMETHOD);
    }
  }
}
