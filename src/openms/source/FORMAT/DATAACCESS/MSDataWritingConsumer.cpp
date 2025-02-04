// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/VALIDATORS/MzMLValidator.h>
// TODO move getVersion to Handler
#include <OpenMS/FORMAT/MzMLFile.h>

#include <utility>

namespace OpenMS
{

  MSDataWritingConsumer::MSDataWritingConsumer(const String& filename) :
    Internal::MzMLHandler(MapType(), filename, MzMLFile().getVersion(), ProgressLogger()),
    started_writing_(false),
    writing_spectra_(false),
    writing_chromatograms_(false),
    spectra_written_(0),
    chromatograms_written_(0),
    spectra_expected_(0),
    chromatograms_expected_(0),
    add_dataprocessing_(false)
  {
    validator_ = new Internal::MzMLValidator(this->mapping_, this->cv_);

    // open file in binary mode to avoid any line ending conversions
    ofs_.open(filename.c_str(), std::ios::out | std::ios::binary);
    ofs_.precision(writtenDigits(double()));
  }

   MSDataWritingConsumer::~MSDataWritingConsumer()
  {
    doCleanup_();
  }

   void MSDataWritingConsumer::setExperimentalSettings(const ExperimentalSettings& exp)
  {
    settings_ = exp;
  }

   void MSDataWritingConsumer::setExpectedSize(Size expectedSpectra, Size expectedChromatograms)
  {
    spectra_expected_ = expectedSpectra;
    chromatograms_expected_ = expectedChromatograms;
  }

   void MSDataWritingConsumer::consumeSpectrum(SpectrumType & s)
  {
    if (writing_chromatograms_)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Cannot write spectra after writing chromatograms.");
    }

    // Process the spectrum 
    SpectrumType scpy = s;
    processSpectrum_(scpy);

    // Add dataprocessing if required
    if (add_dataprocessing_)
    {
      scpy.getDataProcessing().push_back(additional_dataprocessing_);
    }

    if (!started_writing_)
    {
      // This is the first data to be written -> start writing the header
      // We also need to modify the map and add this dummy spectrum in
      // order to write the header correctly
      MapType dummy;
      dummy = settings_;
      dummy.addSpectrum(scpy);

      //--------------------------------------------------------------------
      //header
      //--------------------------------------------------------------------
      Internal::MzMLHandler::writeHeader_(ofs_, dummy, dps_, *validator_);
      started_writing_ = true;
    }
    if (!writing_spectra_)
    {
      // This is the first spectrum, thus write the spectrumList header
      ofs_ << "\t\t<spectrumList count=\"" << spectra_expected_ << "\" defaultDataProcessingRef=\"dp_sp_0\">\n";
      writing_spectra_ = true;
    }
    bool renew_native_ids = false;
    // TODO writeSpectrum assumes that dps_ has at least one value -> assert
    // this here ...
    Internal::MzMLHandler::writeSpectrum_(ofs_, scpy,
            spectra_written_++, *validator_, renew_native_ids, dps_);
  }

   void MSDataWritingConsumer::consumeChromatogram(ChromatogramType & c)
  {
    // make sure to close an open List tag
    if (writing_spectra_)
    {
      ofs_ << "\t\t</spectrumList>\n";
      writing_spectra_ = false;
    }

    // Create copy and add dataprocessing if required
    ChromatogramType ccpy = c;
    processChromatogram_(ccpy);

    if (add_dataprocessing_)
    {
      ccpy.getDataProcessing().push_back(additional_dataprocessing_);
    }

    if (!started_writing_)
    {
      // this is the first data to be written -> start writing the header
      // We also need to modify the map and add this dummy chromatogram in
      // order to write the header correctly
      MapType dummy;
      dummy = settings_;
      dummy.addChromatogram(ccpy);

      //--------------------------------------------------------------------
      //header (fill also dps_ variable)
      //--------------------------------------------------------------------
      Internal::MzMLHandler::writeHeader_(ofs_, dummy, dps_, *validator_);
      started_writing_ = true;
    }
    if (!writing_chromatograms_)
    {
      ofs_ << "\t\t<chromatogramList count=\"" << chromatograms_expected_ << "\" defaultDataProcessingRef=\"dp_sp_0\">\n";
      writing_chromatograms_ = true;
    }
    Internal::MzMLHandler::writeChromatogram_(ofs_, ccpy,
            chromatograms_written_++, *validator_);
  }

   void MSDataWritingConsumer::addDataProcessing(DataProcessing d)
  {
    additional_dataprocessing_ = DataProcessingPtr( new DataProcessing(std::move(d)) );
    add_dataprocessing_ = true;
  }

   Size MSDataWritingConsumer::getNrSpectraWritten() {return spectra_written_;}

   Size MSDataWritingConsumer::getNrChromatogramsWritten() {return chromatograms_written_;}

   void MSDataWritingConsumer::doCleanup_()
  {
    //--------------------------------------------------------------------------------------------
    //cleanup
    //--------------------------------------------------------------------------------------------
    // make sure to close an open List tag
    if (writing_spectra_)
    {
      ofs_ << "\t\t</spectrumList>\n";
    }
    else if (writing_chromatograms_)
    {
      ofs_ << "\t\t</chromatogramList>\n";
    }

    // Only write the footer if we actually did start writing ... 
    if (started_writing_) 
    {
      Internal::MzMLHandlerHelper::writeFooter_(ofs_, options_, spectra_offsets_, chromatograms_offsets_);
    }
    delete validator_;
    ofs_.close();
  }

} // namespace OpenMS
