// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/VALIDATORS/MzMLValidator.h>

namespace OpenMS
{

  MSDataWritingConsumer::MSDataWritingConsumer(String filename) :
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
    additional_dataprocessing_ = DataProcessingPtr( new DataProcessing(d) );
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
      Internal::MzMLHandlerHelper::writeFooter_(ofs_, options_, spectra_offsets_, chromatograms_offsets_);

    delete validator_;
    ofs_.close();
  }

} // namespace OpenMS
