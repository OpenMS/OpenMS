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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/MSDataSqlConsumer.h>

#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h>

namespace OpenMS
{

  MSDataSqlConsumer::MSDataSqlConsumer(String filename, int flush_after, bool full_meta, bool lossy_compression, double linear_mass_acc) :
        filename_(filename),
        handler_(new OpenMS::Internal::MzMLSqliteHandler(filename) ),
        flush_after_(flush_after),
        full_meta_(full_meta)
  {
    spectra_.reserve(flush_after_);
    chromatograms_.reserve(flush_after_);

    handler_->setConfig(full_meta, lossy_compression, linear_mass_acc, flush_after_);
    handler_->createTables();
  }

  MSDataSqlConsumer::~MSDataSqlConsumer()
  {
    flush();

    // Write run level information into the file (e.g. run id, run name and mzML structure)
    bool write_full_meta = full_meta_;
    int run_id = 0;
    peak_meta_.setLoadedFilePath(filename_);
    handler_->writeRunLevelInformation(peak_meta_, write_full_meta, run_id);

    delete handler_;
  }

  void MSDataSqlConsumer::flush()
  {
    if (!spectra_.empty() ) 
    {
      handler_->writeSpectra(spectra_);
      spectra_.clear();
      spectra_.reserve(flush_after_);
    }

    if (!chromatograms_.empty() ) 
    {
      handler_->writeChromatograms(chromatograms_);
      chromatograms_.clear();
      chromatograms_.reserve(flush_after_);
    }
  }

  void MSDataSqlConsumer::consumeSpectrum(SpectrumType & s)
  {
    spectra_.push_back(s);
    s.clear(false);
    if (full_meta_) peak_meta_.addSpectrum(s);

    if (spectra_.size() >= flush_after_) {flush();}
  }

  void MSDataSqlConsumer::consumeChromatogram(ChromatogramType & c)
  {
    chromatograms_.push_back(c);
    c.clear(false);
    if (full_meta_) peak_meta_.addChromatogram(c);

    if (chromatograms_.size() >= flush_after_) {flush();}
  }

  void MSDataSqlConsumer::setExpectedSize(Size /* expectedSpectra */, Size /* expectedChromatograms */) {;}

  void MSDataSqlConsumer::setExperimentalSettings(const ExperimentalSettings& /* exp */) {;}

} // namespace OpenMS

