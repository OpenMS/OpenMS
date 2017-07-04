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

namespace OpenMS
{

  MSDataSqlConsumer::MSDataSqlConsumer(String filename, bool clearData, int flush_after) :
        sql_writer_(filename),
        clearData_(clearData),
        flush_after_(flush_after)
      {
        sql_writer_.createTables();
      }

  MSDataSqlConsumer::~MSDataSqlConsumer()
  {
    flush();
  }

  void MSDataSqlConsumer::flush()
  {
    sql_writer_.writeSpectra(spectra_);
    spectra_.clear();
    sql_writer_.writeChromatograms(chromatograms_);
    chromatograms_.clear();
  }

  void MSDataSqlConsumer::consumeSpectrum(SpectrumType & s)
  {
    spectra_.push_back(s);
    if (spectra_.size() >= flush_after_)
    {
      sql_writer_.writeSpectra(spectra_);
      spectra_.clear();
    }
    if (clearData_) {s.clear(false);}
  }

  void MSDataSqlConsumer::consumeChromatogram(ChromatogramType & c)
  {
    chromatograms_.push_back(c);
    if (chromatograms_.size() >= flush_after_)
    {
      sql_writer_.writeChromatograms(chromatograms_);
      chromatograms_.clear();
    }
    if (clearData_) {c.clear(false);}
  }

  void MSDataSqlConsumer::setExpectedSize(Size /* expectedSpectra */, Size /* expectedChromatograms */) {;}

  void MSDataSqlConsumer::setExperimentalSettings(const ExperimentalSettings& /* exp */) {;}



} // namespace OpenMS

