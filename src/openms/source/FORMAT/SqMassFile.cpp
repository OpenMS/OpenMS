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

#include <OpenMS/FORMAT/SqMassFile.h>

#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h>

namespace OpenMS
{

  SqMassFile::SqMassFile() {}

  SqMassFile::~SqMassFile() {}

  void SqMassFile::load(const String& filename, MapType& map)
  {
    OpenMS::Internal::MzMLSqliteHandler sql_mass(filename);
    sql_mass.setConfig(config_.write_full_meta, config_.use_lossy_numpress, config_.linear_fp_mass_acc);
    sql_mass.readExperiment(map);
  }

  void SqMassFile::store(const String& filename, MapType& map)
  {
    OpenMS::Internal::MzMLSqliteHandler sql_mass(filename);
    sql_mass.setConfig(config_.write_full_meta, config_.use_lossy_numpress, config_.linear_fp_mass_acc);
    sql_mass.createTables();
    sql_mass.writeExperiment(map);
  }

  void SqMassFile::transform(const String& filename_in, Interfaces::IMSDataConsumer* consumer, bool /* skip_full_count */, bool /* skip_first_pass */)
  {
    OpenMS::Internal::MzMLSqliteHandler sql_mass(filename_in);
    sql_mass.setConfig(config_.write_full_meta, config_.use_lossy_numpress, config_.linear_fp_mass_acc);

    // First pass through the file -> get the meta-data and hand it to the consumer
    // if (!skip_first_pass) transformFirstPass_(filename_in, consumer, skip_full_count);
    consumer->setExpectedSize(sql_mass.getNrSpectra(), sql_mass.getNrChromatograms());
    MSExperiment experimental_settings;
    sql_mass.readExperiment(experimental_settings, true);
    consumer->setExperimentalSettings(experimental_settings);

    {
      int batch_size = 500;
      std::vector<int> indices;
      for (size_t batch_idx = 0; batch_idx <= (sql_mass.getNrSpectra() / batch_size); batch_idx++)
      {
        int idx_start, idx_end;
        idx_start = batch_idx * batch_size;
        idx_end = std::max(batch_idx * (batch_size+1), sql_mass.getNrSpectra());

        indices.resize(idx_end - idx_start);
        for (int k = 0; k < idx_end-idx_start; k++)
        {
          indices[k] = idx_start + k;
        }
        std::vector<MSSpectrum> tmp_spectra;
        sql_mass.readSpectra(tmp_spectra, indices, false);
        for (Size k = 0; k < tmp_spectra.size(); k++)
        {
          consumer->consumeSpectrum(tmp_spectra[k]);
        }
      }
    }

    {
      int batch_size = 500;
      std::vector<int> indices;
      for (size_t batch_idx = 0; batch_idx <= (sql_mass.getNrChromatograms() / batch_size); batch_idx++)
      {
        int idx_start, idx_end;
        idx_start = batch_idx * batch_size;
        idx_end = std::max(batch_idx * (batch_size+1), sql_mass.getNrChromatograms());

        indices.resize(idx_end - idx_start);
        for (int k = 0; k < idx_end-idx_start; k++)
        {
          indices[k] = idx_start + k;
        }
        std::vector<MSChromatogram> tmp_chroms;
        sql_mass.readChromatograms(tmp_chroms, indices, false);
        for (Size k = 0; k < tmp_chroms.size(); k++)
        {
          consumer->consumeChromatogram(tmp_chroms[k]);
        }
      }
    }
  }

}
