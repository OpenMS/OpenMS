// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/SqMassFile.h>

#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h>

namespace OpenMS
{

  SqMassFile::SqMassFile() = default;

  SqMassFile::~SqMassFile() = default;

  void SqMassFile::load(const String& filename, MapType& map) const
  {
    OpenMS::Internal::MzMLSqliteHandler sql_mass(filename, 0);
    sql_mass.setConfig(config_.write_full_meta, config_.use_lossy_numpress, config_.linear_fp_mass_acc);
    sql_mass.readExperiment(map);
  }

  void SqMassFile::store(const String& filename, const MapType& map) const
  {
    OpenMS::Internal::MzMLSqliteHandler sql_mass(filename, map.getSqlRunID());
    sql_mass.setConfig(config_.write_full_meta, config_.use_lossy_numpress, config_.linear_fp_mass_acc);
    sql_mass.createTables();
    sql_mass.writeExperiment(map);
  }

  void SqMassFile::transform(const String& filename_in, Interfaces::IMSDataConsumer* consumer, bool /* skip_full_count */, bool /* skip_first_pass */) const
  {
    OpenMS::Internal::MzMLSqliteHandler sql_mass(filename_in, 0);
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
