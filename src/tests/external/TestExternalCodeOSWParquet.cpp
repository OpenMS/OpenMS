// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWParquetWriter.h>
#include <OpenMS/FORMAT/OSWParquetFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/TempFiles.h>

#if __has_include(<OpenMS/FORMAT/ParquetFile.h>)
  #error "The installed consumer can see the private ParquetFile header"
#endif

int main()
{
  OpenMS::TempDir temporary;
  const std::string path = temporary.getPath() + "/workflow";
  OpenMS::File::makeDir(path);
  OpenSwath::LightTargetedExperiment library;
  OpenSwath::LightCompound compound;
  compound.id = "0";
  compound.sequence = "PEPTIDE";
  compound.charge = 2;
  library.compounds.push_back(compound);
  OpenSwath::LightTransition transition;
  transition.transition_name = "3";
  transition.peptide_ref = "0";
  transition.precursor_mz = 500.0;
  transition.product_mz = 200.0;
  library.transitions.push_back(transition);
  OpenMS::Feature feature;
  feature.setUniqueId(1);
  feature.setRT(10.0);
  feature.setIntensity(100.0);
  feature.setMetaValue("PeptideRef", "0");
  OpenMS::FeatureMap features;
  features.push_back(feature);
  OpenMS::OpenSwathOSWParquetWriter().write(path, library, features, 11, "sample.mzML", false);
  OpenMS::OSWParquetFile file(path, library);
  if (file.readRunBasenames().at(11) != "sample") return 1;
  const auto table = file.readOpenSwathFeatureScoreTable(OpenMS::OpenSwathParquetExportConfig {});
  if (table.rows.size() != 1 || table.rows.front().precursor_id != 0) return 2;
  file.writeLevelContextResults({});
  file.commit();
  return 0;
}
