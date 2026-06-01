// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/PSMArrowIO.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/FORMAT/QPXFile.h>
#include <OpenMS/FORMAT/ProteinIdentificationArrowIO.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>

namespace OpenMS
{

namespace
{
  constexpr const char* kPSMs = "psms.parquet";
  constexpr const char* kProteins = "proteins.parquet";
  constexpr const char* kProteinGroups = "protein_groups.parquet";
  constexpr const char* kSearchParams = "search_params.parquet";

  bool ensureDirectory_(const String& dir)
  {
    if (File::exists(dir))
    {
      if (!File::isDirectory(dir))
      {
        OPENMS_LOG_ERROR << "PSMArrowIO: target path exists and is not a directory: " << dir << std::endl;
        return false;
      }
      return true;
    }
    if (!File::makeDir(dir))
    {
      OPENMS_LOG_ERROR << "PSMArrowIO: failed to create directory: " << dir << std::endl;
      return false;
    }
    return true;
  }
}

bool PSMArrowIO::exportToParquet(
  const std::vector<ProteinIdentification>& protein_identifications,
  const PeptideIdentificationList& peptide_identifications,
  const String& dir,
  bool export_all_psms,
  const ParquetWriteConfig& config)
{
  // Mirror XMLHandler::checkUniqueIdentifiers_ — fail before any file is opened
  // so we never leave a partial .idparquet behind.
  ProteinIdentificationArrowIO::checkUniqueIdentifiers(protein_identifications);

  if (!ensureDirectory_(dir)) { return false; }

  // PSMs (PSMSchema, lossless)
  auto psm_table = QPXFile::exportToArrow(
    protein_identifications, peptide_identifications, export_all_psms);
  if (!psm_table)
  {
    OPENMS_LOG_ERROR << "PSMArrowIO: QPXFile::exportToArrow returned null" << std::endl;
    return false;
  }
  if (!ArrowIOHelpers::writeTableToParquet(psm_table, dir + "/" + kPSMs, config))
  {
    OPENMS_LOG_ERROR << "PSMArrowIO: failed to write " << kPSMs << std::endl;
    return false;
  }

  if (!ProteinIdentificationArrowIO::exportProteinsToParquet(
        protein_identifications, dir + "/" + kProteins, config))
  {
    OPENMS_LOG_ERROR << "PSMArrowIO: failed to write " << kProteins << std::endl;
    return false;
  }
  if (!ProteinIdentificationArrowIO::exportProteinGroupsToParquet(
        protein_identifications, dir + "/" + kProteinGroups, config))
  {
    OPENMS_LOG_ERROR << "PSMArrowIO: failed to write " << kProteinGroups << std::endl;
    return false;
  }
  if (!ProteinIdentificationArrowIO::exportSearchParamsToParquet(
        protein_identifications, dir + "/" + kSearchParams, config))
  {
    OPENMS_LOG_ERROR << "PSMArrowIO: failed to write " << kSearchParams << std::endl;
    return false;
  }

  return true;
}

bool PSMArrowIO::importFromParquet(
  const String& dir,
  std::vector<ProteinIdentification>& protein_identifications,
  PeptideIdentificationList& peptide_identifications)
{
  if (!File::exists(dir) || !File::isDirectory(dir))
  {
    OPENMS_LOG_ERROR << "PSMArrowIO: not a directory: " << dir << std::endl;
    return false;
  }
  for (const char* sub : {kPSMs, kProteins, kProteinGroups, kSearchParams})
  {
    if (!File::exists(dir + "/" + sub))
    {
      OPENMS_LOG_ERROR << "PSMArrowIO: missing required file: " << dir + "/" + sub << std::endl;
      return false;
    }
  }

  std::vector<ProteinIdentification> tmp_proteins;
  PeptideIdentificationList tmp_peptides;

  if (!ProteinIdentificationArrowIO::importFromParquet(
        dir + "/" + kProteins,
        dir + "/" + kProteinGroups,
        dir + "/" + kSearchParams,
        tmp_proteins))
  {
    return false;
  }

  auto psm_table = ParquetFile::readTable(dir + "/" + kPSMs);
  if (!psm_table)
  {
    OPENMS_LOG_ERROR << "PSMArrowIO: failed to read " << kPSMs << std::endl;
    return false;
  }

  if (!QPXFile::importFromArrow(psm_table, tmp_proteins, tmp_peptides))
  {
    return false;
  }

  // Mirror IdXMLFile.cpp:530 — synthesize fresh ProtID identifiers on every load
  // and re-stamp pep_ids in lock-step. Synthesis runs only after all 4 tables have
  // loaded with stored identifiers as join keys; from here the in-memory identifier
  // is what downstream tools see.
  auto rename = ProteinIdentificationArrowIO::synthesizeRunIdentifiers(tmp_proteins);
  ProteinIdentificationArrowIO::applyRunIdentifierRename(rename, tmp_peptides);

  protein_identifications.swap(tmp_proteins);
  peptide_identifications.swap(tmp_peptides);
  return true;
}

} // namespace OpenMS
