// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/StreamingTSVWriter.h>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{

  StreamingTSVWriter::StreamingTSVWriter() = default;

  StreamingTSVWriter::~StreamingTSVWriter()
  {
    if (ofs_.is_open())
    {
      close();
    }
  }

  void StreamingTSVWriter::open(const String& filename)
  {
    ofs_.open(filename.c_str());
    if (!ofs_.is_open())
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    header_written_ = false;
    transitions_written_ = 0;
    groups_written_ = 0;
  }

  void StreamingTSVWriter::writeHeader_()
  {
    ofs_ << "PrecursorMz\t"
         << "ProductMz\t"
         << "PrecursorCharge\t"
         << "ProductCharge\t"
         << "LibraryIntensity\t"
         << "NormalizedRetentionTime\t"
         << "PeptideSequence\t"
         << "ModifiedPeptideSequence\t"
         << "ProteinId\t"
         << "TransitionGroupId\t"
         << "TransitionId\t"
         << "Decoy\t"
         << "DetectingTransition\t"
         << "IdentifyingTransition\t"
         << "QuantifyingTransition";

    // Add optional columns for metabolomics
    ofs_ << "\tCompoundName\t"
         << "SumFormula";

    // Add ion mobility
    ofs_ << "\tPrecursorIonMobility";

    ofs_ << "\n";
    header_written_ = true;
  }

  void StreamingTSVWriter::writeTransition_(
    const OpenSwath::LightTransition& transition,
    const OpenSwath::LightCompound& compound,
    const std::vector<std::string>& protein_ids)
  {
    // Build protein string (semicolon-separated)
    std::string protein_str;
    for (size_t i = 0; i < protein_ids.size(); ++i)
    {
      if (i > 0) protein_str += ";";
      protein_str += protein_ids[i];
    }

    ofs_ << transition.precursor_mz << "\t"
         << transition.product_mz << "\t"
         << compound.charge << "\t"
         << transition.fragment_charge << "\t"
         << transition.library_intensity << "\t"
         << compound.rt << "\t"
         << compound.sequence << "\t"
         << compound.sequence << "\t"  // ModifiedPeptideSequence (simplified)
         << protein_str << "\t"
         << transition.peptide_ref << "\t"
         << transition.transition_name << "\t"
         << (transition.decoy ? "1" : "0") << "\t"
         << (transition.detecting_transition ? "1" : "0") << "\t"
         << (transition.identifying_transition ? "1" : "0") << "\t"
         << (transition.quantifying_transition ? "1" : "0");

    // Optional metabolomics columns
    ofs_ << "\t" << compound.compound_name
         << "\t" << compound.sum_formula;

    // Ion mobility
    if (transition.precursor_im > 0)
    {
      ofs_ << "\t" << transition.precursor_im;
    }
    else
    {
      ofs_ << "\t";
    }

    ofs_ << "\n";
    ++transitions_written_;
  }

  void StreamingTSVWriter::writeGroup(
    const std::string& /* group_id */,
    const std::vector<OpenSwath::LightTransition>& transitions,
    const std::vector<OpenSwath::LightTransition>& decoys,
    const OpenSwath::LightCompound& compound,
    const std::vector<std::string>& protein_ids)
  {
    if (!header_written_)
    {
      writeHeader_();
    }

    // Write target transitions
    for (const auto& tr : transitions)
    {
      writeTransition_(tr, compound, protein_ids);
    }

    // Write decoy transitions
    for (const auto& tr : decoys)
    {
      writeTransition_(tr, compound, protein_ids);
    }

    ++groups_written_;
  }

  void StreamingTSVWriter::close()
  {
    if (ofs_.is_open())
    {
      ofs_.flush();
      ofs_.close();
    }
  }

} // namespace OpenMS
