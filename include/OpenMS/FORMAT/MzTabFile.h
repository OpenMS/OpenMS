// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MZTABFILE_H
#define OPENMS_FORMAT_MZTABFILE_H

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/FORMAT/SVOutStream.h>

#include <boost/math/special_functions/fpclassify.hpp>

#include <vector>
#include <algorithm>

namespace OpenMS
{
/**
  @brief File adapter for MzTab files

  @ingroup FileIO
 */
class OPENMS_DLLAPI MzTabFile
{
public:
  ///Default constructor
  MzTabFile();
  ///Destructor
  ~MzTabFile();

  typedef std::map< std::pair< String, String>, std::vector<PeptideHit> > MapAccPepType;

  void store(const String& filename, const std::vector<ProteinIdentification>& protein_ids, const std::vector<PeptideIdentification>& peptide_ids, String in, String document_id) const;

protected:

  static void sortPSM_(std::vector<PeptideIdentification>::iterator begin, std::vector<PeptideIdentification>::iterator end);

  static void keepFirstPSM_(std::vector<PeptideIdentification>::iterator begin, std::vector<PeptideIdentification>::iterator end);

  /// Extract protein and peptide identifications for each run. maps are assumed empty.
  static void partitionIntoRuns(const std::vector<PeptideIdentification>& pep_ids,
                                const std::vector<ProteinIdentification>& pro_ids,
                                std::map<String, std::vector<PeptideIdentification> >& map_run_to_pepids,
                                std::map<String, std::vector<ProteinIdentification> >& map_run_to_proids
                                );


  /// create links from protein to peptides
  static void createProteinToPeptideLinks(const std::map<String, std::vector<PeptideIdentification> >& map_run_to_pepids, MapAccPepType& map_run_accession_to_pephits);

  /// Extracts, if possible a unique protein accession for a peptide hit in mzTab format. Otherwise NA is returned
  static String extractProteinAccession_(const PeptideHit& peptide_hit);

  /// Extracts, modifications and positions of a peptide hit in mzTab format
  static String extractPeptideModifications_(const PeptideHit& peptide_hit);

  /// Map search engine identifier to CV, param etc.
  static String mapSearchEngineToCvParam_(const String& openms_search_engine_name);

  static String mapSearchEngineScoreToCvParam_(const String& openms_search_engine_name, DoubleReal score, String score_type);

  static String extractNumPeptides(const String& common_identifier, const String& protein_accession,
                                   const MapAccPepType& map_run_accesion_to_peptides);

  // mzTab definition of distinct
  static String extractNumPeptidesDistinct(String common_identifier, String protein_accession,
                                           const MapAccPepType& map_run_accesion_to_peptides);

  // same as distinct but additional constraint of uniquenes (=maps to exactly one Protein)
  static String extractNumPeptidesUnambiguous(String common_identifier, String protein_accession,
                                              const MapAccPepType& map_run_accesion_to_peptides);

  static std::map<String, Size> extractNumberOfSubSamples_(const std::map<String, std::vector<ProteinIdentification> >& map_run_to_proids);

  static void writePeptideHeader_( SVOutStream& output, std::map<String, Size> n_sub_samples);

  static void writeProteinHeader_( SVOutStream& output, std::map<String, Size> n_sub_samples);

  static void writeProteinData_(SVOutStream& output,
                                const ProteinIdentification& prot_id,
                                Size run_count,
                                String input_filename,
                                bool has_coverage,
                                const MapAccPepType& map_run_accesion_to_peptides,
                                const std::map<String, Size>& map_run_to_num_sub
                                );

};

} // namespace OpenMS

#endif // OPENMS_FORMAT_MZTABFILE_H
