// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Immanuel Luhn, Leon Kuchenbecker$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDRipper.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <QtCore/QDir>
#include <array>
#include <unordered_set>

using namespace std;

namespace OpenMS
{

  const std::array<std::string, IDRipper::SIZE_OF_ORIGIN_ANNOTATION_FORMAT> IDRipper::names_of_OriginAnnotationFormat = {"file_origin", "map_index", Constants::UserParam::ID_MERGE_INDEX, "unknown"};

  IDRipper::IDRipper() :
    DefaultParamHandler("IDRipper")
  {
  }

  IDRipper::IDRipper(const IDRipper& cp) = default;

  IDRipper::~IDRipper() = default;

  IDRipper& IDRipper::operator=(const IDRipper& rhs) = default;

  IDRipper::IdentificationRuns::IdentificationRuns(const vector<ProteinIdentification>& prot_ids)
  {
    // build index_ map that maps the identifiers in prot_ids to indices 0,1,...
    for (const auto& prot_id : prot_ids)
    {
      String id_run_id = prot_id.getIdentifier();
      if (this->index_map.find(id_run_id) != this->index_map.end())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IdentificationRun IDs are not unique!", id_run_id);
      }
      UInt idx = this->index_map.size();
      this->index_map[id_run_id] = idx;
      const DataValue& mv_spectra_data = prot_id.getMetaValue("spectra_data");
      spectra_data.push_back(mv_spectra_data.isEmpty() ? StringList() : mv_spectra_data.toStringList());
    }
  }

  bool IDRipper::RipFileIdentifierIdxComparator::operator()(const RipFileIdentifier& left, const RipFileIdentifier& right) const
  {
    return std::tie(left.ident_run_idx, left.file_origin_idx)
      < std::tie(right.ident_run_idx, right.file_origin_idx);
  }

  // Identify the output file name features associated via spectra_data or file_origin
  IDRipper::RipFileIdentifier::RipFileIdentifier(const IDRipper::IdentificationRuns& id_runs, const PeptideIdentification& pep_id, const map<String, UInt>& file_origin_map, const IDRipper::OriginAnnotationFormat origin_annotation_fmt, bool split_ident_runs)
  {
      try
      {
          // Numerical identifier of the Identification Run
          this->ident_run_idx   = id_runs.index_map.at(pep_id.getIdentifier());

          // Numerical identifier of the PeptideIdentification origin
          this->file_origin_idx = (origin_annotation_fmt == MAP_INDEX || origin_annotation_fmt == ID_MERGE_INDEX)
              ? pep_id.getMetaValue(names_of_OriginAnnotationFormat[origin_annotation_fmt]).toString().toInt()
              : file_origin_map.at(pep_id.getMetaValue("file_origin").toString());

          // Store the origin full name
          this->origin_fullname = (origin_annotation_fmt == MAP_INDEX || origin_annotation_fmt == ID_MERGE_INDEX)
              ? id_runs.spectra_data.at(this->ident_run_idx).at(this->file_origin_idx)
              : pep_id.getMetaValue("file_origin").toString();

          // Extract the basename, used for output files when --numeric_filenames is not set
          this->out_basename = QFileInfo(this->origin_fullname.toQString()).completeBaseName().toStdString();

          // Drop the identification run identifier if we're not splitting by identification runs
          if (!split_ident_runs)
              this->ident_run_idx = -1u;
      }
      catch (const std::out_of_range& e)
      {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "input file",
                  "Failed to identify corresponding spectra_data element for PeptideIdentification element.");
      }
  }

  UInt IDRipper::RipFileIdentifier::getIdentRunIdx() const
  {
    return ident_run_idx;
  }

  UInt IDRipper::RipFileIdentifier::getFileOriginIdx() const
  {
    return file_origin_idx;
  }

  const String & IDRipper::RipFileIdentifier::getOriginFullname() const
  {
    return origin_fullname;
  }

  const String & IDRipper::RipFileIdentifier::getOutputBasename() const
  {
    return out_basename;
  }

  const std::vector<ProteinIdentification> & IDRipper::RipFileContent::getProteinIdentifications()
  {
    return prot_idents;
  }

  const std::vector<PeptideIdentification> & IDRipper::RipFileContent::getPeptideIdentifications()
  {
    return pep_idents;
  }

  bool IDRipper::registerBasename_(map<String, pair<UInt, UInt> >& basename_to_numeric, const IDRipper::RipFileIdentifier& rfi)
  {
      auto it = basename_to_numeric.find(rfi.out_basename);
      auto p  = make_pair(rfi.ident_run_idx, rfi.file_origin_idx);

      // If we have not seen this basename before, store it in the map
      if (it == basename_to_numeric.end())
      {
          basename_to_numeric[rfi.out_basename] = p;
          return true;
      }
      // Otherwise, check if we save it in the context of the same IdentificationRun and potentially spectra_data position
      return it->second == p;
  }

  void IDRipper::rip(
          RipFileMap& ripped,
          vector<ProteinIdentification>& proteins,
          vector<PeptideIdentification>& peptides,
          bool numeric_filenames,
          bool split_ident_runs)
  {
    // Detect file format w.r.t. origin annotation
    map<String, UInt> file_origin_map;
    IDRipper::OriginAnnotationFormat origin_annotation_fmt = detectOriginAnnotationFormat_(file_origin_map, peptides);

    if (origin_annotation_fmt == UNKNOWN_OAF)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "input file",
        "Unable to detect origin annotation format of provided input file.");
    }

    OPENMS_LOG_DEBUG << "Detected file origin annotation format: " << names_of_OriginAnnotationFormat[origin_annotation_fmt] << std::endl;

    // Build identifier index
    const IdentificationRuns id_runs = IdentificationRuns(proteins);

    // Collect a unique set of representative protein hits. One per accession. Looks at all runs and removes the file origin
    unordered_map<String, const ProteinHit*> acc2protein_hits;
    for (ProteinIdentification& prot : proteins)
    {
      prot.removeMetaValue(names_of_OriginAnnotationFormat[origin_annotation_fmt]);
      const vector<ProteinHit>& protein_hits  = prot.getHits();
      for (const auto& ph : protein_hits)
      {
        acc2protein_hits[ph.getAccession()] = &ph;
      }
    }

    size_t protein_identifier_not_found{};

    map<String, pair<UInt, UInt> > basename_to_numeric;

    // map run identifier to protein accessions that were already added
    map<IDRipper::RipFileIdentifier, unordered_map<String, unordered_set<String>>, RipFileIdentifierIdxComparator> ripped_prot_map;

    //store protein and peptides identifications for each file origin
    for (PeptideIdentification& pep : peptides)
    {
      // Build the output file identifier
      const IDRipper::RipFileIdentifier rfi(id_runs, pep, file_origin_map, origin_annotation_fmt, split_ident_runs);

      // If we are inferring the output file names from the spectra_data or
      // file_origin, make sure they are unique
      if (!numeric_filenames && !registerBasename_(basename_to_numeric, rfi))
      {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Autodetected output file names are not unique. Use -numeric_filenames.");
      }
      
      // remove file origin annotation
      pep.removeMetaValue(names_of_OriginAnnotationFormat[origin_annotation_fmt]);

      // get peptide hits (PSMs) for each peptide identification (spectrum)
      const vector<PeptideHit>& peptide_hits = pep.getHits();
      if (peptide_hits.empty())
      {
        continue;
      }
      // collect all protein accessions that are stored in the peptide hits
      set<String> protein_accessions = getProteinAccessions_(peptide_hits);
      if (protein_accessions.empty())
      {
        OPENMS_LOG_WARN << "Peptide hits with empty protein accession." << std::endl;
        continue;
      }

      // returns all protein hits that are associated with the accessions of current peptide hits
      vector<ProteinHit> proteins_of_accessions;
      getProteinHits_(proteins_of_accessions, acc2protein_hits, protein_accessions);
      if (proteins_of_accessions.empty())
      {
        OPENMS_LOG_WARN << "No proteins found for given accessions." << std::endl;
        continue;
      }

      // search for the protein identification of the peptide identification
      int prot_ident_index = getProteinIdentification_(pep, id_runs);
      if (prot_ident_index == -1)
      {
        ++protein_identifier_not_found;
        OPENMS_LOG_WARN << "Run identifier: " << pep.getIdentifier() << " was not found in protein identification runs." << std::endl;
        continue;
      }

      const ProteinIdentification& merged_protein_id_run = proteins[prot_ident_index]; // protein identification run in the merged file
      const String& merged_prot_identifier = merged_protein_id_run.getIdentifier();        // protein identification run identifier in merged file
    
      if (RipFileMap::iterator it = ripped.find(rfi); it == ripped.end())
      { // file identifier does not exist yet. We need to create it.
        OPENMS_LOG_INFO << "Creating entry for file identifier:\n" 
                        << "File origin: " << rfi.getOriginFullname() << "\n"
                        // << "Identification run index: " << rfi.getIdentRunIdx() << "\n" // not set here?
                        << "Basename: " << rfi.getOutputBasename() << "\n"
                        << "Merged identification file run identifier: " << merged_prot_identifier << "\n"
                        << std::endl;
        // create the protein run but only set the protein hits that are needed for the current peptide identification
        ProteinIdentification p;
        p.copyMetaDataOnly(merged_protein_id_run);
        p.setHits(proteins_of_accessions); // TODO: what about protein groups?
        for (const ProteinHit& prot : proteins_of_accessions)
        { // register protein so we don't add it twice          
          const String& acc = prot.getAccession();
          ripped_prot_map[rfi][merged_prot_identifier].insert(acc);
        }

        vector<ProteinIdentification> protein_idents;
        protein_idents.push_back(std::move(p));

        //create new peptide identification
        vector<PeptideIdentification> peptide_idents;
        peptide_idents.push_back(pep);

        //create and insert new map entry
        ripped.insert(make_pair(rfi, RipFileContent(protein_idents, peptide_idents)));
      }    
      else
      { // if file identifier already exists we attach
        // query all protein identification runs for file identifier in the Ripped data structure
        vector<ProteinIdentification>& ripped_protein_id_runs = it->second.prot_idents;

        bool ripped_protein_identifier_exists{false};

        for (auto& ripped_protein_id_run : ripped_protein_id_runs)
        { // for all protein identification runs associated with the current file identifier...
          const String& ripped_prot_identifier = ripped_protein_id_run.getIdentifier();
          if (merged_prot_identifier == ripped_prot_identifier)
          { // protein identification run already exists in ripped map. just add protein hits if not already present            
            for (const ProteinHit& prot : proteins_of_accessions)
            {               
              // check if protein has already been added              
              const String& acc = prot.getAccession();
              auto& acc_set = ripped_prot_map[rfi][merged_prot_identifier];
              if (auto ri = acc_set.find(acc); ri == acc_set.end())
              { // only add protein once to the run identifier
                ripped_protein_id_run.insertHit(prot);
                acc_set.insert(acc);
                #ifdef DEBUG_IDRIPPER
                  std::cout << "ripped/merged identifier: " << ripped_prot_identifier << " " << prot << std::endl;
                #endif
              }                
            }
            ripped_protein_identifier_exists = true;
            break;
          }
        }

        // file identifier exists but not the protein identification run identifier -  we did not add anything so far to it
        if (!ripped_protein_identifier_exists) 
        {    
          ProteinIdentification p;
          p.copyMetaDataOnly(merged_protein_id_run);

          for (const ProteinHit& prot : proteins_of_accessions)
          {              
            // check if protein has already been added
            const String& acc = prot.getAccession();
             auto& acc_set = ripped_prot_map[rfi][merged_prot_identifier];
            if (auto ri = acc_set.find(acc); ri == acc_set.end())
            { // only add protein once to the run identifier
              p.insertHit(prot);;
              acc_set.insert(acc);
              #ifdef DEBUG_IDRIPPER
                std::cout << "ripped/merged identifier: " << ripped_prot_identifier << " " << prot << std::endl;
              #endif
            }                
          }
          ripped_protein_id_runs.push_back(std::move(p));
        }

        // add current peptide identification
        vector<PeptideIdentification>& ripped_pep = it->second.pep_idents;
        ripped_pep.push_back(pep);
      }
    }

    // Reduce the spectra data string list if that's what we ripped by
    if (origin_annotation_fmt == MAP_INDEX || origin_annotation_fmt == ID_MERGE_INDEX)
    {
      RipFileMap::iterator it;
      for (it = ripped.begin(); it != ripped.end(); ++it)
      {
        const RipFileIdentifier& rfi = it->first;
        RipFileContent& rfc = it->second;

        for (ProteinIdentification& prot_id : rfc.prot_idents)
        {
          StringList old_list;
          prot_id.getPrimaryMSRunPath(old_list);
          StringList new_list;
          new_list.push_back(rfi.origin_fullname);
          prot_id.setPrimaryMSRunPath(new_list);
        }
      }
    }
    if (protein_identifier_not_found > 0)
    {
      OPENMS_LOG_ERROR << "Some protein identification runs referenced in peptide identifications were not found." << std::endl;
    }
  }

  void IDRipper::rip(
            std::vector<RipFileIdentifier> & rfis,
            std::vector<RipFileContent> & rfcs,
            std::vector<ProteinIdentification> & proteins,
            std::vector<PeptideIdentification> & peptides,
            bool numeric_filenames,
            bool split_ident_runs)
  {
      RipFileMap rfm;
      this->rip(rfm, proteins, peptides, numeric_filenames, split_ident_runs);

      rfis.clear();
      rfcs.clear();
      for (RipFileMap::iterator it = rfm.begin(); it != rfm.end(); ++it)
      {
          rfis.push_back(it->first);
          rfcs.push_back(it->second);
      }
  }


bool IDRipper::setOriginAnnotationMode_(short& mode, short const new_value)
{
  if (mode != -1 && mode != new_value)
  {
    return false;
  }
  mode = new_value;
  return true;
}

IDRipper::OriginAnnotationFormat IDRipper::detectOriginAnnotationFormat_(map<String, UInt>& file_origin_map, const std::vector<PeptideIdentification>& peptide_idents)
  {
    // In case we observe 'file_origin' meta values, we assign an index to every unique meta value
    file_origin_map.clear();

    short mode = -1;
    for (vector<PeptideIdentification>::const_iterator it = peptide_idents.begin(); it != peptide_idents.end(); ++it)
    {
      bool mode_identified = false;
      for (size_t i = 0; i < SIZE_OF_ORIGIN_ANNOTATION_FORMAT; ++i)
      {
        if (it->metaValueExists(names_of_OriginAnnotationFormat[i]))
        {
          // Different mode identified for same or different peptide
          if (mode_identified || !setOriginAnnotationMode_(mode, i))
          {
            return UNKNOWN_OAF;
          }
          else
          {
            mode_identified = true;
          }

          if (i == 0) // names_of_OriginAnnotationFormat[0] == "file_origin"
          {
            const String& file_origin = it->getMetaValue("file_origin");
            // Did we already assign an index to this file_origin?
            if (file_origin_map.find(file_origin) == file_origin_map.end())
            {
              // If not, assign a new unique index
              size_t cur_size = file_origin_map.size();
              file_origin_map[file_origin] = cur_size;
            }
          }
        }
      }
      if (!mode_identified)
      {
        return UNKNOWN_OAF;
      }
    }
    if (mode == -1)
    {
      return UNKNOWN_OAF;
    }
    else
    {
      return static_cast<IDRipper::OriginAnnotationFormat>(mode);
    }
  }

  void IDRipper::getProteinHits_(vector<ProteinHit>& result, const unordered_map<String, const ProteinHit*>& acc2protein_hits, const set<String>& protein_accessions)
  {
    for (const String& s : protein_accessions)
    {
      if (auto it = acc2protein_hits.find(s); it != acc2protein_hits.end())
      {
        const ProteinHit* prot_ptr = it->second;
        result.push_back(*prot_ptr);
      }
    }
  }

  std::set<String> IDRipper::getProteinAccessions_(const vector<PeptideHit>& peptide_hits)
  {
    std::set<String> accession_set;
    for (const PeptideHit& it : peptide_hits)
    {
      std::set<String> protein_accessions = it.extractProteinAccessionsSet();
      accession_set.insert(make_move_iterator(protein_accessions.begin()), make_move_iterator(protein_accessions.end()));
    }
    return accession_set;
  }

  int IDRipper::getProteinIdentification_(const PeptideIdentification& pep_ident, const IdentificationRuns& id_runs)
  {
    const String& identifier = pep_ident.getIdentifier();
    if (auto it = id_runs.index_map.find(identifier); it != id_runs.index_map.end())
    {
      return it->second;
    }
    return -1;
  }

} // namespace OpenMS
