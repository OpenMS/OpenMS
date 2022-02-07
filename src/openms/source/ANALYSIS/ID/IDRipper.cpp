// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Immanuel Luhn, Leon Kuchenbecker$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDRipper.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <QDir>

using std::vector;
using std::string;
using std::map;
using std::pair;
using std::make_pair;

namespace OpenMS
{

  const std::array<std::string, IDRipper::SIZE_OF_ORIGIN_ANNOTATION_FORMAT> IDRipper::names_of_OriginAnnotationFormat = {"file_origin", "map_index", "id_merge_index", "unknown"};

  IDRipper::IDRipper() :
    DefaultParamHandler("IDRipper")
  {
  }

  IDRipper::IDRipper(const IDRipper& cp) :
    DefaultParamHandler(cp)
  {
  }

  IDRipper::~IDRipper()
  {
  }

  IDRipper& IDRipper::operator=(const IDRipper& rhs)
  {
    if (this == &rhs)
      return *this;

    DefaultParamHandler::operator=(rhs);
    updateMembers_();

    return *this;
  }

  IDRipper::IdentificationRuns::IdentificationRuns(const vector<ProteinIdentification>& prot_ids)
  {
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

  UInt IDRipper::RipFileIdentifier::getIdentRunIdx()
  {
      return this->ident_run_idx;
  }

  UInt IDRipper::RipFileIdentifier::getFileOriginIdx()
  {
      return this->file_origin_idx;
  }

  const String & IDRipper::RipFileIdentifier::getOriginFullname()
  {
      return this->origin_fullname;
  }

  const String & IDRipper::RipFileIdentifier::getOutputBasename()
  {
      return this->out_basename;
  }

  const std::vector<ProteinIdentification> & IDRipper::RipFileContent::getProteinIdentifications()
  {
      return this->prot_idents;
  }

  const std::vector<PeptideIdentification> & IDRipper::RipFileContent::getPeptideIdentifications()
  {
      return this->pep_idents;
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

    // Collect all protein hits
    vector<ProteinHit> all_protein_hits;
    for (ProteinIdentification& prot : proteins)
    {
      // remove protein identification file origin
      prot.removeMetaValue(names_of_OriginAnnotationFormat[origin_annotation_fmt]);
      vector<ProteinHit>& protein_hits  = prot.getHits();
      all_protein_hits.insert(all_protein_hits.end(), protein_hits.begin(), protein_hits.end());
    }

    map<String, pair<UInt, UInt> > basename_to_numeric;
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

      // try to get peptide hits for peptide identification
      const vector<PeptideHit>& peptide_hits = pep.getHits();
      if (peptide_hits.empty())
      {
        continue;
      }
      // collect all protein accessions that are stored in the peptide hits
      vector<String> protein_accessions;
      getProteinAccessions_(protein_accessions, peptide_hits);

      // returns all protein hits that are associated with the given peptide hits
      vector<ProteinHit> protein2accessions;
      getProteinHits_(protein2accessions, all_protein_hits, protein_accessions);

      // search for the protein identification of the peptide identification
      ProteinIdentification prot_ident;
      getProteinIdentification_(prot_ident, pep, proteins);
      // TODO catch case that ProteinIdentification prot_ident is not found in the for-loop

      RipFileMap::iterator it = ripped.find(rfi);
      // If file_origin already exists
      if (it != ripped.end())
      {
        vector<ProteinIdentification>& prot_tmp = it->second.prot_idents;
        bool flag = true;

        for (vector<ProteinIdentification>::iterator it2 = prot_tmp.begin(); it2 != prot_tmp.end(); ++it2)
        {
          // ProteinIdentification is already there, just add protein hits
          if (prot_ident.getIdentifier().compare(it2->getIdentifier()) == 0)
          {
            for (const ProteinHit& prot : protein2accessions)
            {
              it2->insertHit(prot);
            }
            flag = false;
            break;
          }
        }
        // if it was not found
        if (flag)
        {
          prot_ident.setHits(protein2accessions);
          prot_tmp.push_back(prot_ident);
        }
        vector<PeptideIdentification>& pep_tmp = it->second.pep_idents;
        pep_tmp.push_back(pep);
      }
      else // otherwise create new entry for file_origin
      {
        // create protein identification, TODO parameters
        vector<ProteinIdentification> protein_idents;
        // only use the protein hits that are needed for the peptide identification
        prot_ident.setHits(protein2accessions);
        protein_idents.push_back(prot_ident);

        //create new peptide identification
        vector<PeptideIdentification> peptide_idents;
        peptide_idents.push_back(pep);

        //create and insert new map entry
        ripped.insert(make_pair(rfi, RipFileContent(protein_idents, peptide_idents)));
      }
    }
    // Reduce the spectra data string list if that's what we ripped by
    if (origin_annotation_fmt == MAP_INDEX || origin_annotation_fmt == ID_MERGE_INDEX)
    {
      RipFileMap::iterator it;
      for (it = ripped.begin(); it != ripped.end(); ++it)
      {
        const RipFileIdentifier& rfi = it->first;
        RipFileContent&          rfc = it->second;

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
      for (size_t i = 0; i<SIZE_OF_ORIGIN_ANNOTATION_FORMAT; ++i)
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

  void IDRipper::getProteinHits_(vector<ProteinHit>& result, const vector<ProteinHit>& protein_hits, const vector<String>& protein_accessions)
  {
    for (const String& it : protein_accessions)
    {
      for (const ProteinHit& prot : protein_hits)
      {
        if (prot.getAccession().compare(it) == 0)
        {
          result.push_back(prot);
        }
      }
    }
  }

  void IDRipper::getProteinAccessions_(vector<String>& result, const vector<PeptideHit>& peptide_hits)
  {
    for (const PeptideHit& it : peptide_hits)
    {
      std::set<String> protein_accessions = it.extractProteinAccessionsSet();
      result.insert(result.end(), protein_accessions.begin(), protein_accessions.end());
    }
  }

  void IDRipper::getProteinIdentification_(ProteinIdentification& result, PeptideIdentification pep_ident, std::vector<ProteinIdentification>& prot_idents)
  {
    const String& identifier = pep_ident.getIdentifier();

    for (ProteinIdentification& prot: prot_idents)
    {
      if (identifier.compare(prot.getIdentifier()) == 0)
      {
        result = prot;
        break;
      }
    }
  }

} // namespace OpenMS
