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
// $Maintainer: Timo Sachsenberg $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_IDMerger IDMerger

  @brief Merges several idXML files into one idXML file.

  <center>
  <table>
  <tr>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ IDMerger \f$ \longrightarrow \f$</td>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ConsensusID </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFileConverter </td>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMapper </td>
  </tr>
  </table>
  </center>

  The peptide hits and protein hits of the input files will be written into the single output file. In general, the number of idXML files that can be merged into one file is not limited.

  The combination of search engine and processing date/time should be unique for every identification run over all input files. If this is not the case, the date/time of a conflicting run will be increased in steps of seconds until the combination is unique.

  If an additional file is given through the @p add_to parameter, identifications from the main inputs (@p in) are added to that file, but only for those peptide sequences that were not already present. Only the best peptide hit per identification (MS2 spectrum) is taken into account; peptide identifications and their corresponding protein identifications are transferred.

  Alternatively, with the @p pepxml_protxml option, results from corresponding PeptideProphet and ProteinProphet runs can be combined. In this case, exactly two idXML files are expected as input: one containing data from a pepXML file, and the other containing data from a protXML file that was created based on the pepXML (meaningful results can only be obtained for matching files!). pepXML or protXML can be converted to idXML with the @ref TOPP_IDFileConverter tool.

  @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_IDMerger.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_IDMerger.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDMerger :
  public TOPPBase
{
public:
  TOPPIDMerger() :
    TOPPBase("IDMerger", "Merges several protein/peptide identification files into one file.")
  {

  }

protected:
  void mergePepXMLProtXML_(StringList filenames, vector<ProteinIdentification>&
                           proteins, vector<PeptideIdentification>& peptides)
  {
    IdXMLFile idxml;
    idxml.load(filenames[0], proteins, peptides);
    vector<ProteinIdentification> pepxml_proteins, protxml_proteins;
    vector<PeptideIdentification> pepxml_peptides, protxml_peptides;

    if (proteins[0].getProteinGroups().empty()) // first idXML contains data from the pepXML
    {
      proteins.swap(pepxml_proteins);
      peptides.swap(pepxml_peptides);
      idxml.load(filenames[1], protxml_proteins, protxml_peptides);
      if (protxml_proteins[0].getProteinGroups().empty())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "None of the input files seems to be derived from a protXML file (information about protein groups is missing).");
      }
    }
    else // first idXML contains data from the protXML
    {
      proteins.swap(protxml_proteins);
      peptides.swap(protxml_peptides);
      idxml.load(filenames[1], pepxml_proteins, pepxml_peptides);
    }

    if ((protxml_peptides.size() > 1) || (protxml_proteins.size() > 1))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The idXML derived from a protXML file should contain only one 'ProteinIdentification' and one 'PeptideIdentification' instance.");
    }

    // peptide information comes from the pepXML (additional information in
    // the protXML - adapted peptide hit score, "is_unique", "is_contributing"
    // - is not transferred):
    peptides.swap(pepxml_peptides);

    // prepare scores and coverage values of protein hits from the protXML:
    map<String, pair<double, double> > hit_values;
    ProteinIdentification& protein = protxml_proteins[0];
    for (vector<ProteinHit>::iterator hit_it = protein.getHits().begin();
         hit_it != protein.getHits().end(); ++hit_it)
    {
      hit_values[hit_it->getAccession()] = make_pair(hit_it->getScore(),
                                                     hit_it->getCoverage());
    }

    // merge protein information:
    proteins.swap(pepxml_proteins);
    for (vector<ProteinIdentification>::iterator prot_it = proteins.begin();
         prot_it != proteins.end(); ++prot_it)
    {
      prot_it->getProteinGroups() = protein.getProteinGroups();
      prot_it->getIndistinguishableProteins() =
        protein.getIndistinguishableProteins();
      // TODO: since a protXML file can integrate data from several protein
      // identification runs, the protein groups/indistinguishable proteins
      // that we write to one identification run could contain references to
      // proteins that are not observed in this run, but in others; also, some
      // protein hits without enough evidence may not occur in the protXML
      // (thus also not in the protein groups) - clean this up?

      prot_it->setScoreType(protein.getScoreType());
      prot_it->setHigherScoreBetter(protein.isHigherScoreBetter());
      prot_it->setSignificanceThreshold(protein.getSignificanceThreshold());

      for (vector<ProteinHit>::iterator hit_it = prot_it->getHits().begin();
           hit_it != prot_it->getHits().end(); ++hit_it)
      {
        map<String, pair<double, double> >::const_iterator pos =
          hit_values.find(hit_it->getAccession());
        if (pos == hit_values.end())
        {
          hit_it->setScore(-1);
        }
        else
        {
          hit_it->setScore(pos->second.first);
          hit_it->setCoverage(pos->second.second);
        }
      }
    }
  }

  void generateNewId_(const map<String, ProteinIdentification>& used_ids, const String& search_engine, DateTime& date_time, String& new_id)
  {
    do
    {
      if (date_time.isValid())
      {
        date_time = date_time.addSecs(1);
      }
      else
      {
        date_time = DateTime::now();
      }
      new_id = search_engine + "_" + date_time.toString(Qt::ISODate);
    }
    while (used_ids.find(new_id) != used_ids.end());
  }

  void annotateFileOrigin_(vector<ProteinIdentification>& proteins,
                           vector<PeptideIdentification>& peptides,
                           String filename)
  {
    if (test_mode_) filename = File::basename(filename);
    for (vector<ProteinIdentification>::iterator prot_it = proteins.begin();
         prot_it != proteins.end(); ++prot_it)
    {
      prot_it->setMetaValue("file_origin", DataValue(filename));
    }

    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      pep_it->setMetaValue("file_origin", DataValue(filename));
    }
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blanks");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerInputFile_("add_to", "<file>", "", "Optional input file. IDs from 'in' are added to this file, but only if the (modified) peptide sequences are not present yet (considering only best hits per spectrum).", false);
    setValidFormats_("add_to", ListUtils::create<String>("idXML"));
    registerFlag_("annotate_file_origin", "Store the original filename in each protein/peptide identification (meta value: file_origin).");
    registerFlag_("pepxml_protxml", "Merge idXML files derived from a pepXML and corresponding protXML file.\nExactly two input files are expected in this case. Not compatible with 'add_to'.");
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    StringList file_names = getStringList_("in");
    String out = getStringOption_("out");
    String add_to = getStringOption_("add_to");

    if (file_names.empty())
    {
      // this also allows exactly 1 file, because it might be usefull for
      // a TOPPAS pipeline containing an IDMerger, to run only with one file
      writeLog_("No input filename given. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    bool pepxml_protxml = getFlag_("pepxml_protxml");
    if (pepxml_protxml && (file_names.size() != 2))
    {
      writeLog_("Exactly two input filenames expected for option 'pepxml_protxml'. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }
    if (pepxml_protxml && !add_to.empty())
    {
      // currently not allowed to keep the code simpler and because it doesn't
      // seem useful, but should be possible in principle:
      writeLog_("The options 'add_to' and 'pepxml_protxml' cannot be used together. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    vector<ProteinIdentification> proteins;
    vector<PeptideIdentification> peptides;

    if (pepxml_protxml)
    {
      mergePepXMLProtXML_(file_names, proteins, peptides);
    }
    else
    {
      bool annotate_file_origin = getFlag_("annotate_file_origin");
      map<String, ProteinIdentification> proteins_by_id;
      vector<vector<PeptideIdentification> > peptides_by_file;
      StringList add_to_ids; // IDs from the "add_to" file (if any)

      if (!add_to.empty())
      {
        file_names.erase(std::remove(file_names.begin(), file_names.end(), add_to), file_names.end());
        file_names.insert(file_names.begin(), add_to);
      }

      peptides_by_file.resize(file_names.size());
      for (Size i = 0; i < file_names.size(); ++i)
      {
        const String& file_name = file_names[i];
        vector<ProteinIdentification> additional_proteins;
        IdXMLFile().load(file_name, additional_proteins, peptides_by_file[i]);

        if (annotate_file_origin) // set MetaValue "file_origin" if flag is set
        {
          annotateFileOrigin_(additional_proteins, peptides_by_file[i],
                              file_name);
        }

        for (vector<ProteinIdentification>::iterator prot_it =
               additional_proteins.begin(); prot_it !=
             additional_proteins.end(); ++prot_it)
        {
          String id = prot_it->getIdentifier();
          if (proteins_by_id.find(id) != proteins_by_id.end())
          {
            writeLog_("Warning: The identifier '" + id + "' was used before!");
            // generate a new ID:
            DateTime date_time = prot_it->getDateTime();
            String new_id;
            generateNewId_(proteins_by_id, prot_it->getSearchEngine(),
                           date_time, new_id);
            writeLog_("New identifier '" + new_id +
                      "' generated as replacement.");
            // update fields:
            prot_it->setIdentifier(new_id);
            prot_it->setDateTime(date_time);
            for (vector<PeptideIdentification>::iterator pep_it =
                   peptides_by_file[i].begin(); pep_it !=
                 peptides_by_file[i].end(); ++pep_it)
            {
              if (pep_it->getIdentifier() == id) pep_it->setIdentifier(new_id);
            }
            id = new_id;
          }
          proteins_by_id[id] = *prot_it;
          if (i == 0) add_to_ids.push_back(id);
        }
      }

      if (add_to.empty()) // copy proteins from map into vector for writing
      {
        for (vector<vector<PeptideIdentification> >::iterator pep_it =
               peptides_by_file.begin(); pep_it != peptides_by_file.end();
             ++pep_it)
        {
          peptides.insert(peptides.end(), pep_it->begin(), pep_it->end());
        }
        for (map<String, ProteinIdentification>::iterator map_it =
               proteins_by_id.begin(); map_it != proteins_by_id.end(); ++map_it)
        {
          proteins.push_back(map_it->second);
        }
      }
      else // add only new IDs to an existing file
      {
        // copy over data from reference file ("add_to"):
        map<String, ProteinIdentification> selected_proteins;
        for (StringList::iterator ids_it = add_to_ids.begin();
             ids_it != add_to_ids.end(); ++ids_it)
        {
          selected_proteins[*ids_it] = proteins_by_id[*ids_it];
        }
        // keep track of peptides that shouldn't be duplicated:
        set<AASequence> sequences;
        vector<PeptideIdentification>& base_peptides = peptides_by_file[0];
        for (vector<PeptideIdentification>::iterator pep_it =
               base_peptides.begin(); pep_it != base_peptides.end(); ++pep_it)
        {
          if (pep_it->getHits().empty()) continue;
          pep_it->sort();
          sequences.insert(pep_it->getHits()[0].getSequence());
        }
        peptides.insert(peptides.end(), base_peptides.begin(),
                        base_peptides.end());
        // merge in data from other files:
        for (vector<vector<PeptideIdentification> >::iterator file_it =
               ++peptides_by_file.begin(); file_it != peptides_by_file.end();
             ++file_it)
        {
          set<String> accessions; // keep track to avoid duplicates
          for (vector<PeptideIdentification>::iterator pep_it =
                 file_it->begin(); pep_it != file_it->end(); ++pep_it)
          {
            if (pep_it->getHits().empty()) continue;
            pep_it->sort();
            const PeptideHit& hit = pep_it->getHits()[0];
            LOG_DEBUG << "peptide: " << hit.getSequence().toString() << endl;
            // skip ahead if peptide is not new:
            if (sequences.find(hit.getSequence()) != sequences.end()) continue;
            LOG_DEBUG << "new peptide!" << endl;
            pep_it->getHits().resize(1); // restrict to best hit for simplicity
            peptides.push_back(*pep_it);

            set<String> protein_accessions = hit.extractProteinAccessionsSet();

            // copy over proteins:
            for (set<String>::const_iterator acc_it = protein_accessions.begin(); acc_it != protein_accessions.end(); ++acc_it)
            {
              LOG_DEBUG << "accession: " << *acc_it << endl;
              // skip ahead if accession is not new:
              if (accessions.find(*acc_it) != accessions.end()) continue;
              LOG_DEBUG << "new accession!" << endl;
              // first find the right protein identification:
              const String& id = pep_it->getIdentifier();
              LOG_DEBUG << "identifier: " << id << endl;
              if (proteins_by_id.find(id) == proteins_by_id.end())
              {
                writeLog_("Error: identifier '" + id + "' linking peptides and proteins not found. Skipping.");
                continue;
              }
              ProteinIdentification& protein = proteins_by_id[id];
              // now find the protein hit:
              vector<ProteinHit>::iterator hit_it = protein.findHit(*acc_it);
              if (hit_it == protein.getHits().end())
              {
                writeLog_("Error: accession '" + *acc_it + "' not found in "
                                                           "protein identification '" + id + "'. Skipping.");
                continue;
              }
              // we may need to copy protein ID meta data, if we haven't yet:
              if (selected_proteins.find(id) == selected_proteins.end())
              {
                LOG_DEBUG << "adding protein identification" << endl;
                selected_proteins[id] = protein;
                selected_proteins[id].getHits().clear();
                // remove potentially invalid information:
                selected_proteins[id].getProteinGroups().clear();
                selected_proteins[id].getIndistinguishableProteins().clear();
              }
              selected_proteins[id].insertHit(*hit_it);
              accessions.insert(*acc_it);
              // NOTE: we're only adding the first protein hit for each
              // accession, not taking into account scores or any meta data
            }
          }
        }
        for (map<String, ProteinIdentification>::iterator map_it =
               selected_proteins.begin(); map_it != selected_proteins.end();
             ++map_it)
        {
          proteins.push_back(map_it->second);
        }
      }
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    LOG_DEBUG << "protein IDs: " << proteins.size() << endl
              << "peptide IDs: " << peptides.size() << endl;
    IdXMLFile().store(out, proteins, peptides);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPIDMerger tool;
  return tool.main(argc, argv);
}

/// @endcond
