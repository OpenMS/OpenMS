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
// $Maintainer: Chris Bielow, Hendrik Weisser $
// $Authors: Chris Bielow, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/PepXMLFile.h>

#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/HANDLERS/MascotXMLHandler.h> // for "primary_scan_regex"
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <fstream>
#include <iostream>
#include <boost/regex.hpp>

using namespace std;

namespace OpenMS
{

  PepXMLFile::PepXMLFile() :
    XMLHandler("", "1.12"),
    XMLFile("/SCHEMAS/pepXML_v114.xsd", "1.14"),
    proteins_(nullptr),
    peptides_(nullptr),
    lookup_(nullptr),
    scan_map_(),
    analysis_summary_(false),
    keep_native_name_(false),
    search_score_summary_(false)
  {
    const ElementDB* db = ElementDB::getInstance();
    hydrogen_ = *db->getElement("Hydrogen");
  }

  const double PepXMLFile::mod_tol_ = 0.001;
  const double PepXMLFile::xtandem_artificial_mod_tol_ = 0.0005; // according to cpp in some old version of xtandem somehow very small fixed modification (electron mass?) gets annotated by X!Tandem. Don't add them as they interfere with other modifications.

  PepXMLFile::~PepXMLFile()
  {
  }

  void PepXMLFile::store(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids, const String& mz_file, const String& mz_name, bool peptideprophet_analyzed)
  {
    ofstream f(filename.c_str());
    if (!f)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    String search_engine_name;
    ProteinIdentification::SearchParameters search_params;
    if (!protein_ids.empty())
    {
      if (protein_ids.size() > 1)
      {
        warning(STORE, "More than one protein identification defined; only first search parameters are written into pepXML; search engine must be the same.");
      }
      search_params = protein_ids.begin()->getSearchParameters();
      if (protein_ids.begin()->getSearchEngine() == "XTandem")
      {
        search_engine_name = "X! Tandem";
      }
      else if (protein_ids.begin()->getSearchEngine() == "Mascot")
      {
        search_engine_name = "MASCOT";
      }
      else
      {
        search_engine_name = protein_ids.begin()->getSearchEngine();
        //Comet writes "Comet" in pep.xml, so this is ok
      }
    }

    f.precision(writtenDigits<double>(0.0));
    String raw_data;
    String base_name;
    SpectrumMetaDataLookup lookup;

    // The mz-File (if given)
    if (!mz_file.empty())
    {
      base_name = File::removeExtension(File::basename(mz_file));
      raw_data = FileTypes::typeToName(FileHandler().getTypeByFileName(mz_file));

      PeakMap experiment;
      FileHandler fh;
      fh.loadExperiment(mz_file, experiment, FileTypes::UNKNOWN, ProgressLogger::NONE, false, false);
      lookup.readSpectra(experiment.getSpectra());
    }
    else
    {
      base_name = File::removeExtension(File::basename(filename));
      raw_data = "mzML";
    }
    // mz_name is input from IDFileConverter for 'base_name' attribute, only necessary if different from 'mz_file'.
    if (!mz_name.empty())
    {
      base_name = mz_name;
    }
    if (base_name.hasSubstring(".")) // spectrum query name is split by dot, otherwise correct charge can not be read.
    {
      replace(base_name.begin(), base_name.end(), '.', '_');
    }

    f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << "\n";
    f << "<msms_pipeline_analysis date=\"2007-12-05T17:49:46\" xmlns=\"http://regis-web.systemsbiology.net/pepXML\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v117.xsd\" summary_xml=\".xml\">" << "\n";
    f << "<msms_run_summary base_name=\"" << base_name << "\" raw_data_type=\"raw\" raw_data=\"." << raw_data << "\" search_engine=\"" << search_engine_name << "\">" << "\n";
    String enzyme_name = search_params.digestion_enzyme.getName();
    f << "\t<sample_enzyme name=\"";
    f << enzyme_name.toLower() << "\">" << "\n";
    f << "\t\t<specificity cut=\"";
    if (search_params.digestion_enzyme.getRegEx() != "")
    {
      vector<String> sub_regex;
      search_params.digestion_enzyme.getRegEx().split(")",sub_regex);
      boost::match_results<std::string::const_iterator> results;
      static const boost::regex e("(.*?)([A-Z]+)(.*?)");
      if (boost::regex_match(sub_regex[0], results, e))
      {
        f << results[2];
      }
      if (sub_regex[1].hasSubstring("!P"))
      {
        f << "\" no_cut=\"P";
      }
    }
    f << "\" sense=\"C\"/>" << "\n";
    f << "\t</sample_enzyme>" << "\n";

    f << "\t<search_summary base_name=\"" << base_name;
    f << "\" search_engine=\"" << search_engine_name;
    f << "\" precursor_mass_type=\"";
    if (search_params.mass_type == ProteinIdentification::MONOISOTOPIC)
    {
      f << "monoisotopic";
    }
    else
    {
      f << "average";
    }
    f << "\" fragment_mass_type=\"";
    if (search_params.mass_type == ProteinIdentification::MONOISOTOPIC)
    {
      f << "monoisotopic";
    }
    else
    {
      f << "average";
    }
    f << "\" out_data_type=\"\" out_data=\"\" search_id=\"1\">" << "\n";
    f << "\t\t<search_database local_path=\"" << search_params.db << "\" type=\"AA\"/>" << "\n";


    // register modifications
    set<String> aa_mods;
    set<String> n_term_mods, c_term_mods;
    for (vector<PeptideIdentification>::const_iterator it = peptide_ids.begin();
         it != peptide_ids.end(); ++it)
    {
      if (it->getHits().size() > 0)
      {
        PeptideHit h = *it->getHits().begin();

        if (h.getSequence().isModified())
        {
          AASequence p = h.getSequence();
          if (p.hasNTerminalModification())
          {
            n_term_mods.insert(p.getNTerminalModification()->getFullId());
          }
          if (p.hasCTerminalModification())
          {
            c_term_mods.insert(p.getCTerminalModification()->getFullId());
          }

          for (Size i = 0; i != p.size(); ++i)
          {
            if (p[i].isModified())
            {
              aa_mods.insert(p[i].getModification()->getFullId());
            }
          }
        }
      }
    }

    // write modifications definitions
    // <aminoacid_modification aminoacid="C" massdiff="+58.01" mass="161.014664" variable="Y" binary="N" description="Carboxymethyl (C)"/>
    for (set<String>::const_iterator it = aa_mods.begin();
         it != aa_mods.end(); ++it)
    {
      const ResidueModification& mod = ModificationsDB::getInstance()->getModification(*it, "", ResidueModification::ANYWHERE);

      // compute mass of modified residue
      EmpiricalFormula ef = ResidueDB::getInstance()->getResidue(mod.getOrigin())->getFormula(Residue::Internal);
      ef += mod.getDiffFormula();

      f << "\t\t"
        << "<aminoacid_modification aminoacid=\"" << mod.getOrigin()
        << "\" massdiff=\"" << precisionWrapper(mod.getDiffMonoMass()) << "\" mass=\""
        << precisionWrapper(ef.getMonoWeight())
        << "\" variable=\"Y\" binary=\"N\" description=\"" << *it << "\"/>"
        << "\n";
    }

    for (set<String>::const_iterator it = n_term_mods.begin(); it != n_term_mods.end(); ++it)
    {
      const ResidueModification& mod = ModificationsDB::getInstance()->getModification(*it, "", ResidueModification::N_TERM);
      f << "\t\t"
        << "<terminal_modification terminus=\"n\" massdiff=\""
        << precisionWrapper(mod.getDiffMonoMass()) << "\" mass=\"" << precisionWrapper(mod.getMonoMass())
        << "\" variable=\"Y\" description=\"" << *it
        << "\" protein_terminus=\"\"/>" << "\n";
    }

    for (set<String>::const_iterator it = c_term_mods.begin(); it != c_term_mods.end(); ++it)
    {
      const ResidueModification& mod = ModificationsDB::getInstance()->getModification(*it, "", ResidueModification::C_TERM);
      f << "\t\t"
        << "<terminal_modification terminus=\"c\" massdiff=\""
        << precisionWrapper(mod.getDiffMonoMass()) << "\" mass=\"" << precisionWrapper(mod.getMonoMass())
        << "\" variable=\"Y\" description=\"" << *it
        << "\" protein_terminus=\"\"/>" << "\n";
    }

    f << "\t</search_summary>" << "\n";
    if (peptideprophet_analyzed)
    {
      f << "\t<analysis_timestamp analysis=\"peptideprophet\" time=\"2007-12-05T17:49:52\" id=\"1\"/>" << "\n";
    }

    // Scan index and scan number will be reconstructed if no spectrum lookup is possible to retrieve the values.
    // The scan index is generally zero-based and the scan number generally one-based.
    Int count(0);
    for (vector<PeptideIdentification>::const_iterator it = peptide_ids.begin();
         it != peptide_ids.end(); ++it, ++count)
    {
      if (it->getHits().size() == 0)
      {
        continue;
      }
      for (vector<PeptideHit>::const_iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
      {
        PeptideHit h = *hit;
        AASequence seq = h.getSequence();
        double precursor_neutral_mass = seq.getMonoWeight();

        int scan_index = count;
        int scan_nr = 0;

        if (lookup.empty())
        {
          if (it->metaValueExists("RT_index")) // Setting metaValue "RT_index" in XTandemXMLFile in the case of X! Tandem.
          {
            scan_index = it->getMetaValue("RT_index");
          }
          scan_nr = scan_index + 1;
        }
        else
        {
          scan_index = lookup.findByRT(it->getRT());

          SpectrumMetaDataLookup::SpectrumMetaData meta;
          lookup.getSpectrumMetaData(scan_index, meta);
          scan_nr = meta.scan_number;
        }
        // PeptideProphet requires this format for "spectrum" attribute (otherwise TPP parsing error)
        //  - see also the parser code if iProphet at http://sourceforge.net/p/sashimi/code/HEAD/tree/trunk/trans_proteomic_pipeline/src/Validation/InterProphet/InterProphetParser/InterProphetParser.cxx#l180
        //  strictly required attributes:
        //    - spectrum
        //    - assumed_charge
        //  optional attributes
        //    - retention_time_sec
        //    - swath_assay
        //    - experiment_label

        String spectrum_name = base_name + "." + scan_nr + "." + scan_nr + ".";
        if (it->metaValueExists("pepxml_spectrum_name") && keep_native_name_)
        {
          spectrum_name = it->getMetaValue("pepxml_spectrum_name");
        }

        f << "\t<spectrum_query spectrum=\"" << spectrum_name << h.getCharge() << "\""
          << " start_scan=\"" << scan_nr << "\""
          << " end_scan=\"" << scan_nr << "\""
          << " precursor_neutral_mass=\"" << precisionWrapper(precursor_neutral_mass) << "\""
          << " assumed_charge=\"" << h.getCharge() << "\" index=\"" << scan_index << "\"";

        if (it->hasRT())
        {
          f << " retention_time_sec=\"" << it->getRT() << "\" ";
        }

        if (!it->getExperimentLabel().empty())
        {
          f << " experiment_label=\"" << it->getExperimentLabel() << "\" ";
        }

        // "swath_assay" is an optional parameter used for SWATH-MS mostly and
        // may be set for a PeptideIdentification
        //   note that according to the parsing rules of TPP, this needs to be
        //   "xxx:yyy" where xxx is any string and yyy is probably an integer
        //   indicating the Swath window
        if (it->metaValueExists("swath_assay"))
        {
          f << " swath_assay=\"" << it->getMetaValue("swath_assay") << "\" ";
        }
        // "status" is an attribute that may be target or decoy
        if (it->metaValueExists("status"))
        {
          f << " status=\"" << it->getMetaValue("status") << "\" ";
        }

        f << ">\n";
        f << "\t<search_result>" << "\n";

        vector<PeptideEvidence> pes = h.getPeptideEvidences();

        // select first one if multiple are present as "leader"
        PeptideEvidence pe;
        if (!h.getPeptideEvidences().empty())
        {
          pe = pes[0];
        }

        f << "\t\t<search_hit hit_rank=\"1\" peptide=\""
          << seq.toUnmodifiedString() << "\" peptide_prev_aa=\""
          << pe.getAABefore() << "\" peptide_next_aa=\"" << pe.getAAAfter()
          << "\" protein=\"";

        f << pe.getProteinAccession();

        f << "\" num_tot_proteins=\"1\" num_matched_ions=\"0\" tot_num_ions=\"0\" calc_neutral_pep_mass=\"" << precisionWrapper(precursor_neutral_mass)
          << "\" massdiff=\"0.0\" num_tol_term=\"";
        Int num_tol_term = 1;
        if ((pe.getAABefore() == 'R' || pe.getAABefore() == 'K') && search_params.digestion_enzyme.getName() == "Trypsin")
        {
          num_tol_term = 2;
        }
        f << num_tol_term;
        f << "\" num_missed_cleavages=\"0\" is_rejected=\"0\" protein_descr=\"Protein No. 1\">" << "\n";

        // multiple protein hits: <alternative_protein protein="sp|P0CZ86|GLS24_STRP3" num_tol_term="2" peptide_prev_aa="K" peptide_next_aa="-"/>
        if (pes.size() > 1)
        {
          for (Size k = 1; k != pes.size(); ++k)
          {
            f << "\t\t<alternative_protein protein=\"" << pes[k].getProteinAccession() << "\" num_tol_term=\"" << num_tol_term << "\"";

            if (pes[k].getAABefore() != PeptideEvidence::UNKNOWN_AA)
            {
              f << " peptide_prev_aa=\"" << pes[k].getAABefore() << "\"";
            }

            if (pes[k].getAAAfter() != PeptideEvidence::UNKNOWN_AA)
            {
              f << " peptide_next_aa=\"" << pes[k].getAAAfter() << "\"";
            }
            f << "/>" << "\n";
          }
        }
        if (seq.isModified())
        {
          f << "\t\t\t<modification_info modified_peptide=\""
            << seq.toBracketString() << "\"";

          if (seq.hasNTerminalModification())
          {
            const ResidueModification& mod = *(seq.getNTerminalModification());
            const double mod_nterm_mass = Residue::getInternalToNTerm().getMonoWeight() + mod.getDiffMonoMass();
            f << " mod_nterm_mass=\"" << precisionWrapper(mod_nterm_mass) << "\"";
          }

          if (seq.hasCTerminalModification())
          {
            const ResidueModification& mod = *(seq.getCTerminalModification());
            const double mod_cterm_mass = Residue::getInternalToCTerm().getMonoWeight() + mod.getDiffMonoMass();
            f << " mod_cterm_mass=\"" << precisionWrapper(mod_cterm_mass) << "\"";
          }

          f << ">" << "\n";

          for (Size i = 0; i != seq.size(); ++i)
          {
            if (seq[i].isModified())
            {
              const ResidueModification& mod = *(seq[i].getModification());
              // the modification position is 1-based
              f << "\t\t\t\t<mod_aminoacid_mass position=\"" << (i + 1)
                << "\" mass=\"" <<
                precisionWrapper(mod.getMonoMass() + seq[i].getMonoWeight(Residue::Internal)) << "\"/>" << "\n";
            }
          }

          f << "\t\t\t</modification_info>" << "\n";
        }

        // write out the (optional) search_score_summary that may be associated with peptide prophet results
        bool peptideprophet_written = false;
        if (!h.getAnalysisResults().empty())
        {
          // <analysis_result analysis="peptideprophet">
          //   <peptideprophet_result probability="0.0660" all_ntt_prob="(0.0000,0.0000,0.0660)">
          //     <search_score_summary>
          //       <parameter name="fval" value="0.7114"/>
          //       <parameter name="ntt" value="2"/>
          //       <parameter name="nmc" value="0"/>
          //       <parameter name="massd" value="-0.027"/>
          //       <parameter name="isomassd" value="0"/>
          //     </search_score_summary>
          //   </peptideprophet_result>
          // </analysis_result>

          for (std::vector<PeptideHit::PepXMLAnalysisResult>::const_iterator ar_it = h.getAnalysisResults().begin();
              ar_it != h.getAnalysisResults().end(); ++ar_it)
          {
            f << "\t\t\t<analysis_result analysis=\"" << ar_it->score_type << "\">" << "\n";

            // get name of next tag
            String tagname = "peptideprophet_result";
            if (ar_it->score_type == "peptideprophet")
            {
              peptideprophet_written = true; // remember that we have now already written peptide prophet results
              tagname = "peptideprophet_result";
            }
            else if (ar_it->score_type == "interprophet")
            {
              tagname = "interprophet_result";
            }
            else
            {
              peptideprophet_written = true; // remember that we have now already written peptide prophet results
              warning(STORE, "Analysis type " + ar_it->score_type + " not supported, will use peptideprophet_result.");
            }

            f << "\t\t\t\t<" << tagname <<  " probability=\"" << ar_it->main_score;
            // TODO
            f << "\" all_ntt_prob=\"(" << ar_it->main_score << "," << ar_it->main_score
            << "," << ar_it->main_score << ")\">" << "\n";

            if (!ar_it->sub_scores.empty())
            {
              f << "\t\t\t\t\t<search_score_summary>" << "\n";
              for (std::map<String, double>::const_iterator subscore_it = ar_it->sub_scores.begin();
                  subscore_it != ar_it->sub_scores.end(); ++subscore_it)
              {
                f << "\t\t\t\t\t\t<parameter name=\""<< subscore_it->first << "\" value=\"" << subscore_it->second << "\"/>\n";
              }
              f << "\t\t\t\t\t</search_score_summary>" << "\n";
            }
            f << "\t\t\t\t</" << tagname << ">" << "\n";

            f << "\t\t\t</analysis_result>" << "\n";
          }
        }

        // deprecated way of writing out peptide prophet results (only if
        // requested explicitly and if we have not already written out the
        // peptide prophet results above through AnalysisResults
        if (peptideprophet_analyzed && !peptideprophet_written)
        {
          // if (!h.getAnalysisResults().empty()) { WARNING / }
          f << "\t\t\t<analysis_result analysis=\"peptideprophet\">" << "\n";
          f << "\t\t\t<peptideprophet_result probability=\"" << h.getScore()
            << "\" all_ntt_prob=\"(" << h.getScore() << "," << h.getScore()
            << "," << h.getScore() << ")\">" << "\n";

          f << "\t\t\t</peptideprophet_result>" << "\n";
          f << "\t\t\t</analysis_result>" << "\n";
        }
        else
        {
          if (search_engine_name == "X! Tandem")
          {
            // check if score type is XTandem or qvalue/fdr
            if (it->getScoreType() == "XTandem")
            {
              f << "\t\t\t<search_score" << " name=\"hyperscore\" value=\"" << h.getScore() << "\"" << "/>\n";
              f << "\t\t\t<search_score" << " name=\"nextscore\" value=\"";
              if (h.metaValueExists("nextscore"))
              {
                f << h.getMetaValue("nextscore") << "\"" << "/>\n";
              }
              else
              {
                f << h.getScore() << "\"" << "/>\n";
              }
            }
            else if (h.metaValueExists("XTandem_score"))
            {
              f << "\t\t\t<search_score" << " name=\"hyperscore\" value=\"" << h.getMetaValue("XTandem_score") << "\"" << "/>\n";
              f << "\t\t\t<search_score" << " name=\"nextscore\" value=\"";
              if (h.metaValueExists("nextscore"))
              {
                f << h.getMetaValue("nextscore") << "\"" << "/>\n";
              }
              else
              {
                f << h.getMetaValue("XTandem_score") << "\"" << "/>\n";
              }
            }
            f << "\t\t\t<search_score" << " name=\"expect\" value=\"" << h.getMetaValue("E-Value") << "\"" << "/>\n";
          }
          else if (search_engine_name == "MASCOT")
          {
            f << "\t\t\t<search_score" << " name=\"expect\" value=\"" << h.getMetaValue("EValue") << "\"" << "/>\n";
            f << "\t\t\t<search_score" << " name=\"ionscore\" value=\"" << h.getScore() << "\"" << "/>\n";
          }
          else if (search_engine_name == "OMSSA")
          {
            f << "\t\t\t<search_score" << " name=\"expect\" value=\"" << h.getScore() << "\"" << "/>\n";
          }
          else if (search_engine_name == "MSGFPlus")
          {
            f << "\t\t\t<search_score" << " name=\"expect\" value=\"" << h.getScore() << "\"" << "/>\n";
          }
          else if (search_engine_name == "Percolator")
          {
            double pep_score = static_cast<double>(h.getMetaValue("Percolator_PEP"));
            f << "\t\t\t<search_score" << " name=\"Percolator_score\" value=\"" << h.getMetaValue("Percolator_score") << "\"" << "/>\n";
            f << "\t\t\t<search_score" << " name=\"Percolator_qvalue\" value=\"" << h.getMetaValue("Percolator_qvalue") << "\"" << "/>\n";
            f << "\t\t\t<search_score" << " name=\"Percolator_PEP\" value=\"" << pep_score << "\"" << "/>\n";

            double probability = 1.0 - pep_score;
            f << "\t\t\t<analysis_result" << " analysis=\"peptideprophet\">\n";
            f << "\t\t\t\t<peptideprophet_result" << " probability=\"" << probability << "\"";
            f << " all_ntt_prob=\"(0.0000,0.0000," << probability << ")\"/>\n";
            f << "\t\t\t</analysis_result>" << "\n";
          }
          else
          {
            f << "\t\t\t<search_score" << " name=\"" << it->getScoreType() << "\" value=\"" << h.getScore() << "\"" << "/>\n";
          }
          if (it->getScoreType() == "Posterior Error Probability")
          {
            f << "\t\t\t<search_score" << " name=\"" << it->getScoreType() << "\" value=\"" << h.getScore() << "\"" << "/>\n";
            double probability = 1.0 - h.getScore();
            f << "\t\t\t<analysis_result" << " analysis=\"peptideprophet\">\n";
            f << "\t\t\t\t<peptideprophet_result" << " probability=\"" << probability << "\"";
            f << " all_ntt_prob=\"(0.0000,0.0000," << probability << ")\"/>\n";
            f << "\t\t\t</analysis_result>" << "\n";
          }
        }
        f << "\t\t</search_hit>" << "\n";
        f << "\t</search_result>" << "\n";
        f << "\t</spectrum_query>" << "\n";
      }
    }
    f << "</msms_run_summary>" << "\n";
    f << "</msms_pipeline_analysis>" << "\n";

    f.close();
  }

  void PepXMLFile::matchModification_(const double mass, const String& origin, String& modification_description)
  {
    double mod_mass = mass - ResidueDB::getInstance()->getResidue(origin)->getMonoWeight(Residue::Internal);
    vector<String> mods;
    // try more specific search first:
    ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, mod_mass, mod_tol_, origin, ResidueModification::ANYWHERE);
    if (mods.empty()) ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, mod_mass, mod_tol_, origin);

    // no notification about ambiguities here - that was done when the
    // modification definitions were parsed ("aminoacid_modification" and
    // "terminal_modification" elements)
    if (!mods.empty()) modification_description = mods[0];
  }

  void PepXMLFile::readRTMZCharge_(const xercesc::Attributes& attributes)
  {
    double mass = attributeAsDouble_(attributes, "precursor_neutral_mass");
    charge_ = attributeAsInt_(attributes, "assumed_charge");
    mz_ = (mass + hydrogen_mass_ * charge_) / charge_;
    rt_ = 0;

    bool rt_present = optionalAttributeAsDouble_(rt_, attributes, "retention_time_sec");

    if (!rt_present) // get RT from experiment
    {
      if (lookup_ == nullptr || lookup_->empty())
      {
        // no lookup given, report non-fatal error
        error(LOAD, "Cannot get RT information - no spectra given");
        return;
      }

      // assume only one scan, i.e. ignore "end_scan":
      Size scan = attributeAsInt_(attributes, "start_scan");
      Size index = (scan != 0) ? lookup_->findByScanNumber(scan) :
        lookup_->findByReference(attributeAsString_(attributes, "spectrum"));
      SpectrumMetaDataLookup::SpectrumMetaData meta;
      lookup_->getSpectrumMetaData(index, meta);
      if (meta.ms_level == 2)
      {
        rt_ = meta.rt;
      }
      else
      {
        error(LOAD, "Cannot get RT information - scan mapping is incorrect");
      }
    }
  }

  void PepXMLFile::load(const String& filename, vector<ProteinIdentification>&
                        proteins, vector<PeptideIdentification>& peptides,
                        const String& experiment_name)
  {
    SpectrumMetaDataLookup lookup;
    load(filename, proteins, peptides, experiment_name, lookup);
  }

  void PepXMLFile::load(const String& filename, vector<ProteinIdentification>&
                        proteins, vector<PeptideIdentification>& peptides,
                        const String& experiment_name,
                        const SpectrumMetaDataLookup& lookup)
  {
    // initialize here, since "load" could be called several times:
    exp_name_ = "";
    prot_id_ = "";
    charge_ = 0;
    peptides.clear();
    peptides_ = &peptides;
    proteins.clear();
    proteins_ = &proteins;
    // assume mass type "average" (in case element "search_summary" is missing):
    hydrogen_mass_ = hydrogen_.getAverageWeight();

    file_ = filename; // filename for error messages in XMLHandler

    if (experiment_name != "")
    {
      exp_name_ = File::removeExtension(experiment_name);
      lookup_ = &lookup;
    }

    analysis_summary_ = false;
    wrong_experiment_ = false;
    // without experiment name, don't care about these two:
    seen_experiment_ = exp_name_.empty();
    checked_base_name_ = exp_name_.empty();

    parse_(filename, this);

    if (!seen_experiment_)
    {
      fatalError(LOAD, "Found no experiment with name '" + experiment_name + "'");
    }
    // clean up duplicate ProteinHits in each ProteinIdentification separately:
    // (can't use "sort" and "unique" because no "op<" defined for ProteinHit)
    for (vector<ProteinIdentification>::iterator prot_it = proteins.begin();
         prot_it != proteins.end(); ++prot_it)
    {
      set<String> accessions;
      // modeled after "remove_if" in STL header "algorithm":
      vector<ProteinHit>::iterator first = prot_it->getHits().begin();
      vector<ProteinHit>::iterator result = first;
      for (; first != prot_it->getHits().end(); ++first)
      {
        String accession = first->getAccession();
        bool new_element = accessions.insert(accession).second;
        if (new_element) // don't remove
        {
          *result++ = *first;
        }
      }
      prot_it->getHits().erase(result, first);
    }

    // reset members
    exp_name_.clear();
    prot_id_.clear();
    date_.clear();
    proteins_ = nullptr;
    peptides_ = nullptr;
    lookup_ = nullptr;
    scan_map_.clear();
  }

  /*
    NOTE: numbering schemes for multiple searches
    ---------------------------------------------

    pepXML can contain search results for multiple files and for multiple
    searches of those files. We have seen two different variants (both seem to
    be valid pepXML):

    1. One "msms_run_summary" per searched file, containing multiple
    "search_summary" elements (for each search); counting starts at "1" in each
    "msms_run_summary"; produced e.g. by the TPP:

    <msms_run_summary basename="File1">
      <search_summary search_engine="A" search_id="1">
      <search_summary search_engine="B" search_id="2">
    ...
    <msms_run_summary basename="File2">
      <search_summary search_engine="A" search_id="1">
      <search_summary search_engine="B" search_id="2">
    ...

    2. One "msms_run_summary" per search of a file, containing one
    "search_summary" element; counting is sequential across "msms_run_summary"
    sections; produced e.g. by ProteomeDiscoverer:

    <msms_run_summary basename="File1">
      <search_summary search_engine="A" search_id="1">
    ...
    <msms_run_summary basename="File1">
      <search_summary search_engine="B" search_id="2">
    ...

    The "search_id" numbers are necessary to associate search hits with the
    correct search runs. In the parser, we keep track of the search runs per
    "msms_run_summary" in the "current_proteins_" vector. Importantly, this
    means that for variant 2 of the numbering, the "search_id" number may be
    higher than the number of elements in the vector! However, in this case we
    know that the last - and only - element should be the correct one.
  */

  void PepXMLFile::startElement(const XMLCh* const /*uri*/,
                                const XMLCh* const /*local_name*/,
                                const XMLCh* const qname,
                                const xercesc::Attributes& attributes)
  {
    String element = sm_.convert(qname);

    // cout << "Start: " << element << "\n";

    if (element == "msms_run_summary") // parent: "msms_pipeline_analysis"
    {
      if (!exp_name_.empty())
      {
        String base_name = attributeAsString_(attributes, "base_name");
        if (!base_name.empty())
        {
          wrong_experiment_ = !base_name.hasSuffix(exp_name_);
          seen_experiment_ = seen_experiment_ || !wrong_experiment_;
          checked_base_name_ = true;
        }
        else // really shouldn't happen, but does for Mascot export to pepXML
        {
          error(LOAD, "'base_name' attribute of 'msms_run_summary' element is empty");
          // continue for now, but check later in 'search_summary':
          wrong_experiment_ = false;
          checked_base_name_ = false;
        }
      }
      if (wrong_experiment_) return;

      // create a ProteinIdentification in case "search_summary" is missing:
      ProteinIdentification protein;
      protein.setDateTime(date_);
      prot_id_ = "unknown_" + date_.getDate();
      enzyme_ = "unknown_enzyme";
      // "prot_id_" will be overwritten if elem. "search_summary" is present
      protein.setIdentifier(prot_id_);
      proteins_->push_back(protein);
      current_proteins_.clear();
      current_proteins_.push_back(--proteins_->end());
    }
    else if (element == "analysis_summary") // parent: "msms_pipeline_analysis"
    {
      // this element can contain "search summary" elements, which we only
      // expect as subelements of "msms_run_summary", so skip the whole thing
      analysis_summary_ = true;
    }
    else if (wrong_experiment_ || analysis_summary_)
    {
      // do nothing here (this case exists to prevent parsing of elements for
      // experiments we're not interested in or for analysis summaries)
    }
    // now, elements occurring more frequently are generally closer to the top
    else if (element == "search_score") // parent: "search_hit"
    {
      String name = attributeAsString_(attributes, "name");
      double value;

      // TODO: deal with different scores
      if (name == "expect") // X!Tandem or Mascot E-value
      {
        value = attributeAsDouble_(attributes, "value");
        peptide_hit_.setScore(value);
        current_peptide_.setScoreType(name);
        current_peptide_.setHigherScoreBetter(false);
        if (search_engine_ == "Comet")
        {
          peptide_hit_.setMetaValue("MS:1002257", value); // name: Comet:expectation value
        }
        else if (search_engine_ == "X! Tandem")
        {
          peptide_hit_.setMetaValue("MS:1001330", value); // name: X\!Tandem:expect
        }
        else if (search_engine_ == "Mascot")
        {
          peptide_hit_.setMetaValue("MS:1001172", value); // name: Mascot:expectation value
        }
        //TODO: there is no (generic) umbrella term for expect val in the CV right now
      }
      else if (name == "mvh") // MyriMatch score
      {
        value = attributeAsDouble_(attributes, "value");
        peptide_hit_.setScore(value);
        current_peptide_.setScoreType(name);
        current_peptide_.setHigherScoreBetter(true);
      }
      // if (name == "hyperscore")
      // { // X!Tandem score
      //   value = attributeAsDouble_(attributes, "value");
      //   peptide_hit_.setScore(value);
      //   current_peptide_.setScoreType(name); // add "X!Tandem" to name?
      //   current_peptide_.setHigherScoreBetter(true);
      // }
      else if (name == "xcorr") // Sequest score
      { // and use the mvh
        value = attributeAsDouble_(attributes, "value");
        if (search_engine_ != "MyriMatch") //MyriMatch has also an xcorr, but we want to ignore it
        {
          peptide_hit_.setScore(value);
          current_peptide_.setScoreType(name); // add "Sequest" to name?
          current_peptide_.setHigherScoreBetter(true);
        }
        if (search_engine_ == "Comet")
        {
          peptide_hit_.setMetaValue("MS:1002252", value); // name: Comet:xcorr
        }
        else
        {
          peptide_hit_.setMetaValue("MS:1001155", value); // name: SEQUEST:xcorr
        }
        //TODO: no other xcorr or generic xcorr in the CV right now, use SEQUEST:xcorr
      }
      else if (name == "fval") // SpectraST score
      {
        value = attributeAsDouble_(attributes, "value");
        peptide_hit_.setScore(value);
        current_peptide_.setScoreType(name);
        current_peptide_.setHigherScoreBetter(true);
        peptide_hit_.setMetaValue("MS:1001419", value); // def: "SpectraST spectrum score.
      }
      else
      {
        if (search_engine_ == "Comet")
        {
          if (name == "deltacn")
          {
            value = attributeAsDouble_(attributes, "value");
            peptide_hit_.setMetaValue("MS:1002253", value); // name: Comet:deltacn
          }
          else if (name == "spscore")
          {
            value = attributeAsDouble_(attributes, "value");
            peptide_hit_.setMetaValue("MS:1002255", value); // name: Comet:spscore
          }
          else if (name == "sprank")
          {
            value = attributeAsDouble_(attributes, "value");
            peptide_hit_.setMetaValue("MS:1002256", value); // name: Comet:sprank
          }
        }
      }
    }
    else if (element == "search_hit") // parent: "search_result"
    { // creates a new PeptideHit
      current_sequence_ = attributeAsString_(attributes, "peptide");
      current_modifications_.clear();
      PeptideEvidence pe;
      peptide_hit_ = PeptideHit();
      peptide_hit_.setRank(attributeAsInt_(attributes, "hit_rank"));
      peptide_hit_.setCharge(charge_); // from parent "spectrum_query" tag
      String prev_aa, next_aa;
      if (optionalAttributeAsString_(prev_aa, attributes, "peptide_prev_aa"))
      {
        pe.setAABefore(prev_aa[0]);
      }

      if (optionalAttributeAsString_(next_aa, attributes, "peptide_next_aa"))
      {
        pe.setAAAfter(next_aa[0]);
      }
      if (search_engine_ == "Comet")
      {
        String value;
        if (optionalAttributeAsString_(value, attributes, "num_matched_ions"))
        {
          peptide_hit_.setMetaValue("MS:1002258", value); // name: Comet:matched ions

        }
        if (optionalAttributeAsString_(value, attributes, "tot_num_ions"))
        {
          peptide_hit_.setMetaValue("MS:1002259", value); // name: Comet:total ions
        }
        if (optionalAttributeAsString_(value, attributes, "num_matched_peptides"))
        {
          peptide_hit_.setMetaValue("num_matched_peptides", value);
        }
      }
      String protein = attributeAsString_(attributes, "protein");
      pe.setProteinAccession(protein);
      peptide_hit_.addPeptideEvidence(pe);
      ProteinHit hit;
      hit.setAccession(protein);
      // depending on the numbering scheme used in the pepXML, "search_id_"
      // may appear to be "out of bounds" - see NOTE above:
      current_proteins_[min(UInt(current_proteins_.size()), search_id_) - 1]->insertHit(hit);
    }
    else if (element == "search_result") // parent: "spectrum_query"
    {
      // creates a new PeptideIdentification
      current_peptide_ = PeptideIdentification();
      current_peptide_.setRT(rt_);
      current_peptide_.setMZ(mz_);
      current_peptide_.setBaseName(current_base_name_);

      search_id_ = 1; // default if attr. is missing (ref. to "search_summary")
      optionalAttributeAsUInt_(search_id_, attributes, "search_id");
      // depending on the numbering scheme used in the pepXML, "search_id_"
      // may appear to be "out of bounds" - see NOTE above:
      String identifier = current_proteins_[min(UInt(current_proteins_.size()), search_id_) - 1]->getIdentifier();
      current_peptide_.setIdentifier(identifier);

      // set optional attributes
      if (!native_spectrum_name_.empty() && keep_native_name_)
      {
        current_peptide_.setMetaValue("pepxml_spectrum_name", native_spectrum_name_);
      }
      if (search_engine_ == "Comet")
      {
        current_peptide_.setMetaValue("spectrum_reference", native_spectrum_name_);
        //TODO: we really need something uniform here, like scan number - and not in metainfointerface
      }
      if (!experiment_label_.empty())
      {
        current_peptide_.setExperimentLabel(experiment_label_);
      }
      if (!swath_assay_.empty())
      {
        current_peptide_.setMetaValue("swath_assay", swath_assay_);
      }
      if (!status_.empty())
      {
        current_peptide_.setMetaValue("status", status_);
      }
    }
    else if (element == "spectrum_query") // parent: "msms_run_summary"
    {
      // sample:
      // <spectrum_query spectrum="foobar.02552.02552.2" start_scan="2552"
      // end_scan="2552" precursor_neutral_mass="1168.6176" assumed_charge="2"
      // index="10" retention_time_sec="488.652" experiment_label="urine"
      // swath_assay="EIVLTQSPGTL2:9" status="target">

      readRTMZCharge_(attributes); // sets "rt_", "mz_", "charge_"

      // retrieve optional attributes
      native_spectrum_name_ = "";
      experiment_label_ = "";
      swath_assay_ = "";
      status_ = "";
      optionalAttributeAsString_(native_spectrum_name_, attributes, "spectrum");
      optionalAttributeAsString_(native_spectrum_name_, attributes, "spectrumNativeID"); //some engines write that optional attribute - is preferred to spectrum
      optionalAttributeAsString_(experiment_label_, attributes, "experiment_label");
      optionalAttributeAsString_(swath_assay_, attributes, "swath_assay");
      optionalAttributeAsString_(status_, attributes, "status");


    }
    else if (element == "analysis_result") // parent: "search_hit"
    {
      current_analysis_result_ = PeptideHit::PepXMLAnalysisResult();
      current_analysis_result_.score_type = attributeAsString_(attributes, "analysis");
    }
    else if (element == "search_score_summary")
    {
      search_score_summary_ = true;
    }
    else if (element == "parameter") // parent: "search_score_summary"
    {
      // If we are within a search_score_summary, add the read in values to the current AnalysisResult
      if (search_score_summary_)
      {
        String name = attributeAsString_(attributes, "name");
        double value = attributeAsDouble_(attributes, "value");
        current_analysis_result_.sub_scores[name] = value;
      }
      else if (search_summary_)
      {
        String name = attributeAsString_(attributes, "name");
        if (name == "fragment_bin_tol")
        {
          double value = attributeAsDouble_(attributes, "value");
          params_.fragment_mass_tolerance = value/2.0;
          params_.fragment_mass_tolerance_ppm = false;
        }
        if (name == "peptide_mass_tolerance")
        {
          double value = attributeAsDouble_(attributes, "value");
          params_.precursor_mass_tolerance = value;
        }
        // this is quite comet specific, but so is parameter name peptide_mass_units, see comet configuration file documentation
        if (name == "peptide_mass_units")
        {
          int value = attributeAsDouble_(attributes, "value");
          switch (value) {
          case 0: // comet value 0 type amu
              params_.precursor_mass_tolerance_ppm = false;
              break;
          case 1: // comet value 1 type mmu
              params_.precursor_mass_tolerance_ppm = false;
              break;
          case 2: // comet value 1 type ppm
              params_.precursor_mass_tolerance_ppm = true;
              break;
          default:
              break;
          }
        }
      }
      else
      {
        // currently not handled
      }
    }
    else if (element == "peptideprophet_result") // parent: "analysis_result" (in "search_hit")
    {
      // PeptideProphet probability overwrites original search score
      // maybe TODO: deal with meta data associated with PeptideProphet search
      double value = attributeAsDouble_(attributes, "probability");
      if (current_peptide_.getScoreType() != "InterProphet probability")
      {
        peptide_hit_.setScore(value);
        current_peptide_.setScoreType("PeptideProphet probability");
        current_peptide_.setHigherScoreBetter(true);
      }
      current_analysis_result_.main_score = value;
      current_analysis_result_.higher_is_better = true;
    }
    else if (element == "interprophet_result") // parent: "analysis_result" (in "search_hit")
    {
      // InterProphet probability overwrites PeptideProphet probability and
      // original search score
      double value = attributeAsDouble_(attributes, "probability");
      peptide_hit_.setScore(value);
      current_peptide_.setScoreType("InterProphet probability");
      current_peptide_.setHigherScoreBetter(true);
      current_analysis_result_.main_score = value;
      current_analysis_result_.higher_is_better = true;
    }
    else if (element == "modification_info") // parent: "search_hit" (in "search result")
    {
      // Has N-Term Modification
      double mod_nterm_mass;
      if (optionalAttributeAsDouble_(mod_nterm_mass, attributes, "mod_nterm_mass")) // this specifies a terminal modification
      {
        // look up the modification in the search_summary by mass
        for (vector<AminoAcidModification>::const_iterator it = variable_modifications_.begin(); it != variable_modifications_.end(); ++it)
        {
          if (mod_nterm_mass == it->mass && it->terminus == "n")
          {
            double massdiff = (it->massdiff).toDouble();
            vector<String> mods;
            ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, massdiff, mod_tol_, "", ResidueModification::N_TERM);
            if (!mods.empty())
            {
              current_modifications_.push_back(make_pair(mods[0], 42)); // 42, because position does not matter
            }
            else
            {
              error(LOAD, String("Cannot find N-terminal modification with mass " + String(mod_nterm_mass) + "."));
            }
            break; // only one modification should match, so we can stop the loop here
          }
        }
      }

      // Has C-Term Modification
      double mod_cterm_mass;
      if (optionalAttributeAsDouble_(mod_cterm_mass, attributes, "mod_cterm_mass")) // this specifies a terminal modification
      {
        // look up the modification in the search_summary by mass
        for (vector<AminoAcidModification>::const_iterator it = variable_modifications_.begin(); it != variable_modifications_.end(); ++it)
        {
          if (mod_cterm_mass == it->mass && it->terminus == "c")
          {
            double massdiff = (it->massdiff).toDouble();
            vector<String> mods;
            ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, massdiff, mod_tol_, "", ResidueModification::C_TERM);
            if (!mods.empty())
            {
              current_modifications_.push_back(make_pair(mods[0], 42)); // 42, because position does not matter
            }
            else
            {
              error(LOAD, String("Cannot find C-terminal modification with mass " + String(mod_cterm_mass) + "."));
            }
            break; // only one modification should match, so we can stop the loop here
          }
        }
      }
    }
    else if (element == "alternative_protein") // parent: "search_hit"
    {
      String protein = attributeAsString_(attributes, "protein");
      PeptideEvidence pe;
      pe.setProteinAccession(protein);
      peptide_hit_.addPeptideEvidence(pe);
      ProteinHit hit;
      hit.setAccession(protein);
      // depending on the numbering scheme used in the pepXML, "search_id_"
      // may appear to be "out of bounds" - see NOTE above:
      current_proteins_[min(UInt(current_proteins_.size()), search_id_) - 1]->
      insertHit(hit);
    }
    else if (element == "mod_aminoacid_mass") // parent: "modification_info" (in "search_hit")
    {
      double modification_mass = attributeAsDouble_(attributes, "mass");
      Size modification_position = attributeAsInt_(attributes, "position");
      String origin = String(current_sequence_[modification_position - 1]);
      String temp_description = "";

      matchModification_(modification_mass, origin, temp_description);

      if (!temp_description.empty())
      {
        // the modification position is 1-based
        current_modifications_.push_back(make_pair(temp_description, modification_position));
      }
      else
      {
        error(LOAD, String("Cannot find modification '") + String(modification_mass) + "' of residue " + String(origin) + " at position " + String(modification_position) + " in '" + current_sequence_ + "'");
      }
    }
    else if (element == "aminoacid_modification" || element == "terminal_modification") // parent: "search_summary"
    {
      AminoAcidModification aa_mod;
      optionalAttributeAsString_(aa_mod.description, attributes, "description");
      aa_mod.massdiff = attributeAsString_(attributes, "massdiff");
      aa_mod.mass = attributeAsDouble_(attributes, "mass");
      String is_variable = attributeAsString_(attributes, "variable");
      ResidueModification::TermSpecificity term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY;
      if (element == "aminoacid_modification")
      {
        aa_mod.aminoacid = attributeAsString_(attributes, "aminoacid");
        // can't set term_spec to "ANYWHERE", because terminal mods. may not be registered as such!
      }
      else
      {
        // Somehow very small fixed modifications (electron mass?) get annotated by X!Tandem. Don't add them as they interfere with other mods.
        if (fabs(aa_mod.massdiff.toDouble()) < xtandem_artificial_mod_tol_)
        {
          return;
        }

        optionalAttributeAsString_(aa_mod.aminoacid, attributes, "aminoacid");
        aa_mod.terminus = String(attributeAsString_(attributes, "terminus")).toLower();
        if (aa_mod.terminus == "n")
        {
          term_spec = ResidueModification::N_TERM;
        }
        else if (aa_mod.terminus == "c")
        {
          term_spec = ResidueModification::C_TERM;
        }
      }
      String desc = "";
      // check if the modification is uniquely defined:
      if (!aa_mod.description.empty())
      {
        try
        {
          desc = ModificationsDB::getInstance()->getModification(aa_mod.description, aa_mod.aminoacid, term_spec).getFullId();
        }
        catch (Exception::BaseException)
        {
          error(LOAD, "Modification '" + aa_mod.description + "' of residue '" + aa_mod.aminoacid + "' could not be matched. Trying by modification mass.");
        }
      }
      else
      {
        error(LOAD, "No modification description given. Trying to define by modification mass.");
      }
      if (desc.empty())
      {
        vector<String> mods;
        if (term_spec == ResidueModification::NUMBER_OF_TERM_SPECIFICITY) // try more specific search first
        {
          ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(
            mods, aa_mod.massdiff.toDouble(), mod_tol_, aa_mod.aminoacid, ResidueModification::ANYWHERE);
        }
        if (mods.empty())
        {
          ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(
            mods, aa_mod.massdiff.toDouble(), mod_tol_, aa_mod.aminoacid, term_spec);
        }
        if (mods.empty() && aa_mod.massdiff.toDouble() != 0)
        {
          desc = aa_mod.aminoacid;
          if (aa_mod.massdiff.toDouble() >= 0)
          {
            desc += "+" + String(aa_mod.massdiff);
          }
          else
          {
            desc += String(aa_mod.massdiff);
          }
          //Modification unknown, but trying to continue as we want to be able to read the rest despite of the modifications but warning this will fail downstream
          error(LOAD, "Modification '" + String(aa_mod.mass) + "/delta " + String(aa_mod.massdiff) + "' is unknown. Resuming with '" + desc +  "', which could lead to failures using the data downstream.");
        }
        else if (!mods.empty())
        {
          desc = mods[0];
          if (mods.size() > 1)
          {
            String mod_str = mods[0];
            for (vector<String>::const_iterator mit = ++mods.begin(); mit != mods.end(); ++mit)
            {
              mod_str += ", " + *mit;
            }
            error(LOAD, "Modification '" + String(aa_mod.mass) + "' is not uniquely defined by the given data. Using '" + mods[0] +  "' to represent any of '" + mod_str + "'.");
          }
        }
      }
      if (!desc.empty())
      {
        if (is_variable == "Y")
        {
          variable_modifications_.push_back(aa_mod);
          params_.variable_modifications.push_back(desc);
        }
        else
        {
          fixed_modifications_.push_back(aa_mod);
          params_.fixed_modifications.push_back(desc);
        }
      }
    }
    else if (element == "search_summary") // parent: "msms_run_summary"
    { // creates a new ProteinIdentification (usually)
      search_summary_ = true;
      current_base_name_ = "";
      optionalAttributeAsString_(current_base_name_, attributes, "base_name");
      if (!checked_base_name_) // work-around for files exported by Mascot
      {
        if (current_base_name_.hasSuffix(exp_name_))
        {
          seen_experiment_ = true;
        }
        else // wrong experiment after all - roll back changes that were made
        {
          proteins_->pop_back();
          current_proteins_.clear();
          wrong_experiment_ = true;
          return;
        }
      }

      fixed_modifications_.clear();
      variable_modifications_.clear();
      params_ = ProteinIdentification::SearchParameters();
      params_.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(enzyme_));
      String mass_type = attributeAsString_(attributes, "precursor_mass_type");
      if (mass_type == "monoisotopic")
      {
        hydrogen_mass_ = hydrogen_.getMonoWeight();
      }
      else
      {
        hydrogen_mass_ = hydrogen_.getAverageWeight();
        if (mass_type != "average")
        {
          error(LOAD, "'precursor_mass_type' attribute of 'search_summary' tag should be 'monoisotopic' or 'average', not '" + mass_type + "' (assuming 'average')");
        }
      }
      // assuming "SearchParameters::mass_type" refers to the fragment mass
      mass_type = attributeAsString_(attributes, "fragment_mass_type");
      if (mass_type == "monoisotopic")
      {
        params_.mass_type = ProteinIdentification::MONOISOTOPIC;
      }
      else
      {
        if (mass_type == "average")
        {
          params_.mass_type = ProteinIdentification::AVERAGE;
        }
        else
        {
          error(LOAD, "'fragment_mass_type' attribute of 'search_summary' tag should be 'monoisotopic' or 'average', not '" + mass_type + "'");
        }
      }

      search_engine_ = attributeAsString_(attributes, "search_engine");

      // generate a unique identifier for every search engine run.
      prot_id_ = search_engine_ + "_" + date_.getDate() + "_" + date_.getTime();

      search_id_ = 1;
      optionalAttributeAsUInt_(search_id_, attributes, "search_id");
      vector<ProteinIdentification>::iterator prot_it;
      if (search_id_ <= proteins_->size()) // ProteinIdent. was already created for "msms_run_summary" -> add to it
      {
        prot_it = current_proteins_.back();
      }
      else // create a new ProteinIdentification
      {
        ProteinIdentification protein;
        protein.setDateTime(date_);
        proteins_->push_back(protein);
        prot_it = --proteins_->end();
        prot_id_ = prot_id_ + "_" + search_id_; // make sure the ID is unique
      }
      prot_it->setSearchEngine(search_engine_);
      prot_it->setIdentifier(prot_id_);
    }
    else if (element == "sample_enzyme") // parent: "msms_run_summary"
    { // special case: search parameter that occurs *before* "search_summary"!
      enzyme_ = attributeAsString_(attributes, "name");
      if (enzyme_ == "nonspecific") enzyme_ = "unspecific cleavage";
      if (ProteaseDB::getInstance()->hasEnzyme(enzyme_.toLower()))
      {
        params_.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(enzyme_));
      }
    }
    else if (element == "enzymatic_search_constraint") // parent: "search_summary"
    {
      ///<enzymatic_search_constraint enzyme="nonspecific" max_num_internal_cleavages="1" min_number_termini="2"/>
      enzyme_ = attributeAsString_(attributes, "enzyme");
      if (enzyme_ == "nonspecific") enzyme_ = "unspecific cleavage";
      if (ProteaseDB::getInstance()->hasEnzyme(enzyme_))
      {
        params_.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(enzyme_.toLower()));
      }

      int mc = attributeAsInt_(attributes, "max_num_internal_cleavages");
      params_.missed_cleavages = mc;
    }
    else if (element == "search_database") // parent: "search_summary"
    {
      params_.db = attributeAsString_(attributes, "local_path");
      if (params_.db.empty())
      {
        optionalAttributeAsString_(params_.db, attributes, "database_name");
      }
    }
    else if (element == "msms_pipeline_analysis") // root
    {
      String date = attributeAsString_(attributes, "date");
      // fix for corrupted xs:dateTime format:
      if ((date[4] == ':') && (date[7] == ':') && (date[10] == ':'))
      {
        error(LOAD, "Format of attribute 'date' in tag 'msms_pipeline_analysis' does not comply with standard 'xs:dateTime'");
        date[4] = '-';
        date[7] = '-';
        date[10] = 'T';
      }
      date_ = asDateTime_(date);
    }
  }

  void PepXMLFile::endElement(const XMLCh* const /*uri*/,
                              const XMLCh* const /*local_name*/,
                              const XMLCh* const qname)
  {
    String element = sm_.convert(qname);

    // cout << "End: " << element << "\n";

    if (element == "analysis_summary")
    {
      analysis_summary_ = false;
    }
    else if (element == "search_score_summary")
    {
      search_score_summary_ = false;
    }
    else if (element == "analysis_result") // parent: "search_hit"
    {
      peptide_hit_.addAnalysisResults(current_analysis_result_);
    }
    else if (wrong_experiment_ || analysis_summary_)
    {
      // do nothing here (skip all elements that belong to the wrong experiment
      // or to an analysis summary)
    }
    else if (element == "spectrum_query")
    {
      // clear optional attributes
      native_spectrum_name_ = "";
      experiment_label_ = "";
      swath_assay_ = "";
      status_ = "";
    }
    else if (element == "search_hit")
    {
      AASequence temp_aa_sequence = AASequence::fromString(current_sequence_);

      // modification position is 1-based
      for (vector<pair<String, Size> >::const_iterator it = current_modifications_.begin(); it != current_modifications_.end(); ++it)
      {
        // e.g. Carboxymethyl (C)
        vector<String> mod_split;
        it->first.split(' ', mod_split);
        if (it->first.hasSubstring("C-term"))
        {
          temp_aa_sequence.setCTerminalModification(it->first);
        }
        else if (it->first.hasSubstring("N-term"))
        {
          temp_aa_sequence.setNTerminalModification(it->first);
        }
        else if (mod_split.size() == 2)
        {
          temp_aa_sequence.setModification(it->second - 1, mod_split[0]);
        }
        else
        {
          error(LOAD, String("Cannot parse modification '") + it->first + "@" + it->second + "'");
        }
      }

      // fixed modifications
      for (vector<AminoAcidModification>::const_iterator it = fixed_modifications_.begin(); it != fixed_modifications_.end(); ++it)
      {
        const Residue* residue = ResidueDB::getInstance()->getResidue(it->aminoacid);

        if (residue == nullptr)
        {
          double new_mass = it->massdiff.toDouble();
          if (it->aminoacid == "" && it->terminus =="n")
          {
            vector<String> mods;
            ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, new_mass, mod_tol_, "", ResidueModification::N_TERM);
            if (!mods.empty())
            {
              if (!temp_aa_sequence.hasNTerminalModification())
              {
                temp_aa_sequence.setNTerminalModification(mods[0]);
              }
              else
              {
                error(LOAD, String("Trying to add modification to modified terminal '") + it->aminoacid + "', delta mass: " +
                      it->massdiff + " mass: " + String(it->mass) + " variable: " + String(it->variable) + " terminus: " + it->terminus + " description: " + it->description);
              }
            }
            else
            {
              error(LOAD, String("Cannot find terminal modification '") + it->aminoacid + "', delta mass: " +
                    it->massdiff + " mass: " + String(it->mass) + " variable: " + String(it->variable) + " terminus: " + it->terminus + " description: " + it->description);
            }
          }
          else if (it->aminoacid == "" && it->terminus =="c")
          {
            vector<String> mods;
            ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, new_mass, mod_tol_, "", ResidueModification::C_TERM);
            if (!mods.empty())
            {
              if (!temp_aa_sequence.hasCTerminalModification())
              {
                temp_aa_sequence.setCTerminalModification(mods[0]);
              }
              else
              {
                error(LOAD, String("Trying to add modification to modified terminal '") + it->aminoacid + "', delta mass: " +
                      it->massdiff + " mass: " + String(it->mass) + " variable: " + String(it->variable) + " terminus: " + it->terminus + " description: " + it->description);
              }
            }
            else
            {
              error(LOAD, String("Cannot find terminal modification '") + it->aminoacid + "', delta mass: " +
                    it->massdiff + " mass: " + String(it->mass) + " variable: " + String(it->variable) + " terminus: " + it->terminus + " description: " + it->description);
            }
          }
          else
          {
            error(LOAD, String("Cannot parse modification of unknown amino acid '") + it->aminoacid + "', delta mass: " +
                  it->massdiff + " mass: " + String(it->mass) + " variable: " + String(it->variable) + " terminus: " + it->terminus + " description: " + it->description);
          }
        }
        else
        {
          double new_mass = it->mass - residue->getMonoWeight(Residue::Internal);
          vector<String> mods;
          ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, new_mass, mod_tol_, it->aminoacid, ResidueModification::ANYWHERE);
          if (mods.empty()) ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, new_mass, mod_tol_, it->aminoacid);
          if (!mods.empty())
          {
            for (Size i = 0; i < temp_aa_sequence.size(); ++i)
            {
              if (it->aminoacid.hasSubstring(temp_aa_sequence[i].getOneLetterCode()))
              {
                temp_aa_sequence.setModification(i, mods[0]);
              }
            }
          }
          else
          {
            error(LOAD, String("Cannot parse modification of amino acid '") + it->aminoacid + "'");
          }
        }
      }

      peptide_hit_.setSequence(temp_aa_sequence);
      current_peptide_.insertHit(peptide_hit_);
    }
    else if (element == "search_result")
    {
      peptides_->push_back(current_peptide_);
    }
    else if (element == "search_summary")
    {
      // In idXML we only store search engine and date as identifier, but to distinguish two identification runs these values must be unique.
      // As a workaround to support multiple runs, we make the date unique by adding one second for every additional identification run.
      UInt hour, minute, second;
      date_.getTime(hour, minute, second);
      hour = (hour + (minute + (second + 1) / 60) / 60) % 24;
      minute = (minute + (second + 1) / 60) % 60;
      second = (second + 1) % 60;
      date_.setTime(hour, minute, second);

      current_proteins_.back()->setSearchParameters(params_);
      search_summary_ = false;
    }
  }

} // namespace OpenMS

