// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/MascotXMLHandler.h> // for "primary_scan_regex"
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/SYSTEM/File.h>
#include <fstream>
#include <iostream>
#include <boost/regex.hpp>

using namespace std;

namespace OpenMS
{

  PepXMLFile::PepXMLFile() :
    XMLHandler("", "1.12"),
    XMLFile("/SCHEMAS/pepXML_v114.xsd", "1.14"),
    proteins_(0),
    peptides_(0),
    experiment_(0),
    scan_map_(),
    rt_tol_(10.0),
    mz_tol_(10.0)
  {
    const ElementDB* db = ElementDB::getInstance();
    hydrogen_ = *db->getElement("Hydrogen");
  }

  PepXMLFile::~PepXMLFile()
  {
  }

  void PepXMLFile::store(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids, const String& mz_file, const String& mz_name, bool peptideprophet_analyzed)
  {
    ofstream f(filename.c_str());
    if (!f)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
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
      }
    }

    f.precision(writtenDigits<double>(0.0));
    String raw_data;
    String base_name;
    // The mz-File (if given)
    if (!mz_file.empty())
    {
      base_name = File::basename(mz_file);
      raw_data = FileTypes::typeToName(FileHandler().getTypeByFileName(mz_file));
    }
    else
    {
      base_name = File::basename(filename);
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
    // If enzyme is not trypsin, skip it here and specify in TPP parser.
    if (search_params.enzyme == ProteinIdentification::TRYPSIN || search_params.digestion_enzyme.getName() == "Trypsin")
    {
      f << "\t<sample_enzyme name=\"" << "trypsin" << "\">" << "\n";
      f << "\t\t<specificity cut=\"KR\" no_cut=\"P\" sense=\"C\"/>" << "\n";
      f << "\t</sample_enzyme>" << "\n";
    }
    else if (search_params.digestion_enzyme.getRegEx() != "")
    {
      f << "\t<sample_enzyme name=\"" << search_params.digestion_enzyme.getName() << "\">" << "\n";
      f << "\t\t<specificity cut=\"";
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
      f << "\" sense=\"C\"/>" << "\n";
      f << "\t</sample_enzyme>" << "\n";
    }

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
            n_term_mods.insert(ModificationsDB::getInstance()->getTerminalModification(p.getNTerminalModification(), ResidueModification::N_TERM).getFullId());
          }
          if (p.hasCTerminalModification())
          {
            c_term_mods.insert(ModificationsDB::getInstance()->getTerminalModification(p.getCTerminalModification(), ResidueModification::C_TERM).getFullId());
          }

          for (Size i = 0; i != p.size(); ++i)
          {
            if (p[i].isModified())
            {
              aa_mods.insert(ModificationsDB::getInstance()->getModification(p[i].getOneLetterCode(), p[i].getModification(), ResidueModification::ANYWHERE).getFullId());
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
      const ResidueModification& mod = ModificationsDB::getInstance()->getModification(*it);

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
      const ResidueModification& mod = ModificationsDB::getInstance()->
                                       getModification(*it);
      f << "\t\t"
        << "<terminal_modification terminus=\"n\" massdiff=\""
        << precisionWrapper(mod.getDiffMonoMass()) << "\" mass=\"" << precisionWrapper(mod.getMonoMass())
        << "\" variable=\"Y\" description=\"" << *it
        << "\" protein_terminus=\"\"/>" << "\n";
    }

    for (set<String>::const_iterator it = c_term_mods.begin(); it != c_term_mods.end(); ++it)
    {
      const ResidueModification& mod = ModificationsDB::getInstance()->getModification(*it);
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

    Int count(1);
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

        Int scan_index;
        if (it->metaValueExists("RT_index")) // Setting metaValue "RT_index" in XTandemXMLFile in the case of X! Tandem.
        {
          scan_index = it->getMetaValue("RT_index");
        }
        else
        {
          scan_index = count;
        }
        // PeptideProphet requires this format for "spectrum" attribute (otherwise TPP parsing error)
        f << "\t<spectrum_query spectrum=\"" << base_name << ".00000.00000." << h.getCharge() << "\""
          << " start_scan=\"" << scan_index << "\""
          << " end_scan=\"" << scan_index << "\""
          << " precursor_neutral_mass=\"" << precisionWrapper(precursor_neutral_mass) << "\""
          << " assumed_charge=\"" << h.getCharge() << "\" index=\"" << count << "\"";

        if (it->hasRT())
        {
          f << " retention_time_sec=\"" << it->getRT() << "\" ";
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
        if ((pe.getAABefore() == 'R' || pe.getAABefore() == 'K') && search_params.enzyme == ProteinIdentification::TRYPSIN)
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
            << seq << "\"";

          if (seq.hasNTerminalModification())
          {
            const ResidueModification& mod = ModificationsDB::getInstance()->getTerminalModification(seq.getNTerminalModification(), ResidueModification::N_TERM);
            f << " mod_nterm_mass=\"" <<
              precisionWrapper(mod.getMonoMass() + seq[(Size)0].getMonoWeight(Residue::Internal)) << "\"";
          }

          if (seq.hasCTerminalModification())
          {
            const ResidueModification& mod = ModificationsDB::getInstance()->getTerminalModification(seq.getCTerminalModification(), ResidueModification::C_TERM);
            f << "mod_cterm_mass=\"" <<
              precisionWrapper(mod.getMonoMass() + seq[seq.size() - 1].getMonoWeight(Residue::Internal)) << "\"";
          }

          f << ">" << "\n";

          for (Size i = 0; i != seq.size(); ++i)
          {
            if (seq.isModified(i))
            {
              const ResidueModification& mod = ModificationsDB::getInstance()->getModification(seq[i].getOneLetterCode(), seq[i].getModification(), ResidueModification::ANYWHERE);
              // the modification position is 1-based
              f << "\t\t\t\t<mod_aminoacid_mass position=\"" << (i + 1)
                << "\" mass=\"" <<
                precisionWrapper(mod.getMonoMass() + seq[i].getMonoWeight(Residue::Internal)) << "\"/>" << "\n";
            }
          }

          f << "\t\t\t</modification_info>" << "\n";
        }
        if (peptideprophet_analyzed)
        {
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
              if (it->metaValueExists("nextscore"))
              {
                f << h.getMetaValue("nextscore") << "\"" << "/>\n";
              }
              else
              {
                f << h.getScore() << "\"" << "/>\n";
              }
            }
            else if (it->metaValueExists("XTandem_score"))
            {
              f << "\t\t\t<search_score" << " name=\"hyperscore\" value=\"" << h.getMetaValue("XTandem_score") << "\"" << "/>\n";
              f << "\t\t\t<search_score" << " name=\"nextscore\" value=\"";
              if (it->metaValueExists("nextscore"))
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
            f << "\t\t\t<search_score" << " name=\"expect\" value=\"" << h.getMetaValue("E-Value") << "\"" << "/>\n";
          }
          else if (search_engine_name == "OMSSA")
          {
            f << "\t\t\t<search_score" << " name=\"expect\" value=\"" << h.getScore() << "\"" << "/>\n";
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
    ModificationsDB::getInstance()->getModificationsByDiffMonoMass(mods, origin, mod_mass, 0.001);

    // no notification about ambiguities here - that was done when the
    // modification definitions were parsed ("aminoacid_modification" and
    // "terminal_modification" elements)
    if (!mods.empty()) modification_description = mods[0];
  }

  void PepXMLFile::makeScanMap_()
  {
    scan_map_.clear();
    Size scan = 0;
    for (MSExperiment<>::ConstIterator e_it = experiment_->begin(); e_it != experiment_->end(); ++e_it, ++scan)
    {
      String id = e_it->getNativeID();
      bool failed = false;
      try
      {
        // expected format: "spectrum=#" (mzData) or "scan=#" (mzXML)
        Int num_id = id.suffix('=').toInt();
        if (num_id >= 0)
        {
          scan_map_.insert(scan_map_.end(), pair<Size, Size>(num_id, scan));
        }
        else
        {
          failed = true;
        }
      }
      catch (Exception::ConversionError)
      {
        failed = true;
      }
      if (failed)
      {
        scan_map_.clear();
        error(LOAD, "Could not construct mapping of native scan numbers to indexes");
      }
    }
  }

  void PepXMLFile::readRTMZCharge_(const xercesc::Attributes& attributes)
  {
    double mass = attributeAsDouble_(attributes, "precursor_neutral_mass");
    charge_ = attributeAsInt_(attributes, "assumed_charge");
    mz_ = (mass + hydrogen_mass_ * charge_) / charge_;
    rt_ = 0;

    bool rt_present = optionalAttributeAsDouble_(rt_, attributes, "retention_time_sec");

    if (!rt_present || use_precursor_data_) // get RT from experiment
    {
      if (!experiment_)
      {
        error(LOAD, "Cannot get precursor information - no experiment given");
        return;
      }

      // assume only one scan, i.e. ignore "end_scan":
      Size scan = attributeAsInt_(attributes, "start_scan");
      if (scan == 0) // work-around for pepXMLs exported from Mascot
      {
        String spectrum = attributeAsString_(attributes, "spectrum");
        boost::regex re(Internal::MascotXMLHandler::primary_scan_regex,
                        boost::regex::perl | boost::regex::icase);
        boost::smatch match;
        if (boost::regex_search(spectrum, match, re))
        {
          scan = String(match["SCAN"].str()).toInt();
        }
      }

      if (!scan_map_.empty()) scan = scan_map_[scan];
      const MSSpectrum<>& spec = (*experiment_)[scan];
      bool success = false;
      if (spec.getMSLevel() == 2)
      {
        if (!use_precursor_data_)
        {
          rt_ = spec.getRT();
          success = true;
        }
        else if (!rt_present || Math::approximatelyEqual(spec.getRT(), rt_, 0.001))
        {
          double prec_mz = 0, prec_rt = 0;
          vector<Precursor> precursors = spec.getPrecursors();
          if (!precursors.empty())
          {
            prec_mz = precursors[0].getMZ(); // assume only one precursor
          }
          MSExperiment<>::ConstIterator it = experiment_->getPrecursorSpectrum(experiment_->begin() + scan);
          if (it != experiment_->end())
          {
            prec_rt = it->getRT();
          }

          // check if "rt"/"mz" are similar to "prec_rt"/"prec_mz"
          // (otherwise, precursor mapping is wrong)
          if ((prec_mz > 0) && Math::approximatelyEqual(prec_mz, mz_, mz_tol_)    && (prec_rt > 0) && (!rt_present || Math::approximatelyEqual(prec_rt, rt_, rt_tol_)))
          {
            // double diff;
            // diff = mz_ - prec_mz;
            // cout << "m/z difference: " << diff << " ("
            //       << diff / max(mz_, prec_mz) * 100 << "%)\n";
            // diff = rt_ - prec_rt;
            // cout << "RT difference: " << diff << " ("
            //       << diff / max(rt_, prec_rt) * 100 << "%)\n" << "\n";
            mz_ = prec_mz;
            rt_ = prec_rt;
            success = true;
          }
        }
      }
      if (!success)
      {
        error(LOAD, "Cannot get precursor information - scan mapping is incorrect");
      }
    }
  }

  void PepXMLFile::load(const String& filename, vector<ProteinIdentification>&
                        proteins, vector<PeptideIdentification>& peptides,
                        const String& experiment_name)
  {
    MSExperiment<> exp;
    load(filename, proteins, peptides, experiment_name, exp, false);
  }

  void PepXMLFile::load(const String& filename, vector<ProteinIdentification>&
                        proteins, vector<PeptideIdentification>& peptides,
                        const String& experiment_name, const MSExperiment<>&
                        experiment, bool use_precursor_data)
  {
    // initialize here, since "load" could be called several times:
    exp_name_ = "";
    experiment_ = 0;
    use_precursor_data_ = use_precursor_data;
    prot_id_ = "";
    charge_ = 0;
    rt_tol_ = 10.0;
    mz_tol_ = 10.0;
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

      if (!experiment.empty()) // use experiment only if we know the name
      {
        experiment_ = &experiment;
        MSExperiment<>::AreaType area = experiment_->getDataRange();
        // set tolerance to 1% of data range (if above a sensible minimum):
        rt_tol_ = max((area.maxX() - area.minX()) * 0.01, rt_tol_);
        mz_tol_ = max((area.maxY() - area.minY()) * 0.01, mz_tol_);
        makeScanMap_();
      }
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
    proteins_ = 0;
    peptides_ = 0;
    experiment_ = 0;
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
      enzyme_ = ProteinIdentification::UNKNOWN_ENZYME;
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
      else if (name == "xcorr" && search_engine_ != "MyriMatch") // Sequest score; MyriMatch has also an xcorr, but we want to ignore it
      { // and use the mvh
        value = attributeAsDouble_(attributes, "value");
        peptide_hit_.setScore(value);
        current_peptide_.setScoreType(name); // add "Sequest" to name?
        current_peptide_.setHigherScoreBetter(true);
      }
      else if (name == "fval") // SpectraST score
      {
        value = attributeAsDouble_(attributes, "value");
        peptide_hit_.setScore(value);
        current_peptide_.setScoreType(name);
        current_peptide_.setHigherScoreBetter(true);
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
    { // creates a new PeptideIdentification
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
    }
    else if (element == "spectrum_query") // parent: "msms_run_summary"
    {
      readRTMZCharge_(attributes); // sets "rt_", "mz_", "charge_"
    }
    else if (element == "peptideprophet_result") // parent: "analysis_result" (in "search_hit")
    {
      // PeptideProphet probability overwrites original search score
      // maybe TODO: deal with meta data associated with PeptideProphet search
      if (current_peptide_.getScoreType() != "InterProphet probability")
      {
        double value = attributeAsDouble_(attributes, "probability");
        peptide_hit_.setScore(value);
        current_peptide_.setScoreType("PeptideProphet probability");
        current_peptide_.setHigherScoreBetter(true);
      }
    }
    else if (element == "interprophet_result") // parent: "analysis_result" (in "search_hit")
    {
      // InterProphet probability overwrites PeptideProphet probability and
      // original search score
      double value = attributeAsDouble_(attributes, "probability");
      peptide_hit_.setScore(value);
      current_peptide_.setScoreType("InterProphet probability");
      current_peptide_.setHigherScoreBetter(true);
    }
    else if (element == "modification_info") // parent: "search_hit" (in "search result")
    {
      // Has N-Term Modification
      double mod_nterm_mass;
      if (optionalAttributeAsDouble_(mod_nterm_mass, attributes, "mod_nterm_mass")) // this specifies a terminal modification
      {
        // lookup the modification in the search_summary by mass
        for (vector<AminoAcidModification>::const_iterator it = variable_modifications_.begin(); it != variable_modifications_.end(); ++it)
        {
          if (mod_nterm_mass == it->mass && it->terminus == "n")
          {
            double massdiff = (it->massdiff).toDouble();
            vector<String> mods;
            ModificationsDB::getInstance()->getTerminalModificationsByDiffMonoMass(mods, massdiff, 0.001, ResidueModification::N_TERM);
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
        // lookup the modification in the search_summary by mass
        for (vector<AminoAcidModification>::const_iterator it = variable_modifications_.begin(); it != variable_modifications_.end(); ++it)
        {
          if (mod_cterm_mass == it->mass && it->terminus == "c")
          {
            double massdiff = (it->massdiff).toDouble();
            vector<String> mods;
            ModificationsDB::getInstance()->getTerminalModificationsByDiffMonoMass(mods, massdiff, 0.001, ResidueModification::C_TERM);
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
    else if (element == "aminoacid_modification") // parent: "search_summary"
    {
      AminoAcidModification aa_mod;
      optionalAttributeAsString_(aa_mod.description, attributes, "description");
      aa_mod.massdiff = attributeAsString_(attributes, "massdiff");
      aa_mod.aminoacid = attributeAsString_(attributes, "aminoacid");
      aa_mod.mass = attributeAsDouble_(attributes, "mass");
      String is_variable = attributeAsString_(attributes, "variable");
      if (is_variable == "Y")
      {
        if (aa_mod.description != "")
        {
          params_.variable_modifications.push_back(aa_mod.description); // TODO
        }
        else
        {
          String desc = aa_mod.aminoacid;
          if (aa_mod.massdiff.toDouble() >= 0)
          {
            desc += "+" + String(aa_mod.massdiff);
          }
          else
          {
            desc += String(aa_mod.massdiff);
          }
          params_.variable_modifications.push_back(desc);
        }
      }
      else
      {
        fixed_modifications_.push_back(aa_mod);
        if (aa_mod.description != "")
        {
          params_.fixed_modifications.push_back(aa_mod.description); // TODO
        }
        else
        {
          String desc = aa_mod.aminoacid;
          if (aa_mod.massdiff.toDouble() >= 0)
          {
            desc += "+" + String(aa_mod.massdiff);
          }
          else
          {
            desc += String(aa_mod.massdiff);
          }
          params_.fixed_modifications.push_back(desc);
        }
      }
      // check if the modification is uniquely defined:
      vector<String> mods;
      ModificationsDB::getInstance()->getModificationsByDiffMonoMass(
        mods, aa_mod.aminoacid, aa_mod.massdiff.toDouble(), 0.001);
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
    else if (element == "terminal_modification") // parent: "search_summary"
    {
      // <terminal_modification terminus="n" massdiff="+108.05" mass="109.06" variable="N" protein_terminus="" description="dNIC (N-term)"/>
      AminoAcidModification aa_mod;
      optionalAttributeAsString_(aa_mod.description, attributes, "description");
      aa_mod.massdiff = attributeAsString_(attributes, "massdiff");

      // somehow very small fixed modification (electron mass?) gets annotated by X!Tandem. Don't add them as they interfere with other modifications.
      if (fabs(aa_mod.massdiff.toDouble()) < 0.0005)
      {
        return;
      }
      optionalAttributeAsString_(aa_mod.aminoacid, attributes, "aminoacid");
      aa_mod.mass = attributeAsDouble_(attributes, "mass");
      aa_mod.terminus = attributeAsString_(attributes, "terminus");
      String is_variable = attributeAsString_(attributes, "variable");

      if (is_variable == "Y")
      {
        variable_modifications_.push_back(aa_mod);
        if (aa_mod.description != "")
        {
          params_.variable_modifications.push_back(aa_mod.description); // TODO
        }
        else
        {
          String desc = aa_mod.aminoacid;
          if (aa_mod.massdiff.toDouble() > 0)
          {
            desc += "+" + String(aa_mod.massdiff);
          }
          else
          {
            desc += String(aa_mod.massdiff);
          }
          params_.variable_modifications.push_back(desc);
        }
      }
      else
      {
        fixed_modifications_.push_back(aa_mod);
        if (aa_mod.description != "")
        {
          params_.fixed_modifications.push_back(aa_mod.description); // TODO
        }
        else
        {
          String desc = aa_mod.aminoacid;
          if (aa_mod.massdiff.toDouble() > 0)
          {
            desc += "+" + String(aa_mod.massdiff);
          }
          else
          {
            desc += String(aa_mod.massdiff);
          }
          params_.fixed_modifications.push_back(desc);
        }
      }
      // check if the modification is uniquely defined:
      vector<String> mods;
      ModificationsDB::getInstance()->getModificationsByDiffMonoMass(
        mods, aa_mod.aminoacid, aa_mod.massdiff.toDouble(), 0.001);
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
    else if (element == "search_summary") // parent: "msms_run_summary"
    { // creates a new ProteinIdentification (usually)
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
      params_.enzyme = enzyme_;
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
      String name = attributeAsString_(attributes, "name");
      name.toLower();
      // spelling of enzyme names in pepXML?
      if (name.hasPrefix("trypsin"))
        enzyme_ = ProteinIdentification::TRYPSIN;
      else if (name.hasPrefix("pepsin"))
        enzyme_ = ProteinIdentification::PEPSIN_A;
      else if (name.hasPrefix("protease"))
        enzyme_ = ProteinIdentification::PROTEASE_K;
      else if (name.hasPrefix("chymotrypsin"))
        enzyme_ = ProteinIdentification::CHYMOTRYPSIN;
      else
        enzyme_ = ProteinIdentification::UNKNOWN_ENZYME;

      ProteinIdentification::SearchParameters params =
        current_proteins_.front()->getSearchParameters();
      params.enzyme = enzyme_;
      current_proteins_.front()->setSearchParameters(params);
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
    else if (wrong_experiment_ || analysis_summary_)
    {
      // do nothing here (skip all elements that belong to the wrong experiment
      // or to an analysis summary)
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

        if (residue == 0)
        {
          double new_mass = it->massdiff.toDouble();
          if (it->aminoacid == "" && it->terminus =="n")
          {
            vector<String> mods;
            ModificationsDB::getInstance()->getTerminalModificationsByDiffMonoMass(mods, new_mass, 0.001, ResidueModification::N_TERM);
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
            ModificationsDB::getInstance()->getTerminalModificationsByDiffMonoMass(mods, new_mass, 0.001, ResidueModification::C_TERM);
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
          ModificationsDB::getInstance()->getModificationsByDiffMonoMass(mods, it->aminoacid, new_mass, 0.001);
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
    }
  }

} // namespace OpenMS
