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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlModificationsGenerator.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlReport.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlMarkerIonExtractor.h>

#include <QtCore/QStringList>
#include <QtCore/QProcess>
#include <QtCore/QObject>
#include <QFileInfo>

#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace OpenMS;

#define RT_FACTOR 10000000
#define RT_FACTOR_PRECISION 1000
#define RT_MODULO_FACTOR 10000 // last 4 digits is the index



/**
    @page UTILS_RNPxl RNPxl

    @brief Perform preotein-RNA cross-linking experiments.

    <CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ RNPxl \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> RNPxlXICFilter </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
        </tr>
    </table>
</CENTER>

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_RNPxl.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_RNPxl.html
*/

class TOPPRNPxl :
  public TOPPBase
{
public:
  TOPPRNPxl() :
    TOPPBase("RNPxl", "Tool for RNP cross linking experiment analysis.", false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    // input files
    registerInputFile_("in_mzML", "<file>", "", "Input file");
    setValidFormats_("in_mzML", ListUtils::create<String>("mzML"));

    registerIntOption_("length", "", 4, "Oligonucleotide maximum length.", false);

    registerStringOption_("sequence", "", "", "Sequence to restrict the generation of oligonucleotide chains. (disabled for empty sequence)", false);

    StringList target_nucleotides;
    target_nucleotides << "A=C10H14N5O7P" << "C=C9H14N3O8P" << "G=C10H14N5O8P" << "U=C9H13N2O9P";
    registerStringList_("target_nucleotides", "", target_nucleotides, "format:  target nucleotide=empirical formula of nucleoside monophosphate \n e.g. A=C10H14N5O7P, ..., U=C10H14N5O7P, X=C9H13N2O8PS  where X represents e.g. tU \n or e.g. Y=C10H14N5O7PS where Y represents tG", false, false);

    StringList mapping;
    mapping << "A->A" << "C->C" << "G->G" << "U->U";
    registerStringList_("mapping", "", mapping, "format: source->target e.g. A->A, ..., U->U, U->X", false, false);

    StringList restrictions;
    restrictions << "A=0" << "C=0" << "U=0" << "G=0";
    registerStringList_("restrictions", "", restrictions, "format: target nucleotide=min_count: e.g U=1 if at least one U must be in the generated sequence.", false, false);

    StringList modifications;
    modifications << "-H2O" << "" << "-H2O-HPO3" << "-HPO3" << "-H2O+HPO3" << "+HPO3";
    registerStringList_("modifications", "", modifications, "format: empirical formula e.g -H2O, ..., H2O+PO3", false, false);

    registerDoubleOption_("peptide_mass_threshold", "<threshold>", 600, "Lower peptide mass (Da) threshold.", false);
    registerDoubleOption_("precursor_variant_mz_threshold", "<threshold>", 260, "Lower m/z (Th) threshold for precursor variant.", false);

    registerFlag_("CysteineAdduct", "Use this flag if the +152 adduct is expected.");

    registerFlag_("continue", "Do not recreate intermediate files to continue after unexpected crash.", true);
    // search
    registerInputFile_("in_OMSSA_ini", "<file>", "", "Ini file for the OMSSA search engine\n");
    setValidFormats_("in_OMSSA_ini", ListUtils::create<String>("ini"));

    // indexing
    registerInputFile_("in_fasta", "<file>", "", "Fasta file for search result annotation\n");
    setValidFormats_("in_fasta", ListUtils::create<String>("fasta"));

    // reporting
    registerDoubleOption_("marker_ions_tolerance", "<tolerance>", 0.05, "mz tolerance used to determine marker ions.", false);
    registerOutputFile_("out_idXML", "<file>", "", "idXML output file\n");
    setValidFormats_("out_idXML", ListUtils::create<String>("idXML"));
    registerOutputFile_("out_csv", "<file>", "", "csv output file\n");
    setValidFormats_("out_csv", ListUtils::create<String>("csv"));
  }

  ExitCodes main_(int, const char**) override
  {
    const string in_mzml(getStringOption_("in_mzML"));

    // string format:  target,formula e.g. "A=C10H14N5O7P", ..., "U=C10H14N5O7P", "X=C9H13N2O8PS"  where X represents tU
    StringList target_nucleotides = getStringList_("target_nucleotides");

    // string format:  source->target e.g. "A->A", ..., "U->U", "U->X"
    StringList mappings = getStringList_("mapping");

    // string format: target,min_count: e.g "X=1" if at least one tU must be in the generated sequence.
    // All target nucleotides must be included. X=0 -> disable restriction
    StringList restrictions = getStringList_("restrictions");

    StringList modifications = getStringList_("modifications");

    String sequence_restriction = getStringOption_("sequence");

    Int max_length = getIntOption_("length");

    bool cysteine_adduct = getFlag_("CysteineAdduct");

    Size debug_level = (Size)getIntOption_("debug");

    double small_peptide_mass_filter_threshold = getDoubleOption_("peptide_mass_threshold");

    double precursor_variant_mz_threshold = getDoubleOption_("precursor_variant_mz_threshold");

    const String in_fasta_file(getStringOption_("in_fasta"));

    const String out_idXML = getStringOption_("out_idXML");
    const string out_csv = getStringOption_("out_csv");

    RNPxlModificationMassesResult mm = RNPxlModificationsGenerator::initModificationMassesRNA(target_nucleotides, mappings, restrictions, modifications, sequence_restriction, cysteine_adduct, max_length);

    String base_name = QFileInfo(QString::fromStdString(in_mzml)).baseName();

    vector<String> file_list_variants_mzML;

    PeakMap exp;
    MzMLFile().load(in_mzml, exp);

    String tmp_path = File::getTempDirectory();
    tmp_path.substitute('\\', '/');

    //REPORT:
    cout << "Theoretical precursor variants: " << mm.mod_masses.size() << endl;
    Size count_MS2 = 0;
    for (PeakMap::ConstIterator it = exp.begin(); it != exp.end(); ++it)
    {
      if (it->getMSLevel() == 2)
      {
        count_MS2++;
      }
    }
    cout << "Tandem spectra: " << count_MS2 << endl;

    Size fractional_mass_filtered = 0;
    Size small_peptide_weight_filtered = 0;
    Size precursor_variant_mz_filtered = 0;

    Size count(0);
    for (PeakMap::ConstIterator it = exp.begin(); it != exp.end(); ++it)
    {
      ++count;
      if (count % 100 == 0)
      {
        cout << (double) count / exp.size() * 100.0 << "%" << endl;
      }

      PeakMap new_exp;
      if (it->getMSLevel() != 2)
      {
        continue;
      }
      if (it->getPrecursors().size() == 0 || it->getPrecursors().begin()->getPosition()[0] == 0)
      {
        cerr << "Warning: no precursors found or no precursors with m/z > 0 found, skipping spectrum!" << endl;
        continue;
      }

      double prec_pos = it->getPrecursors().begin()->getPosition()[0];
      Int prec_charge = it->getPrecursors().begin()->getCharge();
      if (prec_charge == 0)
      {
        cerr << "Warning: precursor charge of spectrum RT=" << it->getRT() << " is zero, assuming double charged!" << endl;
        prec_charge = 2;
      }

      // add spec without any modifications
      double orig_rt = it->getRT();

      Size orig_rt_mul = (Size)(orig_rt * RT_FACTOR_PRECISION + 0.5) * RT_FACTOR / RT_FACTOR_PRECISION;

      Size mod_count(0);
      PeakSpectrum new_spec = *it;
      //cout << "add: " << setprecision(15) << (double)(orig_rt_mul + mod_count) << endl;
      new_spec.setRT(orig_rt_mul + mod_count++);
      Precursor new_prec;
      new_prec.setMZ(prec_pos);
      new_prec.setCharge(prec_charge);
      vector<Precursor> new_precs;
      new_precs.push_back(new_prec);
      new_spec.setPrecursors(new_precs);
      new_spec.setName("no_name");
      new_spec.setComment("no_comment");

      // Filter: peptide mass < 1750 and first deciamal place < 0.2
      double peptide_weight = prec_pos * prec_charge - prec_charge * Constants::PROTON_MASS_U;
      if (peptide_weight < 1750 && peptide_weight - floor(peptide_weight) < 0.2)
      {
        fractional_mass_filtered++;

        if (debug_level >= 1)
        {
          cout << orig_rt << "\t" << prec_pos << "\t" << "peptide weight < 1750 Da and first decimal place < 0.2" << endl;
        }
        continue;
      }

      // Filter: peptide mass < small_peptide_mass_filter_threshold (usually 600)
      if (peptide_weight < small_peptide_mass_filter_threshold)
      {
        small_peptide_weight_filtered++;
        if (debug_level >= 1)
        {
          cout << orig_rt << "\t" << prec_pos << "\t" << "peptide weight < 600 Da" << endl;
        }
        continue;
      }

      if (debug_level >= 1)
      {
        cout << orig_rt << "\t" << prec_pos << "\t" << "added with ";
      }
      new_exp.addSpectrum(new_spec);

      // add a new spec with each of the modifications
      int valid_mod_count = 0;
      for (Map<String, double>::ConstIterator mit = mm.mod_masses.begin(); mit != mm.mod_masses.end(); ++mit)
      {
        PeakSpectrum spec = *it;
        Precursor p;
        double prec_variant_mz = prec_pos - mit->second / (double)prec_charge;
        p.setMZ(prec_variant_mz);
        p.setCharge(prec_charge);
        spec.setName(String(prec_pos) + mit->first);
        spec.setComment(String(prec_pos) + mit->first);
        vector<Precursor> precs;
        precs.push_back(p);
        spec.setPrecursors(precs);
        spec.setRT(orig_rt_mul + mod_count++);

        //  Filter: check for minimum mass of the precursor variant => skip modification

        if (prec_variant_mz < precursor_variant_mz_threshold)
        {
          precursor_variant_mz_filtered++;
          if (debug_level > 2)
          {
            cout << spec.getRT() << "\t" << prec_variant_mz << "\t" << "m/z < " << precursor_variant_mz_threshold << endl;
          }
          continue;
        }
        valid_mod_count++;
        new_exp.addSpectrum(spec);
      }

      if (debug_level >= 1)
      {
        cout << valid_mod_count << " modifications." << endl;
      }
      String rt_string = String(it->getRT());
      String mz_string = String(prec_pos);

      if (!rt_string.has('.'))
      {
        rt_string += ".000";
      }

      if (!mz_string.has('.'))
      {
        mz_string += ".000";
      }

      String file_name_variant = tmp_path + "/" + base_name + "_"  + rt_string + "_" + mz_string + "_variant.mzML";
      file_list_variants_mzML.push_back(file_name_variant);

      if (getFlag_("continue") && File::exists(file_name_variant))
      {
        continue;
      }

      MzMLFile().store(file_name_variant, new_exp);
    }

    cout << base_name << ": " << "Spectra filtered by fractional mass: " << fractional_mass_filtered << endl;
    cout << base_name << ": " << "Spectra filtered by peptide weight: " << small_peptide_weight_filtered << endl;
    cout << base_name << ": " << "Precursor variants filtered by m/z: " << precursor_variant_mz_filtered << endl;

    Size sum_before = count_MS2 * mm.mod_masses.size();
    Size sum_after = (count_MS2 - fractional_mass_filtered - small_peptide_weight_filtered) * mm.mod_masses.size() - precursor_variant_mz_filtered;
    cout << base_name << ": " << "Before filtering: " << sum_before << " theoretical precursor variants." << endl;
    cout << base_name << ": " << "After filtering:  " << sum_after << " theoretical precursor variants." << endl;

    // perform OMSSA search
    {
      // get names of precursor variants mzml files
      const String in_OMSSA_ini(getStringOption_("in_OMSSA_ini"));
      for (vector<String>::const_iterator it = file_list_variants_mzML.begin(); it != file_list_variants_mzML.end(); ++it)
      {
        String in_string = *it;
        String out_string = in_string;
        out_string.substitute(".mzML", ".idXML");

        if (getFlag_("continue") && File::exists(out_string))
        {
          continue;
        }

        // Compose argument list and run OMSSA with new ini
        QStringList args;
        args << "-ini" << in_OMSSA_ini.toQString() << "-in" << in_string.toQString() << "-out" << out_string.toQString() << "-database" << in_fasta_file.toQString() << "-no_progress";

        // forward debug level to adapter
        if (getIntOption_("debug") != 0)
        {
          args << "-debug" << String(getIntOption_("debug")).toQString();
        }

        QProcess* p = new QProcess();
        p->setProcessChannelMode(QProcess::MergedChannels);
        p->start("OMSSAAdapter", args);
        p->waitForFinished(999999999);
        QString std_output = QString(p->readAllStandardOutput());
        if (getIntOption_("debug") != 0)
        {
          cout << std_output.toStdString() << endl;
        }
        delete(p);
      }
    }

    // create report
    vector<RNPxlReportRow> csv_rows;

    const double marker_tolerance = getDoubleOption_("marker_ions_tolerance");


    // protein and peptide identifications for all spectra
    vector<PeptideIdentification> whole_experiment_filtered_peptide_ids;
    vector<ProteinIdentification> whole_experiment_filtered_protein_ids;

    Size counter(0);
    for (vector<String>::const_iterator it = file_list_variants_mzML.begin(); it != file_list_variants_mzML.end(); ++it, ++counter)
    {
      if (*it == "")
      {
        continue;
      }

      String mzml_string = *it;
      String idxml_string = *it;
      idxml_string.substitute(".mzML", ".idXML");

      vector<ProteinIdentification> prot_ids;
      vector<PeptideIdentification> pep_ids;
      IdXMLFile().load(idxml_string, prot_ids, pep_ids);

      // copy protein identifications as is - they are not really needed in the later output
      whole_experiment_filtered_protein_ids.insert(whole_experiment_filtered_protein_ids.end(), prot_ids.begin(), prot_ids.end());

      // load map with all precursor variations (originating from one single precursor) that corresponds to this identification run
      PeakMap exp;
      MzMLFile().load(mzml_string, exp);

      // find marker ions
      RNPxlMarkerIonExtractor::MarkerIonsType marker_ions = RNPxlMarkerIonExtractor::extractMarkerIons(*exp.begin(), marker_tolerance);

      // case 1: no peptide identification
      RNPxlReportRow row;
      if (pep_ids.size() == 0)
      {
        row.no_id = true;
        row.rt = exp.begin()->getRT() / (double)RT_FACTOR;
        row.original_mz = exp.begin()->getPrecursors().begin()->getMZ();
        row.marker_ions = marker_ions;
        csv_rows.push_back(row);
        continue;
      }
      // end: case 1

      // case 2: peptide identifications

      // For every precursor variant spectrum of the map can produce a peptide identification with a single peptide hit TODO: check why there are sometimes more than 1
      // The best peptide hit of the whole variant map is retained and stored

      // copy all peptide hits (should be only one) from all peptide identifications (there can be one for every variant)
      vector<PeptideHit> pep_hits;
      for (vector<PeptideIdentification>::const_iterator pit = pep_ids.begin(); pit != pep_ids.end(); ++pit)
      {
        for (vector<PeptideHit>::const_iterator hit = pit->getHits().begin(); hit != pit->getHits().end(); ++hit)
        {
          pep_hits.push_back(*hit);
          pep_hits.back().setMetaValue("RT", pit->getRT());
          pep_hits.back().setMetaValue("MZ", pit->getMZ());
        }
      }

      // create new peptide identification and reassign all hits
      PeptideIdentification new_pep_id = *pep_ids.begin();
      new_pep_id.setHigherScoreBetter(false);
      new_pep_id.setHits(pep_hits);
      new_pep_id.assignRanks(); //sort by score and assign ranks
      pep_hits = new_pep_id.getHits();

      // only retain top hit if multiple peptide hits are present
      if (pep_hits.size() > 1)
      {
        pep_hits.resize(1);
      }
      new_pep_id.setHits(pep_hits); // assign top hit

      // store best peptide identification
      whole_experiment_filtered_peptide_ids.push_back(new_pep_id);

      for (vector<PeptideHit>::const_iterator hit = pep_hits.begin(); hit != pep_hits.end(); ++hit)
      {
        Size orig_rt = (double)hit->getMetaValue("RT");
        double orig_mz = (double)hit->getMetaValue("MZ");
        Size xlink_idx = orig_rt % RT_MODULO_FACTOR;
        String xlink_name = "";

        if (xlink_idx != 0) // non-modified
        {
          xlink_name = *(mm.mod_combinations[mm.mod_formula_idx[xlink_idx]].begin()); // take first one (if ambiguous)
        }

        double rt = (double)orig_rt / (double)RT_FACTOR;
        double pep_weight = hit->getSequence().getMonoWeight();
        double rna_weight = EmpiricalFormula(mm.mod_formula_idx[xlink_idx]).getMonoWeight();

        // xlink weight for different charge states
        double weight_z1 = (pep_weight + rna_weight + 1.0 * Constants::PROTON_MASS_U);
        double weight_z2 = (pep_weight + rna_weight + 2.0 * Constants::PROTON_MASS_U) / 2.0;
        double weight_z3 = (pep_weight + rna_weight + 3.0 * Constants::PROTON_MASS_U) / 3.0;
        double weight_z4 = (pep_weight + rna_weight + 4.0 * Constants::PROTON_MASS_U) / 4.0;

        Size charge = hit->getCharge();
        double ppm_difference(0), absolute_difference(0);
        double exp_mz = orig_mz + rna_weight / (double)charge;
        double theo_mz = (pep_weight + rna_weight + (double)charge * Constants::PROTON_MASS_U) / (double)charge;
        absolute_difference = theo_mz - exp_mz;
        ppm_difference = absolute_difference / theo_mz * 1000000;

        String protein_accessions;
        set<String> accs = hit->extractProteinAccessionsSet();

        for (set<String>::const_iterator a_it = accs.begin(); a_it != accs.end(); ++a_it)
        {
          if (a_it != accs.begin())
          {
            protein_accessions += ",";
          }
          protein_accessions += *a_it;
        }

        row.no_id = false;
        row.rt = rt;
        row.original_mz = exp.begin()->getPrecursors().begin()->getMZ();
        row.accessions = protein_accessions;
        row.RNA = xlink_name;
        row.peptide = hit->getSequence().toString();
        row.charge = hit->getCharge();
        row.score = hit->getScore();
        row.peptide_weight = pep_weight;
        row.RNA_weight = rna_weight;
        row.xl_weight = pep_weight + rna_weight;

        whole_experiment_filtered_peptide_ids.back().setMZ(exp_mz);
        whole_experiment_filtered_peptide_ids.back().setMetaValue("cross link id", DataValue(xlink_idx));
        whole_experiment_filtered_peptide_ids.back().setMetaValue("RNA", DataValue(xlink_name));
        whole_experiment_filtered_peptide_ids.back().setMetaValue("peptide mass", DataValue(pep_weight));
        whole_experiment_filtered_peptide_ids.back().setMetaValue("RNA mass", DataValue(rna_weight));
        whole_experiment_filtered_peptide_ids.back().setMetaValue("cross link mass", (double)pep_weight + rna_weight);

        for (Map<String, vector<pair<double, double> > >::const_iterator it = marker_ions.begin(); it != marker_ions.end(); ++it)
        {
          for (Size i = 0; i != it->second.size(); ++i)
          {
            whole_experiment_filtered_peptide_ids.back().setMetaValue(it->first + "_" + it->second[i].first, (double)it->second[i].second * 100.0);
          }
        }

        row.marker_ions = marker_ions;
        row.abs_prec_error = absolute_difference;
        row.rel_prec_error = ppm_difference;
        row.m_H = weight_z1;
        row.m_2H = weight_z2;
        row.m_3H = weight_z3;
        row.m_4H = weight_z4;

        csv_rows.push_back(row);

        whole_experiment_filtered_peptide_ids.back().setMetaValue("Da difference", (double)absolute_difference);
        whole_experiment_filtered_peptide_ids.back().setMetaValue("ppm difference", (double)ppm_difference);
        whole_experiment_filtered_peptide_ids.back().setMetaValue("z1 mass", (double)weight_z1);
        whole_experiment_filtered_peptide_ids.back().setMetaValue("z2 mass", (double)weight_z2);
        whole_experiment_filtered_peptide_ids.back().setMetaValue("z3 mass", (double)weight_z3);
        whole_experiment_filtered_peptide_ids.back().setMetaValue("z4 mass", (double)weight_z4);
      }
    }

    vector<ProteinIdentification> pr_tmp;
    pr_tmp.push_back(ProteinIdentification());
    for (size_t k = 0; k != whole_experiment_filtered_protein_ids.size(); ++k)
    {
      vector<ProteinHit> ph_tmp = whole_experiment_filtered_protein_ids[k].getHits();
      for (vector<ProteinHit>::iterator it = ph_tmp.begin(); it != ph_tmp.end(); ++it)
      {
        pr_tmp[0].insertHit(*it);
        //cout << it->getAccession() << endl;
      }
    }

    // create new peptide identifications and copy over data
    vector<PeptideIdentification> pt_tmp;
    for (size_t k = 0; k != whole_experiment_filtered_peptide_ids.size(); ++k)
    {
      for (vector<PeptideHit>::const_iterator hit = whole_experiment_filtered_peptide_ids[k].getHits().begin(); hit != whole_experiment_filtered_peptide_ids[k].getHits().end(); ++hit)
      {
        PeptideIdentification np;
        double rt = whole_experiment_filtered_peptide_ids[k].getRT();
        double orig_rt = rt / (double)RT_FACTOR;
        np.setRT(orig_rt);
        np.setMZ(whole_experiment_filtered_peptide_ids[k].getMZ());

        vector<PeptideHit> phs;
        PeptideHit ph = *hit;
        std::vector<String> keys;
        whole_experiment_filtered_peptide_ids[k].getKeys(keys);
        for (size_t i = 0; i != keys.size(); ++i)
        {
          DataValue dv = whole_experiment_filtered_peptide_ids[k].getMetaValue(keys[i]);
          if (dv.valueType() == DataValue::DOUBLE_VALUE)
          {
            ph.setMetaValue(keys[i], (double)dv);
          }
          else
          {
            ph.setMetaValue(keys[i], dv);
          }
        }
        phs.push_back(ph);
        np.setHits(phs);
        np.assignRanks(); //sort by score and assign ranks
        pt_tmp.push_back(np);
      }
    }

    if (!pr_tmp.empty())
    {
      ProteinIdentification &p_tmp = pr_tmp[0];
      StringList ms_runs;
      exp.getPrimaryMSRunPath(ms_runs);
      p_tmp.setPrimaryMSRunPath(ms_runs);
    }
    IdXMLFile().store(out_idXML, pr_tmp, pt_tmp, "summary");


    QProcess* p;

    // index final result
    QStringList args;
    args << "-missing_decoy_action" << "warn" << "-fasta" << in_fasta_file.toQString() << "-in" << out_idXML.toQString() << "-out" << out_idXML.toQString() << "-no_progress";
    p = new QProcess();
    p->setProcessChannelMode(QProcess::MergedChannels);
    p->start("PeptideIndexer", args);
    p->waitForFinished(-1);
    QString peptide_indexer_stdout = QString(p->readAllStandardOutput());
    if (getIntOption_("debug") > 0)
    {
      cout << peptide_indexer_stdout.toStdString() << endl;
    }
    delete(p);

    // load indexed idXML
    pr_tmp.clear();
    pt_tmp.clear();
    IdXMLFile().load(out_idXML, pr_tmp, pt_tmp);

    // reindex tabular data to contain protein ids
    map<double, String> map_rt_2_accession;
    for (vector<PeptideIdentification>::const_iterator pit = pt_tmp.begin(); pit != pt_tmp.end(); ++pit)
    {
      for (vector<PeptideHit>::const_iterator hit = pit->getHits().begin(); hit != pit->getHits().end(); ++hit)
      {
        double rt = pit->getRT();

        set<String> accessions = hit->extractProteinAccessionsSet();

        String accession_string;
        Size j = 0;
        for (set<String>::const_iterator a_it = accessions.begin(); a_it != accessions.end(); ++a_it, ++j)
        {
          if (j < 3)
          {
            accession_string += *a_it + " ";
          }
          else
          {
            accession_string += "...";
            break;
          }
        }
        map_rt_2_accession[rt] = accession_string;
      }
    }

    for (vector<RNPxlReportRow>::iterator rit = csv_rows.begin(); rit != csv_rows.end(); ++rit)
    {
      double current_rt = rit->rt;
      map<double, String>::iterator before = map_rt_2_accession.lower_bound(current_rt);
      map<double, String>::iterator min_distance_it;

      if (before == map_rt_2_accession.begin())
      {
        min_distance_it = before;
      }
      else if (before == map_rt_2_accession.end())
      {
        min_distance_it = --before;
      }
      else
      {
        map<double, String>::iterator after = before;
        --before;
        if ((after->first - current_rt) < (current_rt - before->first))
        {
          min_distance_it = after;
        }
        else
        {
          min_distance_it = before;
        }
      }

      cout << min_distance_it->first << " " << current_rt << endl;
      if (fabs(min_distance_it->first - current_rt) < 0.01)
      {
        rit->accessions = min_distance_it->second;
      }
    }

    // write csv
    ofstream csv_file(out_csv.c_str());

    // header of table
    csv_file << RNPxlReportRowHeader().getString("\t") << endl;

    for (Size i = 0; i != csv_rows.size(); ++i)
    {
      csv_file << csv_rows[i].getString("\t") << endl;
    }

    csv_file.close();

    // cleanup
    if (debug_level < 1)
    {
      Size mzml_removed = 0;
      Size idxml_removed = 0;
      for (vector<String>::const_iterator it = file_list_variants_mzML.begin(); it != file_list_variants_mzML.end(); ++it, ++counter)
      {
        if (*it == "")
        {
          continue;
        }
        String mzml_string = *it;
        String idxml_string = *it;
        idxml_string.substitute(".mzML", ".idXML");
        if (QFile(mzml_string.toQString()).remove())
        {
          mzml_removed++;
        }

        if (QFile(idxml_string.toQString()).remove())
        {
          idxml_removed++;
        }
      }
      cout << "Cleaning up. Removed " << mzml_removed << " temporary mzML files and " << idxml_removed << " temporary idXML files." << endl;
    }
    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPRNPxl tool;
  return tool.main(argc, argv);
}
