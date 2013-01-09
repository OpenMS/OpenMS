// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Mathias Walzer $
// $Author: Mathias Walzer, Sven Nahnsen $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/QcMLFile.h>

#include <QFileInfo>

#include <vector>
#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_QCCalculator QCCalculator

    @brief This application is used to provide data export from raw, id and feature data files generated via TOPP pipelines. It is intended to provide tables that can be read into R where QC metrics will be calculated.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_QCCalculator.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_QCCalculator.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPQCCalculator :
  public TOPPBase
{
public:
  TOPPQCCalculator() :
    TOPPBase("QCCalculator", "produces table data dedicted for R import. These data is produced based on mzML, featureXMl and/ or idXML files", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "raw data input file (this is relevant if you want to look at MS1, MS2 and precursor peak information)");
    setValidFormats_("in", StringList::create("mzML"));
    registerOutputFile_("out", "<file>", "", "Your qcML file.");
    setValidFormats_("out", StringList::create("qcML"));
    registerInputFile_("id", "<file>", "", "Input idXML file containing the identifications. Your identifications will be exported in an easy-to-read format", false);
    setValidFormats_("id", StringList::create("idXML"));
    registerInputFile_("feature", "<file>", "", "feature input file (this is relevant for most QC issues)", false);
    setValidFormats_("feature", StringList::create("featureXML"));
    registerInputFile_("consensus", "<file>", "", "consensus input file (this is only used for charge state deconvoluted output. Use the consensusXML output form the DeCharger)", false);
    setValidFormats_("consensus", StringList::create("consensusXML"));
    registerFlag_("AllHits", "If set, the exporter takes all peptide candidates per spectrum into account.");
    registerFlag_("remove_duplicate_features", "This flag should be set, if you work with a set of merged features.");
    registerFlag_("MS1", "This flag should be set, if you want to work with MS1 stats.");
    registerFlag_("MS2", "This flag should be set, if you want to work with MS2 stats.");
  }

  ExitCodes main_(int, const char **)
  {
    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;
    ProteinHit temp_protein_hit;

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String inputfile_id         = getStringOption_("id");
    String inputfile_feature    = getStringOption_("feature");
    String inputfile_consensus  = getStringOption_("consensus");
    String inputfile_raw        = getStringOption_("in");
    String outputfile_name      = getStringOption_("out");

    bool AllHits(getFlag_("AllHits"));
    bool Ms1(getFlag_("MS1"));
    bool Ms2(getFlag_("MS2"));
    bool remove_duplicate_features(getFlag_("remove_duplicate_features"));
    //-------------------------------------------------------------
    // reading input
    //------------------------------------------------------------

    QcMLFile qcmlfile;
    String base_name = QFileInfo(QString::fromStdString(inputfile_raw)).baseName();

    cout << "Reading mzML file..." << endl;
    MzMLFile mz_data_file;
    MSExperiment<Peak1D> exp;
    MzMLFile().load(inputfile_raw, exp);

    //---file origin qp
    QcMLFile::QualityParameter qp;
    qp.name = "raw data file"; ///< Name
    qp.id = base_name + "_origin_rawfile_name"; ///< Identifier
    qp.cvRef = "PSI"; ///< cv reference
    qp.cvAcc = "MS:1000577"; ///< cv accession
    qp.value = base_name;
    qcmlfile.addQualityParameter(base_name, qp);

    //---mzml stats qp
    qp = QcMLFile::QualityParameter() ;
    qp.name = "mzmlstats"; ///< Name
    qp.id = base_name + "_mzmlstats" ; ///< Identifier
    qp.cvRef = "QC"; ///< cv reference
    qp.cvAcc = "QC:xxxxxxxx"; ///< cv accession
    qp.colTypes.push_back("name");
    qp.colTypes.push_back("value");
    std::vector<String> row;
    row.push_back("#spectra");
    row.push_back(exp.size());
    qp.tableRows.push_back(row);
    row.clear();
    map<Size, UInt> counts;
    for (MSExperiment<Peak1D>::iterator it = exp.begin(); it != exp.end(); ++it)
    {
      counts[it->getMSLevel()]++;
    }
    for (map<Size, UInt>::iterator it = counts.begin(); it != counts.end(); ++it)
    {
      row.push_back("in_ms" + String(it->first));
      row.push_back(String(it->second));
      qp.tableRows.push_back(row);
      row.clear();
    }
    row.push_back("#peaks");
    row.push_back(exp.getSize());
    qp.tableRows.push_back(row);
    row.clear();
    row.push_back("#chromatograms");
    row.push_back(exp.getChromatograms().size());
    qp.tableRows.push_back(row);
    row.clear();
    qcmlfile.addQualityParameter(base_name, qp);

    //---tic qp
    qp = QcMLFile::QualityParameter() ;
    qp.name = "ticpoints"; ///< Name
    qp.id = base_name + "_TIC" ; ///< Identifier
    qp.cvRef = "QC"; ///< cv reference
    qp.cvAcc = "QC:xxxxxxxx"; ///< cv accession

    //~ qp.colTypes.push_back("Native_ID");
    qp.colTypes.push_back("RT_(sec)");
    qp.colTypes.push_back("TIC");
    for (Size i = 0; i < exp.size(); ++i)
    {
      UInt sum = 0;
      for (Size j = 0; j < exp[i].size(); ++j)
      {
        sum += exp[i][j].getIntensity();
      }
      std::vector<String> row;
      //~ String nid = exp[i].getNativeID();
      //~ row.push_back(nid.removeWhitespaces());
      row.push_back(exp[i].getRT());
      row.push_back(sum);
      qp.tableRows.push_back(row);
    }
    qcmlfile.addQualityParameter(base_name, qp);

    //---precursor qp
    qp = QcMLFile::QualityParameter() ;
    qp.name = "precursorpoints"; ///< Name
    qp.id = base_name + "_precursor" ; ///< Identifier
    qp.cvRef = "QC"; ///< cv reference
    qp.cvAcc = "QC:xxxxxxxx"; ///< cv accession

    //~ qp.colTypes.push_back("Native_ID");
    qp.colTypes.push_back("RT_(sec)");
    qp.colTypes.push_back("Precursor");
    for (Size i = 0; i < exp.size(); ++i)
    {
      if (exp[i].getMSLevel() == 2)
      {
        std::vector<String> row;
        //~ String nid = exp[i].getNativeID();
        //~ row.push_back(nid.removeWhitespaces());
        row.push_back(exp[i].getRT());
        row.push_back(exp[i].getPrecursors()[0].getMZ());
        qp.tableRows.push_back(row);
      }
    }
    qcmlfile.addQualityParameter(base_name, qp);


    if (inputfile_id != "")
    {
      IdXMLFile().load(inputfile_id, prot_ids, pep_ids);
      cerr << "idXML read ended. Found " << pep_ids.size() << " peptide identifications." << endl;

      ProteinIdentification::SearchParameters params = prot_ids[0].getSearchParameters();
      vector<String> var_mods = params.variable_modifications;

      //---idxml stats qp
      qp = QcMLFile::QualityParameter() ;
      qp.name = "idxmlstats"; ///< Name
      qp.id = base_name + "_idxmlstats" ; ///< Identifier
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:xxxxxxxx"; ///< cv accession
      qp.colTypes.push_back("name");
      qp.colTypes.push_back("value");
      std::vector<String> row;
      row.push_back("database");
      row.push_back(prot_ids.at(0).getSearchParameters().db);
      qp.tableRows.push_back(row);
      row.clear();
      row.push_back("db_version");
      row.push_back(prot_ids.at(0).getSearchParameters().db_version);
      qp.tableRows.push_back(row);
      row.clear();
      row.push_back("taxonomy");
      row.push_back(prot_ids.at(0).getSearchParameters().taxonomy);
      qp.tableRows.push_back(row);
      row.clear();

      UInt spectrum_count = 0;
      Size peptide_hit_count = 0;
      UInt runs_count = 0;
      Size protein_hit_count = 0;
      set<String> peptides;
      set<String> proteins;
      for (Size i = 0; i < pep_ids.size(); ++i)
      {
        if (!pep_ids[i].empty())
        {
          ++spectrum_count;
          peptide_hit_count += pep_ids[i].getHits().size();
          const vector<PeptideHit> & temp_hits = pep_ids[i].getHits();
          for (Size j = 0; j < temp_hits.size(); ++j)
          {
            peptides.insert(temp_hits[j].getSequence().toString());
          }
        }
      }
      for (Size i = 0; i < prot_ids.size(); ++i)
      {
        ++runs_count;
        protein_hit_count += prot_ids[i].getHits().size();
        const vector<ProteinHit> & temp_hits = prot_ids[i].getHits();
        for (Size j = 0; j < temp_hits.size(); ++j)
        {
          proteins.insert(temp_hits[j].getAccession());
        }
      }

      row.push_back("#id_runs");
      row.push_back(runs_count);
      qp.tableRows.push_back(row);
      row.clear();
      row.push_back("#protein_hits");
      row.push_back(runs_count);
      qp.tableRows.push_back(row);
      row.clear();
      row.push_back("#unique_proteins");
      row.push_back(proteins.size());
      qp.tableRows.push_back(row);
      row.clear();
      row.push_back("#sim");
      row.push_back(spectrum_count);
      qp.tableRows.push_back(row);
      row.clear();
      row.push_back("#peptide_hits");
      row.push_back(peptide_hit_count);
      qp.tableRows.push_back(row);
      row.clear();
      row.push_back("#unique_peptides");
      row.push_back(peptides.size());
      qp.tableRows.push_back(row);
      row.clear();

      qcmlfile.addQualityParameter(base_name, qp);

      //---id accuracy stats qp
      qp = QcMLFile::QualityParameter() ;
      qp.name = "pepidpoints";
      qp.id = base_name + "_peptide_identifications" ;
      qp.cvRef = "QC";
      qp.cvAcc = "QC:xxxxxxxx";

      qp.colTypes.push_back("RT");
      qp.colTypes.push_back("MZ");
      qp.colTypes.push_back("uniqueness");
      qp.colTypes.push_back("ProteinID");
      qp.colTypes.push_back("target/decoy");
      qp.colTypes.push_back("Score");
      qp.colTypes.push_back("PeptideSequence");
      qp.colTypes.push_back("Annots");
      qp.colTypes.push_back("Similarity");
      qp.colTypes.push_back("Charge");
      qp.colTypes.push_back("TheoreticalWeight");
      for (UInt w = 0; w < var_mods.size(); ++w)
      {
        qp.colTypes.push_back(String(var_mods[w]).substitute(' ','_'));
      }

      prot_ids[0].getSearchParameters();
      for (vector<PeptideIdentification>::iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
      {
        if (it->getHits().size() > 0)
        {
          std::vector<String> row;
          row.push_back(it->getMetaValue("RT"));
          row.push_back(it->getMetaValue("MZ"));
          PeptideHit tmp;
          if (!AllHits)
          {
            tmp = *it->getHits().begin();
            String logo, logo2;
            String logo3 = "NA";
            vector<String> logo1;
            logo = tmp.getMetaValue("protein_references");
            logo1 = tmp.getProteinAccessions();
            if (tmp.metaValueExists("target_decoy"))
            {
              logo3 = tmp.getMetaValue("target_decoy");
            }
            if (logo1.size() > 0)
            {
              logo2 = logo1[0];
              for (UInt ii = 1; ii < logo1.size(); ++ii)
              {
                logo2 += logo1[ii];
              }
            }
            vector<UInt> pep_mods;
            for (UInt w = 0; w < var_mods.size(); ++w)
            {
              pep_mods.push_back(0);
            }
            //cout<<var_mods[0]<<SEP<<var_mods[1];
            for (AASequence::ConstIterator z =  tmp.getSequence().begin(); z != tmp.getSequence().end(); ++z)
            {
              Residue res = *z;
              String temp;
              if (res.getModification().size() > 0 && res.getModification() != "Carbamidomethyl")
              {
                temp = res.getModification() + " (" + res.getOneLetterCode()  + ")";
                //cout<<res.getModification()<<endl;
                for (UInt w = 0; w < var_mods.size(); ++w)
                {
                  if (temp == var_mods[w])
                  {
                    //cout<<temp;
                    pep_mods[w] += 1;
                  }
                }
              }
            }

            row.push_back(logo);
            row.push_back(logo2);
            row.push_back(logo3);
            row.push_back(tmp.getScore());
            row.push_back(tmp.getSequence().toString());
            row.push_back(tmp.getMetaValue("Number of annotations"));
            row.push_back(tmp.getMetaValue("similarity"));
            row.push_back(tmp.getCharge());
            row.push_back(String((tmp.getSequence().getMonoWeight() + tmp.getCharge() * Constants::PROTON_MASS_U) / tmp.getCharge()));

            for (UInt w = 0; w < var_mods.size(); ++w)
            {
              row.push_back(pep_mods[w]);
            }

            qp.tableRows.push_back(row);
            //exchange with this line, if you want to integrate consensusID parameters
            //out << logo  << SEP << logo2 << SEP << logo3 << SEP << tmp.getScore() << SEP << tmp.getSequence() << SEP << tmp.getCharge() << SEP << tmp.getMetaValue("Number of annotations") << SEP << tmp.getMetaValue("similarity") << endl;;
          }

          else
          {
            for (vector<PeptideHit>::const_iterator tit = it->getHits().begin(); tit != it->getHits().end(); ++tit)
            {
              std::vector<String> row_allhits = row;
              tmp = *tit;
              String logo, logo2;
              String logo3 = "NA";
              vector<String> logo1;
              logo = tmp.getMetaValue("protein_references");
              logo1 = tmp.getProteinAccessions();
              if (tmp.metaValueExists("target_decoy"))
              {
                logo3 = tmp.getMetaValue("target_decoy");
              }
              if (logo1.size() > 0)
              {
                logo2 = logo1[0];
              }
              row_allhits.push_back(logo);
              row_allhits.push_back(logo2);
              row_allhits.push_back(logo3);
              row_allhits.push_back(tmp.getMetaValue("predicted_PT"));
              row_allhits.push_back(tmp.getSequence().toString());
              row_allhits.push_back(tmp.getCharge());
              row_allhits.push_back(tmp.getMetaValue("Number of annotations"));
              row_allhits.push_back(tmp.getMetaValue("similarity"));
              qp.tableRows.push_back(row_allhits);
            }
          }
        }
      }
      qcmlfile.addQualityParameter(base_name, qp);
    }


    if (inputfile_feature != "")
    {
      cout << "Reading featureXML file..." << endl;

      FeatureMap<> map;
      FeatureXMLFile f;
      f.load(inputfile_feature, map);
      //~ UInt fiter = 0;
      map.sortByRT();
      map.updateRanges();

      //---fxml stats qp
      qp = QcMLFile::QualityParameter() ;
      qp.name = "fxmlstats"; ///< Name
      qp.id = base_name + "_fxmlstats" ; ///< Identifier
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:xxxxxxxx"; ///< cv accession
      qp.colTypes.push_back("name");
      qp.colTypes.push_back("value");
      std::vector<String> row;
      row.push_back("#features");
      row.push_back(map.size());
      qp.tableRows.push_back(row);
      row.clear();

      // Charge distribution and TIC
      Map<UInt, UInt> charges;
      DoubleReal tic = 0.0;
      for (Size i = 0; i < map.size(); ++i)
      {
        charges[map[i].getCharge()]++;
        tic += map[i].getIntensity();
      }
      row.push_back("#feature_tic");
      row.push_back(tic);
      qp.tableRows.push_back(row);
      row.clear();

      for (Map<UInt, UInt>::const_iterator it = charges.begin(); it != charges.end(); ++it)
      {
        row.push_back(String("charge") + String(it->first) + String("_distribution"));
        row.push_back(String(it->second));
        qp.tableRows.push_back(row);
        row.clear();
      }

      qcmlfile.addQualityParameter(base_name, qp);
    }

    if (inputfile_feature != "" && !remove_duplicate_features)
    {
      cout << "Reading featureXML file..." << endl;

      qp = QcMLFile::QualityParameter() ;
      qp.name = "featurepoints"; ///< Name
      qp.id = base_name + "_features" ; ///< Identifier
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:xxxxxxxx"; ///< cv accession

      qp.colTypes.push_back("MZ");
      qp.colTypes.push_back("RT");
      qp.colTypes.push_back("Intensity");
      qp.colTypes.push_back("Charge");

      FeatureMap<> map;
      FeatureXMLFile f;
      f.load(inputfile_feature, map);
      UInt fiter = 0;
      map.sortByRT();
      //ofstream out(outputfile_name.c_str());

      while (fiter < map.size())
      {
        std::vector<String> row;
        row.push_back(map[fiter].getMZ());
        row.push_back(map[fiter].getRT());
        row.push_back(map[fiter].getIntensity());
        row.push_back(map[fiter].getCharge());
        fiter++;
        qp.tableRows.push_back(row);
      }
      qcmlfile.addQualityParameter(base_name, qp);
    }
    else if (inputfile_feature != "" && remove_duplicate_features)
    {
      cout << "Reading featureXML file..." << endl;

      qp = QcMLFile::QualityParameter() ;
      qp.name = "featurepoints"; ///< Name
      qp.id = base_name + "_features" ; ///< Identifier
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:xxxxxxxx"; ///< cv accession

      qp.colTypes.push_back("MZ");
      qp.colTypes.push_back("RT");
      qp.colTypes.push_back("Intensity");
      qp.colTypes.push_back("Charge");

      FeatureMap<> map, map_out;
      FeatureXMLFile f;
      f.load(inputfile_feature, map);
      UInt fiter = 0;
      map.sortByRT();
      while (fiter < map.size())
      {
        FeatureMap<> map_tmp;
        for (UInt k = fiter; k <= map.size(); ++k)
        {
          if (abs(map[fiter].getRT() - map[k].getRT()) < 0.1)
          {
            cout << fiter << endl;
            map_tmp.push_back(map[k]);
          }
          else
          {
            fiter = k;
            break;
          }
        }
        map_tmp.sortByMZ();
        UInt retif = 1;
        map_out.push_back(map_tmp[0]);
        while (retif < map_tmp.size())
        {
          if (abs(map_tmp[retif].getMZ() - map_tmp[retif - 1].getMZ()) > 0.01)
          {
            cout << "equal RT, but mass different" << endl;
            map_out.push_back(map_tmp[retif]);
          }
          retif++;
        }
      }

      //~ FeatureXMLFile().store(out_feature, map_out); //TODO into vec of vecs ...
      qcmlfile.addQualityParameter(base_name, qp);

    }
    if (inputfile_consensus != "")
    {
      cout << "Reading consensusXML file..." << endl;
      ConsensusXMLFile f;
      ConsensusMap map;
      f.load(inputfile_consensus, map);
      //~ String CONSENSUS_NAME = "_consensus.tsv";
      //~ String combined_out = outputfile_name + CONSENSUS_NAME;
      //~ ofstream out(combined_out.c_str());

      qp = QcMLFile::QualityParameter() ;
      qp.name = "consensuspoints"; ///< Name
      qp.id = base_name + "_consensuses" ; ///< Identifier
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:xxxxxxxx"; ///< cv accession

      qp.colTypes.push_back("Native_spectrum_ID");
      qp.colTypes.push_back("DECON_RT_(sec)");
      qp.colTypes.push_back("DECON_MZ_(Th)");
      qp.colTypes.push_back("DECON_Intensity");
      qp.colTypes.push_back("Feature_RT_(sec)");
      qp.colTypes.push_back("Feature_MZ_(Th)");
      qp.colTypes.push_back("Feature_Intensity");
      qp.colTypes.push_back("Feature_Charge");
      for (ConsensusMap::const_iterator cmit = map.begin(); cmit != map.end(); ++cmit)
      {
        ConsensusFeature CF = *cmit;
        for (ConsensusFeature::const_iterator cfit = cmit->begin(); cfit != cmit->end(); ++cfit)
        {
          std::vector<String> row;
          FeatureHandle FH = *cfit;
          row.push_back(CF.getMetaValue("spectrum_native_id"));
          row.push_back(CF.getRT());row.push_back(CF.getMZ());
          row.push_back(CF.getIntensity());
          row.push_back(FH.getRT());
          row.push_back(FH.getMZ());
          row.push_back(FH.getCharge());
          qp.tableRows.push_back(row);
        }
      }
      qcmlfile.addQualityParameter(base_name, qp);
    }
    if (Ms1)
    {
      qp = QcMLFile::QualityParameter() ;
      qp.name = "ms1stats"; ///< Name
      qp.id = base_name + "_ms1" ; ///< Identifier
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:xxxxxxxx"; ///< cv accession

      qp.colTypes.push_back("Native_ID");
      qp.colTypes.push_back("RT_(sec)");
      qp.colTypes.push_back("MZ_(Th)");
      qp.colTypes.push_back("Intensity");
      for (Size i = 0; i < exp.size(); ++i)
      {
        if (exp[i].getMSLevel() == 1)
        {
          for (Size j = 0; j < exp[i].size(); ++j)
          {
            std::vector<String> row;
            String nid = exp[i].getNativeID();
            row.push_back(nid.removeWhitespaces());
            row.push_back(exp[i].getRT());
            row.push_back(exp[i][j].getMZ());
            row.push_back(exp[i][j].getIntensity());
            qp.tableRows.push_back(row);
          }
        }
      }
      qcmlfile.addQualityParameter(base_name, qp);
    }
    if (Ms2)
    {
      qp = QcMLFile::QualityParameter() ;
      qp.name = "ms2stats"; ///< Name
      qp.id = base_name + "_ms2" ; ///< Identifier
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:xxxxxxxx"; ///< cv accession

      qp.colTypes.push_back("Native_ID");
      qp.colTypes.push_back("RT_(sec)");
      qp.colTypes.push_back("MZ_(Th)");
      qp.colTypes.push_back("Intensity");
      qp.colTypes.push_back("Precursor");
      for (Size i = 0; i < exp.size(); ++i)
      {
        if (exp[i].getMSLevel() == 2)
        {
          for (Size j = 0; j < exp[i].size(); ++j)
          {
            std::vector<String> row;
            String nid = exp[i].getNativeID();
            row.push_back(nid.removeWhitespaces());
            row.push_back(exp[i].getRT());
            row.push_back(exp[i][j].getMZ());
            row.push_back(exp[i][j].getIntensity());
            row.push_back(exp[i].getPrecursors()[0].getMZ());
            qp.tableRows.push_back(row);
          }
        }
      }
      qcmlfile.addQualityParameter(base_name, qp);
    }
    qcmlfile.store(outputfile_name);
    return EXECUTION_OK;
  }


};
int main(int argc, const char ** argv)
{
  TOPPQCCalculator tool;
  return tool.main(argc, argv);
}

/// @endcond
