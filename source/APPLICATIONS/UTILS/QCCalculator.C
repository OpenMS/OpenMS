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
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <QFileInfo>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/median.hpp>

#include <vector>
#include <map>

using namespace OpenMS;
using namespace std;
using namespace boost::accumulators;

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

  DoubleReal getMassDifference(DoubleReal theo_mz, DoubleReal exp_mz, bool use_ppm)
  {
    DoubleReal error(exp_mz - theo_mz);
    if (use_ppm)
    {
      error = error / theo_mz * (DoubleReal)1e6;
      //~ error = (1-exp_mz/theo_mz) * (DoubleReal)1e6;
    }
    return error;
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
    qp.name = "rmzML file"; ///< Name
    qp.id = base_name + "_origin_rawfile_name"; ///< Identifier
    qp.cvRef = "PSI"; ///< cv reference
    qp.cvAcc = "MS:1000584";
    qp.value = base_name;
    qcmlfile.addQualityParameter(base_name, qp);

    //---MS2 distribution qp
    exp.sortSpectra();
    UInt min_mz = std::numeric_limits<UInt>::max();
    UInt max_mz = 0;
    std::map<Size, UInt> mslevelcounts;

    qp = QcMLFile::QualityParameter() ;
    qp.name = "precursor distribution"; ///< Name
    qp.id = base_name + "_precursors" ; ///< Identifier
    qp.cvRef = "QC"; ///< cv reference
    qp.cvAcc = "QC:0000005";

    qp.colTypes.push_back("RT_(sec)");
    qp.colTypes.push_back("Precursor");
    for (Size i = 0; i < exp.size(); ++i)
    {
      mslevelcounts[exp[i].getMSLevel()]++;
      if (exp[i].getMSLevel() == 2)
      {
        if (exp[i].getPrecursors().front().getMZ() < min_mz)
        {
          min_mz = exp[i].getPrecursors().front().getMZ();
        }
        if (exp[i].getPrecursors().front().getMZ() > max_mz)
        {
          max_mz = exp[i].getPrecursors().front().getMZ();
        }
        std::vector<String> row;
        row.push_back(exp[i].getRT());
        row.push_back(exp[i].getPrecursors().front().getMZ());
        qp.tableRows.push_back(row);
      }
    }
    qcmlfile.addQualityParameter(base_name, qp);

    //---aquisition results qp
    qp = QcMLFile::QualityParameter() ;
    qp.name = "aquisition ranges"; ///< Name
    qp.id = base_name + "_aquisition" ; ///< Identifier
    qp.cvRef = "QC"; ///< cv reference
    qp.cvAcc = "QC:0000003"; ///< cv accession for "aquisition results"
    qp.colTypes.push_back("name");
    qp.colTypes.push_back("value");
    std::vector<String> row;
    row.push_back("#spectra");
    row.push_back(exp.size());
    qp.tableRows.push_back(row);
    row.clear();
    for (map<Size, UInt>::iterator it = mslevelcounts.begin(); it != mslevelcounts.end(); ++it)
    {
      row.push_back("in_ms" + String(it->first));
      row.push_back(String(it->second));
      qp.tableRows.push_back(row);
      row.clear();
    }
    row.push_back("#peaks");
    row.push_back(exp.getSize()); //TODO this is bug-ish
    qp.tableRows.push_back(row);
    row.clear();
    row.push_back("#chromatograms");
    row.push_back(exp.getChromatograms().size());
    qp.tableRows.push_back(row);
    row.clear();
    row.push_back("RT_min");
    row.push_back(exp.begin()->getRT());
    qp.tableRows.push_back(row);
    row.clear();
    row.push_back("RT_max");
    row.push_back(exp.back().getRT());
    qp.tableRows.push_back(row);
    row.clear();
    row.push_back("mz_min");
    row.push_back(min_mz);
    qp.tableRows.push_back(row);
    row.clear();
    row.push_back("mz_max");
    row.push_back(max_mz);
    qp.tableRows.push_back(row);
    row.clear();
    qcmlfile.addQualityParameter(base_name, qp);

    //---ion current stability ( & tic ) qp
    qp = QcMLFile::QualityParameter() ;
    qp.name = "total ion current distribution"; ///< Name
    qp.id = base_name + "_tic" ; ///< Identifier
    qp.cvRef = "QC"; ///< cv reference
    qp.cvAcc = "QC:0000007";
    //~ qp.colTypes.push_back("Native_ID");
    qp.colTypes.push_back("RT_(sec)");
    qp.colTypes.push_back("TIC");
    UInt max = 0;
    Size below_10k = 0;
    for (Size i = 0; i < exp.size(); ++i)
    {
      UInt sum = 0;
      for (Size j = 0; j < exp[i].size(); ++j)
      {
        sum += exp[i][j].getIntensity();
      }
      if (sum > max)
      {
        max = sum;
      }
      if (sum < 10000)
      {
        ++below_10k;
      }
      std::vector<String> row;
      row.push_back(exp[i].getRT());
      row.push_back(sum);
      qp.tableRows.push_back(row);
    }
    qcmlfile.addQualityParameter(base_name, qp);


    qp = QcMLFile::QualityParameter() ;
    qp.name = "ion current stability"; ///< Name
    qp.id = base_name + "_ics" ; ///< Identifier
    qp.cvRef = "QC"; ///< cv reference
    qp.cvAcc = "QC:0000017";
    qp.colTypes.push_back("name");
    qp.colTypes.push_back("value");
    row.clear();
    row.push_back("max_ic");
    row.push_back(max);
    row.push_back("%ions<10k");
    row.push_back((100/exp.size())*below_10k);

    if (inputfile_id != "")
    {
      IdXMLFile().load(inputfile_id, prot_ids, pep_ids);
      cerr << "idXML read ended. Found " << pep_ids.size() << " peptide identifications." << endl;

      ProteinIdentification::SearchParameters params = prot_ids[0].getSearchParameters();
      vector<String> var_mods = params.variable_modifications;

      //---idxml stats qp
      qp = QcMLFile::QualityParameter() ;
      qp.name = "search input details"; ///< Name
      qp.id = base_name + "_search_input" ; ///< Identifier
      qp.cvRef = "MS"; ///< cv reference
      qp.cvAcc = "MS:1001249"; ///< cv accession "basic identification results"
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
      qcmlfile.addQualityParameter(base_name, qp);

      qp = QcMLFile::QualityParameter() ;
      qp.name = "Spectrum identification result details"; ///< Name
      qp.id = base_name + "_search_result" ; ///< Identifier
      qp.cvRef = "MS"; ///< cv reference
      qp.cvAcc = "MS:1001405"; ///< cv accession "basic identification results"
      qp.colTypes.push_back("name");
      qp.colTypes.push_back("value");
      row.push_back("#protein_hits");
      row.push_back(runs_count);
      qp.tableRows.push_back(row);
      row.clear();
      row.push_back("#unique_proteins");
      row.push_back(proteins.size());
      qp.tableRows.push_back(row);
      row.clear();
      row.push_back("#specta_identified");
      row.push_back(spectrum_count);
      qp.tableRows.push_back(row);
      row.clear();
      row.push_back("#PSM");
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
      qp.name = "delta ppm distribution";
      qp.id = base_name + "_delta_ppm" ;
      qp.cvRef = "QC";
      qp.cvAcc = "QC:0000008";

      qp.colTypes.push_back("RT");
      qp.colTypes.push_back("MZ");
      qp.colTypes.push_back("Score");
      qp.colTypes.push_back("PeptideSequence");
      qp.colTypes.push_back("Charge");
      qp.colTypes.push_back("TheoreticalWeight");
      qp.colTypes.push_back("delta_ppm");
      for (UInt w = 0; w < var_mods.size(); ++w)
      {
        qp.colTypes.push_back(String(var_mods[w]).substitute(' ','_'));
      }

      std::vector<DoubleReal> deltas;
      //~ prot_ids[0].getSearchParameters();
      for (vector<PeptideIdentification>::iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
      {
        if (it->getHits().size() > 0)
        {
          std::vector<String> row;
          row.push_back(it->getMetaValue("RT"));
          row.push_back(it->getMetaValue("MZ"));
          PeptideHit tmp;
          vector<UInt> pep_mods;
          for (UInt w = 0; w < var_mods.size(); ++w)
          {
            pep_mods.push_back(0);
          }
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
          row.push_back(tmp.getScore());
          row.push_back(tmp.getSequence().toString().removeWhitespaces());
          row.push_back(tmp.getCharge());
          row.push_back(String((tmp.getSequence().getMonoWeight() + tmp.getCharge() * Constants::PROTON_MASS_U) / tmp.getCharge()));
          DoubleReal dppm = std::abs(getMassDifference(((tmp.getSequence().getMonoWeight() + tmp.getCharge() * Constants::PROTON_MASS_U) / tmp.getCharge()), double(it->getMetaValue("MZ")), true));
          row.push_back(String(dppm));
          deltas.push_back(dppm);
          for (UInt w = 0; w < var_mods.size(); ++w)
          {
            row.push_back(pep_mods[w]);
          }
          qp.tableRows.push_back(row);
        }
      }
      qcmlfile.addQualityParameter(base_name, qp);

      //---mass accuracy stats qp
      qp = QcMLFile::QualityParameter() ;
      qp.name = "mass accuracy"; ///< Name
      qp.id = base_name + "_mass_accuracy" ; ///< Identifier
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000018";
      qp.colTypes.push_back("name");
      qp.colTypes.push_back("value");
      row.clear();
      row.push_back("mean");
      //~ row.push_back(mean(acc));
      row.push_back(OpenMS::Math::mean(deltas.begin(),deltas.end()));
      qp.tableRows.push_back(row);
      row.clear();
      row.push_back("median");
      //~ row.push_back(median(acc));
      row.push_back(OpenMS::Math::median(deltas.begin(),deltas.end()));
      qp.tableRows.push_back(row);
      row.clear();
      qcmlfile.addQualityParameter(base_name, qp);

      //---mass accuracy stats qp
      qp = QcMLFile::QualityParameter() ;
      qp.name = "id ratio"; ///< Name
      qp.id = base_name + "_ratio_id" ; ///< Identifier
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000019";
      qp.colTypes.push_back("name");
      qp.colTypes.push_back("value");
      row.clear();
      row.push_back("ratio");
      row.push_back(String(DoubleReal(pep_ids.size())/DoubleReal(mslevelcounts[2])));
      qp.tableRows.push_back(row);
      row.clear();
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
      qp.name = "featurefinding result details"; ///< Name
      qp.id = base_name + "_featurefinding" ; ///< Identifier
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000020";
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
      qp.name = "feature distribution"; ///< Name
      qp.id = base_name + "_features" ; ///< Identifier
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000009";

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
      qp.name = "feature distribution"; ///< Name
      qp.id = base_name + "_features" ; ///< Identifier
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000009";

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
      qp.cvAcc = "QC:xxxxxxxx"; ///< cv accession "featuremapper results"

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
      qp.cvAcc = "QC:xxxxxxxx"; ///< cv accession for "aquisition results"

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
      qp.cvAcc = "QC:xxxxxxxx"; ///< cv accession for "aquisition results"

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
