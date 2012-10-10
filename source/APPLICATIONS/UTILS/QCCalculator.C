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
// $Maintainer: Sven Nahnsen $
// $Author: Sven Nahnsen $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
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
#include <vector>
#include <fstream>

#include <map>

using namespace OpenMS;
using namespace std;

#define SEP "\t"


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
    registerInputFile_("in_id", "<file>", "", "Input idXML file containing the identifications. Your identifications will be exported in an easy-to-read format", false);
    setValidFormats_("in_id", StringList::create("idXML"));
    registerInputFile_("feature", "<file>", "", "feature input file (this is relevant for most QC issues)", false);
    setValidFormats_("feature", StringList::create("featureXML"));
    registerInputFile_("consensus", "<file>", "", "consensus input file (this is only used for charge state deconvoluted output. Use the consensusXML output form the DeCharger)", false);
    setValidFormats_("consensus", StringList::create("consensusXML"));
    registerInputFile_("raw", "<file>", "", "raw data input file (this is relevant if you want to look at MS1, MS2 and precursor peak information)", false);
    setValidFormats_("raw", StringList::create("mzML"));
    registerOutputFile_("out_id", "<file>", "", "Your export file for id stats.", false, true);
    setValidFormats_("out_id", StringList::create("_id.tsv"));
    registerOutputFile_("out_consensus", "<file>", "", "Your export file for consensus stats.", false, true);
    setValidFormats_("out_consensus", StringList::create("_consensus.tsv"));
    registerOutputFile_("out_feature", "<file>", "", "Your export file for feature stats.", false, true);
    setValidFormats_("out_feature", StringList::create("tsv"));
    registerOutputFile_("out_precursor", "<file>", "", "Your export file for precursor stats.", false, true);
    setValidFormats_("out_precursor", StringList::create("_PRECURSOR.tsv"));
    registerOutputFile_("out_MS2", "<file>", "", "Your export file for MS2 stats.", false, true);
    setValidFormats_("out_MS2", StringList::create("_MS2.tsv"));
    registerOutputFile_("out_MS1", "<file>", "", "Your export file for MS1 stats.", false, true);
    setValidFormats_("out_MS1", StringList::create("_MS1.tsv"));
    registerOutputFile_("out_TIC", "<file>", "", "Your export file for TIC stats,if you are exporting MS1 and MS2 data, but you want to have TIC intensities instead of the MS1 table. Note that TICs, in this case are the sums of all intensities per spectrum.", false, true);
    setValidFormats_("out_TIC", StringList::create("_TIC.tsv"));
    registerFlag_("AllHits", "If set, the exporter takes all peptide candidates per spectrum into account.");
    registerFlag_("remove_duplicate_features", "This flag should be set, if you work with a set of merged features.");
  }

  ExitCodes main_(int, const char **)
  {
    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;
    ProteinHit temp_protein_hit;
    String inputfile_name = "";
    String outputfile_name = "";
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    inputfile_name              = getStringOption_("in_id");
    String inputfile_feature    = getStringOption_("feature");
    String inputfile_consensus  = getStringOption_("consensus");
    String inputfile_raw        = getStringOption_("raw");
    //~ outputfile_name             = getStringOption_("out");
		
		String out_id = getStringOption_("out_id");
    String out_feature = getStringOption_("out_feature");
    String out_consensus = getStringOption_("out_consensus");
    String out_ms1 = getStringOption_("out_MS1");
    String out_ms2 = getStringOption_("out_MS2");
    String out_precursor = getStringOption_("out_precursor");
    String out_tic = getStringOption_("out_TIC");

    bool AllHits(getFlag_("AllHits"));
    bool remove_duplicate_features(getFlag_("remove_duplicate_features"));
    //-------------------------------------------------------------
    // reading input
    //------------------------------------------------------------
    if (inputfile_name != "" && out_id != "")      // -InclusionList was given
    {
      //~ String ID_NAME = "_id.tsv";
      IdXMLFile().load(inputfile_name, prot_ids, pep_ids);
      cerr << "ended. Found " << pep_ids.size() << " peptide identifications." << endl;
      cerr << "Writing text file..." << endl;
      ProteinIdentification::SearchParameters params = prot_ids[0].getSearchParameters();
      vector<String> var_mods = params.variable_modifications;
      //~ String combined_out = outputfile_name + ID_NAME;
      //~ ofstream out(combined_out.c_str());
      ofstream out(out_id.c_str());
      out << "RT" << SEP << "MZ" << SEP << "uniqueness" << SEP << "ProteinID" << SEP << "target/decoy" << SEP << "Score" << SEP << "PeptideSequence" << SEP << "Annots" << SEP << "Similarity" << SEP << "Charge" << SEP << "TheoreticalWeight";
      for (UInt w = 0; w < var_mods.size(); ++w)
      {
        out << SEP << var_mods[w];
      }
      out << endl;
      prot_ids[0].getSearchParameters();
      for (vector<PeptideIdentification>::iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
      {
        if (it->getHits().size() > 0)
        {
          out << it->getMetaValue("RT") << SEP << it->getMetaValue("MZ") << SEP;
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

            out << logo  << SEP << logo2 << SEP << logo3 << SEP << tmp.getScore() << SEP << tmp.getSequence() << SEP << tmp.getMetaValue("Number of annotations") << SEP << tmp.getMetaValue("similarity") << SEP << tmp.getCharge() << SEP << precisionWrapper((tmp.getSequence().getMonoWeight() + tmp.getCharge() * Constants::PROTON_MASS_U) / tmp.getCharge());
            for (UInt w = 0; w < var_mods.size(); ++w)
            {
              out << SEP << pep_mods[w];
            }
            out << endl;
            //exchange with this line, if you want to integrate consensusID parameters
            //out << logo  << SEP << logo2 << SEP << logo3 << SEP << tmp.getScore() << SEP << tmp.getSequence() << SEP << tmp.getCharge() << SEP << tmp.getMetaValue("Number of annotations") << SEP << tmp.getMetaValue("similarity") << endl;;
          }
          else
          {
            for (vector<PeptideHit>::const_iterator tit = it->getHits().begin(); tit != it->getHits().end(); ++tit)
            {
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
              out << logo  << SEP << logo2 << SEP << logo3 << SEP << tmp.getMetaValue("predicted_PT") << SEP << tmp.getSequence() << SEP << tmp.getCharge() << SEP << tmp.getMetaValue("Number of annotations") << SEP << tmp.getMetaValue("similarity");
              out << endl;
            }
          }
        }

      }
      out.close();
    }
    if (inputfile_feature != "" && !remove_duplicate_features && out_feature	!= "")      //
    {
      //~ String FEATURE_NAME_NO_REMOVE = "_features.tsv";
      //~ String combined_out = outputfile_name + FEATURE_NAME_NO_REMOVE;
      //~ ofstream out(combined_out.c_str());
      ofstream out(out_feature.c_str());
      FeatureMap<> map;
      FeatureXMLFile f;
      f.load(inputfile_feature, map);
      UInt fiter = 0;
      map.sortByRT();
      //ofstream out(outputfile_name.c_str());
      out << "MZ" << SEP << "RT" << SEP << "Intensity" << SEP << "Charge" << endl;
      while (fiter < map.size())
      {
        out << map[fiter].getMZ() << SEP << map[fiter].getRT() << SEP << map[fiter].getIntensity() << SEP << map[fiter].getCharge() << endl;
        fiter++;
      }
      out.close();
    }
    else if (inputfile_feature != "" && remove_duplicate_features && out_feature != "")      //
    {
      FeatureMap<> map, map_out;
      FeatureXMLFile f;
      f.load(inputfile_feature, map);
      UInt fiter = 0;
      map.sortByRT();
      //ofstream out(outputfile_name.c_str());
      //map_tmp.push_back(map[0]);
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

      //~ String FEATURE_NAME_REMOVE = "_features.tsv";
      //~ String combined_out = outputfile_name + FEATURE_NAME_REMOVE;
      //~ ofstream out(combined_out.c_str());
      ofstream out(out_feature.c_str());
      FeatureXMLFile().store(out_feature, map_out);
    }
    if (inputfile_consensus != "" && out_consensus != "")
    {
      cout << "Reading consensusXML file..." << endl;
      ConsensusXMLFile f;
      ConsensusMap map;
      f.load(inputfile_consensus, map);
      //~ String CONSENSUS_NAME = "_consensus.tsv";
      //~ String combined_out = outputfile_name + CONSENSUS_NAME;
      //~ ofstream out(combined_out.c_str());
      ofstream out(out_consensus.c_str());
      cout << "Writing text file..." << endl;
      out << "Native spectrum ID" << SEP << "DECON RT (sec)" << SEP << "DECON MZ (Th)" << SEP << "DECON Intensity" << SEP << " Feature RT (sec)" << SEP << " Feature MZ (Th)" << SEP << "Feature Intensity" << SEP << "Feature Charge" << endl;
      for (ConsensusMap::const_iterator cmit = map.begin(); cmit != map.end(); ++cmit)
      {
        ConsensusFeature CF = *cmit;
        for (ConsensusFeature::const_iterator cfit = cmit->begin(); cfit != cmit->end(); ++cfit)
        {
          FeatureHandle FH = *cfit;
          out << CF.getMetaValue("spectrum_native_id") << SEP << CF.getRT() << SEP << CF.getMZ() << SEP << CF.getIntensity() << SEP << FH.getRT() << SEP << FH.getMZ() << SEP << FH.getCharge() << endl;
        }
      }
      out.close();
    }
		if (out_ms1 != "")
		{
      cout << "Reading mzML file..." << endl;
      MzMLFile mz_data_file;
      MSExperiment<Peak1D> exp;
      MzMLFile().load(inputfile_raw, exp);
      cout << "Writing text file..." << endl;
      //~ String MS1_NAME = "_MS1.tsv";
      //~ String combined_out = outputfile_name + MS1_NAME;
      //~ ofstream out(combined_out.c_str());
      ofstream out(out_ms1.c_str());
      //three different output formats can be created depending on the choice of the ouput type

      out << "Native ID" << SEP << "RT (sec)" << SEP << "MZ (Th)" << SEP << "Intensity" << endl;
      for (Size i = 0; i < exp.size(); ++i)
      {
        if (exp[i].getMSLevel() == 1)
        {
          for (Size j = 0; j < exp[i].size(); ++j)
          {
            out << exp[i].getNativeID() << SEP << exp[i].getRT() << SEP << exp[i][j].getMZ() << SEP << exp[i][j].getIntensity() << endl;
          }
        }
      }
      out.close();
    }
		if (out_ms2 != "")
		{
      cout << "Reading mzML file..." << endl;
      MzMLFile mz_data_file;
      MSExperiment<Peak1D> exp;
      MzMLFile().load(inputfile_raw, exp);
      cout << "Writing text file..." << endl;
      //~ String MS2_NAME = "_MS2.tsv";
      //~ String combined_out = outputfile_name + MS2_NAME;
      //~ ofstream out(combined_out.c_str());
      ofstream out(out_ms2.c_str());
      out << "Native ID" << SEP << "RT (sec)" << SEP << "MZ (Th)" << SEP << "Intensity" << SEP << "Precursor" << endl;
      for (Size i = 0; i < exp.size(); ++i)
      {
        if (exp[i].getMSLevel() == 2)
        {
          for (Size j = 0; j < exp[i].size(); ++j)
          {
            out << exp[i].getNativeID() << SEP << exp[i].getRT() << SEP << exp[i][j].getMZ() << SEP << exp[i][j].getIntensity() << SEP << exp[i].getPrecursors()[0].getMZ() << endl;
          }
        }
      }
      out.close();
    }
		if (out_precursor != "")
		{
      cout << "Reading mzML file..." << endl;
      MzMLFile mz_data_file;
      MSExperiment<Peak1D> exp;
      MzMLFile().load(inputfile_raw, exp);
      cout << "Writing text file..." << endl;
      //~ String PREC_NAME = "_PRECURSOR.tsv";
      //~ String combined_out = outputfile_name + PREC_NAME;
      //~ ofstream out(combined_out.c_str());
      ofstream out(out_precursor.c_str());
      out << "Native ID" << SEP << "RT (sec)" << SEP << "Precursor" << endl;
      for (Size i = 0; i < exp.size(); ++i)
      {
        if (exp[i].getMSLevel() == 2)
        {
          out << exp[i].getNativeID() << SEP << exp[i].getRT() << SEP << exp[i].getPrecursors()[0].getMZ() << endl;
        }
      }
      out.close();
    }

		if (out_tic != "")
		{
      cout << "Reading mzML file..." << endl;
      MzMLFile mz_data_file;
      MSExperiment<Peak1D> exp;
      MzMLFile().load(inputfile_raw, exp);
      cout << "Writing text file..." << endl;
      //~ String MS1_NAME = "_TIC.tsv";
      //~ String combined_out = outputfile_name + MS1_NAME;
      //~ ofstream out(combined_out.c_str());
      ofstream out(out_tic.c_str());
      //three different output formats can be created depending on the choice of the ouput type

      out << "Native ID" << SEP << "RT (sec)" << SEP << "TIC" << endl;
      for (Size i = 0; i < exp.size(); ++i)
      {
        UInt sum = 0;
        for (Size j = 0; j < exp[i].size(); ++j)
        {
          sum += exp[i][j].getIntensity();
        }
        out << exp[i].getNativeID() << SEP << exp[i].getRT() << SEP << sum << endl;
      }
      out.close();
    }


    return EXECUTION_OK;
  }

};
int main(int argc, const char ** argv)
{
  TOPPQCCalculator tool;
  return tool.main(argc, argv);
}

/// @endcond
