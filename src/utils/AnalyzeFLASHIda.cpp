//
// Created by Kyowon Jeong on 2/22/21.
//

// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <QFile>
#include <iomanip>

#include <OpenMS/METADATA/SpectrumLookup.h>
using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page
*/

/// @cond

class TOPPAnalyzeFLASHIda :
    public TOPPBase
{
public:
  TOPPAnalyzeFLASHIda() :
      TOPPBase("AnalyzeFLASHIda", ".", false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    //registerOutputFileList_("out_spec", "<file for MS1, file for MS2, ...>", {""},
    //                            "output files (tsv) - spectrum level deconvoluted masses per ms level", false);
    //

    registerInputFileList_("in", "<files>", {}, "Input files");
    registerOutputFile_("out", "<file>", "", "");
    //setValidFormats_("i1n", ListUtils::create<String>("tsv"));
    //setValidFormats_("i1n", ListUtils::create<String>("tsv"));
  }

  struct TopPicItem
  {
  public:
    TopPicItem() = default;

    TopPicItem(String in)
    {
      str_ = in;
      vector<String> results;
      stringstream tmp_stream(in);
      String str;

      while (getline(tmp_stream, str, ','))
      {
        results.push_back(str);
      }
      //cout<<str<<endl;
      protein_acc_ = results[0];
      proteform_id_ = stoi(results[1]);
      rt_ = stod(results[2]);
      precursor_mass_ = stod(results[4]);
      precursor_mz_ = stod(results[5]);
      precursor_intensity = stod(results[6]);
        mass_intensity_ = stod(results[7]);
      charge_ = stoi(results[8]);
      unexp_mod_ = stod(results[9]);
      q_value_ = stod(results[15]);
      e_value_ = stod(results[16]);
    }

    //const String header="Findex,";// = "Findex,ACC,ProID,RT,PrecursorMonoMass,PrecursorAvgMass,PrecursorMz,PrecursorIntensity,PrecursorCharge,PTM,ChargeCos,ChargeSNR,Cos,SNR,ChargeScore,Qscore,Evalue,Class";

    String str_;
    double rt_;
    int charge_;
    double precursor_intensity;
    double mass_intensity_;
    double precursor_mass_;
    double precursor_mz_;
    int proteform_id_;
    String protein_acc_;
    double unexp_mod_;
    double e_value_;
    double q_value_;

  };

  ExitCodes main_(int, const char **) override
  {
    auto infiles = getStringList_("in");
    auto outfile = getStringOption_("out");
    double tol = 20.0;
    fstream outstream;
    outstream.open(outfile, fstream::out); //

    MSExperiment map;
    MzMLFile mzml;

    mzml.setLogType(log_type_);
    mzml.load(infiles[0], map);
    auto iso_mzs = DoubleList{126.127725, 127.124760, 128.134433, 129.131468, 130.141141, 131.138176};

    auto mz2report = std::map<double, std::vector<double>>();
    auto mz2reportcount = std::map<double, int>();
    auto msscan2mz = std::map<int, double>();

    for (auto &it: map)
    {
      if (it.getMSLevel() == 1)
      {
        continue;
      }

      String am;
      for (auto &activation_method: it.getPrecursors()[0].getActivationMethods())
      {
        am = Precursor::NamesOfActivationMethodShort[activation_method];
        //std::cout<<am<<" ";
      }
      //std::cout<<std::endl;

      int scan = SpectrumLookup::extractScanNumber(it.getNativeID(),
                                                   map.getSourceFiles()[0].getNativeIDTypeAccession());
      double pmz = it.getPrecursors()[0].getMZ();

      int current_ch = 0;
      int nonzero = 0;
      auto iso_intensities = std::vector<double>(iso_mzs.size(), 0.0);
      for (auto &peak: it)
      {
        if (current_ch >= iso_mzs.size())
        {
          break;
        }
        if (peak.getMZ() < iso_mzs[current_ch] - iso_mzs[current_ch] * tol * 1e-6)
        {
          continue;
        }
        if (peak.getMZ() > iso_mzs[current_ch] + iso_mzs[current_ch] * tol * 1e-6)
        {
          if (iso_intensities[current_ch] > 0)
          {
            nonzero++;
          }
          current_ch++;
          continue;
        }
        iso_intensities[current_ch] += peak.getIntensity();
      }

      msscan2mz[scan] = pmz;

      if (mz2reportcount.find(pmz) == mz2reportcount.end())
      {
        mz2reportcount[pmz] = 0;
      }
      int nonzeromax = mz2reportcount[pmz];
      if (nonzeromax < nonzero)
      {
        mz2reportcount[pmz] = nonzero;
      }
      else
      {
        continue;
      }

      double sum = .0;
      for (auto v: iso_intensities)
      {
        sum += v;
      }
      if (sum > 0)
      {
        for (int i = 0; i < 6; i++)
        {
          iso_intensities[i] /= sum;
        }
      }

      mz2report[pmz] = iso_intensities;
    }

    auto in = infiles[1];
    std::ifstream in_tsv(in);
    String line;
    bool start = false;
    while (std::getline(in_tsv, line))
    {

      line = line.substr(0, line.length() - 1);

      if (line.hasPrefix("Data"))
      {
        start = true;
        outstream << line << "\tch1\tch2\tch3\tch4\tch5\tch6\n";
        continue;
      }
      if (!start)
      {
        continue;
      }

      if (!line.hasSubstring("[IodoTMTsixplex]"))
      {
        continue;
      }
      vector<String> tokens;
      stringstream tmp_stream(line);
      String str;

      while (getline(tmp_stream, str, '\t'))
      {
        tokens.push_back(str);
      }
      int scan = stoi(tokens[4]);

      double pmz = msscan2mz[scan];

      if (mz2report.find(pmz) == mz2report.end())
      {
        continue;
      }
      auto rions = mz2report[pmz];
      outstream << line;
      for (int i = 0; i < 6; i++)
      {
        outstream << "\t" << rions[i];
      }
      outstream << "\n";
    }

    in_tsv.close();
    outstream.close();
    return EXECUTION_OK;
  }


  ExitCodes main_1(int, const char **) //override
  {
    auto infiles = getStringList_("in");
    auto outfile = getStringOption_("out");

    fstream outstream;
    String header = "";
    outstream.open(outfile, fstream::out); //
    vector<vector<TopPicItem>> results;
    map<String, map<int, vector<TopPicItem>>> acc_map; // acc, input index, matches

    for (int i = 0; i < infiles.size(); i++)
    {
      auto in = infiles[i];
      std::ifstream in_trainstream(in);
      String line;
      //bool start = false;
      results.emplace_back();

      while (std::getline(in_trainstream, line))
      {
        if (line.hasPrefix("ACC"))
        {
          header = "FileIndex," + line;
          continue;
        }
        if (line.hasSuffix("F"))
        {
          //outstream << "FileIndex,"<<line << "\n";
          continue;
        }
        TopPicItem item(line);

        //if (item.e_value > 1e-2)
        //{ //
        //  continue;
        //}
        results[i].push_back(item);
      }
      cout << "# entries: " << results[i].size() << endl;
      in_trainstream.close();
    }
    outstream << header << "\n";

    //map<String, map<int, vector<TopPicItem>>> acc_map; // acc, input index, matches
    for (int i = 0; i < infiles.size(); i++)
    {
      auto r = results[i];
      for (auto &item : r)
      {
        if (item.unexp_mod_ == .0)
        { // TODO
          continue;
        }
        if (acc_map.find(item.protein_acc_) == acc_map.end())
        {
          acc_map[item.protein_acc_] = map<int, vector<TopPicItem>>();
          for (int j = 0; j < infiles.size(); j++)
          {
            acc_map[item.protein_acc_][j] = vector<TopPicItem>();
          }
        }
        acc_map[item.protein_acc_][i].push_back(item);
      }
    }

    for (auto &item: acc_map)
    {
      map<int, set<double>> s_map; // input index, set of  masses for this acc
      for (auto &s_item : item.second)
      {
        s_map[s_item.first] = set<double>();
        for (auto &id : s_item.second)
        {
          s_map[s_item.first].insert(id.precursor_mass_);
        }
      }
      // cout<<item.first<<endl;

      set<double> intersect;

      for (int j = 0; j < infiles.size(); j++)
      {
        auto s1 = s_map[j];
        for (int k = j + 1; k < infiles.size(); k++)
        {
          auto s2 = s_map[k];
          if (s1.empty() || s2.empty())
          {
            continue; // output
          }
          int s1o = 0, s2o = 0;
          for (auto n1: s1)
          {
            bool intersected = false;
            for (auto n2: s2)
            {
              if (abs(n1 - n2) < 3.0)
              {
                intersect.insert(n1);
                intersect.insert(n2);
                intersected = true;
              }
            }
            if (intersected)
            {
              s1o++;
            }
          }

          for (auto n2: s2)
          {
            bool intersected = false;
            for (auto n1: s1)
            {
              if (abs(n1 - n2) < 3.0)
              {
                intersect.insert(n1);
                intersect.insert(n2);
                intersected = true;
              }
            }
            if (intersected)
            {
              s2o++;
            }
          }
          //cout<< s1.size() << " " << s2.size() << " " << s1o<< " " << s2o << endl;

          //cout<< item.first << " : " << s1.size() << " and " << s2.size() << " : intersect " << intersect.size() << endl;
          //cout<< s1.size() << " " << s2.size() << " " << s1o<< endl;
        }
      }
      // acc, input index, matches
      for (auto &s_item : item.second)
      {
        for (auto &id : s_item.second)
        {
          auto nm = id.precursor_mass_;
          bool intersected = false;
          for (auto is:intersect)
          {
            if (abs(is - nm) > .1)
            {
              continue;
            }
            intersected = true;
            break;
          }
          if (intersected)
          {
            continue;
          }

          //if(s_map[0].empty() ||s_map[1].empty() ){ // TODO only contain overlapping proteins
          //  continue;
          //}
          //if(s_item.first<3)
          {
            outstream << s_item.first << "," << id.str_ << "\n";
          }
          // write

        }
      }
    }
    outstream.close();
    return EXECUTION_OK;
  }
};

int main(int argc, const char **argv)
{
  TOPPAnalyzeFLASHIda tool;
  return tool.main(argc, argv);
}

/// @endcond
