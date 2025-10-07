// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Author: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/QcMLFile.h>

#include <QByteArray>
#include <QFile>
#include <QString>
#include <QFileInfo>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/MATH/StatisticFunctions.h>

//~ #include <QIODevice>
#include <algorithm>
#include <fstream>
#include <vector>
#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_QCMerger QCMerger

@brief Merges two qcml files together.

<CENTER>
  <table>
    <tr>
    <th ALIGN = "center"> pot. predecessor tools </td>
    <td VALIGN="middle" ROWSPAN=3> &rarr; QCCalculator &rarr;</td>
    <th ALIGN = "center"> pot. successor tools </td>
    </tr>
    <tr>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_QCCalculator </td>
    </tr>
    <tr>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_QCShrinker </td>
    </tr>
  </table>
</CENTER>

The two or more given files (see parameter @p in) are merged. If a run/set exisits in several files, the quality parameters of these are merged as well.
Several runs from qcml files can be comprised in a set.

- @p setname If the runs of the given input files are to be comprised in a set, this will be the name of the set.

Output is in qcML format (see parameter @p out) which can be viewed directly in a modern browser (chromium, firefox, safari).

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_QCMerger.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_QCMerger.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPQCMerger :
  public TOPPBase
{
public:
  TOPPQCMerger() :
    TOPPBase("QCMerger", "Merges two qcml files together.", 
    true, {{ "Walzer M, Pernas LE, Nasso S, Bittremieux W, Nahnsen S, Kelchtermans P,  Martens, L", "qcML: An Exchange Format for Quality Control Metrics from Mass Spectrometry Experiments", "Molecular & Cellular Proteomics 2014; 13(8)" , "10.1074/mcp.M113.035907"}})
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "List of qcml files to be merged.");
    setValidFormats_("in", ListUtils::create<String>("qcML"));
    registerOutputFile_("out", "<file>", "", "Output extended/reduced qcML file");
    setValidFormats_("out",ListUtils::create<String>("qcML"));
    registerStringOption_("setname", "<string>", "", "Use only when all given qcml files belong to one set, which will be held under the given name.", false);
  }

  void addBoxPlotQPs(std::map<String,String> nums, std::map<String,String> nams, String set, QcMLFile& qcmlfile)
  {
    for (std::map<String, String >::const_iterator it = nums.begin(); it != nums.end(); ++it)
    {
      QcMLFile::QualityParameter qp;
      qp.name = nams[it->first]; ///< Name
      qp.id = set + it->first; ///< Identifier
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = it->first;
      qp.value = it->second;
      qcmlfile.addSetQualityParameter(set, qp);
    }
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    StringList in_files     = getStringList_("in");
    String out              = getStringOption_("out");
    String setname          = getStringOption_("setname");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    QcMLFile qcmlfile;
    if (!setname.empty())
    {
      qcmlfile.registerSet(setname,setname,std::set< String >());
    }
    for (Size i = 0; i < in_files.size(); ++i)
    {
      QcMLFile tmpfile;
      tmpfile.load(in_files[i]);
      qcmlfile.merge(tmpfile,setname);
    }

    if (!setname.empty())
    {
//        // make #ms2 set stats
//        std::vector<String> ms2nums_strings;
//        qcmlfile.collectSetParameter(setname,"QC:0000007", ms2nums_strings);
//        std::vector<Int> ms2nums;
//        for (std::vector<String>::iterator it = ms2nums_strings.begin(); it != ms2nums_strings.end(); ++it) //transform is too ugly and errorprone
//        {
//          ms2nums.push_back(it->toInt());
//        }

//        std::sort(ms2nums.begin(), ms2nums.end());

//        if (ms2nums.size()>0)
//        {
//          std::map<String,String> nums;
//          std::map<String,String> nams;
//          //~ min,q1,q2,q3,max
//          nums["QC:0000043"] = String(ms2nums.front());
//          nams["QC:0000043"] = "min ms2 number";
//          nums["QC:0000044"] = String(OpenMS::Math::quantile1st(ms2nums.begin(), ms2nums.end(),true));
//          nams["QC:0000044"] = "Q1 ms2 number";
//          nums["QC:0000045"] = String(OpenMS::Math::median(ms2nums.begin(), ms2nums.end(), true));
//          nams["QC:0000045"] = "Q2 ms2 number";
//          nums["QC:0000046"] = String(OpenMS::Math::quantile3rd(ms2nums.begin(), ms2nums.end(),true));
//          nams["QC:0000046"] = "Q3 ms2 number";
//          nums["QC:0000047"] = String(ms2nums.back());
//          nams["QC:0000047"] = "max ms2 number";

//          addBoxPlotQPs(nums, nams, setname, qcmlfile);
//        }

//        // make #id-psm set stats
//        std::vector<String> idnums_strings;
//        qcmlfile.collectSetParameter(setname,"QC:0000029", idnums_strings);
//        std::vector<Int> idnums;
//        for (std::vector<String>::iterator it = idnums_strings.begin(); it != idnums_strings.end(); ++it) //transform is too ugly and errorprone
//        {
//          idnums.push_back(it->toInt());
//        }

//        std::sort(idnums.begin(), idnums.end());

//        if (idnums.size()>0)
//        {
//          std::map<String,String> nums;
//          std::map<String,String> nams;
//          //~ min,q1,q2,q3,max

//          nums["QC:0000053"] = String(idnums.front());
//          nams["QC:0000053"] = "min id numbers";
//          nums["QC:0000054"] = String(OpenMS::Math::quantile1st(idnums.begin(), idnums.end()));
//          nams["QC:0000054"] = "Q1 id numbers";
//          nums["QC:0000055"] = String(OpenMS::Math::median(idnums.begin(), idnums.end()));
//          nams["QC:0000055"] = "Q2 id numbers";
//          nums["QC:0000056"] = String(OpenMS::Math::quantile3rd(idnums.begin(), idnums.end()));
//          nams["QC:0000056"] = "Q3 id numbers";
//          nums["QC:0000057"] = String(idnums.back());
//          nams["QC:0000057"] = "max id number";

//          addBoxPlotQPs(nums, nams, setname, qcmlfile);
//        }
    }

    qcmlfile.store(out);
    return EXECUTION_OK;
  }

};
int main(int argc, const char ** argv)
{
  TOPPQCMerger tool;
  return tool.main(argc, argv);
}

/// @endcond
