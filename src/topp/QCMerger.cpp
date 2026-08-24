// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Author: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/QcMLFile.h>

#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/MATH/StatisticFunctions.h>

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
    setValidFormats_("in", ListUtils::create<std::string>("qcML"));
    registerOutputFile_("out", "<file>", "", "Output extended/reduced qcML file");
    setValidFormats_("out",ListUtils::create<std::string>("qcML"));
    registerStringOption_("setname", "<string>", "", "Use only when all given qcml files belong to one set, which will be held under the given name.", false);
  }

  void addBoxPlotQPs(std::map<std::string, std::string> nums, std::map<std::string, std::string> nams, std::string set, QcMLFile& qcmlfile)
  {
    for (std::map<std::string, std::string>::const_iterator it = nums.begin(); it != nums.end(); ++it)
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
    std::string out              = getStringOption_("out");
    std::string setname          = getStringOption_("setname");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    QcMLFile qcmlfile;
    if (!setname.empty())
    {
      qcmlfile.registerSet(setname,setname,std::set<std::string>());
    }
    for (Size i = 0; i < in_files.size(); ++i)
    {
      QcMLFile tmpfile;
      tmpfile.load(in_files[i]);
      qcmlfile.merge(tmpfile,setname);
    }

    if (!setname.empty())
    {
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
