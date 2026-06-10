// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker, Chris Bielow, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/APPLICATIONS/TOPPASWorkflow.h>
#include <OpenMS/APPLICATIONS/PipelineExecutor.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <iostream>
#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_ExecutePipeline ExecutePipeline

@brief Executes workflows created by TOPPAS.

This tool is the non-GUI, i.e. command line version for non-interactive execution of TOPPAS pipelines.
In order to really use this tool in batch-mode, you can provide a TOPPAS resource file (.trf) which specifies the
input files for the input nodes in your pipeline.

<B> *.trf files </B>

A TOPPAS resource file (<TT>*.trf</TT>) specifies the locations of input files for a pipeline.
It is an XML file following the normal TOPP INI file schema, i.e. it can be edited using the INIFileEditor or filled using a script (we do NOT provide one - sorry).
It can be exported from TOPPAS (<TT>File -> Save TOPPAS resource file</TT>). For two input nodes 1 and 2 with files (<TT>dataA.mzML</TT>, <TT>dataB.mzML</TT>) and (<TT>dataC.mzML</TT>) respectively it has the following format.

\code
<?xml version="1.0" encoding="ISO-8859-1"?>
<PARAMETERS version="1.3" xsi:noNamespaceSchemaLocation="http://open-ms.sourceforge.net/schemas/Param_1_3.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <NODE name="1" description="">
    <ITEMLIST name="url_list" type="string" description="">
      <LISTITEM value="file:///Users/jeff/dataA.mzML"/>
      <LISTITEM value="file:///Users/jeff/dataB.mzML"/>
    </ITEMLIST>
  </NODE>
  <NODE name="2" description="">
    <ITEMLIST name="url_list" type="string" description="">
      <LISTITEM value="file:///Users/jeff/dataC.mzML"/>
    </ITEMLIST>
  </NODE>
</PARAMETERS>
\endcode

The node names refer to the topological numbers of the input nodes in the workflow.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_ExecutePipeline.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_ExecutePipeline.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPExecutePipeline :
  public TOPPBase
{
public:
  TOPPExecutePipeline() :
    TOPPBase("ExecutePipeline",
             "Executes workflows created by TOPPAS.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "The workflow to be executed.");
    setValidFormats_("in", ListUtils::create<std::string>("toppas"));
    registerStringOption_("out_dir", "<directory>", "", "Directory for output files (default: user's home directory)", false);
    registerStringOption_("resource_file", "<file>", "", "A TOPPAS resource file (*.trf) specifying the files this workflow is to be applied to", false);
    registerIntOption_("num_jobs", "<integer>", 1, "Maximum number of jobs running in parallel", false, false);
    setMinInt_("num_jobs", 1);
  }

  ExitCodes main_(int, const char**) override
  {
    std::string toppas_file = getStringOption_("in");
    std::string out_dir = getStringOption_("out_dir");
    std::string resource_file = getStringOption_("resource_file");
    int num_jobs = getIntOption_("num_jobs");

    if (num_jobs > 1)
    {
      writeLogWarn_("Note: parallel execution is not implemented yet; running the pipeline serially (num_jobs is ignored).");
    }

    // create a fresh temporary working directory (deleted at the end)
    std::string tmp_path = File::getTempDirectory() + "/" + File::getUniqueName();
    if (!File::makeDir(tmp_path))
    {
      writeLogError_("Could not create temporary directory '" + tmp_path + "'.");
      return CANNOT_WRITE_OUTPUT_FILE;
    }

    // determine and prepare the output directory
    if (out_dir.empty())
    {
      out_dir = File::getUserDirectory() + "/" + File::stemName(toppas_file);
      std::cout << "No output directory specified. Using '" << out_dir << "'." << std::endl;
    }
    else
    {
      out_dir = File::absolutePath(out_dir);
    }
    if (!((File::exists(out_dir) && File::isDirectory(out_dir)) || File::makeDir(out_dir)))
    {
      writeLogError_("Cannot create or access output directory '" + out_dir + "'.");
      return CANNOT_WRITE_OUTPUT_FILE;
    }

    // load the workflow and build the executable pipeline
    TOPPASWorkflow wf;
    wf.load(toppas_file);

    PipelineExecutor executor;
    executor.fromWorkflow(wf, toppas_file);

    std::map<std::string, std::vector<std::string>> overrides;
    if (!resource_file.empty())
    {
      overrides = PipelineExecutor::loadResourceFile(resource_file);
    }

    PipelineExecutor::Options opts;
    opts.out_dir = out_dir;
    opts.tmp_dir = tmp_path;
    opts.num_jobs = num_jobs;

    std::string error_msg;
    PipelineExecutor::Result result = executor.execute(opts, overrides, error_msg);

    // clean up the temporary files (safety: only if the path is below the system temp directory)
    if (StringUtils::hasPrefix(File::absolutePath(tmp_path), File::getTempDirectory()))
    {
      File::removeDirRecursively(tmp_path);
    }

    switch (result)
    {
      case PipelineExecutor::OK:
        return EXECUTION_OK;
      case PipelineExecutor::TOOL_NOT_FOUND:
        writeLogError_(error_msg);
        return EXTERNAL_PROGRAM_NOTFOUND;
      default:
        writeLogError_(error_msg);
        return UNKNOWN_ERROR;
    }
  }

};


int main(int argc, const char ** argv)
{
  TOPPExecutePipeline tool;
  return tool.main(argc, argv);
}

/// @endcond
