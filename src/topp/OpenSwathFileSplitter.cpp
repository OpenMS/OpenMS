// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

// Files
#include <OpenMS/ANALYSIS/OPENSWATH/SwathQC.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/SYSTEM/File.h>
#include <QDir>

using namespace OpenMS;

// OpenMS base classes
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_OpenSwathFileSplitter OpenSwathFileSplitter

@brief A tool for splitting a single SWATH / DIA file into a set of files, each containing one SWATH window (plus one file for the MS1 survey scans).

Will use the input SWATH / DIA file to generate one output file containing
the MS1 survey scans and \a n individual files for each SWATH / DIA window in
the input file. The number of windows is read from the input file itself.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_OpenSwathFileSplitter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_OpenSwathFileSplitter.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPOpenSwathFileSplitter : public TOPPBase
{
public:
  TOPPOpenSwathFileSplitter() : TOPPBase("OpenSwathFileSplitter", "Splits SWATH files into n files, each containing one window.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<files>", "", "Input file (SWATH/DIA file)");
    setValidFormats_("in", ListUtils::create<String>("mzML,mzXML"));
    registerOutputPrefix_("outputDirectory", "<output>", "./", "Output file prefix", false, true);
    // additional QC data
    registerOutputFile_("out_qc", "<file>", "", "Optional QC meta data (charge distribution in MS1). Only works with mzML input files.", false, true);
    setValidFormats_("out_qc", ListUtils::create<String>("json"));
  }

  void loadSwathFiles(const String& file_in, const String& tmp, const String& readoptions, boost::shared_ptr<ExperimentalSettings>& exp_meta, std::vector<OpenSwath::SwathMap>& swath_maps,
                      Interfaces::IMSDataConsumer* plugin_consumer = nullptr)
  {
    SwathFile swath_file;
    swath_file.setLogType(log_type_);

    FileTypes::Type in_file_type = FileHandler::getTypeByFileName(file_in);
    if (in_file_type == FileTypes::MZML)
    {
      swath_maps = swath_file.loadMzML(file_in, tmp, exp_meta, readoptions, plugin_consumer);
    }
    else if (in_file_type == FileTypes::MZXML)
    {
      swath_maps = swath_file.loadMzXML(file_in, tmp, exp_meta, readoptions);
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Input file needs to have ending .mzML(.gz) or .mzXML(.gz)");
    }
  }

  ExitCodes main_(int, const char**) override
  {
    ///////////////////////////////////
    // Prepare Parameters
    ///////////////////////////////////
    String file_in = getStringOption_("in");

    // make sure tmp is a directory with proper separator at the end (downstream methods simply do path + filename)
    // (do not use QDir::separator(), since its platform specific (/ or \) while absolutePath() will always use '/')
    String tmp_dir = String(QDir(getStringOption_("outputDirectory").c_str()).absolutePath()).ensureLastChar('/');

    QFileInfo fi(file_in.toQString());
    String tmp = tmp_dir + String(fi.baseName());

    String out_qc = getStringOption_("out_qc");

    ///////////////////////////////////
    // Load the SWATH files
    ///////////////////////////////////
    boost::shared_ptr<ExperimentalSettings> exp_meta(new ExperimentalSettings);
    std::vector<OpenSwath::SwathMap> swath_maps;

    // collect some QC data
    if (out_qc.empty())
    {
      loadSwathFiles(file_in, tmp, "split", exp_meta, swath_maps);
    }
    else
    {
      OpenSwath::SwathQC qc(30, 0.04);
      MSDataTransformingConsumer qc_consumer; // apply some transformation
      qc_consumer.setSpectraProcessingFunc(qc.getSpectraProcessingFunc());
      qc_consumer.setExperimentalSettingsFunc(qc.getExpSettingsFunc());
      loadSwathFiles(file_in, tmp, "split", exp_meta, swath_maps, &qc_consumer);
      qc.storeJSON(out_qc);
    }

    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPOpenSwathFileSplitter tool;
  return tool.main(argc, argv);
}

/// @endcond
