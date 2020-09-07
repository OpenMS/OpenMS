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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

// Files
#include <OpenMS/ANALYSIS/OPENSWATH/SwathQC.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>

#include <OpenMS/METADATA/ExperimentalSettings.h>

#include <OpenMS/SYSTEM/File.h>


#include <QDir>

using namespace OpenMS;

// OpenMS base classes
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_OpenSwathFileSplitter OpenSwathFileSplitter

  @brief A tool for splitting a single SWATH / DIA file into a set of files, each containing one SWATH window (plus one file for the MS1 survey scans).

  Will use the input SWATH / DIA file to generate one output file containing
  the MS1 survey scans and \a n individual files for each SWATH / DIA window in
  the input file. The number of windows is read from the input file itself.

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_OpenSwathFileSplitter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_OpenSwathFileSplitter.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPOpenSwathFileSplitter
  : public TOPPBase
{
public:

  TOPPOpenSwathFileSplitter()
    : TOPPBase("OpenSwathFileSplitter", "Splits SWATH files into n files, each containing one window.", false)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<files>", "", "Input file (SWATH/DIA file)");
    setValidFormats_("in", ListUtils::create<String>("mzML,mzXML"));

    registerStringOption_("outputDirectory", "<output>", "./", "Output path to store the split files", false, true);
    
    // additional QC data
    registerOutputFile_("out_qc", "<file>", "", "Optional QC meta data (charge distribution in MS1). Only works with mzML input files.", false, true);
    setValidFormats_("out_qc", ListUtils::create<String>("json"));
  }

  void loadSwathFiles(const String& file_in,
                      const String& tmp,
                      const String& readoptions,
                      boost::shared_ptr<ExperimentalSettings >& exp_meta,
                      std::vector< OpenSwath::SwathMap >& swath_maps,
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
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Input file needs to have ending .mzML(.gz) or .mzXML(.gz)");
    }
  }

  ExitCodes main_(int, const char **) override
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
    std::vector< OpenSwath::SwathMap > swath_maps;

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

int main(int argc, const char ** argv)
{
  TOPPOpenSwathFileSplitter tool;
  return tool.main(argc, argv);
}

/// @endcond
