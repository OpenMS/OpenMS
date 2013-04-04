// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#include "OpenMS/ANALYSIS/OPENSWATH/CachedmzML.h"

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <fstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_OpenSwathMzMLFileCacher OpenSwathMzMLFileCacher

  @brief Serialize a spectra and/or chromatogram mzML file

  This class will serialize a spectra and/or chromatogram mzML file and store
  it in a binary format that contains ONLY the spectra and chromatogram data
  (no metadata).
 
  This is implemented using the write_memdump and read_memdump functions.
  For reading there are 2 options
  - read the whole file into the OpenMS datastructures
  - read only an index (read_memdump_idx) of the spectra and chromatograms and then use
    random-access to retrieve a specific spectra from the disk (read_memdump_spectra)

  @note This tool is experimental!

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_OpenSwathMzMLFileCacher.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_OpenSwathMzMLFileCacher.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPOpenSwathMzMLFileCacher
  : public TOPPBase,
    public ProgressLogger
{
 public:

  TOPPOpenSwathMzMLFileCacher()
    : TOPPBase("OpenSwathMzMLFileCacher","This tool caches the spectra and chromatogram data of an mzML to disk.", false)
  {
  }

 typedef MSExperiment<Peak1D> MapType;

 protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in","<file>","","transition file ('csv')");
    setValidFormats_("in", StringList::create("mzML"));

    registerOutputFile_("out","<file>","","output file");

    //registerStringOption_("out_meta","<file>","","output file", false);
    //setValidFormats_("out_meta",StringList::create("mzML"));

    registerFlag_("convert_back", "Convert back to mzML");

  }

  ExitCodes main_(int , const char**)
  {
    String in = getStringOption_("in");
    String out_meta = getStringOption_("out");
    String in_cached = in + ".cached";
    String out_cached = out_meta + ".cached";
    bool convert_back =  getFlag_("convert_back");

    if (!convert_back)
    {
      MapType exp;
      CachedmzML cacher;
      MzMLFile f;

      cacher.setLogType(log_type_);
      f.setLogType(log_type_);

      f.load(in,exp);
      cacher.writeMemdump(exp, out_cached);

      DataProcessing dp;
      std::set<DataProcessing::ProcessingAction> actions;
      actions.insert(DataProcessing::FORMAT_CONVERSION);
      dp.setProcessingActions(actions);
      dp.setMetaValue("cached_data", "true");
      for (Size i=0; i<exp.size(); ++i)
      {
        exp[i].getDataProcessing().push_back(dp);
      }
      std::vector<MSChromatogram<ChromatogramPeak> > chromatograms = exp.getChromatograms();
      for (Size i=0; i<chromatograms.size(); ++i)
      {
        chromatograms[i].getDataProcessing().push_back(dp);
      }
      exp.setChromatograms(chromatograms);
      cacher.writeMetadata(exp, out_meta);
    }
    else
    {
      MzMLFile f;
      MapType meta_exp;
      CachedmzML cacher;
      MapType exp_reading;

      cacher.setLogType(log_type_);
      f.setLogType(log_type_);

      f.load(in,meta_exp);
      cacher.readMemdump(exp_reading, in_cached);

      std::cout << " read back, got " << exp_reading.size() << " spectra " << exp_reading.getChromatograms().size() << " chromats " << std::endl;

      {
      for (Size i=0; i<meta_exp.size(); ++i)
      {
        for (Size j = 0; j < meta_exp[i].getDataProcessing().size(); j++)
        {
          DataProcessing& dp = meta_exp[i].getDataProcessing()[j];
          if (dp.metaValueExists("cached_data"))
          {
            dp.removeMetaValue("cached_data");
          }
        }
      }

      std::vector<MSChromatogram<ChromatogramPeak> > chromatograms = meta_exp.getChromatograms();
      for (Size i=0; i<chromatograms.size(); ++i)
      {
        for (Size j = 0; j < chromatograms[i].getDataProcessing().size(); j++)
        {
          DataProcessing& dp = chromatograms[i].getDataProcessing()[j];
          if (dp.metaValueExists("cached_data"))
          {
            dp.removeMetaValue("cached_data");
          }
        }
      }
      }

      if (meta_exp.size() != exp_reading.size())
      {
        std::cerr << " Both experiments need to have the same size!";
      }

      for (Size i=0; i<exp_reading.size(); ++i)
      {
        for (Size j = 0; j < exp_reading[i].size(); j++)
        {
          meta_exp[i].push_back(exp_reading[i][j]);
        }
      }
      std::vector<MSChromatogram<ChromatogramPeak> > chromatograms = exp_reading.getChromatograms();
      std::vector<MSChromatogram<ChromatogramPeak> > old_chromatograms = meta_exp.getChromatograms();
      for (Size i=0; i<chromatograms.size(); ++i)
      {
        for (Size j = 0; j < chromatograms[i].size(); j++)
        {
          old_chromatograms[i].push_back(chromatograms[i][j]);
        }
      }
      meta_exp.setChromatograms(old_chromatograms);


      f.store(out_meta,meta_exp);
    }

    return EXECUTION_OK;
  }
};

int main( int argc, const char** argv )
{

  TOPPOpenSwathMzMLFileCacher tool;
  return tool.main(argc,argv);
}

/// @endcond
