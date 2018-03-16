// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <OpenMS/FORMAT/CachedMzML.h>
#include <OpenMS/FORMAT/SqMassFile.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <fstream>

#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataCachedConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataSqlConsumer.h>

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

 typedef PeakMap MapType;

 protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in","<file>","","Input mzML file");
    registerStringOption_("in_type", "<type>", "", "input file type -- default: determined from file extension or content\n", false);
    String formats("mzML,sqMass");
    setValidFormats_("in", ListUtils::create<String>(formats));
    setValidStrings_("in_type", ListUtils::create<String>(formats));

    formats = "mzML,sqMass";
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>(formats));
    registerStringOption_("out_type", "<type>", "", "Output file type -- default: determined from file extension or content\nNote: that not all conversion paths work or make sense.", false);
    setValidStrings_("out_type", ListUtils::create<String>(formats));

    //registerStringOption_("out_meta","<file>","","output file", false);
    //setValidFormats_("out_meta",ListUtils::create<String>("mzML"));

    registerFlag_("convert_back", "Convert back to mzML");

    registerStringOption_("lossy_compression", "<type>", "true", "Use numpress compression to achieve optimally small file size (attention: may cause small loss of precision; only for mzML data).", false);
    setValidStrings_("lossy_compression", ListUtils::create<String>("true,false"));
    registerStringOption_("full_meta", "<type>", "true", "Write full meta information into sqMass file (may require large amounts of memory)", false);
    setValidStrings_("full_meta", ListUtils::create<String>("true,false"));

    registerDoubleOption_("lossy_mass_accuracy", "<error>", -1.0, "Desired (absolute) m/z accuracy for lossy compression (e.g. use 0.0001 for a mass accuracy of 0.2 ppm at 500 m/z, default uses -1.0 for maximal accuracy).", false, true);

    registerFlag_("process_lowmemory", "Whether to process the file on the fly without loading the whole file into memory first (only for conversions of mzXML/mzML to mzML).\nNote: this flag will prevent conversion from spectra to chromatograms.", true);
    registerIntOption_("lowmem_batchsize", "<number>", 500, "The batch size of the low memory conversion", false, true);
    setMinInt_("lowmem_batchsize", 0);
  }

  ExitCodes main_(int , const char**) override
  {
    String out_meta = getStringOption_("out");
    String out_cached = out_meta + ".cached";
    bool convert_back =  getFlag_("convert_back");
    bool process_lowmemory = getFlag_("process_lowmemory");
    int batchSize = (int)getIntOption_("lowmem_batchsize");

    bool full_meta = (getStringOption_("full_meta") == "true");
    bool lossy_compression = (getStringOption_("lossy_compression") == "true");
    double mass_acc = getDoubleOption_("lossy_mass_accuracy");

    FileHandler fh;

    //input file type
    String in = getStringOption_("in");
    String in_cached = in + ".cached";
    FileTypes::Type in_type = FileTypes::nameToType(getStringOption_("in_type"));

    if (in_type == FileTypes::UNKNOWN)
    {
      in_type = fh.getType(in);
      writeDebug_(String("Input file type: ") + FileTypes::typeToName(in_type), 2);
    }

    if (in_type == FileTypes::UNKNOWN)
    {
      writeLog_("Error: Could not determine input file type!");
      return PARSE_ERROR;
    }

    //output file names and types
    String out = getStringOption_("out");
    FileTypes::Type out_type = FileTypes::nameToType(getStringOption_("out_type"));

    if (out_type == FileTypes::UNKNOWN)
    {
      out_type = fh.getTypeByFileName(out);
    }

    if (out_type == FileTypes::UNKNOWN)
    {
      writeLog_("Error: Could not determine output file type!");
      return PARSE_ERROR;
    }

    if (in_type == FileTypes::SQMASS && out_type == FileTypes::MZML)
    {
      MapType exp;
      SqMassFile sqfile;
      MzMLFile f;
      sqfile.load(in, exp);
      f.store(out, exp);
      return EXECUTION_OK;
    }
    else if (in_type == FileTypes::MZML && out_type == FileTypes::SQMASS && process_lowmemory)
    {
      MSDataSqlConsumer consumer(out, batchSize, full_meta, lossy_compression, mass_acc);
      MzMLFile f;
      PeakFileOptions opt = f.getOptions();
      opt.setMaxDataPoolSize(batchSize); 
      f.setOptions(opt);
      f.transform(in, &consumer, true, true);
      return EXECUTION_OK;
    }
    else if (in_type == FileTypes::SQMASS && out_type == FileTypes::SQMASS && process_lowmemory)
    {
      PlainMSDataWritingConsumer consumer(out);
      consumer.getOptions().setWriteIndex(true);
      SqMassFile f;
      f.transform(in, &consumer, true, true);
      return EXECUTION_OK;
    }
    else if (in_type == FileTypes::MZML && out_type == FileTypes::SQMASS)
    {
      MzMLFile f;

      SqMassFile::SqMassConfig config;
      config.write_full_meta = full_meta;
      config.use_lossy_numpress = lossy_compression;
      config.linear_fp_mass_acc = mass_acc;

      SqMassFile sqfile;
      sqfile.setConfig(config);

      MapType exp;
      f.load(in, exp);
      sqfile.store(out, exp);
      return EXECUTION_OK;
    }


    if (!convert_back)
    {
      MapType exp;
      CachedmzML cacher;
      MzMLFile f;

      cacher.setLogType(log_type_);
      f.setLogType(log_type_);

      f.load(in,exp);
      cacher.writeMemdump(exp, out_cached);
      cacher.writeMetadata(exp, out_meta, true);
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
          if (meta_exp[i].getDataProcessing()[j]->metaValueExists("cached_data"))
          {
            meta_exp[i].getDataProcessing()[j]->removeMetaValue("cached_data");
          }
        }
      }

      for (Size i=0; i < meta_exp.getNrChromatograms(); ++i)
      {
        for (Size j = 0; j < meta_exp.getChromatogram(i).getDataProcessing().size(); j++)
        {
          if (meta_exp.getChromatogram(i).getDataProcessing()[j]->metaValueExists("cached_data"))
          {
            meta_exp.getChromatogram(i).getDataProcessing()[j]->removeMetaValue("cached_data");
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
      std::vector<MSChromatogram > chromatograms = exp_reading.getChromatograms();
      std::vector<MSChromatogram > old_chromatograms = meta_exp.getChromatograms();
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
