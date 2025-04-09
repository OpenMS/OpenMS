// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/HANDLERS/CachedMzMLHandler.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/FORMAT/CachedMzML.h>
#include <OpenMS/FORMAT/FileTypes.h>
//TODO add support for indexed mzml to handler
#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/SYSTEM/SysInfo.h>

#include <numeric>

using namespace OpenMS;
using namespace std;
using namespace OpenSwath;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_TICCalculator TICCalculator

@brief Calculates the TIC of a raw mass spectrometric file. 

This class was developed to benchmark multiple methods inside OpenMS for
reading raw mass spectrometric data. It can be used for benchmarking these
different methods as well as benchmarking external tools. Of course you can
also calculate the TIC with this tool.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_TICCalculator.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_TICCalculator.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TICConsumer : 
    public Interfaces::IMSDataConsumer
{

    typedef PeakMap MapType;
    typedef MapType::SpectrumType SpectrumType;
    typedef MapType::ChromatogramType ChromatogramType;

public:
  double TIC;
  int nr_spectra;
  long int nr_peaks;

  // Create new consumer, set TIC to zero
  TICConsumer() :
    TIC(0.0),
    nr_spectra(0.0),
    nr_peaks(0)
    {}

  void consumeSpectrum(SpectrumType & s) override
  {
    for (Size i = 0; i < s.size(); i++) 
    { 
      TIC += s[i].getIntensity(); 
    }
    nr_peaks += s.size();
    nr_spectra++;
  }

  void consumeChromatogram(ChromatogramType& /* c */) override {}
  void setExpectedSize(Size /* expectedSpectra */, Size /* expectedChromatograms */) override {}
  void setExperimentalSettings(const ExperimentalSettings& /* exp */) override {}
};

/**
  @brief Abstraction of a std::ifstream 

  Useful for parallel access to the file when each thread is given its own
  instance of this class. Each thread will then have its own file stream and
  access the file independently.

*/
class FileAbstraction
{

public:
  // constructor
  explicit FileAbstraction(std::string filename) :
    filename_(filename)
  {
    ifs_.open(filename.c_str(), std::ios::binary);
  }

  // copy constructor
  FileAbstraction(const FileAbstraction& source) :
    filename_(source.filename_),
    ifs_(source.filename_.c_str())
  {}

  // access to underlying stream
  std::ifstream & getStream()
  {
    return ifs_;
  }

private:
  std::string filename_;
  std::ifstream ifs_;
};

class TOPPTICCalculator :
  public TOPPBase
{
public:
  TOPPTICCalculator() :
    TOPPBase("TICCalculator", 
    "Calculates the TIC from a mass spectrometric raw file (useful for benchmarking).", 
    true)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file to convert.");
    registerStringOption_("in_type", "<type>", "", "Input file type -- default: determined from file extension or content\n", false);
    String formats("mzData,mzXML,mzML,cachedMzML,dta,dta2d,mgf,featureXML,consensusXML,ms2,fid,tsv,peplist,kroenik,edta");
    setValidFormats_("in", ListUtils::create<String>(formats));
    setValidStrings_("in_type", ListUtils::create<String>(formats));
    
    registerStringOption_("read_method", "<method>", "regular", "Method to read the file", false);
    String method("regular,indexed,indexed_parallel,streaming,cached,cached_parallel");
    setValidStrings_("read_method", ListUtils::create<String>(method));

    registerStringOption_("loadData", "<method>", "true", "Whether to actually load and decode the binary data (or whether to skip decoding the binary data)", false);
    String loadData("true,false");
    setValidStrings_("loadData", ListUtils::create<String>(loadData));
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    //input file names
    String in = getStringOption_("in");
    String read_method = getStringOption_("read_method");
    bool load_data = getStringOption_("loadData") == "true";

    if (read_method == "streaming")
    {
      std::cout << "Read method: streaming" << std::endl;

      // Create the consumer, set output file name, transform
      TICConsumer consumer;
      MzMLFile mzml;
      mzml.setLogType(log_type_);

      PeakFileOptions opt = mzml.getOptions();
      opt.setFillData(load_data); // whether to actually load any data
      opt.setSkipXMLChecks(true); // save time by not checking base64 strings for whitespaces 
      opt.setMaxDataPoolSize(100);
      opt.setAlwaysAppendData(false);
      mzml.setOptions(opt);
      mzml.transform(in, &consumer, true, true);

      std::cout << "There are " << consumer.nr_spectra << " spectra and " << consumer.nr_peaks << " peaks in the input file." << std::endl;
      std::cout << "The total ion current is " << consumer.TIC << std::endl;
      size_t after;
      SysInfo::getProcessMemoryConsumption(after);
      std::cout << " Memory consumption after " << after << std::endl;
    }
    else if (read_method == "regular")
    {
      std::cout << "Read method: regular" << std::endl;

      MzMLFile mzml;
      mzml.setLogType(log_type_);
      PeakFileOptions opt = mzml.getOptions();
      opt.setFillData(load_data); // whether to actually load any data
      opt.setSkipXMLChecks(true); // save time by not checking base64 strings for whitespaces 
      mzml.setOptions(opt);
      PeakMap map;
      mzml.load(in, map);
      double TIC = 0.0;
      long int nr_peaks = 0;
      for (Size i =0; i < map.size(); i++)
      {
        nr_peaks += map[i].size();
        for (Size j = 0; j < map[i].size(); j++)
        {
          TIC += map[i][j].getIntensity();
        }
      }

      std::cout << "There are " << map.size() << " spectra and " << nr_peaks << " peaks in the input file." << std::endl;
      std::cout << "The total ion current is " << TIC << std::endl;
      size_t after;
      SysInfo::getProcessMemoryConsumption(after);
      std::cout << " Memory consumption after " << after << std::endl;
    }
    else if (read_method == "indexed")
    {
      std::cout << "Read method: indexed" << std::endl;
      
      IndexedMzMLFileLoader imzml;
      // load data from an indexed MzML file
      OnDiscPeakMap map;
      imzml.load(in, map);
      double TIC = 0.0;
      long int nr_peaks = 0;
      if (load_data)
      {
        for (Size i =0; i < map.getNrSpectra(); i++)
        {
          OpenMS::Interfaces::SpectrumPtr sptr = map.getSpectrumById(i);

          nr_peaks += sptr->getIntensityArray()->data.size();
          TIC += std::accumulate(sptr->getIntensityArray()->data.begin(), sptr->getIntensityArray()->data.end(), 0.0);
        }
      }

      std::cout << "There are " << map.getNrSpectra() << " spectra and " << nr_peaks << " peaks in the input file." << std::endl;
      std::cout << "The total ion current is " << TIC << std::endl;
      size_t after;
      SysInfo::getProcessMemoryConsumption(after);
      std::cout << " Memory consumption after " << after << std::endl;
    }
    else if (read_method == "indexed_parallel")
    {
      std::cout << "Read method: indexed (parallel)" << std::endl;
      
      IndexedMzMLFileLoader imzml;
      PeakFileOptions opt = imzml.getOptions();
      opt.setFillData(load_data); // whether to actually load any data
      imzml.setOptions(opt);

      // load data from an indexed MzML file
      OnDiscPeakMap map;
      map.openFile(in, true);
      map.setSkipXMLChecks(true);

      double TIC = 0.0;
      long nr_peaks = 0;

      if (load_data)
      {
        // firstprivate means that each thread has its own instance of the
        // variable, each copy initialized with the initial value 
#pragma omp parallel for firstprivate(map) reduction(+: TIC, nr_peaks)
        for (SignedSize i =0; i < (SignedSize)map.getNrSpectra(); i++)
        {
          OpenMS::Interfaces::SpectrumPtr sptr = map.getSpectrumById(i);
          nr_peaks += sptr->getIntensityArray()->data.size();
          TIC += std::accumulate(sptr->getIntensityArray()->data.begin(), sptr->getIntensityArray()->data.end(), 0.0);
        }
      }

      std::cout << "There are " << map.getNrSpectra() << " spectra and " << nr_peaks << " peaks in the input file." << std::endl;
      std::cout << "The total ion current is " << TIC << std::endl;
      size_t after;
      SysInfo::getProcessMemoryConsumption(after);
      std::cout << " Memory consumption after " << after << std::endl;
    }
    else if (read_method == "cached")
    {
      std::cout << "Read method: cached" << std::endl;

      // Special handling of cached mzML as input types: 
      // we expect two paired input files which we should read into exp
      std::vector<String> split_out;
      in.split(".cachedMzML", split_out);
      if (split_out.size() != 2)
      {
        OPENMS_LOG_ERROR << "Cannot deduce base path from input '" << in << 
          "' (note that '.cachedMzML' should only occur once as the final ending)" << std::endl;
        return ILLEGAL_PARAMETERS;
      }
      String in_meta = split_out[0] + ".mzML";

      MzMLFile f;
      f.setLogType(log_type_);

      Internal::CachedMzMLHandler cache;
      cache.createMemdumpIndex(in);
      const std::vector<std::streampos> spectra_index = cache.getSpectraIndex();

      std::ifstream ifs_;
      ifs_.open(in.c_str(), std::ios::binary);

      double TIC = 0.0;
      long int nr_peaks = 0;
      for (Size i=0; i < spectra_index.size(); ++i)
      {

        BinaryDataArrayPtr mz_array(new BinaryDataArray);
        BinaryDataArrayPtr intensity_array(new BinaryDataArray);
        int ms_level = -1;
        double rt = -1.0;
        ifs_.seekg(spectra_index[i]);
        Internal::CachedMzMLHandler::readSpectrumFast(mz_array, intensity_array, ifs_, ms_level, rt);

        nr_peaks += intensity_array->data.size();
        for (Size j = 0; j < intensity_array->data.size(); j++)
        {
          TIC += intensity_array->data[j];
        }
      }

      std::cout << "There are " << spectra_index.size() << " spectra and " << nr_peaks << " peaks in the input file." << std::endl;
      std::cout << "The total ion current is " << TIC << std::endl;
      size_t after;
      SysInfo::getProcessMemoryConsumption(after);
      std::cout << " Memory consumption after " << after << std::endl;
    }
    else if (read_method == "cached_parallel")
    {
      std::cout << "Read method: cached parallel" << std::endl;

      // Special handling of cached mzML as input types: 
      // we expect two paired input files which we should read into exp
      std::vector<String> split_out;
      in.split(".cachedMzML", split_out);
      if (split_out.size() != 2)
      {
        OPENMS_LOG_ERROR << "Cannot deduce base path from input '" << in << 
          "' (note that '.cachedMzML' should only occur once as the final ending)" << std::endl;
        return ILLEGAL_PARAMETERS;
      }
      String in_meta = split_out[0] + ".mzML";

      MzMLFile f;
      f.setLogType(log_type_);

      Internal::CachedMzMLHandler cache;
      cache.createMemdumpIndex(in);
      const std::vector<std::streampos> spectra_index = cache.getSpectraIndex();

      FileAbstraction filestream(in);

      double TIC = 0.0;
      long nr_peaks = 0;
#pragma omp parallel for firstprivate(filestream) reduction(+: TIC, nr_peaks)
      for (SignedSize i=0; i < (SignedSize)spectra_index.size(); ++i)
      {
        BinaryDataArrayPtr mz_array(new BinaryDataArray);
        BinaryDataArrayPtr intensity_array(new BinaryDataArray);
        int ms_level = -1;
        double rt = -1.0;
        // we only change the position of the thread-local filestream
        filestream.getStream().seekg(spectra_index[i]);
        Internal::CachedMzMLHandler::readSpectrumFast(mz_array, intensity_array, filestream.getStream(), ms_level, rt);

        nr_peaks += intensity_array->data.size();
        TIC += std::accumulate(intensity_array->data.begin(), intensity_array->data.end(), 0.0);
      }

      std::cout << "There are " << spectra_index.size() << " spectra and " << nr_peaks << " peaks in the input file." << std::endl;
      std::cout << "The total ion current is " << TIC << std::endl;
      size_t after;
      SysInfo::getProcessMemoryConsumption(after);
      std::cout << " Memory consumption after " << after << std::endl;
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{

  // Print usage if used without arguments
  if (argc == 1)
  {
    TOPPTICCalculator tool;
    tool.main(argc, argv);
    return 0;
  }

  // Add -test at the end of the arguments in order to avoid calling the OpenMS
  // server for usage statistics (and thus making the benchmark slower)
  char testflag[] = "-test";
  std::vector<const char *> newArgs(argc+1); // vector containing one element more than required
  for (int arg = 0; arg < argc; ++arg)
  {
    newArgs[arg] = argv[arg];
  }
  newArgs[argc] = testflag;

  TOPPTICCalculator tool;
  size_t after, before;
  SysInfo::getProcessMemoryConsumption(before);
  std::cout << " Memory consumption before " << before << std::endl;
  tool.main(argc+1, (const char **)&newArgs[0]);
  SysInfo::getProcessMemoryConsumption(after);
  std::cout << " Memory consumption final " << after << std::endl;

  return 0;
}

/// @endcond
