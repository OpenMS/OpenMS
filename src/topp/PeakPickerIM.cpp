// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Author: Mohammed Alhigaylan $
// $Maintainer: Timo Sachsenberg $
// -------------------------------------------------------------------------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>

using namespace OpenMS;
using namespace std;

class TOPPPeakPickerIM : public TOPPBase
{
public:
  TOPPPeakPickerIM() : TOPPBase("PeakPickerIM", "Applies PeakPickerIM to an mzML file", false) {}

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input mzML file");
    setValidFormats_("in", { "mzML" });

    registerOutputFile_("out", "<file>", "", "Output mzML file");
    setValidFormats_("out", { "mzML" });

    registerStringOption_("processOption", "<name>", "inmemory", "Whether to load all data and process them in-memory or whether to process the data on the fly (lowmemory) without loading the whole file into memory first", false, true);
    setValidStrings_("processOption", ListUtils::create<String>("inmemory,lowmemory"));

  }

    /**
    @brief Helper class for the Low Memory peak-picking
  */
  class Consumer : public MSDataWritingConsumer
  {
    public:

    Consumer(String filename, const PeakPickerIM& pp) :
      MSDataWritingConsumer(std::move(filename))
    {
      pp_ = pp;
    }

    void processSpectrum_(MapType::SpectrumType& spectrum) override
    {
      pp_.pickIMTraces(spectrum);
    }

    void processChromatogram_(MapType::ChromatogramType & c) override
    {
    }

  private:
    PeakPickerIM pp_;
  };


  ExitCodes doLowMemAlgorithm(const PeakPickerIM& pp, const String& input_file, const String& output_file)
  {
    ///////////////////////////////////
    // Create the consumer object, add data processing
    ///////////////////////////////////
    Consumer pp_consumer(output_file, pp);
    pp_consumer.addDataProcessing(getProcessingInfo_(DataProcessing::PEAK_PICKING));

    ///////////////////////////////////
    // Create new MSDataReader and set our consumer
    ///////////////////////////////////
    MzMLFile mz_data_file;
    mz_data_file.setLogType(log_type_);
    mz_data_file.transform(input_file, &pp_consumer);

    return EXECUTION_OK;
  }

  ExitCodes main_(int, const char**) override
  {
    // Get input and output file paths
    String input_file = getStringOption_("in");
    String output_file = getStringOption_("out");
    String process_option = getStringOption_("processOption");

    PeakPickerIM picker;
    if (process_option == "lowmemory")
    {
      return doLowMemAlgorithm(picker, input_file, output_file);
    }
    else
    {
      // Load input mzML file
      PeakMap exp;
      MzMLFile mzml;
      mzml.load(input_file, exp);

      // Process each spectrum with PeakPickerIM
      for (MSSpectrum& spectrum : exp)
      {
        OPENMS_LOG_DEBUG << "Processing MS" << spectrum.getMSLevel() << " spectrum with " 
          << spectrum.size() << " peaks in the IM frame." << std::endl;
        picker.pickIMTraces(spectrum);
        OPENMS_LOG_DEBUG << "Processed spectrum has " << spectrum.size() << " centroided IM peaks." << std::endl;
      }
      // Save output mzML file
      OPENMS_LOG_DEBUG << "Saving output mzML file: " << output_file << std::endl;
      mzml.store(output_file, exp);

      return EXECUTION_OK;
    }
  }
};

int main(int argc, const char** argv)
{
  TOPPPeakPickerIM tool;
  return tool.main(argc, argv);
}

