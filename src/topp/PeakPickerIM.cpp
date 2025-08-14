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
    setValidStrings_("processOption", { "inmemory", "lowmemory" } );

    registerStringOption_("method", "<name>", "", "Method to pick peaks in IM dimension", false, true);
    setValidStrings_("method", { "mobilogram", "cluster", "traces" } );

    // --- pickIMTraces parameters ---
    registerDoubleOption_("pickIMTraces:sum_tolerance_mz", "<num>", 0.1,
                          "Tolerance for summing adjacent m/z peaks (ppm)", false);
    registerDoubleOption_("pickIMTraces:gauss_ppm_tolerance", "<num>", 5.0,
                          "Gaussian smoothing m/z tolerance in ppm", false);
    registerDoubleOption_("pickIMTraces:sum_tolerance_im", "<num>", 0.0006,
                          "Tolerance for summing adjacent ion mobility peaks (1/k0)", false);
    registerIntOption_("pickIMTraces:sgolay_frame_length", "<num>", 5,
                       "Savitzky-Golay smoothing frame length", false);
    registerIntOption_("pickIMTraces:sgolay_polynomial_order", "<num>", 3,
                       "Savitzky-Golay smoothing polynomial order", false);
    registerIntOption_("pickIMTraces:padding_points", "<num>", 15,
                       "Padding points on each side for mobilograms", false);

    // --- pickIMCluster parameters ---
    registerDoubleOption_("pickIMCluster:ppm_tolerance", "<num>", 10.0,
                          "m/z tolerance in ppm for clustering", false);
    registerDoubleOption_("pickIMCluster:im_tolerance", "<num>", 0.1,
                          "Ion mobility tolerance in 1/k for clustering", false);

    // --- pickIMElutionProfiles  ---
    registerDoubleOption_("pickIMElutionProfiles:ppm_tolerance", "<ppm>", 10.0,
                          "Mass trace m/z tolerance in ppm.", false);
  }

    /**
    @brief Helper class for the Low Memory peak-picking
  */
  class Consumer : public MSDataWritingConsumer
  {
    public:

    Consumer(String filename, const String& method, const PeakPickerIM& pp) :
      MSDataWritingConsumer(std::move(filename)), pp_(pp), method_(method) {}

    void processSpectrum_(MapType::SpectrumType& spectrum) override
    {
      if (method_ == "mobilogram")
      {
        pp_.pickIMTraces(spectrum);
      }
      else if (method_ == "cluster")
      {
        PeakPickerIM::pickIMCluster(spectrum, 50.0, 0.08);
      }
      else if (method_ == "traces")
      {
        PeakPickerIM::pickIMElutionProfiles(spectrum, 50.0);
      }
    }

    void processChromatogram_(MapType::ChromatogramType & ) override
    {
    }

  private:
    PeakPickerIM pp_;
    String method_ = "mobilogram";
  };


  ExitCodes doLowMemAlgorithm(const String method, const PeakPickerIM& pp, const String& input_file, const String& output_file)
  {
    ///////////////////////////////////
    // Create the consumer object, add data processing
    ///////////////////////////////////
    Consumer pp_consumer(output_file, method, pp);
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
    String method = getStringOption_("method");

    double cluster_ppm = getDoubleOption_("pickIMCluster:ppm_tolerance");
    double cluster_im  = getDoubleOption_("pickIMCluster:im_tolerance");
    double traces_ppm  = getDoubleOption_("pickIMElutionProfiles:ppm_tolerance");

    // pull PeakIMTraces parameters
    PeakPickerIM picker;
    Param pickimtr_params = getParam_().copy("pickIMTraces:", /*remove_prefix=*/true);
    if (process_option == "lowmemory")
    {
      return doLowMemAlgorithm(method, picker, input_file, output_file); // TODO: needs parallelization
    }
    else
    {
      // Load input mzML file
      PeakMap exp;
      MzMLFile mzml;
      mzml.load(input_file, exp);

      // Process each spectrum with PeakPickerIM
      #pragma omp parallel for
      for (Int64 i = 0; i != exp.size(); i++)
      {
        MSSpectrum& spectrum = exp[i];
        OPENMS_LOG_DEBUG << "Processing MS" << spectrum.getMSLevel() << " spectrum with " 
          << spectrum.size() << " peaks in the IM frame." << std::endl;
        if (method == "mobilogram")
        {
          picker.pickIMTraces(spectrum);
        }
        else if (method == "cluster")
        {
          const double ppm_tol = getDoubleOption_("pickIMCluster:ppm_tolerance");
          const double im_tol  = getDoubleOption_("pickIMCluster:im_tolerance");
          PeakPickerIM::pickIMCluster(spectrum, cluster_ppm, cluster_im);
        }
        else if (method == "traces")
        {
          double ppm_tol = getDoubleOption_("pickIMElutionProfiles:ppm_tolerance");
          PeakPickerIM::pickIMElutionProfiles(spectrum, ppm_tol);
        }
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

