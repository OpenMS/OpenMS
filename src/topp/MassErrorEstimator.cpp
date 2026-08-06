#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/PROCESSING/CALIBRATION/MassErrorEstimator.h>

using namespace OpenMS;

class TOPPMassErrorEstimator :
  public TOPPBase
{
public:
  TOPPMassErrorEstimator() :
    TOPPBase("MassErrorEstimator",
             "Estimates mass error from a polysiloxane contaminant marker.", 
             false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input mzML file");
    setValidFormats_("in", {"mzML"});
  }

  ExitCodes main_(int, const char**) override
  {
    String in = getStringOption_("in");

    MSExperiment exp;
    MzMLFile().load(in, exp);

    // detect whether MS1 spectra are profile or centroided and if needed apply peak picking
    bool profile = false;
    for (const auto& spec : exp.getSpectra())
    {
      if (spec.getMSLevel() != 1) continue;

      if (spec.getType() == MSSpectrum::SpectrumType::PROFILE)
      {
        profile = true;
        break;
      }
    }

    MSExperiment exp_centroided;

    if (profile)
    {
      writeLogInfo_("Profile data detected -> applying peak picking");
      PeakPickerHiRes picker;
      picker.pickExperiment(exp, exp_centroided);
    }
    else
    {
      writeLogInfo_("Input data already centroided");
      exp_centroided = exp;
    }

    MassErrorEstimator estimator;
    estimator.estimate(exp_centroided);

    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPMassErrorEstimator tool;
  return tool.main(argc, argv);
}