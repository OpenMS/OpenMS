
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <boost/shared_ptr.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

using namespace OpenMS;


class IMSDriftConverter :
  public TOPPBase
{
private:
  //base functor
  class BaseFunctor
  {
  public:
    BaseFunctor() : temp0_(273.15), pressure0_(1013.25) {}
    virtual double convert (double x) = 0;
    virtual double revert (double x) = 0;
  protected:
    const double temp0_ ; //0 celsius temperature in kelvin
    const double pressure0_; //atmospheric pressure in millibars
  };
  typedef boost::shared_ptr<BaseFunctor> BaseFunctorPtr;

  //functor to do no conversion
  class UnitFunctor : public BaseFunctor
  {
  public:
    double convert (double x) { return x; }
    double revert (double x) { return x; }
  };

  //functor to convert from raw to inverse mobility and vice-versa
  class InverseMobilityFunctor : public BaseFunctor
  {
  public:
    InverseMobilityFunctor(double a, double b) : current_(a), tube_length_(b) {}
    double convert (double x)
    {
      return x * current_ / (tube_length_ * tube_length_);
    }
    double revert (double x)
    {
      return x * (tube_length_ * tube_length_) / current_;
    }
  private:
    double current_;
    double tube_length_;
  };

  //functor to convert from raw to inverse reduced mobility and vice-versa
  class InverseReducedMobilityFunctor : public BaseFunctor
  {
  public:
    InverseReducedMobilityFunctor(double a, double b, double c, double d)
    : current_(a), tube_length_(b), pressure_(c), temp_(d) {}
    double convert (double x)
    {
      return (x * current_ / (tube_length_ * tube_length_)) * (pressure0_ / pressure_) * (temp_ / temp0_);
    }
    double revert (double x)
    {
      return (x * (tube_length_ * tube_length_) / current_) * (pressure_ / pressure0_) * (temp0_ / temp_ );
    }
  private:
    double current_;
    double tube_length_;
    double pressure_;
    double temp_;
  };

  //functor to convert from raw to corrected inverse reduced mobility and vice-versa
  class CorrectedInverseReducedMobilityFunctor : public BaseFunctor
  {
  public:
    CorrectedInverseReducedMobilityFunctor(double a, double b, double c, double d, double e)
    : current_(a), tube_length_(b), pressure_(c), temp_(d), unit_correction_factor_(e) {}
    double convert (double x)
    {
      return (x * current_ / (tube_length_ * tube_length_)) * (pressure0_ / pressure_) * (temp_ / temp0_) * unit_correction_factor_;
    }
    double revert (double x)
    {
      return (x * (tube_length_ * tube_length_) / current_) * (pressure_ / pressure0_) * (temp0_ / temp_ ) * (1 / unit_correction_factor_);
    }
  private:
    double current_;
    double tube_length_;
    double pressure_;
    double temp_;
    double unit_correction_factor_;
  };

  //determine the source drift intensity unit and chose a functor
  BaseFunctorPtr calcRevertFunctor_(MSExperiment<> & exp)
  {
    double current = std::fabs((double)exp.getMetaValue("U"));//fabs
    double tube_length = exp.getMetaValue("L");
    double pressure = exp.getMetaValue("P");
    double temp = exp.getMetaValue("T");
    double unit_correction_factor = exp.getMetaValue("Z");

    BaseFunctorPtr ptr;
    if(exp.getMetaValue("drift_unit") == "raw")
      ptr = BaseFunctorPtr( new UnitFunctor() );
    else if(exp.getMetaValue("drift_unit") == "inverse_mobility")
      ptr = BaseFunctorPtr( new InverseMobilityFunctor(current, tube_length) );
    else if(exp.getMetaValue("drift_unit") == "inverse_reduced_mobility")
      ptr = BaseFunctorPtr( new InverseReducedMobilityFunctor(current, tube_length, pressure, temp) );
    else if(exp.getMetaValue("drift_unit") == "corrected_inverse_reduced_mobility")
      ptr = BaseFunctorPtr( new CorrectedInverseReducedMobilityFunctor(current, tube_length, pressure, temp, unit_correction_factor) );
    else
      throw;

    return ptr;
  }

  //determine the target drift intensity unit and chose a functor
  BaseFunctorPtr calcConvertFunctor_(MSExperiment<> & exp)
  {
    double current = std::fabs((double)exp.getMetaValue("U"));//fabs
    double tube_length = exp.getMetaValue("L");
    double pressure = exp.getMetaValue("P");
    double temp = exp.getMetaValue("T");
    double unit_correction_factor = exp.getMetaValue("Z");

    BaseFunctorPtr ptr;
    if(getStringOption_("drift_unit") == "raw")
      ptr = BaseFunctorPtr( new UnitFunctor() );
    else if(getStringOption_("drift_unit") == "inverse_mobility")
      ptr = BaseFunctorPtr( new InverseMobilityFunctor(current, tube_length) );
    else if(getStringOption_("drift_unit") == "inverse_reduced_mobility")
      ptr = BaseFunctorPtr( new InverseReducedMobilityFunctor(current, tube_length, pressure, temp) );
    else if(getStringOption_("drift_unit") == "corrected_inverse_reduced_mobility")
      ptr = BaseFunctorPtr( new CorrectedInverseReducedMobilityFunctor(current, tube_length, pressure, temp, unit_correction_factor) );
    else
      throw;

    return ptr;
  }

  //transform in raw drift intensities
  void transformToRaw_(MSExperiment<> & exp, BaseFunctorPtr functor)
  {
    MSExperiment<>::iterator expIter;
    for(expIter=exp.begin(); expIter!=exp.end(); ++expIter)
    {
      MSSpectrum<>::ContainerType::iterator specIter;
      for(specIter=expIter->begin(); specIter!=expIter->end(); ++specIter)
      {
	//std::cout << "orig: " << specIter->getMZ() << " revert: " << functor->revert(specIter->getMZ()) << std::endl;//DEBUG
        specIter->setMZ( functor->revert(specIter->getMZ()) );
      }
    }
  }

  //transform to target drift intensities, e.g. inverse mobility
  void transformToTarget_(MSExperiment<> & exp, BaseFunctorPtr functor)
  {
    MSExperiment<>::iterator expIter;
    for(expIter=exp.begin(); expIter!=exp.end(); ++expIter)
    {
      MSSpectrum<>::ContainerType::iterator specIter;
      for(specIter=expIter->begin(); specIter!=expIter->end(); ++specIter)
      {
        std::cout << "orig: " << specIter->getMZ() << "convert: " << functor->convert(specIter->getMZ()) << std::endl;//DEBUG
        specIter->setMZ( functor->convert(specIter->getMZ()) );
      }
    }
  }

  //check if the data will not change by the transformation
  bool checkUnitTransformation_(MSExperiment<> & exp)
  {
    return exp.getMetaValue("drift_unit") == getStringOption_("drift_unit");
  }
  //transform drift intensities
  void transform_(MSExperiment<> & exp)
  {
    writeLog_("Transforming IMS data from '" + String(exp.getMetaValue("drift_unit")) + "' to '" + getStringOption_("drift_unit") + "'");
    if(checkUnitTransformation_(exp)) return;

    //convert to raw drift intensities
    transformToRaw_(exp, calcRevertFunctor_(exp));

    //convert to target drift intensities, e.g. inverse mobility
    transformToTarget_(exp, calcConvertFunctor_(exp));

    exp.setMetaValue("drift_unit", getStringOption_("drift_unit"));
  }

public:
  IMSDriftConverter() :
    TOPPBase("IMSDriftConverter", "Reads a binary file from an IMS machine and exports the data as mzML.", false, true)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<input-file>", "", "The binary input files.",true);
    registerOutputFile_("out", "<output-file>", "", "The mzML file containing the IMS spectra.",true);
    setValidFormats_("out", ListUtils::create<String>("mzML"));

    registerStringOption_("drift_unit", "<unit>", "inverse_mobility", "Unit of IMS drift parameter.", false);
    setValidStrings_("drift_unit", ListUtils::create<String>("raw,inverse_mobility,inverse_reduced_mobility,corrected_inverse_reduced_mobility"));
  }

  ExitCodes main_(int, const char **)
  {
    String input_file = getStringOption_("in");
    String output_file = getStringOption_("out");

    MSExperiment<> exp;
    MzMLFile f;
    f.setLogType(log_type_);
    f.load(input_file, exp);

    if (exp.empty())
    {
      LOG_WARN << "The loaded file is empty.";
      return INCOMPATIBLE_INPUT_DATA;
    }

    transform_(exp);

    MzMLFile().store(output_file, exp);

    return EXECUTION_OK;
  }
};

int main(int argc, const char ** argv)
{
  IMSDriftConverter tool;
  return tool.main(argc, argv);
}
