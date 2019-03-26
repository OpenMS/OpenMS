#pragma once
#include <OpenMS/QC/QCBase.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Types.h>

class MZcalibration : QCBase
{
	/**
   * @brief This class calculate the difference between m/z before and after correction

   */
public:

	MZcalibration() = default;
	void calculate (FeatureMap& features, const MSExperiment& exp);
	void clearResult();

private:
	double getMZraw(double rt, const MSExperiment& exp);
	double EPSILON{ 0.05 };
	//ProteinIdentification::PeakMassType mass_type;
	Status requires() const override;
};