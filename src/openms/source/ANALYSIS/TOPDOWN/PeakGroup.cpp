//
// Created by Kyowon Jeong on 4/22/20.
//
#include "include/OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h"
namespace OpenMS
{

  PeakGroup::~PeakGroup()
  {
    std::vector<LogMzPeak>().swap(peaks);
    clearChargeInfo();
  }

  void PeakGroup::push_back(FLASHDeconvHelperStructs::LogMzPeak &p)
  {
    peaks.push_back(p);
  }

  void PeakGroup::reserve(Size n)
  {
    peaks.reserve(n);
  }

  void PeakGroup::clearChargeInfo()
  {
    for (auto& item : perChargeInfo)
    {
      std::vector<float>().swap(item.second);
    }
    std::unordered_map<int, std::vector<float>>().swap(perChargeInfo);
  }

  bool PeakGroup::operator<(const PeakGroup &a) const
  {
    //if(this->spec->getRT() == a.spec->getRT()){
    if(this->monoisotopicMass == a.monoisotopicMass){
      return this->intensity < a.intensity;
    }
    return this->monoisotopicMass < a.monoisotopicMass;
    //}
    //return this->spec->getRT() < a.spec->getRT();
  }

  bool PeakGroup::operator>(const PeakGroup &a) const
  {
    //if(this->spec->getRT() == a.spec->getRT()){
    if(this->monoisotopicMass == a.monoisotopicMass){
      return this->intensity > a.intensity;
    }
    return this->monoisotopicMass > a.monoisotopicMass;
    //}
    //return this->spec->getRT() > a.spec->getRT();
  }

  bool PeakGroup::operator==(const PeakGroup &a) const
  {
    return// this->spec->getRT() == a.spec->getRT() &&
        this->monoisotopicMass == a.monoisotopicMass
        && this->intensity == a.intensity;
  }

  void PeakGroup::updateMassesAndIntensity(FLASHDeconvHelperStructs::PrecalculatedAveragine &averagines,
                                           double chargeMass,
                                           int offset,
                                           int maxIsoIndex)
  {
    //

    if (offset != 0)
    {
      std::vector<LogMzPeak> tmpPeaks;
      tmpPeaks.swap(peaks);
      peaks.reserve(tmpPeaks.size());

      for (auto &p : tmpPeaks)
      {
        p.isotopeIndex -= offset;
        if (p.isotopeIndex < 0 || p.isotopeIndex >= maxIsoIndex)
        {
          continue;
        }
        peaks.push_back(p);
      }
    }

    intensity = .0;
    double nominator = .0;

    for (auto &p : peaks)
    {
      double pi = p.intensity;
      intensity += pi;
      nominator += pi * (p.getUnchargedMass(chargeMass) - p.isotopeIndex * Constants::C13C12_MASSDIFF_U);

    }
    monoisotopicMass = nominator / intensity;
    auto massDelta = averagines.getAverageMassDelta(monoisotopicMass);
    avgMass = monoisotopicMass + massDelta;

  }
}