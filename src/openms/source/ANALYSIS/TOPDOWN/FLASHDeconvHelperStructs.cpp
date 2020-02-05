// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>


namespace OpenMS
{
  FLASHDeconvHelperStructs::PrecalcularedAveragine::PrecalcularedAveragine(double m,
                                                                           double M,
                                                                           double delta,
                                                                           CoarseIsotopePatternGenerator *gen)
      :
      massInterval(delta), minMass(m)
  {
    int i = 0;
    while (true)
    {
      double a = i * massInterval;
      i++;
      if (a < m)
      {
        continue;
      }
      if (a > M)
      {
        break;
      }

      //estimateForFragmentFromPeptideWeight //
      generator = gen;
      auto iso = generator->estimateFromPeptideWeight(a);
      //iso.trimIntensities()

      //std::cout<< a << " "<<iso[10].getMZ() - iso[9].getMZ()<<std::endl;

      auto factor = .001;
      iso.trimRight(factor * iso.getMostAbundant().getIntensity());

      // double norm = .0;
      // double mean = .0;
      Size mostAbundantIndex = 0;
      double mostAbundantInt = 0;

      for (Size k = 0; k < iso.size(); k++)
      {
        if (mostAbundantInt >= iso[k].getIntensity())
        {
          continue;
        }
        mostAbundantInt = iso[k].getIntensity();
        mostAbundantIndex = k;
      }

      Size leftIndex = mostAbundantIndex;
      Size rightIndex = iso.size() - 1 - mostAbundantIndex;

      //iso.averageMass()
      //iso[0].getMZ()

      rightIndices.push_back(rightIndex);
      leftIndices.push_back(leftIndex);
      isotopes.push_back(iso);

      //std::cout<< a << " " << mostAbundantIndex <<" " << endIndex << std::endl;
      //auto mostAbundant = iso.getMostAbundant();
    }
  }

  IsotopeDistribution FLASHDeconvHelperStructs::PrecalcularedAveragine::get(double mass)
  {
    Size i = (Size) (.5 + (mass - minMass) / massInterval);
    i = i >= isotopes.size() ? isotopes.size() - 1 : i;
    return isotopes[i];
  }

  IsotopeDistribution FLASHDeconvHelperStructs::PrecalcularedAveragine::get(double precursorMass, double mass, std::set<UInt>& precursors)
  {
    //return get(mass);
    if(precursors.size() == 0 || precursorMass <= mass){
      return get(mass);
    }
    auto id_fragment = get(mass);
    auto id_comp_fragment = get(precursorMass - mass);
    return generator->calcFragmentIsotopeDist(id_fragment, id_comp_fragment, precursors, mass);
  }

  Size FLASHDeconvHelperStructs::PrecalcularedAveragine::getLeftIndex(double mass)
  {
    Size i = (Size) (.5 + (mass - minMass) / massInterval);
    i = i >= isotopes.size() ? isotopes.size() - 1 : i;
    return leftIndices[i];
  }

  double FLASHDeconvHelperStructs::PrecalcularedAveragine::getAverageMassDelta(IsotopeDistribution dist)
  {
    return dist.averageMass() - dist[0].getMZ();
  }

  Size FLASHDeconvHelperStructs::PrecalcularedAveragine::getRightIndex(double mass)
  {
    Size i = (Size) (.5 + (mass - minMass) / massInterval);
    i = i >= isotopes.size() ? isotopes.size() - 1 : i;
    return rightIndices[i];
  }

  FLASHDeconvHelperStructs::LogMzPeak::LogMzPeak() :
      mz(0),
      intensity(0),
      logMz(-1000),
      charge(0),
      isotopeIndex(0)
  {
  }

  FLASHDeconvHelperStructs::LogMzPeak::LogMzPeak(Peak1D &peak) :
      mz(peak.getMZ()),
      intensity(peak.getIntensity()),
      logMz(getLogMz(peak.getMZ())),
      charge(0),
      isotopeIndex(0)
  {
  }

  FLASHDeconvHelperStructs::LogMzPeak::LogMzPeak(LogMzPeak &peak, int c, int i) :
      mz(peak.mz),
      intensity(peak.intensity),
      logMz(peak.logMz),
      charge(c),
      isotopeIndex(i)
  {
  }

  FLASHDeconvHelperStructs::LogMzPeak::~LogMzPeak()
  {
  }

  double FLASHDeconvHelperStructs::LogMzPeak::getUnchargedMass()
  {
    if (mass <= 0)
    {
      mass = (mz - Constants::PROTON_MASS_U) * charge;
    }
    return mass;
  }

  bool FLASHDeconvHelperStructs::LogMzPeak::operator<(const LogMzPeak &a) const
  {
    return this->logMz < a.logMz;
  }

  bool FLASHDeconvHelperStructs::LogMzPeak::operator>(const LogMzPeak &a) const
  {
    return this->logMz >= a.logMz;
  }

  FLASHDeconvHelperStructs::PeakGroup::~PeakGroup()
  {
    std::vector<LogMzPeak>().swap(peaks);
  }

  void FLASHDeconvHelperStructs::PeakGroup::push_back(FLASHDeconvHelperStructs::LogMzPeak &p)
  {
    peaks.push_back(p);
  }

  void FLASHDeconvHelperStructs::PeakGroup::reserve(Size n)
  {
    peaks.reserve(n);
  }

  bool FLASHDeconvHelperStructs::PeakGroup::operator<(const FLASHDeconvHelperStructs::PeakGroup &a) const
  {
    //if(this->spec->getRT() == a.spec->getRT()){
    return  this->monoisotopicMass < a.monoisotopicMass;
    //}
    //return this->spec->getRT() < a.spec->getRT();
  }

  bool FLASHDeconvHelperStructs::PeakGroup::operator>(const FLASHDeconvHelperStructs::PeakGroup &a) const
  {
    //if(this->spec->getRT() == a.spec->getRT()){
    return  this->monoisotopicMass > a.monoisotopicMass;
    //}
    //return this->spec->getRT() > a.spec->getRT();
  }

  bool FLASHDeconvHelperStructs::PeakGroup::operator==(const PeakGroup &a) const
  {
    return// this->spec->getRT() == a.spec->getRT() &&
        this->monoisotopicMass == a.monoisotopicMass
        && this->intensity == a.intensity;
  }

  void FLASHDeconvHelperStructs::PeakGroup::updateMassesAndIntensity(FLASHDeconvHelperStructs::PrecalcularedAveragine &averagines,
                                                                     unsigned int &msLevel,
                                                                     Parameter& param,
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
    //double nominator = .0;
    //double denominator = .0;
    //auto masses = std::vector<double>();
    double maxIntensityForMonoIsotopeMass = -1;
    for (auto &p : peaks)
    {
      double pi = p.intensity;
      intensity += pi;
      // auto w = std::sqrt(pi);
      // denominator += w;
      // nominator += w * (p.getMass() - p.isotopeIndex * Constants::C13C12_MASSDIFF_U);
      //masses.push_back(p.getUnchargedMass() - p.isotopeIndex * Constants::ISOTOPE_MASSDIFF_55K_U);

      if (maxIntensityForMonoIsotopeMass > pi)
      {
        continue;
      }
      maxIntensityForMonoIsotopeMass = pi;
      monoisotopicMass = p.getUnchargedMass() - p.isotopeIndex * Constants::ISOTOPE_MASSDIFF_55K_U;
      //  int mostAbundantIndex = averagines.getMostAbundantIndex(monoisotopicMass);
      // avgMass = p.getMass() + (mostAbundantIndex - p.isotopeIndex) * Constants::C13C12_MASSDIFF_U;

    }
    //sort(masses.begin(), masses.end());
    //monoisotopicMass = masses[masses.size()/2];
    //monoisotopicMass = nominator / denominator;
    if(maxIsoIndex>0)
    {
      auto dist = msLevel == 1 ?
                  averagines.get(monoisotopicMass) :
                  averagines.get(param.currentMaxMass, monoisotopicMass, param.precursorIsotopes);
      auto massDelta = averagines.getAverageMassDelta(dist);
      avgMass = monoisotopicMass + massDelta;
    }
  }

  double FLASHDeconvHelperStructs::getLogMz(double mz)
  {
    return std::log(mz - Constants::PROTON_MASS_U);
  }
}
