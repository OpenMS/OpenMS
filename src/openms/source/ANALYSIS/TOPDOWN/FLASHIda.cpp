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
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{
  // constructor
  FLASHIda::FLASHIda(char *arg)
  {
    std::unordered_map<std::string, std::vector<double>> inputs;
    char *token = std::strtok(arg, " ");
    std::string key;
    while (token != nullptr)
    {
      auto tokenString = std::string(token);
      auto num = atof(tokenString.c_str());

      if (num == 0.0)
      {
        key = tokenString;
        inputs[key] = DoubleList();
      }
      else
      {
        inputs[key].push_back(num);
      }
      token = std::strtok(nullptr, " ");
    }

    Param fd_defaults = FLASHDeconvAlgorithm().getDefaults();
    // overwrite algorithm default so we export everything (important for copying back MSstats results)
    fd_defaults.setValue("min_charge", (int) inputs["min_charge"][0]);
    fd_defaults.setValue("max_charge", (int) inputs["max_charge"][0]);
    fd_defaults.setValue("min_mass", inputs["min_mass"][0]);
    fd_defaults.setValue("max_mass", inputs["max_mass"][0]);
   
    fd_defaults.setValue("tol", inputs["tol"]);
    fd_defaults.setValue("RT_window", 20.0, "");

    auto massCountd = inputs["max_mass_count"];

    for (double j : massCountd)
    {
      massCount.push_back((int) j);
    }

    fd_defaults.setValue("min_mass_count", massCount);

    fd.setParameters(fd_defaults);
    fd.calculateAveragine(false);

    RTwindow = inputs["RT_window"][0];
    qScoreThreshold = inputs["score_threshold"][0];

    avg = fd.getAveragine();
    selected = std::map<int, std::vector<double>>(); // int mass, rt, qscore
    std::cout << fd_defaults << std::endl;
  }

  int FLASHIda::getPeakGroups(double *mzs,
                              double *ints,
                              int length,
                              double rt,
                              int msLevel,
                              char *name)
  {

    auto spec = makeMSSpectrum(mzs, ints, length, rt, msLevel, name);
    deconvolutedSpectrum = DeconvolutedSpectrum(spec, 1);
    if (msLevel == 1)
    {
      //currentMaxMass = maxMass;
      //currentChargeRange = chargeRange;
    }
    else
    {
      //TODO precursor infor here
    }
   
    // per spec deconvolution
    //    int specIndex = 0, massIndex = 0; // ..
    fd.fillPeakGroupsInDeconvolutedSpectrum(deconvolutedSpectrum, 0);
    FLASHIda::filterPeakGroupsUsingMassExclusion(spec, msLevel);

   // spec.clear(true);

    return deconvolutedSpectrum.size();
  }

  void FLASHIda::filterPeakGroupsUsingMassExclusion(MSSpectrum &spec, int msLevel)
  {
    double rt = spec.getRT();

    std::map<int, std::vector<double>> nselected;

    for (auto &item : selected)
    {
      if (item.first < rt - RTwindow)
      {
        continue;
      }
      nselected[item.first] = item.second;
    }

    //auto shorRTwindow = RTwindow / 10.0;
    std::vector<PeakGroup> filtered;
    filtered.reserve(massCount.size());
    std::sort(deconvolutedSpectrum.begin(), deconvolutedSpectrum.end(), QscoreComparator);
    std::set<int> cselected; // current selected masses
    int mc =  massCount[msLevel - 1];

    for (auto &pg : deconvolutedSpectrum)
    {
      if (filtered.size() >= mc)
      {
        break;
      }

      auto massDelta = avg.getAverageMassDelta(pg.getMonoMass());
      auto m = FLASHDeconvAlgorithm::getNominalMass(pg.getMonoMass() + massDelta);
      if (cselected.find(m) != cselected.end()){
        continue;
      }
      auto qScore = pg.getQScore();
      if (nselected.find(m) != nselected.end())
      {
        continue;
      }

      nselected[m] = std::vector<double>(2);
      cselected.insert(m);
      nselected[m][0] = rt;
      nselected[m][1] = qScore;
      filtered.push_back(pg);
    }

    if (filtered.size() < mc)
    {
      for (auto &pg : deconvolutedSpectrum)
      {
        if (filtered.size() >= mc)
        {
          break;
        }

        auto massDelta = avg.getAverageMassDelta(pg.getMonoMass());
        auto m = FLASHDeconvAlgorithm::getNominalMass(pg.getMonoMass() + massDelta);
        if (cselected.find(m) != cselected.end()){
          continue;
        }
        auto qScore = pg.getQScore();
        if (nselected.find(m) != nselected.end())
        {
          if (qScore < nselected[m][1])
          {
            continue;
          }
        }
        else
        {
          nselected[m] = std::vector<double>(2);
        }
        cselected.insert(m);
        nselected[m][0] = rt;
        nselected[m][1] = qScore;
        filtered.push_back(pg);
      }
    }


    if (filtered.size() < mc)
    {
      for (auto &pg : deconvolutedSpectrum)
      {
        if (filtered.size() >= mc)
        {
          break;
        }

        auto massDelta = avg.getAverageMassDelta(pg.getMonoMass());
        auto m = FLASHDeconvAlgorithm::getNominalMass(pg.getMonoMass() + massDelta);
        if (cselected.find(m) != cselected.end()){
          continue;
        }
        auto qScore = pg.getQScore();
        if (nselected.find(m) == nselected.end())
        {
          nselected[m] = std::vector<double>(2);
        }
        cselected.insert(m);
        nselected[m][0] = rt;
        nselected[m][1] = qScore;
        filtered.push_back(pg);
      }
    }

    nselected.swap(selected);
    std::map<int, std::vector<double>>().swap(nselected);

    deconvolutedSpectrum.swap(filtered);
    std::vector<PeakGroup>().swap(filtered);
  }

  void FLASHIda::getIsolationWindows(double *wstart, double *wend, double *qScores, int *charges, double *avgMasses)
  {
    std::sort(deconvolutedSpectrum.begin(), deconvolutedSpectrum.end(), QscoreComparator);

    for (auto i = 0; i < deconvolutedSpectrum.size(); i++)
    {
      auto qrange = deconvolutedSpectrum[i].getMzxQScoreMzRange();
      wstart[i] = std::get<0>(qrange) - 1.2;
      wend[i] = std::get<1>(qrange) + 1.2;

      qScores[i] = deconvolutedSpectrum[i].getQScore();
      charges[i] = deconvolutedSpectrum[i].getRepCharge();
      auto massDelta = avg.getAverageMassDelta(deconvolutedSpectrum[i].getMonoMass());
      avgMasses[i] = massDelta + deconvolutedSpectrum[i].getMonoMass();
    }
    std::vector<PeakGroup>().swap(deconvolutedSpectrum);
  }

  MSSpectrum FLASHIda::makeMSSpectrum(double *mzs, double *ints, int length, double rt, int msLevel, char *name)
  {
    auto spec = MSSpectrum();
    for (auto i = 0; i < length; i++)
    {
      if (ints[i] <= 0){
        continue;
      }
      spec.push_back(Peak1D(mzs[i], ints[i]));
    }
    spec.setMSLevel(msLevel);
    spec.setName(name);
    spec.setRT(rt);//
    return spec;
  }
}