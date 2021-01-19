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
  FLASHIda::FLASHIda(char* arg)
  {

    std::unordered_map<std::string, std::vector<double>> inputs;
    char* token = std::strtok(arg, " ");
    std::string key;
    while (token != nullptr)
    {
      String tokenString = std::string(token);
      double num = atof(tokenString.c_str());

      if (num == 0 && !isdigit(tokenString[tokenString.size() - 1]))
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
    fd_defaults.setValue("min_charge", (int)inputs["min_charge"][0]);
    fd_defaults.setValue("max_charge", (int)inputs["max_charge"][0]);
    fd_defaults.setValue("min_mass", inputs["min_mass"][0]);
    fd_defaults.setValue("max_mass", inputs["max_mass"][0]);

    fd_defaults.setValue("tol", inputs["tol"]);
    fd_defaults.setValue("RT_window", 20.0, "");

    auto massCountd = inputs["max_mass_count"];

    for (double j : massCountd)
    {
      mass_count.push_back((int)j);
    }

    fd_defaults.setValue("min_mass_count", mass_count);

    fd.setParameters(fd_defaults);
    fd.calculateAveragine(false);

    rt_window = inputs["RT_window"][0];
    qscore_threshold = inputs["score_threshold"][0];

    avg = fd.getAveragine();
    std::cout << "QScore threshold: " << qscore_threshold << std::endl;
    std::cout << fd_defaults << std::endl;
  }

  int FLASHIda::getPeakGroups(const double* mzs,
                              const double* ints,
                              const int length,
                              const double rt,
                              const int ms_level,
                              const char* name)
  {

    auto spec = makeMSSpectrum(mzs, ints, length, rt, ms_level, name);
    deconvoluted_spectrum = DeconvolutedSpectrum(spec, 1);
    if (ms_level == 1)
    {
      //current_max_mass = max_mass;
      //currentChargeRange = chargeRange;
    }
    else
    {
      //TODO precursor infor here
    }

    // per spec deconvolution
    //    int specIndex = 0, massIndex = 0; // ..
    fd.fillPeakGroupsInDeconvolutedSpectrum(deconvoluted_spectrum, 0);
    FLASHIda::filterPeakGroupsUsingMassExclusion(spec, ms_level);

    // spec.clear(true);

    return deconvoluted_spectrum.size();
  }


  void FLASHIda::filterPeakGroupsUsingMassExclusion(const MSSpectrum& spec, const int ms_level)
  {
    double rt = spec.getRT();
    std::sort(deconvoluted_spectrum.begin(), deconvoluted_spectrum.end(), QscoreComparator);
    int mc = mass_count[ms_level - 1];

    const auto order = std::vector<char>({'B', 'R', 'G', 'b', 'r'});

    std::unordered_map<int, std::vector<double>> nall;;
    std::unordered_map<int, char> ncolor;

    for (auto& item : all)
    {
      if (item.second[0] < rt - rt_window)
      {
        continue;
      }
      nall[item.first] = item.second;
    }

    for (auto& item : color)
    {
      if (nall.find(item.first) == nall.end())
      {
        continue;
      }
      ncolor[item.first] = item.second;
    }

    for (auto& pg : deconvoluted_spectrum) // now update color and nall
    {
      int m = FLASHDeconvAlgorithm::getNominalMass(pg.getMonoMass());
      double qScore = pg.getQScore();

      if (nall.find(m) == nall.end()){ // new mass
        nall[m] = std::vector<double>(2);
        nall[m][1] = qScore;
        ncolor[m] = 'G';
      }else if (nall[m][1] < qScore) { // increasing mass
        nall[m][1] = qScore;
        if(ncolor.find(m) == ncolor.end() || ncolor[m] == 'G' || ncolor[m] == 'R'){
          ncolor[m] = 'B';
        }else if(ncolor[m] == 'r'){
          ncolor[m] = 'b';
        }
      }else{ // decreasing mass
        if(ncolor.find(m) == ncolor.end() || ncolor[m] == 'G' || ncolor[m] == 'B'){
          ncolor[m] = 'R';
        }else if(ncolor[m] == 'b'){
          ncolor[m] = 'r';
        }
      }
      nall[m][0] = rt;
    }

    nall.swap(all);
    std::unordered_map<int, std::vector<double>>().swap(nall);

    ncolor.swap(color);
    std::unordered_map<int, char>().swap(ncolor);

    std::vector<PeakGroup> filtered;
    filtered.reserve(mass_count.size());
    std::set<int> cselected; // current selected masses
    std::set<int> cselectedMz; // current selected mzs

    for (auto& c: order)
    {
      for (auto& pg : deconvoluted_spectrum)
      {
        if (filtered.size() >= mc){
          break;
        }
        if (pg.getQScore() < qscore_threshold)
        {
          break;
        }
        int mz = (int)round((std::get<0>(pg.getMzxQScoreMzRange()) + std::get<1>(pg.getMzxQScoreMzRange())) / 2.0);
        if (cselectedMz.find(mz) != cselectedMz.end()) {
          continue;
        }

        int m = FLASHDeconvAlgorithm::getNominalMass(pg.getMonoMass());
        if (cselected.find(m) != cselected.end()) {
          continue;
        }
        if(color.find(m) == color.end() || color[m] != c){
          continue;
        }

        char pc = color[m];
        if(pc == 'B' && pg.getQScore() > 0){
          color[m] = 'b';
        }else if (pc == 'R'&& pg.getQScore() > 0){
          color[m] = 'r';
        }

        filtered.push_back(pg);
        cselected.insert(m);
        cselectedMz.insert(mz);
      }
    }

    deconvoluted_spectrum.swap(filtered);
    std::vector<PeakGroup>().swap(filtered);
  }


  /*
  void FLASHIda::filterPeakGroupsUsingMassExclusion(MSSpectrum& spec, int msLevel)
  {
    double rt = spec.getRT();
    std::sort(deconvoluted_spectrum.begin(), deconvoluted_spectrum.end(), QscoreComparator);

    std::map<int, std::vector<double>> nall;
    std::map<int, std::vector<double>> nallMz;

    std::map<int, std::vector<double>> nselected;
    std::map<int, std::vector<double>> nselectedMz; // m/z * 20
    std::vector<PeakGroup> toConsider;
    toConsider.reserve(deconvoluted_spectrum.size());
    std::vector<PeakGroup> toConsider2;
    toConsider2.reserve(deconvoluted_spectrum.size());
    std::vector<PeakGroup> toConsider3;
    toConsider3.reserve(deconvoluted_spectrum.size());

    for (auto& item : all)
    {
      if (item.second[0] < rt - rt_window)
      {
        continue;
      }
      nall[item.first] = item.second;
      if(selected.find(item.first) == selected.end()){
        continue;
      }
      nselected[item.first] = selected[item.first];
    }

    for (auto& item : allMz)
    {
      if (item.second[0] < rt - rt_window)
      {
        continue;
      }
      nallMz[item.first] = item.second;
      if(selectedMz.find(item.first) == selectedMz.end()){
        continue;
      }
      nselectedMz[item.first] = selectedMz[item.first];
    }

    for (auto& pg : deconvoluted_spectrum)
    {
      if (pg.getQScore() < qscore_threshold)
      {
        break;
      }
      auto massDelta = avg.getAverageMassDelta(pg.getMonoMass());
      auto m = FLASHDeconvAlgorithm::getNominalMass(pg.getMonoMass() + massDelta);
      auto mz = (int)round((std::get<0>(pg.getMzxQScoreMzRange()) + std::get<1>(pg.getMzxQScoreMzRange())) / 2.0);

      auto qscore = pg.getQScore();
      if ((nall.find(m) != nall.end() && nall[m][1] < qscore)
          || (nallMz.find(mz) != nallMz.end() && nallMz[mz][1] < qscore))
      {
        toConsider.push_back(pg);
        continue;
      }

      if (nall.find(m) != nall.end()
          || nallMz.find(mz) != nallMz.end())
      {
        toConsider2.push_back(pg);
        continue;
      }
      toConsider3.push_back(pg);
    }

    for (auto& pg : deconvoluted_spectrum)
    {
      auto massDelta = avg.getAverageMassDelta(pg.getMonoMass());
      auto m = FLASHDeconvAlgorithm::getNominalMass(pg.getMonoMass() + massDelta);
      auto mz = (int)round((std::get<0>(pg.getMzxQScoreMzRange()) + std::get<1>(pg.getMzxQScoreMzRange())) / 2.0);

      auto qscore = pg.getQScore();
      if (nall.find(m) == nall.end()) {
        nall[m] = std::vector<double>(2);
        nall[m][1] = qscore;
      }
      else if (nall[m][1] < qscore) {
        nall[m][1] = qscore;
      }

      nall[m][0] = rt;

      if (nallMz.find(mz) == nallMz.end()) {
        nallMz[mz] = std::vector<double>(2);
        nallMz[mz][1] = qscore;
      }
      else if (nallMz[mz][1] < qscore) {
        nallMz[mz][1] = qscore;
      }

      nallMz[mz][0] = rt;
    }

    nall.swap(all);
    std::map<int, std::vector<double>>().swap(nall);

    nallMz.swap(allMz);
    std::map<int, std::vector<double>>().swap(nallMz);

    std::vector<PeakGroup> filtered;
    filtered.reserve(mass_count.size());
    std::set<int> cselected; // current selected masses
    std::set<int> cselectedMz; // current selected mzs

    int mc = mass_count[msLevel - 1];

    if (filtered.size() < mc)
    {
      for (auto& pg : toConsider) // increasing and not selected
      {
        if (filtered.size() >= mc)
        {
          break;
        }

        auto mz = (int)round((std::get<0>(pg.getMzxQScoreMzRange()) + std::get<1>(pg.getMzxQScoreMzRange())) / 2.0);
        if (cselectedMz.find(mz) != cselectedMz.end()) {
          continue;
        }
        auto massDelta = avg.getAverageMassDelta(pg.getMonoMass());
        auto m = FLASHDeconvAlgorithm::getNominalMass(pg.getMonoMass() + massDelta);
        if (cselected.find(m) != cselected.end()) {
          continue;
        }
        auto qscore = pg.getQScore();
        if (nselectedMz.find(mz) != nselectedMz.end())
        {
          continue;
        }

        if (nselected.find(m) != nselected.end())
        {
          continue;
        }

        cselected.insert(m);
        cselectedMz.insert(mz);

        if (nselectedMz.find(mz) == nselectedMz.end()) {
          nselectedMz[mz] = std::vector<double>(2);
          nselectedMz[mz][1] = qscore;
        }
        else if (nselectedMz[mz][1] < qscore) {
          nselectedMz[mz][1] = qscore;
        }

        nselectedMz[mz][0] = rt;

        if (nselected.find(m) == nselected.end()) {
          nselected[m] = std::vector<double>(2);
          nselected[m][1] = qscore;
        }
        else if (nselected[m][1] < qscore) {
          nselected[m][1] = qscore;
        }

        nselected[m][0] = rt;

        filtered.push_back(pg);
      }
    }

    if (filtered.size() < mc)
    {
      for (auto& pg : toConsider2) // decreasing but not selected
      {
        if (filtered.size() >= mc)
        {
          break;
        }

        auto mz = (int)round((std::get<0>(pg.getMzxQScoreMzRange()) + std::get<1>(pg.getMzxQScoreMzRange())) / 2.0);
        if (cselectedMz.find(mz) != cselectedMz.end()) {
          continue;
        }
        auto massDelta = avg.getAverageMassDelta(pg.getMonoMass());
        auto m = FLASHDeconvAlgorithm::getNominalMass(pg.getMonoMass() + massDelta);
        if (cselected.find(m) != cselected.end()) {
          continue;
        }
        auto qscore = pg.getQScore();
        if (nselectedMz.find(mz) != nselectedMz.end())
        {
          continue;
        }

        if (nselected.find(m) != nselected.end())
        {
          continue;
        }

        cselected.insert(m);
        cselectedMz.insert(mz);

        if (nselectedMz.find(mz) == nselectedMz.end()) {
          nselectedMz[mz] = std::vector<double>(2);
          nselectedMz[mz][1] = qscore;
        }
        else if (nselectedMz[mz][1] < qscore) {
          nselectedMz[mz][1] = qscore;
        }

        nselectedMz[mz][0] = rt;

        if (nselected.find(m) == nselected.end()) {
          nselected[m] = std::vector<double>(2);
          nselected[m][1] = qscore;
        }
        else if (nselected[m][1] < qscore) {
          nselected[m][1] = qscore;
        }

        nselected[m][0] = rt;

        filtered.push_back(pg);
      }
    }

    if (filtered.size() < mc)
    {
      for (auto& pg : toConsider) // increasing but preselected
      {
        if (filtered.size() >= mc)
        {
          break;
        }

        auto mz = (int)round((std::get<0>(pg.getMzxQScoreMzRange()) + std::get<1>(pg.getMzxQScoreMzRange())) / 2.0);
        if (cselectedMz.find(mz) != cselectedMz.end()) {
          continue;
        }
        auto massDelta = avg.getAverageMassDelta(pg.getMonoMass());
        auto m = FLASHDeconvAlgorithm::getNominalMass(pg.getMonoMass() + massDelta);
        if (cselected.find(m) != cselected.end()) {
          continue;
        }
        auto qscore = pg.getQScore();

        cselected.insert(m);
        cselectedMz.insert(mz);
        if (nselectedMz.find(mz) == nselectedMz.end()) {
          nselectedMz[mz] = std::vector<double>(2);
          nselectedMz[mz][1] = qscore;
        }
        else if (nselectedMz[mz][1] < qscore) {
          nselectedMz[mz][1] = qscore;
        }

        nselectedMz[mz][0] = rt;

        if (nselected.find(m) == nselected.end()) {
          nselected[m] = std::vector<double>(2);
          nselected[m][1] = qscore;
        }
        else if (nselected[m][1] < qscore) {
          nselected[m][1] = qscore;
        }

        nselected[m][0] = rt;

        filtered.push_back(pg);
      }
    }

    if (filtered.size() < mc)
    {
      for (auto& pg : toConsider3) // newly emerging
      {
        if (filtered.size() >= mc)
        {
          break;
        }

        auto mz = (int)round((std::get<0>(pg.getMzxQScoreMzRange()) + std::get<1>(pg.getMzxQScoreMzRange())) / 2.0);
        if (cselectedMz.find(mz) != cselectedMz.end()) {
          continue;
        }
        auto massDelta = avg.getAverageMassDelta(pg.getMonoMass());
        auto m = FLASHDeconvAlgorithm::getNominalMass(pg.getMonoMass() + massDelta);
        if (cselected.find(m) != cselected.end()) {
          continue;
        }
        auto qscore = pg.getQScore();

        cselected.insert(m);
        cselectedMz.insert(mz);

        if (nselectedMz.find(mz) == nselectedMz.end()) {
          nselectedMz[mz] = std::vector<double>(2);
          nselectedMz[mz][1] = qscore;
        }
        else if (nselectedMz[mz][1] < qscore) {
          nselectedMz[mz][1] = qscore;
        }

        nselectedMz[mz][0] = rt;

        if (nselected.find(m) == nselected.end()) {
          nselected[m] = std::vector<double>(2);
          nselected[m][1] = qscore;
        }
        else if (nselected[m][1] < qscore) {
          nselected[m][1] = qscore;
        }

        nselected[m][0] = rt;

        filtered.push_back(pg);
      }
    }

    if (filtered.size() < mc)
    {
      for (auto& pg : toConsider2) // decreasing and preselected
      {
        if (filtered.size() >= mc)
        {
          break;
        }

        auto mz = (int)round((std::get<0>(pg.getMzxQScoreMzRange()) + std::get<1>(pg.getMzxQScoreMzRange())) / 2.0);
        if (cselectedMz.find(mz) != cselectedMz.end()) {
          continue;
        }
        auto massDelta = avg.getAverageMassDelta(pg.getMonoMass());
        auto m = FLASHDeconvAlgorithm::getNominalMass(pg.getMonoMass() + massDelta);
        if (cselected.find(m) != cselected.end()) {
          continue;
        }
        auto qscore = pg.getQScore();

        cselected.insert(m);
        cselectedMz.insert(mz);
        if (nselectedMz.find(mz) == nselectedMz.end()) {
          nselectedMz[mz] = std::vector<double>(2);
          nselectedMz[mz][1] = qscore;
        }
        else if (nselectedMz[mz][1] < qscore) {
          nselectedMz[mz][1] = qscore;
        }

        nselectedMz[mz][0] = rt;

        if (nselected.find(m) == nselected.end()) {
          nselected[m] = std::vector<double>(2);
          nselected[m][1] = qscore;
        }
        else if (nselected[m][1] < qscore) {
          nselected[m][1] = qscore;
        }

        nselected[m][0] = rt;

        filtered.push_back(pg);
      }
    }

    nselected.swap(selected);
    std::map<int, std::vector<double>>().swap(nselected);

    nselectedMz.swap(selectedMz);
    std::map<int, std::vector<double>>().swap(nselectedMz);

    deconvoluted_spectrum.swap(filtered);
    std::vector<PeakGroup>().swap(filtered);
  }
*/
  void FLASHIda::getIsolationWindows(double* wstart, double* wend, double* qscores, int* charges, double* avg_masses)
  {
    std::sort(deconvoluted_spectrum.begin(), deconvoluted_spectrum.end(), QscoreComparator);

    for (int i = 0; i < deconvoluted_spectrum.size(); i++)
    {
      auto qrange = deconvoluted_spectrum[i].getMzxQScoreMzRange();
      wstart[i] = std::get<0>(qrange) - 1.2;
      wend[i] = std::get<1>(qrange) + 1.2;

      qscores[i] = deconvoluted_spectrum[i].getQScore();
      charges[i] = deconvoluted_spectrum[i].getRepCharge();
      double massDelta = avg.getAverageMassDelta(deconvoluted_spectrum[i].getMonoMass());
      avg_masses[i] = massDelta + deconvoluted_spectrum[i].getMonoMass();
    }
    std::vector<PeakGroup> empty;
    deconvoluted_spectrum.swap(empty);
  }

  MSSpectrum FLASHIda::makeMSSpectrum(const double* mzs, const double* ints, const int length, const double rt, const int ms_level, const char* name)
  {
    auto spec = MSSpectrum();
    for (int i = 0; i < length; i++)
    {
      if (ints[i] <= 0) {
        continue;
      }
      spec.push_back(Peak1D(mzs[i], ints[i]));
    }
    spec.setMSLevel(ms_level);
    spec.setName(name);
    spec.setRT(rt);//
    return spec;
  }
}