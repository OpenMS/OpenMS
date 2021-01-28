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
      String token_string = std::string(token);
      double num = atof(token_string.c_str());

      if (num == 0 && !isdigit(token_string[token_string.size() - 1]))
      {
        key = token_string;
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
    fd_defaults.setValue("min_peaks", IntList{2, 1}); // more sensitive

    auto mass_count_double = inputs["max_mass_count"];

    for (double j : mass_count_double)
    {
      mass_count_.push_back((int)j);
    }

    fd_defaults.setValue("min_mass_count", mass_count_);

    fd_.setParameters(fd_defaults);
    fd_.calculateAveragine(false);

    rt_window_ = inputs["RT_window"][0];
    qscore_threshold_ = inputs["score_threshold"][0];

    averagine_ = fd_.getAveragine();
    std::cout << "QScore threshold: " << qscore_threshold_ << std::endl;
    std::cout << fd_defaults << std::endl;
  }

  int FLASHIda::getPeakGroups(const double* mzs,
                              const double* ints,
                              const int length,
                              const double rt,
                              const int ms_level,
                              const char* name)
  {

    auto spec = makeMSSpectrum_(mzs, ints, length, rt, ms_level, name);
    //deconvoluted_spectrum_ = DeconvolutedSpectrum(spec, 1);
    if (ms_level == 1)
    {
      //current_max_mass_ = max_mass;
      //currentChargeRange = chargeRange;
    }
    else
    {
      //TODO precursor infor here
    }

    std::vector<DeconvolutedSpectrum> tmp;
    deconvoluted_spectrum_ = fd_.getDeconvolutedSpectrum(spec, tmp, 0);

    // per spec deconvolution
    //    int specIndex = 0, massIndex = 0; // ..
    //fd_.getDeconvolutedSpectrum(deconvoluted_spectrum_, 0);
    FLASHIda::filterPeakGroupsUsingMassExclusion_(spec, ms_level);

    // spec.clear(true);

    return deconvoluted_spectrum_.size();
  }


  void FLASHIda::filterPeakGroupsUsingMassExclusion_(const MSSpectrum& spec, const int ms_level)
  {
    double rt = spec.getRT();
    std::sort(deconvoluted_spectrum_.begin(), deconvoluted_spectrum_.end(), QscoreComparator_);
    int mass_count = mass_count_[ms_level - 1];

    const auto color_order = std::vector<char>({'B', 'R', 'G', 'b', 'r'});

    std::unordered_map<int, std::vector<double>> new_mass_rt_qscore_map;
    std::unordered_map<int, char> new_color_map;

    for (auto& item : mass_rt_qscore_map_)
    {
      if (item.second[0] < rt - rt_window_)
      {
        continue;
      }
      new_mass_rt_qscore_map[item.first] = item.second;
    }

    for (auto& item : mass_color_map_)
    {
      if (new_mass_rt_qscore_map.find(item.first) == new_mass_rt_qscore_map.end())
      {
        continue;
      }
      new_color_map[item.first] = item.second;
    }

    for (auto& pg : deconvoluted_spectrum_) // now update color and new_mass_rt_qscore_map
    {
      int m = FLASHDeconvAlgorithm::getNominalMass(pg.getMonoMass());
      double qscore = pg.getQScore();

      if (new_mass_rt_qscore_map.find(m) == new_mass_rt_qscore_map.end()){ // new mass
        new_mass_rt_qscore_map[m] = std::vector<double>(2);
        new_mass_rt_qscore_map[m][1] = qscore;
        new_color_map[m] = 'G';
      }else if (new_mass_rt_qscore_map[m][1] < qscore) { // increasing mass
        new_mass_rt_qscore_map[m][1] = qscore;
        if(new_color_map.find(m) == new_color_map.end() || new_color_map[m] == 'G' || new_color_map[m] == 'R'){
          new_color_map[m] = 'B';
        }else if(new_color_map[m] == 'r'){
          new_color_map[m] = 'b';
        }
      }else{ // decreasing mass
        if(new_color_map.find(m) == new_color_map.end() || new_color_map[m] == 'G' || new_color_map[m] == 'B'){
          new_color_map[m] = 'R';
        }else if(new_color_map[m] == 'b'){
          new_color_map[m] = 'r';
        }
      }
      new_mass_rt_qscore_map[m][0] = rt;
    }

    new_mass_rt_qscore_map.swap(mass_rt_qscore_map_);
    std::unordered_map<int, std::vector<double>>().swap(new_mass_rt_qscore_map);

    new_color_map.swap(mass_color_map_);
    std::unordered_map<int, char>().swap(new_color_map);

    std::vector<PeakGroup> filtered_peakgroups;
    filtered_peakgroups.reserve(mass_count_.size());
    std::set<int> current_selected_masses; // current selected masses
    std::set<int> current_selected_mzs; // current selected mzs

    for (auto& c: color_order)
    {
      for (auto& pg : deconvoluted_spectrum_)
      {
        if (filtered_peakgroups.size() >= mass_count){
          break;
        }
        if (pg.getQScore() < qscore_threshold_)
        {
          break;
        }
        int mz = (int)round((std::get<0>(pg.getMaxQScoreMzRange()) + std::get<1>(pg.getMaxQScoreMzRange())) / 2.0);
        if (current_selected_mzs.find(mz) != current_selected_mzs.end()) {
          continue;
        }

        int nominal_mass = FLASHDeconvAlgorithm::getNominalMass(pg.getMonoMass());
        if (current_selected_masses.find(nominal_mass) != current_selected_masses.end()) {
          continue;
        }
        if(mass_color_map_.find(nominal_mass) == mass_color_map_.end() || mass_color_map_[nominal_mass] != c){
          continue;
        }

        char prev_color = mass_color_map_[nominal_mass];
        if(prev_color == 'B' && pg.getQScore() > 0){
          mass_color_map_[nominal_mass] = 'b';
        }else if (prev_color == 'R' && pg.getQScore() > 0){
          mass_color_map_[nominal_mass] = 'r';
        }

        filtered_peakgroups.push_back(pg);
        current_selected_masses.insert(nominal_mass);
        current_selected_mzs.insert(mz);
      }
    }

    deconvoluted_spectrum_.swap(filtered_peakgroups);
    std::vector<PeakGroup>().swap(filtered_peakgroups);
  }


  /*
  void FLASHIda::filterPeakGroupsUsingMassExclusion_(MSSpectrum& spec, int msLevel)
  {
    double rt = spec.getRT();
    std::sort(deconvoluted_spectrum_.begin(), deconvoluted_spectrum_.end(), QscoreComparator_);

    std::map<int, std::vector<double>> nall;
    std::map<int, std::vector<double>> nallMz;

    std::map<int, std::vector<double>> nselected;
    std::map<int, std::vector<double>> nselectedMz; // m/z * 20
    std::vector<PeakGroup> toConsider;
    toConsider.reserve(deconvoluted_spectrum_.size());
    std::vector<PeakGroup> toConsider2;
    toConsider2.reserve(deconvoluted_spectrum_.size());
    std::vector<PeakGroup> toConsider3;
    toConsider3.reserve(deconvoluted_spectrum_.size());

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

    for (auto& pg : deconvoluted_spectrum_)
    {
      if (pg.getQScore() < qscore_threshold_)
      {
        break;
      }
      auto massDelta = avg.getAverageMassDelta(pg.getMonoMass());
      auto m = FLASHDeconvAlgorithm::getNominalMass(pg.getMonoMass() + massDelta);
      auto mz = (int)round((std::get<0>(pg.getMzxQScoreMzRange()) + std::get<1>(pg.getMzxQScoreMzRange())) / 2.0);

      auto qscore_ = pg.getQScore();
      if ((nall.find(m) != nall.end() && nall[m][1] < qscore_)
          || (nallMz.find(mz) != nallMz.end() && nallMz[mz][1] < qscore_))
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

    for (auto& pg : deconvoluted_spectrum_)
    {
      auto massDelta = avg.getAverageMassDelta(pg.getMonoMass());
      auto m = FLASHDeconvAlgorithm::getNominalMass(pg.getMonoMass() + massDelta);
      auto mz = (int)round((std::get<0>(pg.getMzxQScoreMzRange()) + std::get<1>(pg.getMzxQScoreMzRange())) / 2.0);

      auto qscore_ = pg.getQScore();
      if (nall.find(m) == nall.end()) {
        nall[m] = std::vector<double>(2);
        nall[m][1] = qscore_;
      }
      else if (nall[m][1] < qscore_) {
        nall[m][1] = qscore_;
      }

      nall[m][0] = rt;

      if (nallMz.find(mz) == nallMz.end()) {
        nallMz[mz] = std::vector<double>(2);
        nallMz[mz][1] = qscore_;
      }
      else if (nallMz[mz][1] < qscore_) {
        nallMz[mz][1] = qscore_;
      }

      nallMz[mz][0] = rt;
    }

    nall.swap(all);
    std::map<int, std::vector<double>>().swap(nall);

    nallMz.swap(allMz);
    std::map<int, std::vector<double>>().swap(nallMz);

    std::vector<PeakGroup> filtered;
    filtered.reserve(mass_count_.size());
    std::set<int> cselected; // current selected masses
    std::set<int> cselectedMz; // current selected mzs

    int mc = mass_count_[msLevel - 1];

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
        auto qscore_ = pg.getQScore();
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
          nselectedMz[mz][1] = qscore_;
        }
        else if (nselectedMz[mz][1] < qscore_) {
          nselectedMz[mz][1] = qscore_;
        }

        nselectedMz[mz][0] = rt;

        if (nselected.find(m) == nselected.end()) {
          nselected[m] = std::vector<double>(2);
          nselected[m][1] = qscore_;
        }
        else if (nselected[m][1] < qscore_) {
          nselected[m][1] = qscore_;
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
        auto qscore_ = pg.getQScore();
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
          nselectedMz[mz][1] = qscore_;
        }
        else if (nselectedMz[mz][1] < qscore_) {
          nselectedMz[mz][1] = qscore_;
        }

        nselectedMz[mz][0] = rt;

        if (nselected.find(m) == nselected.end()) {
          nselected[m] = std::vector<double>(2);
          nselected[m][1] = qscore_;
        }
        else if (nselected[m][1] < qscore_) {
          nselected[m][1] = qscore_;
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
        auto qscore_ = pg.getQScore();

        cselected.insert(m);
        cselectedMz.insert(mz);
        if (nselectedMz.find(mz) == nselectedMz.end()) {
          nselectedMz[mz] = std::vector<double>(2);
          nselectedMz[mz][1] = qscore_;
        }
        else if (nselectedMz[mz][1] < qscore_) {
          nselectedMz[mz][1] = qscore_;
        }

        nselectedMz[mz][0] = rt;

        if (nselected.find(m) == nselected.end()) {
          nselected[m] = std::vector<double>(2);
          nselected[m][1] = qscore_;
        }
        else if (nselected[m][1] < qscore_) {
          nselected[m][1] = qscore_;
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
        auto qscore_ = pg.getQScore();

        cselected.insert(m);
        cselectedMz.insert(mz);

        if (nselectedMz.find(mz) == nselectedMz.end()) {
          nselectedMz[mz] = std::vector<double>(2);
          nselectedMz[mz][1] = qscore_;
        }
        else if (nselectedMz[mz][1] < qscore_) {
          nselectedMz[mz][1] = qscore_;
        }

        nselectedMz[mz][0] = rt;

        if (nselected.find(m) == nselected.end()) {
          nselected[m] = std::vector<double>(2);
          nselected[m][1] = qscore_;
        }
        else if (nselected[m][1] < qscore_) {
          nselected[m][1] = qscore_;
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

        auto mz = (int)round((std::get<0>(pg.getMzxQScoreMzRange()) + std::get<1>(pg.getMaxQScoreMzRange())) / 2.0);
        if (cselectedMz.find(mz) != cselectedMz.end()) {
          continue;
        }
        auto massDelta = avg.getAverageMassDelta(pg.getMonoMass());
        auto m = FLASHDeconvAlgorithm::getNominalMass(pg.getMonoMass() + massDelta);
        if (cselected.find(m) != cselected.end()) {
          continue;
        }
        auto qscore_ = pg.getQScore();

        cselected.insert(m);
        cselectedMz.insert(mz);
        if (nselectedMz.find(mz) == nselectedMz.end()) {
          nselectedMz[mz] = std::vector<double>(2);
          nselectedMz[mz][1] = qscore_;
        }
        else if (nselectedMz[mz][1] < qscore_) {
          nselectedMz[mz][1] = qscore_;
        }

        nselectedMz[mz][0] = rt;

        if (nselected.find(m) == nselected.end()) {
          nselected[m] = std::vector<double>(2);
          nselected[m][1] = qscore_;
        }
        else if (nselected[m][1] < qscore_) {
          nselected[m][1] = qscore_;
        }

        nselected[m][0] = rt;

        filtered.push_back(pg);
      }
    }

    nselected.swap(selected);
    std::map<int, std::vector<double>>().swap(nselected);

    nselectedMz.swap(selectedMz);
    std::map<int, std::vector<double>>().swap(nselectedMz);

    deconvoluted_spectrum_.swap(filtered);
    std::vector<PeakGroup>().swap(filtered);
  }
*/
  void FLASHIda::getIsolationWindows(double* window_start, double* window_end, double* qscores, int* charges, double* avg_masses)
  {
    std::sort(deconvoluted_spectrum_.begin(), deconvoluted_spectrum_.end(), QscoreComparator_);

    for (int i = 0; i < deconvoluted_spectrum_.size(); i++)
    {
      auto mz_range = deconvoluted_spectrum_[i].getMaxQScoreMzRange();
      window_start[i] = std::get<0>(mz_range) - min_isolation_window_half_;
      window_end[i] = std::get<1>(mz_range) + min_isolation_window_half_;

      qscores[i] = deconvoluted_spectrum_[i].getQScore();
      charges[i] = deconvoluted_spectrum_[i].getRepCharge();
      double mass_diff = averagine_.getAverageMassDelta(deconvoluted_spectrum_[i].getMonoMass());
      avg_masses[i] = mass_diff + deconvoluted_spectrum_[i].getMonoMass();
    }
    std::vector<PeakGroup> empty;
    deconvoluted_spectrum_.swap(empty);
  }

  MSSpectrum FLASHIda::makeMSSpectrum_(const double* mzs, const double* ints, const int length, const double rt, const int ms_level, const char* name)
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
