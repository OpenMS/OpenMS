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
//
// Created by Kyowon Jeong on 4/22/20.
//

#include "include/OpenMS/ANALYSIS/TOPDOWN/DeconvolutedSpectrum.h"

namespace OpenMS
{
  DeconvolutedSpectrum::DeconvolutedSpectrum()
  {
  }

  DeconvolutedSpectrum::DeconvolutedSpectrum(MSSpectrum &s, int n) :
      spec(&s), scanNumber(n)
  {
  }

  DeconvolutedSpectrum::~DeconvolutedSpectrum()
  {
    //delete peakGroups;
    //std::vector<LogMzPeak>().swap(peaks);
    std::vector<PeakGroup>().swap(peakGroups);
    //std::unordered_map<int, PeakGroup>().swap(peakGroupMap);
  }

  bool DeconvolutedSpectrum::empty() const
  {
    return peakGroups.empty();
  }

  //void DeconvolutedSpectrum::clearChargeSNRMap(){
  //  for(auto &pg : peakGroups){
  //    pg.clearChargeInfo();
  //  }
  //}
  /*void DeconvolutedSpectrum::updatePeakGroupMap()
  {
    if (!peakGroupMap.empty())
    {
      return;
    }
    //std::cout<< peaks.size()<<std::endl;

    for (auto &pg : peakGroups)
    {
      for (auto &p : pg.peaks)
      {
        auto position = p.index;
        if (peakGroupMap.find(position) == peakGroupMap.end())
        {
          peakGroupMap[position] = pg;
        }
        else
        {
          if (peakGroupMap[position].isotopeCosineScore < pg.isotopeCosineScore)
          {
            peakGroupMap[position] = pg;
          }
          continue;
        }
      }
    }
    //std::cout<<"2"<<std::endl;
  }*/

  MSSpectrum DeconvolutedSpectrum::toSpectrum(){
    auto outSpec = MSSpectrum(*spec);
    outSpec.clear(false);
    for (auto &pg : peakGroups)
    {
      if (pg.peaks.empty())
      {
        continue;
      }
      auto p = Peak1D(pg.monoisotopicMass, pg.intensity);
      outSpec.push_back(p);
    }

    return outSpec;
  }


  void DeconvolutedSpectrum::writeDeconvolutedMasses(std::fstream &fs,
                                                     FLASHDeconvHelperStructs::Parameter &param)//, fstream &fsm, fstream &fsp)
  {
    if (empty())
    {
      return;
    }

    for (auto &pg : peakGroups)
    {
      if (pg.peaks.empty())
      {
        continue;
      }
      double &m = pg.monoisotopicMass;
      double &am = pg.avgMass;
      double &intensity = pg.intensity;
      int minCharge = param.chargeRange + param.minCharge;
      int maxCharge = -1;
      for (auto &p : pg.peaks)
      {
        minCharge = minCharge < p.charge ? minCharge : p.charge;
        maxCharge = maxCharge > p.charge ? maxCharge : p.charge;
      }

      //"MassIndex\tSpecIndex\tScanNum\tRetentionTime\tMassCountInSpec\tAvgMass\tMonoMass\t"
      //               "SumIntensity\tMinCharge\tMaxCharge\t"
      //               "PeakCount
      //               \tPeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\tPeakIntensities\t"

      //               "IsotopeCosine\tChargeIntensityCosine\tMassSNR\tMaxQScoreCharge\tMaxQScoreMzStart\tMaxQScoreMzEnd\tQScore\n";
      //
      // "MassIndex\tSpecIndex\tScanNum\tRetentionTime\tMassCountInSpec\tAvgMass\tMonoMass\t"
      //               "SumIntensity\tMinCharge\tMaxCharge\t"
      //               "PeakCount
      //               \tPeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\tPeakIntensities\t"

      //               "PrecursorSpecIndex\tPrecursorMz\tPrecursorCharge\tPrecursorMonoMass\tPrecursorIntensity\tPrecursorQScore\t"
      //               "IsotopeCosine\tMassSNR\tMaxQScoreCharge\tMaxQScoreMzStart\tMaxQScoreMzEnd\tQScore\n";
      //

      fs << pg.massIndex << "\t" << pg.specIndex << "\t" << param.fileName << "\t" << pg.scanNumber << "\t"
         << std::to_string(spec->getRT())<< "\t"
         << pg.peaks.size() << "\t"
         << std::to_string(am) << "\t" << std::to_string(m) << "\t" << intensity << "\t"
         << minCharge << "\t" << maxCharge << "\t"
         << pg.peaks.size() << "\t";

      if (param.writeDetail)
      {
        fs << std::fixed << std::setprecision(2);
        for (auto &p : pg.peaks)
        {
          fs << p.mz << ";";
        }
        fs << "\t";

        fs << std::fixed << std::setprecision(1);
        for (auto &p : pg.peaks)
        {
          fs << p.intensity << ";";
        }
        fs << "\t";
        fs <<std::setprecision(-1);

        for (auto &p : pg.peaks)
        {
          fs << p.charge << ";";
        }
        fs << "\t";
        for (auto &p : pg.peaks)
        {
          fs << p.getUnchargedMass(param.chargeMass) << ";";
        }
        fs << "\t";
        for (auto &p : pg.peaks)
        {
          fs << p.isotopeIndex << ";";
        }
        fs << "\t";

        for (auto &p : pg.peaks)
        {
          auto tm = pg.monoisotopicMass + p.isotopeIndex * Constants::ISOTOPE_MASSDIFF_55K_U;
          auto diff = (tm / p.charge + param.chargeMass - p.mz) / p.mz;

          fs << 1e6 * diff << ";";
        }
        fs << "\t";


      }

      //               "IsotopeCosine\tChargeIntensityCosine\tMassSNR\tMaxQScoreCharge\tMaxQScoreMzStart\tMaxQScoreMzEnd\tQScore\n";
      //               "PrecursorSpecIndex\tPrecursorMz\tPrecursorCharge\tPrecursorMonoMass\tPrecursorIntensity\tPrecursorQScore\t"
      //               "IsotopeCosine\tMassSNR\tMaxQScoreCharge\tMaxQScoreMzStart\tMaxQScoreMzEnd\tQScore\n";

      if (spec->getMSLevel() > 1)
      {
        //PrecursorScanNum	PrecursorMz	PrecursorIntensity PrecursorCharge	PrecursorMonoMass		PrecursorQScore
        fs << precursorScanNumber<< "\t" << std::to_string(precursorPeak.getMZ()) << "\t" << precursorPeak.getIntensity() << "\t"
           << precursorPeak.getCharge()
           << "\t";

        if (precursorPeakGroup == nullptr){
          fs <<"N/A\tN/A\t";
        }else{
          fs << std::to_string(precursorPeakGroup->monoisotopicMass) << "\t"
              << precursorPeakGroup->qScore<<"\t";
        }
      }
      fs << pg.isotopeCosineScore<< "\t" ;

      if (spec->getMSLevel() == 1)
      {
        fs << pg.chargeCosineScore << "\t";
      }
      fs << pg.totalSNR << "\t"
       << pg.maxQScoreCharge << "\t" << std::to_string(pg.maxQScoreMzStart) << "\t" << std::to_string(pg.maxQScoreMzEnd) << "\t"
       << pg.qScore << "\n" << std::setprecision(-1); //
    }
  }


  void DeconvolutedSpectrum::writeDeconvolutedMassesHeader(std::fstream &fs, int &n, bool detail)
  {
    if (detail)
    {
      if (n == 1)
      {
        fs
            << "MassIndex\tSpecIndex\tFileName\tScanNum\tRetentionTime\tMassCountInSpec\tAvgMass\tMonoMass\t"
               "SumIntensity\tMinCharge\tMaxCharge\t"
               "PeakCount\tPeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
               "IsotopeCosine\tChargeIntensityCosine\tMassSNR\tMaxQScoreCharge\tMaxQScoreMzStart\tMaxQScoreMzEnd\tQScore\n";
      }
      else
      {
        fs
            << "MassIndex\tSpecIndex\tFileName\tScanNum\tRetentionTime\tMassCountInSpec\tAvgMass\tMonoMass\t"
               "SumIntensity\tMinCharge\tMaxCharge\t"
               "PeakCount\tPeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
               "PrecursorScanNum\tPrecursorMz\tPrecursorIntensity\tPrecursorCharge\tPrecursorMonoMass\tPrecursorQScore\t"
               "IsotopeCosine\tMassSNR\tMaxQScoreCharge\tMaxQScoreMzStart\tMaxQScoreMzEnd\tQScore\n";
      }
    }
    else
    {
      if (n == 1)
      {
        fs
            << "MassIndex\tSpecIndex\tFileName\tScanNum\tRetentionTime\tMassCountInSpec\tAvgMass\tMonoMass\t"
               "SumIntensity\tMinCharge\tMaxCharge\t"
               "PeakCount\t"
               //"PeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\tPeakIntensities\t"
               "IsotopeCosine\tChargeIntensityCosine\tMassSNR\tMaxQScoreCharge\tMaxQScoreMzStart\tMaxQScoreMzEnd\tQScore\n";

      }
      else
      {
        fs
            << "MassIndex\tSpecIndex\tFileName\tScanNum\tRetentionTime\tMassCountInSpec\tAvgMass\tMonoMass\t"
               "SumIntensity\tMinCharge\tMaxCharge\t"
               "PeakCount\t"
               //"PeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\tPeakIntensities\t"
               "PrecursorScanNum\tPrecursorMz\tPrecursorIntensity\tPrecursorCharge\tPrecursorMonoMass\tPrecursorQScore\t"
               "IsotopeCosine\tMassSNR\tMaxQScoreCharge\tMaxQScoreMzStart\tMaxQScoreMzEnd\tQScore\n";


      }

    }
    //pg.maxSNRcharge << "\t" << pg.maxSNR << "\t" << pg.maxSNRminMz << "\t" << pg.maxSNRmaxMz
    //MinScan	MaxScan	 RepScan	RepCharge	RepMz ApexScanNum    Envelope

  }

  void DeconvolutedSpectrum::writeAttCsvHeader(std::fstream &fs){
    fs<<"ScanNumber,RetentionTime,PrecursorScanNumber,Charge,ChargeSNR,PeakIntensity,EnvIntensity,EnvIsotopeCosine,PeakMz,"
        "MonoMass,MassSNR,IsotopeCosine,ChargeIntensityCosine,MassIntensity,QScore,Class\n";
  }

  void DeconvolutedSpectrum::writeThermoInclusionHeader(std::fstream &fs){
    fs<<"Compound,Formula,Adduct,m/z,z,t start (min),t stop (min),Isolation Window (m/z),Normalized AGC Target (%)\n";
  }


  void DeconvolutedSpectrum::clearPeakGroupsChargeInfo(){
    for (auto& pg: peakGroups)
    {
      pg.clearChargeInfo();
    }
  }

  void DeconvolutedSpectrum::writeAttCsv(std::fstream &fs, int msLevel, double qScoreThreshold = -1000, int numMaxMS2 = -1){

    if (msLevel>1)
    {
      if(precursorPeakGroup != nullptr)
      {
        /*auto avgDiff = .0;

        for (auto &p : precursorPeakGroup->peaks)
        {
          auto tm = precursorPeakGroup->monoisotopicMass + p.isotopeIndex * Constants::ISOTOPE_MASSDIFF_55K_U;
          avgDiff += (tm / p.charge + Constants::PROTON_MASS_U - p.mz) / p.mz;
        }
        avgDiff /= precursorPeakGroup->peaks.size();
        avgDiff *= 1e6;*/

        fs << scanNumber << "," << precursorPeakGroup->spec->getRT() << "," << precursorPeakGroup->scanNumber << "," << precursorPeak.getCharge()
           << "," << log10(precursorPeakGroup->perChargeInfo[precursorPeak.getCharge()][0] + 1e-3) << ","
           << log10(precursorPeak.getIntensity() + 1)
            << "," <<log10(precursorPeakGroup->perChargeInfo[precursorPeak.getCharge()][2] + 1)
            << "," << precursorPeakGroup->perChargeInfo[precursorPeak.getCharge()][1]
            << "," << precursorPeak.getMZ()
           << "," << precursorPeakGroup->monoisotopicMass << "," << log10(precursorPeakGroup->totalSNR + 1e-3) << ","
           << precursorPeakGroup->isotopeCosineScore
           << "," << precursorPeakGroup->chargeCosineScore << "," << log10(precursorPeakGroup->intensity + 1)
           << "," << precursorPeakGroup->qScore
           << ",f\n";
      }else{
          fs << scanNumber << ","<<spec->getRT()<<"," << 0 << "," << 0
             << ",?,"
             << log10(precursorPeak.getIntensity() + 1)
             << ",?,?," << (precursorPeak.getMZ())
             << "," << 0 << ",?,?,?,?"
             << "," << precursorPeakGroup->qScore
             << ",f\n";
      }
    }else{
      //double scoreThreshold = 0;
      std::vector<double> scores;

      if(numMaxMS2>0 && peakGroups.size() > (Size)numMaxMS2)// max peak count for TopPic
      {
        scores.reserve(peakGroups.size());
        for (auto &pg : peakGroups)
        {
          scores.push_back(pg.qScore);
        }
        std::sort(scores.begin(), scores.end());
        qScoreThreshold = std::max(qScoreThreshold, scores[scores.size() - numMaxMS2]);
        std::vector<double>().swap(scores);
      }

      int size = 0;
      for (auto &pg : peakGroups)
      {
        if(pg.qScore < qScoreThreshold){
          continue;
        }
        if (size > numMaxMS2){
          break;
        }

        LogMzPeak *peak = nullptr;
        double maxIntensity = 0;
        for (auto &p : pg.peaks)
        {
          if (p.charge != pg.maxQScoreCharge){
            continue;
          }

          if (p.mz < pg.maxQScoreMzStart){
            continue;
          }
          if (p.mz > pg.maxQScoreMzEnd){
            break;
          }

          if(p.intensity < maxIntensity){
            continue;
          }
          maxIntensity = p.intensity;
          peak = &p;
        }
        if (peak == nullptr){
          continue;
        }

        size++;
        fs << scanNumber << "," << spec->getRT() << ",0," << pg.maxQScoreCharge
           << "," << log10(pg.perChargeInfo[pg.maxQScoreCharge][0] + 1e-3) //<< "," << log10(peak->intensity + 1)
           << "," << log10(peak->intensity + 1)
            << "," <<log10(pg.perChargeInfo[pg.maxQScoreCharge][2] + 1)
            << "," << pg.perChargeInfo[pg.maxQScoreCharge][1]
            << "," << peak->mz
          << "," << pg.monoisotopicMass << "," << log10(pg.totalSNR+1e-3) << "," << pg.isotopeCosineScore
          << "," << pg.chargeCosineScore << "," << log10(pg.intensity+1)
          << "," << pg.qScore
          <<  ",f\n";
      }
    }
  }


  void DeconvolutedSpectrum::writeMassList(std::fstream &fs, double retDelta, double qScoreThreshold = -1000, int numMaxMS2 = -1){

    static std::map<double, std::map<int, double>> selected; // rt, int mass, qscore
    static double prevRT = -1;
    const double retWindow1 = 10;
    const double retWindow2 = 60;

    double rt = spec->getRT();
    if(rt < prevRT){
      std::map<double, std::map<int, double>>().swap(selected);
    }
    prevRT = rt;

    std::map<double, std::map<int, double>> nselected;
    std::set<int> toExclude;
    std::map<int, double> allMasses;
    for(auto& item : selected){
      if(item.first < rt - retWindow2){
        continue;
      }
      nselected[item.first] = item.second;

      for (auto &i : item.second)
      {
        allMasses[i.first] = i.second;
        if(item.first >= rt - retWindow1)
        {
          toExclude.insert(i.first);
        }
      }
    }
    //std::cout<<toExclude.size()<<std::endl;
    nselected.swap(selected);
    std::map<double, std::map<int, double>>().swap(nselected);

    //double scoreThreshold = 0;
    std::vector<double> scores;

    if(numMaxMS2>0 && peakGroups.size() > (Size)numMaxMS2)// max peak count
    {
      scores.reserve(peakGroups.size());
      for (auto &pg : peakGroups)
      {
        int nmass= (int) (pg.monoisotopicMass * 0.999497 + .5); // TODO optimize!!
        if(toExclude.find(nmass) != toExclude.end()){
          continue;
        }

        if(allMasses.find(nmass) != allMasses.end()){
          if(allMasses[nmass] > pg.intensity){
            continue;
          }
        }
        scores.push_back(pg.qScore);
      }
      std::sort(scores.begin(), scores.end());
      qScoreThreshold = std::max(qScoreThreshold, scores[scores.size() - numMaxMS2]);
      std::vector<double>().swap(scores);
    }

    int size = 0;
    std::map<int, double> tselected;

    for (auto &pg : peakGroups)
    {
      if(pg.qScore < qScoreThreshold){
        continue;
      }
      if (size > numMaxMS2){
        break;
      }
      int nmass= (int) (pg.monoisotopicMass * 0.999497 + .5);

      if(numMaxMS2>0 && peakGroups.size() > (Size)numMaxMS2){
        if(toExclude.find(nmass) != toExclude.end()){
          continue;
        }

        if(allMasses.find(nmass) != allMasses.end()){
          if(allMasses[nmass] > pg.intensity){
            continue;
          }
        }
      }

      LogMzPeak *peak = nullptr;
      double maxIntensity = 0;
      for (auto &p : pg.peaks)
      {
        if (p.charge != pg.maxQScoreCharge){
          continue;
        }

        if (p.mz < pg.maxQScoreMzStart){
          continue;
        }
        if (p.mz > pg.maxQScoreMzEnd){
          break;
        }

        if(p.intensity < maxIntensity){
          continue;
        }
        maxIntensity = p.intensity;
        peak = &p;
      }
      if (peak == nullptr){
        continue;
      }

      tselected[nmass] = pg.intensity;
      size++;
      auto tmz = (pg.maxQScoreMzStart + pg.maxQScoreMzEnd)/2.0;
      auto rtinMin = pg.spec->getRT()/60.0;
      auto z = pg.maxQScoreCharge;
      auto iso = (pg.maxQScoreMzEnd - pg.maxQScoreMzStart)/2.0 + .2;
      iso = iso < .4? .4 : iso; // TODO

      fs << pg.qScore  //"spec"<<scanNumber<<"_mz"<<tmz<<"_z"<<z<<"_m"<<pg.monoisotopicMass
      <<","<<pg.monoisotopicMass<<"," << maxIntensity << ","<<std::to_string(tmz)<<","<<z<<","<<rtinMin<<","<<(rtinMin + retDelta/60.0)<< ","<<iso<<",400\n";
    }
    selected[rt] = tselected;
  }

  void DeconvolutedSpectrum::writeTopFD(std::fstream &fs, int id)//, fstream &fsm, fstream &fsp)
  {
    if (spec->getMSLevel() == 1)
    {
      return;
    }
    fs << std::fixed << std::setprecision(2);
    fs << "BEGIN IONS\n"
       << "ID=" << id << "\n"
       << "SCANS=" << scanNumber << "\n"
       << "RETENTION_TIME=" << spec->getRT() << "\n";

      fs << "ACTIVATION=" << activationMethod << "\n";

      if(precursorPeakGroup != nullptr)
      {
        fs << "MS_ONE_ID=" << precursorPeakGroup->specIndex << "\n"
           << "MS_ONE_SCAN=" << precursorPeakGroup->scanNumber << "\n"
           //<< "PRECURSOR_MZ_ori="
           //<< std::to_string(spec->getPrecursors()[0].getMZ()) << "\n"
           //<< "PRECURSOR_CHARGE_ori=" << spec->getPrecursors()[0].getCharge() << "\n"
           << "PRECURSOR_MZ="
           << std::to_string(precursorPeak.getMZ()) << "\n"
           << "PRECURSOR_CHARGE=" << precursorPeak.getCharge() << "\n"
           << "PRECURSOR_MASS=" << std::to_string(precursorPeakGroup->monoisotopicMass) << "\n"
           << "PRECURSOR_INTENSITY=" << precursorPeak.getIntensity() << "\n";
      }else{
        fs << "MS_ONE_ID=" << 0 << "\n"
           << "MS_ONE_SCAN=" << 0 << "\n"
           << "PRECURSOR_MZ="
           << std::to_string(precursorPeak.getMZ()) << "\n"
           << "PRECURSOR_CHARGE=" << precursorPeak.getCharge() << "\n"
           //<< "PRECURSOR_MZ="
           //<< std::to_string(precursorPeak->mz) << "\n"
           //<< "PRECURSOR_CHARGE=" << precursorPeak->charge << "\n"
           << "PRECURSOR_MASS=" << std::to_string(precursorPeak.getMZ() *  precursorPeak.getCharge()) << "\n"
           << "PRECURSOR_INTENSITY=" << precursorPeak.getIntensity() << "\n";
      }
    fs << std::setprecision(-1);


    double scoreThreshold = 0;
    std::vector<double> scores;

    if(peakGroups.size()>500)// max peak count for TopPic
    {
      scores.reserve(peakGroups.size());
      for (auto &pg : peakGroups)
      {
        scores.push_back(pg.isotopeCosineScore);
      }
      std::sort(scores.begin(), scores.end());
      scoreThreshold = scores[scores.size() - 500];
      std::vector<double>().swap(scores);
    }

    int size = 0;
    for (auto &pg : peakGroups)
    {
      if (pg.isotopeCosineScore < scoreThreshold){
        continue;
      }
      if (size >= 500){
        break;
      }

      size++;
      fs << std::fixed << std::setprecision(2);
      fs << std::to_string(pg.monoisotopicMass) << "\t" << pg.intensity << "\t" << pg.maxQScoreCharge
         //  << "\t" << log10(pg.precursorSNR+1e-10) << "\t" << log10(pg.precursorTotalSNR+1e-10)
         //  << "\t" << log10(pg.maxSNR + 1e-10) << "\t" << log10(pg.totalSNR + 1e-10)
         << "\n";
      fs << std::setprecision(-1);
    }

    fs << "END IONS\n\n";
  }

  bool DeconvolutedSpectrum::registerPrecursor(DeconvolutedSpectrum &precursorSpectrum)
  {
    //precursorSpectrum.updatePeakGroupMap();
    for (auto &p: spec->getPrecursors())
    {
      for (auto &act :  p.getActivationMethods())
      {
        activationMethod = Precursor::NamesOfActivationMethodShort[act];
        break;
      }
      precursorPeak = p;
      precursorScanNumber = precursorSpectrum.scanNumber;
      auto startMz = p.getIsolationWindowLowerOffset() > 100.0 ?
                     p.getIsolationWindowLowerOffset() :
                     -p.getIsolationWindowLowerOffset() + p.getMZ();
      auto endMz = p.getIsolationWindowUpperOffset() > 100.0 ?
                   p.getIsolationWindowUpperOffset() :
                   p.getIsolationWindowUpperOffset() + p.getMZ();

      double maxSumIntensity = 0.0;
      for (auto &pg: (precursorSpectrum.peakGroups))
      {
        if (pg.peaks[0].mz > endMz || pg.peaks[pg.peaks.size() - 1].mz < startMz)
        {
          continue;
        }

        double sumIntensity = .0;
        double maxIntensity = .0;
        LogMzPeak *tmp = nullptr;
        for (auto &pt:pg.peaks)
        {
          if (pt.mz < startMz)
          {
            continue;
          }
          if (pt.mz > endMz)
          {
            break;
          }
          sumIntensity += pt.intensity;

          if(pt.intensity < maxIntensity){
            continue;
          }
          maxIntensity = pt.intensity;
          precursorPeak.setCharge(pt.charge);
          tmp = &pt;
        }

        if (sumIntensity <= maxSumIntensity || tmp == nullptr)
        {
          continue;
        }

        maxSumIntensity = sumIntensity;
        precursorPeakGroup = &pg;
      }

      /*if (precursorPeak == nullptr)
      {
        double minDistance = 10000.0;

        auto position =
            std::lower_bound(precursorSpectrum.spec->begin(), precursorSpectrum.spec->end(), LogMzPeak(startMz)) -
                precursorSpectrum.spec->begin();

        for (; position < (int) precursorSpectrum.spec->size(); position++)
        {
          auto mz = (*precursorSpectrum.spec)[position].getMZ();
          if (mz < startMz)
          {
            continue;
          }
          if (mz > endMz)
          {
            break;
          }

          auto distance = abs(p.getMZ() - mz);
          if (distance > minDistance)
          {
            continue;
          }

          minDistance = distance;
          precursorPeak = &((*precursorSpectrum.spec)[position]);
        }
      }*/
    }

     return precursorPeakGroup != nullptr;
   }



 }