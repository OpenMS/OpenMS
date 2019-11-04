//
// Created by Kyowon Jeong on 9/25/19.
//

#include "OpenMS/ANALYSIS/TOPDOWN/MassFeatureTrace.h"
namespace OpenMS
{

  void MassFeatureTrace::findFeatures(std::vector<PeakGroup> &peakGroups,
                    int &featureCntr,
                    std::fstream &fsf,
                    PrecalcularedAveragine &averagines,
                                      Param &mtd_param,
                                      Parameter &param)
  {

    MSExperiment map;
    boost::unordered_map<float, PeakGroup> *peakGroupMap;
    // boost::unordered_map<float, MSSpectrum> rtOrignalSpecMap;
    boost::unordered_map<float, int> rtSpecMap;
    std::map<int, MSSpectrum> indexSpecMap;

    int maxSpecIndex = 0;

    for (auto &pg : peakGroups)
    {
      auto &spec = pg.spec;
      if(spec->getMSLevel() != 1){
        continue;
      }

      if (indexSpecMap.find(pg.specIndex) == indexSpecMap.end()){
        indexSpecMap[pg.specIndex] = MSSpectrum();
      }
      auto &deconvSpec = indexSpecMap[pg.specIndex];
      rtSpecMap[spec->getRT()] = pg.specIndex;
      maxSpecIndex = maxSpecIndex > pg.specIndex ? maxSpecIndex : pg.specIndex ;

      deconvSpec.setRT(spec->getRT());
      Peak1D tp(pg.monoisotopicMass, (float) pg.intensity);
      deconvSpec.push_back(tp);
    }

    //int tmp = 0;
    for(auto iter = indexSpecMap.begin(); iter != indexSpecMap.end(); ++iter)
    {
      //tmp+=(iter->second).size();
      map.addSpectrum(iter->second);
    }

    if(map.size() < 3){
      return;
    }
    //std::cout<<map.size()<< " " <<tmp<< std::endl;
    peakGroupMap = new boost::unordered_map<float, PeakGroup>[maxSpecIndex + 1];

    for (auto &pg : peakGroups)
    {
      auto &spec = pg.spec;
      if(spec->getMSLevel() != 1){
        continue;
      }
      auto &pgMap = peakGroupMap[pg.specIndex];

      pgMap[pg.monoisotopicMass] = pg;
      //std::cout<<pg.monoisotopicMass<< " " << pg.specIndex << std::endl;

    }

    for (auto it = map.begin(); it != map.end(); ++it)
    {
      it->sortByPosition();
      // cout<<it->size()<<endl;
    }

    MassTraceDetection mtdet;

    //mtd_param.setValue("mass_error_da", .3,// * (param.chargeRange+ param.minCharge),
    //                   "Allowed mass deviation (in da).");
    mtd_param.setValue("mass_error_ppm", param.tolerance * 1e6 * 2, "");
    mtd_param.setValue("trace_termination_criterion", "outlier", "");

    mtd_param.setValue("reestimate_mt_sd", "true", "");
    mtd_param.setValue("quant_method", "area", "");
    mtd_param.setValue("noise_threshold_int", .0, "");

    //double rtDuration = (map[map.size() - 1].getRT() - map[0].getRT()) / ms1Cntr;
    mtd_param.setValue("min_sample_rate", 0.01, "");
    mtd_param.setValue("trace_termination_outliers", param.numOverlappedScans, "");
    mtd_param.setValue("min_trace_length", param.minRTSpan, "");
    //mtd_param.setValue("max_trace_length", 1000.0, "");
    mtdet.setParameters(mtd_param);

    std::vector<MassTrace> m_traces;

    mtdet.run(map, m_traces);  // m_traces : output of this function

    auto *perChargeIntensity = new double[param.chargeRange + param.minCharge + 1];
    auto *perChargeMaxIntensity = new double[param.chargeRange + param.minCharge + 1];
    auto *perChargeMz = new double[param.chargeRange + param.minCharge + 1];
    auto *perIsotopeIntensity = new double[param.maxIsotopeCount];

    for (auto &mt : m_traces)
    {
      int minCharge = param.chargeRange + param.minCharge + 1;
      int maxCharge = 0;
      boost::dynamic_bitset<> charges(param.chargeRange + param.minCharge + 1);
      std::fill_n(perChargeIntensity, param.chargeRange + param.minCharge + 1, 0);
      std::fill_n(perChargeMaxIntensity, param.chargeRange + param.minCharge + 1, 0);
      std::fill_n(perChargeMz, param.chargeRange + param.minCharge + 1, 0);
      std::fill_n(perIsotopeIntensity, param.maxIsotopeCount, 0);
      double massDiff = 0;
      double max_intensity = -1;

      for (auto &p2 : mt)
      {
       // std::cout << p2.getRT() << " " << p2.getMZ() << std::endl;
        int specIndex = rtSpecMap[(float) p2.getRT()];
        auto &pgMap = peakGroupMap[specIndex];
        auto &pg = pgMap[(float) p2.getMZ()];
        minCharge = minCharge < pg.minCharge? minCharge : pg.minCharge;
        maxCharge = maxCharge > pg.maxCharge? maxCharge : pg.maxCharge;

        //std::cout<<1<<std::endl;
        if (pg.intensity > max_intensity)
        {
          max_intensity = pg.intensity;
          massDiff = pg.avgMass - pg.monoisotopicMass;
        }
        //std::cout<<2<<std::endl;
        for (auto &p : pg.peaks)
        {
          //std::cout<<1<<std::endl;

          if (p.isotopeIndex < 0 || p.isotopeIndex >= param.maxIsotopeCount || p.charge < 0 ||
              p.charge >= param.chargeRange + param.minCharge + 1)
          {
            continue;
          }

          charges[p.charge] = true;
          perChargeIntensity[p.charge] += p.intensity;
          perIsotopeIntensity[p.isotopeIndex] += p.intensity;
          if (perChargeMaxIntensity[p.charge] > p.intensity)
          {
            continue;
          }
          perChargeMaxIntensity[p.charge] = p.intensity;
          perChargeMz[p.charge] = p.mz;
         // std::cout<<3<<std::endl;
        }
        //std::cout<<3<<std::endl;
      }

      if (massDiff <= 0)
      {
        continue;
      }

      double chargeScore = PeakGroupScoring::getChargeFitScore(perChargeIntensity, param.minCharge + param.chargeRange + 1);
      if (chargeScore < param.minChargeCosine) //
      {
        continue;
      }

      int offset = 0;
      double mass = mt.getCentroidMZ();
      double isoScore = PeakGroupScoring::getIsotopeCosineAndDetermineIsotopeIndex(mass,
                                                                                       perIsotopeIntensity,
                                                                                       param.maxIsotopeCount,
                                                                                        offset,averagines);
      if (isoScore < param.minIsotopeCosine)
      {
        continue;
      }

      if (offset != 0)
      {
        mass += offset * Constants::C13C12_MASSDIFF_U;
        //avgMass += offset * Constants::C13C12_MASSDIFF_U;
        //p.isotopeIndex -= offset;
      }
      //auto mass = mt.getCentroidMZ();
      fsf << ++featureCntr << "\t" << param.fileName << "\t" << std::to_string(mass) << "\t"
          << std::to_string(mass + massDiff) << "\t"
          << mt.getSize() << "\t"
          //fsf << ++featureCntr << "\t" << param.fileName << "\t" << mass << "\t"
          //<< getNominalMass(mass) << "\t"
          << mt.begin()->getRT() << "\t"
          << mt.rbegin()->getRT() << "\t"
          << mt.getTraceLength() << "\t"
          << mt[mt.findMaxByIntPeak()].getRT() << "\t"
          << mt.getMaxIntensity(false) << "\t"
          // << mt.computePeakArea() << "\t"
          << minCharge << "\t"
          << maxCharge << "\t"
          << charges.count() << "\t"
          << isoScore << "\t"
          << chargeScore << "\n";
    }
    delete[] perIsotopeIntensity;
    delete[] perChargeMz;
    delete[] perChargeMaxIntensity;
    delete[] perChargeIntensity;
    delete[] peakGroupMap;
  }
}
