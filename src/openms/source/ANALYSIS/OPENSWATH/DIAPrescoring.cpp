// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DIAPrescoring.h>

//#include <OpenMS/OPENSWATHALGO/DATAACCESS/SpectrumHelpers.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionHelper.h>
#include <OpenMS/OPENSWATHALGO/ALGO/StatsHelpers.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <boost/lexical_cast.hpp>

#include <iostream>
#include <algorithm>
#include <utility>

namespace OpenMS
{

  void getNormalizedLibraryIntensities(
    const std::vector<OpenSwath::LightTransition>& transitions,
    std::vector<double>& normalizedLibraryIntensities //normalized intensities
    )
  {
    double totalInt = 0.;
    for (std::size_t i = 0; i < transitions.size(); ++i)
    {
      double libInt = transitions[i].getLibraryIntensity();
      if (libInt < 0.)
        libInt = 0.;
      totalInt += libInt;
      normalizedLibraryIntensities.push_back(libInt);
    }
    std::transform(normalizedLibraryIntensities.begin(),
                   normalizedLibraryIntensities.end(),
                   normalizedLibraryIntensities.begin(),
                   [totalInt](auto && PH1) { return std::divides<double>()(std::forward<decltype(PH1)>(PH1), totalInt); });
  }

  void getMZIntensityFromTransition(const std::vector<OpenSwath::LightTransition>& trans,
                                    std::vector<std::pair<double, double> >& res)
  {
    for (std::size_t i = 0; i < trans.size(); ++i)
    {
      res.emplace_back(trans[i].product_mz, trans[i].library_intensity);
    }
  }

  void DiaPrescore::operator()(const OpenSwath::SpectrumAccessPtr& swath_ptr,
                               OpenSwath::LightTargetedExperiment& transition_exp_used, const RangeMobility& im_range,
                               OpenSwath::IDataFrameWriter* ivw) const
  {
    //getParams();
    typedef std::map<std::string, std::vector<OpenSwath::LightTransition> > Mmap;
    Mmap transmap;
    OpenSwath::TransitionHelper::convert(transition_exp_used, transmap);
    // std::cout << "nr peptides : " << transmap.size() << std::endl;

    Mmap::iterator m_begin = transmap.begin();
    Mmap::iterator m_end = transmap.end();
    std::vector<std::string> transitionsNames;

    for (; m_begin != m_end; ++m_begin)
    {
      transitionsNames.push_back(m_begin->first);
    }

    ivw->colnames(transitionsNames);
    //iterate over spectra

    for (UInt i = 0; i < swath_ptr->getNrSpectra(); ++i)
    {
      OpenSwath::SpectrumPtr s = swath_ptr->getSpectrumById(i);
      SpectrumSequence spec;
      spec.push_back(s);
      OpenSwath::SpectrumMeta specmeta = swath_ptr->getSpectrumMetaById(i);
      std::cout << "Processing Spectrum  " << i << "RT " << specmeta.RT << std::endl;

      //iterate over spectra
      size_t xx = 0;
      Mmap::iterator beg = transmap.begin();
      Mmap::iterator end = transmap.end();
      std::vector<double> score1v;
      std::vector<double> score2v;
      for (; beg != end; ++beg, ++xx)
      {
        //std::cout << "analysing transition" << xx << beg->second.size()
        //    << " " << beg->first << std::endl;
        double score1;
        double score2;
        //OpenSwath::LightPeptide pep;
        score(spec, beg->second, im_range, score1, score2);

        score1v.push_back(score1);
        score2v.push_back(score2);
      } //end of for loop over transitions

      //std::string ispectrum = boost::lexical_cast<std::string>(i);
      std::string specRT = boost::lexical_cast<std::string>(specmeta.RT);
      ivw->store("score1_" + specRT, score1v);
      ivw->store("score2_" + specRT, score2v);
    } //end of for loop over spectra
  }

  void DiaPrescore::score(const SpectrumSequence& spec,
                          const std::vector<OpenSwath::LightTransition>& lt,
                          const RangeMobility& im_range,
                          double& dotprod,
                          double& manhattan) const
  {
    std::vector<std::pair<double, double> > res;
    std::vector<std::pair<double, double> > spectrumWIso, spectrumWIsoNegPreIso;
    int chg;
    // add expected isotope intensities for every transition productMZ based on averagine
    //TODO allow usage of annotated formulas from transition.compound.sum_formula
    for (const auto& transition : lt)
    {
      chg = 1;
      if (transition.fragment_charge != 0) chg = transition.fragment_charge;
      DIAHelpers::addSinglePeakIsotopes2Spec(transition.getProductMZ(),
                                             transition.getLibraryIntensity(),
                                             spectrumWIso,
                                             nr_isotopes_,
                                             chg);
    }
    // duplicate since we will add differently weighted preIsotope intensities
    spectrumWIsoNegPreIso.reserve(spectrumWIso.size());
    std::copy(spectrumWIso.begin(), spectrumWIso.end(), back_inserter(spectrumWIsoNegPreIso));
    UInt nrNegPeaks = 2;
    double avgTheorTransitionInt = std::accumulate(lt.begin(),lt.end(),0.,[](double val, const OpenSwath::LightTransition& lt){return val + lt.getLibraryIntensity();});
    avgTheorTransitionInt /= lt.size();
    double negWeight = 0.5 * avgTheorTransitionInt; // how much of ONE transition should be negatively weighted at the prePeaks (distributed equally on them)
    // for every transition add either zero weighted (for manhattan) or negatively weighted (for dotprod) preIsotope intensities
    for (const auto& transition : lt)
    {
      chg = 1.;
      if (transition.fragment_charge != 0) chg = transition.fragment_charge;
      DIAHelpers::addPreisotopeWeights(transition.getProductMZ(), spectrumWIso, nrNegPeaks, 0.0,
                                       Constants::C13C12_MASSDIFF_U,
                                       chg);
      DIAHelpers::addPreisotopeWeights(transition.getProductMZ(),
                                       spectrumWIsoNegPreIso,
                                       nrNegPeaks,
                                       -negWeight,
                                       Constants::C13C12_MASSDIFF_U,
                                       chg);
    }
    //sort by mz
    DIAHelpers::sortByFirst(spectrumWIso);
    DIAHelpers::sortByFirst(spectrumWIsoNegPreIso);

    // compare against the spectrum with 0 weight preIsotope peaks
    std::vector<double> mzTheor, intTheor;
    DIAHelpers::extractFirst(spectrumWIso, mzTheor);
    DIAHelpers::extractSecond(spectrumWIso, intTheor);
    std::vector<double> intExp, mzExp, imExp;
    DIAHelpers::integrateWindows(spec, mzTheor, dia_extract_window_, intExp, mzExp, imExp, im_range);
    std::transform(intExp.begin(), intExp.end(), intExp.begin(), [](double val){return std::sqrt(val);});
    std::transform(intTheor.begin(), intTheor.end(), intTheor.begin(), [](double val){return std::sqrt(val);});

    // get sum for normalization. All entries in both should be positive
    double intExpTotal = std::accumulate(intExp.begin(), intExp.end(), 0.0);
    double intTheorTotal = std::accumulate(intTheor.begin(), intTheor.end(), 0.0);

    OpenSwath::normalize(intExp, intExpTotal, intExp);
    OpenSwath::normalize(intTheor, intTheorTotal, intTheor);

    //TODO think about normalizing the distance by dividing by the max value 2.
    // Generally I think a combined manhattan distance is not the best feature here, since because of normalization,
    // different transitions affect each other (e.g. if one transition is missing, the other(s) get a much higher
    // normalized value and the whole distance is "penalized twice")
    // Maybe we could use two features, one for the average manhattan distance and one for matching of the total intensities to the
    // library intensities. Also maybe normalising by the max-value or the monoisotope (instead of the total sum) helps?
    manhattan = OpenSwath::manhattanDist(intExp.begin(), intExp.end(), intTheor.begin());

    // compare against the spectrum with negative weight preIsotope peaks
    std::vector<double> intTheorNeg;
    // WARNING: This was spectrumWIso and therefore with 0 preIso weights in earlier versions! Was this a bug?
    // Otherwise, we don't need the second spectrum at all.
    DIAHelpers::extractSecond(spectrumWIsoNegPreIso, intTheorNeg);
    // Sqrt does not work if we actually have negative values
    //std::transform(intTheorNeg.begin(), intTheorNeg.end(), intTheorNeg.begin(), OpenSwath::mySqrt());
    double intTheorNegEuclidNorm = OpenSwath::norm(intTheorNeg.begin(), intTheorNeg.end()); // use Euclidean norm since we have negative values
    OpenSwath::normalize(intTheorNeg, intTheorNegEuclidNorm, intTheorNeg);

    // intExp is normalized already, but we can normalize again with euclidean norm to have the same norm (not sure if it makes much of a difference)
    double intExpEuclidNorm = OpenSwath::norm(intExp.begin(), intExp.end());
    double intTheorEuclidNorm = OpenSwath::norm(intTheor.begin(), intTheor.end());
    OpenSwath::normalize(intExp, intExpEuclidNorm, intExp);
    OpenSwath::normalize(intTheor, intTheorEuclidNorm, intTheor);

    //calculate maximum possible value and maximum negative value to rescale
    // depends on the amount of relative weight is negative
    // TODO check if it is the same amount for every spectrum, then we could leave it out.
    double negVal = (-negWeight/intTheorNegEuclidNorm) * sqrt(nrNegPeaks*lt.size());
    std::vector<double> intTheorNegBest;
    intTheorNegBest.resize(intTheorNeg.size());
    std::transform(intTheorNeg.begin(), intTheorNeg.end(), intTheorNegBest.begin(),
                   [&](double val){
                   if (val >= 0)
                   {
                     return val * nrNegPeaks * lt.size() * negWeight/intTheorNegEuclidNorm;
                   }
                   else
                   {
                     return 0.;
                   }
    });
    double intTheorNegBestEuclidNorm = OpenSwath::norm(intTheorNegBest.begin(), intTheorNegBest.end());
    OpenSwath::normalize(intTheorNegBest, intTheorNegBestEuclidNorm, intTheorNegBest);
    double posVal = OpenSwath::dotProd(intTheorNegBest.begin(), intTheorNegBest.end(), intTheorNeg.begin());

    dotprod = OpenSwath::dotProd(intExp.begin(), intExp.end(), intTheorNeg.begin());
    //simplified: dotprod = (((dotprod - negVal) * (1. - -1.)) / (posVal - negVal)) + -1.;
    dotprod = (((dotprod - negVal) * 2.) / (posVal - negVal)) - 1.;
  }

  void DiaPrescore::updateMembers_()
  {
    dia_extract_window_ = (double) param_.getValue(
      "dia_extraction_window");
    nr_isotopes_ = (int) param_.getValue("nr_isotopes");
    //TODO nr_charges_ is never used???
    nr_charges_ = (int) param_.getValue("nr_charges");
  }

  void DiaPrescore::defineDefaults()
  {
    defaults_.setValue("dia_extraction_window", 0.1,
                       "DIA extraction window in Th.");
    defaults_.setMinFloat("dia_extraction_window", 0.0); //done
    defaults_.setValue("nr_isotopes", 4, "nr of istopes");
    defaults_.setValue("nr_charges", 4, "nr charges");
    defaultsToParam_();
  }

  DiaPrescore::DiaPrescore(double dia_extract_window, int nr_isotopes, int nr_charges) :
    DefaultParamHandler("DIAPrescore"),
    dia_extract_window_(dia_extract_window),
    nr_isotopes_(nr_isotopes),
    nr_charges_(nr_charges)
  {
  }

  DiaPrescore::DiaPrescore() :
    DefaultParamHandler("DIAPrescore")
  {
    defineDefaults();
  }

}
