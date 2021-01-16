// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DIAPrescoring.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SpectrumHelpers.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionHelper.h>
#include <OpenMS/OPENSWATHALGO/ALGO/StatsHelpers.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>

#include <iostream>

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
                   boost::bind(std::divides<double>(), _1, totalInt));
  }

  void getMZIntensityFromTransition(const std::vector<OpenSwath::LightTransition>& trans,
                                    std::vector<std::pair<double, double> >& res)
  {
    for (std::size_t i = 0; i < trans.size(); ++i)
    {
      res.emplace_back(trans[i].product_mz, trans[i].library_intensity);
    }
  }

  void DiaPrescore::operator()(OpenSwath::SpectrumAccessPtr swath_ptr,
                               OpenSwath::LightTargetedExperiment& transition_exp_used,
                               OpenSwath::IDataFrameWriter* ivw)
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

      OpenSwath::SpectrumPtr spec = swath_ptr->getSpectrumById(i);
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
        score(spec, beg->second, score1, score2);

        score1v.push_back(score1);
        score2v.push_back(score2);
      } //end of forloop over transitions

      //std::string ispectrum = boost::lexical_cast<std::string>(i);
      std::string specRT = boost::lexical_cast<std::string>(specmeta.RT);
      ivw->store("score1_" + specRT, score1v);
      ivw->store("score2_" + specRT, score2v);
    } //end of forloop over spectra
  }


  void DiaPrescore::score(OpenSwath::SpectrumPtr spec,
                          const std::vector<OpenSwath::LightTransition>& lt,
                          double& dotprod,
                          double& manhattan) const
  {
    std::vector<std::pair<double, double> > res;
    std::vector<std::pair<double, double> > spectrumWIso, spectrumWIsoNegPreIso;
    int chg = 1;
    for (const auto& transition : lt)
    {
      if (transition.fragment_charge != 0) chg = transition.fragment_charge;
      DIAHelpers::addIsotopes2Spec(transition.getProductMZ(), transition.getLibraryIntensity(), spectrumWIso, nr_isotopes_, chg);
    }
    spectrumWIsoNegPreIso.resize(spectrumWIso.size());
    std::copy(spectrumWIso.begin(), spectrumWIso.end(), spectrumWIsoNegPreIso.begin());
    for (const auto& transition : lt)
    {
      if (transition.fragment_charge != 0) chg = transition.fragment_charge;
      DIAHelpers::addPreisotopeWeights(transition.getProductMZ(), spectrumWIso, 2, 0.0);
      DIAHelpers::addPreisotopeWeights(transition.getProductMZ(), spectrumWIsoNegPreIso, 2, -0.5);
    }

    // compare against the spectrum with 0 weight preIsotope peaks
    std::vector<double> mzTheor, intTheor;
    DIAHelpers::extractFirst(spectrumWIso, mzTheor);
    DIAHelpers::extractSecond(spectrumWIso, intTheor);
    std::vector<double> intExp, mzExp;
    DIAHelpers::integrateWindows(spec, mzTheor, dia_extract_window_, intExp, mzExp);
    std::transform(intExp.begin(), intExp.end(), intExp.begin(), OpenSwath::mySqrt());
    std::transform(intTheor.begin(), intTheor.end(), intTheor.begin(), OpenSwath::mySqrt());

    double intExpTotal = std::accumulate(intExp.begin(), intExp.end(), 0.0);
    double intTheorTotal = std::accumulate(intTheor.begin(), intTheor.end(), 0.0);

    OpenSwath::normalize(intExp, intExpTotal, intExp);
    OpenSwath::normalize(intTheor, intTheorTotal, intTheor);

    manhattan = OpenSwath::manhattanDist(intExp.begin(), intExp.end(), intTheor.begin());

    // compare against the spectrum with negative weight preIsotope peaks
    std::vector<double> intTheor2;
    // WARNING: This was spectrumWIso and therefore with 0 preIso weights before! Was this a bug?
    // Otherwise we dont need the second spectrum at all.
    DIAHelpers::extractSecond(spectrumWIsoNegPreIso, intTheor2);
    std::transform(intTheor2.begin(), intTheor2.end(), intTheor2.begin(), OpenSwath::mySqrt());

    intExpTotal = OpenSwath::norm(intExp.begin(), intExp.end());
    intTheorTotal = OpenSwath::norm(intTheor2.begin(), intTheor2.end());

    OpenSwath::normalize(intExp, intExpTotal, intExp);
    OpenSwath::normalize(intTheor2, intTheorTotal, intTheor2);

    dotprod = OpenSwath::dotProd(intExp.begin(), intExp.end(), intTheor2.begin());
  }

  /*void DiaPrescore::scoreDDA(OpenSwath::SpectrumPtr spec,
                          const std::vector<OpenSwath::LightTransition>& lt,
                          double& dotprod,
                          double& manhattan) const
  {
    std::vector<std::pair<double, double> > res;
    getMZIntensityFromTransition(lt, res);
    double monoIsoMZ = lt.begin()->getPrecursorMZ();
    std::vector<double> transitionmzs, theomasses;
    DIAHelpers::extractFirst(res, transitionmzs);
    std::vector<std::pair<double, double> > spectrumWIso, spectrumWIsoNegPreIso;
    DIAHelpers::addIsotopes2Spec(res, spectrumWIso, nr_charges_);
    spectrumWIsoNegPreIso.resize(spectrumWIso.size());
    std::copy(spectrumWIso.begin(), spectrumWIso.end(), spectrumWIsoNegPreIso.begin());
    //std::cout << spectrumWIso.size() << std::endl;
    DIAHelpers::addPreisotopeWeights(transitionmzs, spectrumWIso, 2, 0.0);
    //extracts masses from spectrum
    DIAHelpers::extractFirst(spectrumWIso, theomasses);
    std::vector<double>  intTheor;
    DIAHelpers::extractSecond(spectrumWIso, intTheor);
    std::vector<double> intExp, mzExp;
    DIAHelpers::integrateWindows(spec, theomasses, dia_extract_window_, intExp, mzExp);
    std::transform(intExp.begin(), intExp.end(), intExp.begin(), OpenSwath::mySqrt());
    std::transform(intTheor.begin(), intTheor.end(), intTheor.begin(), OpenSwath::mySqrt());

    double intExpTotal = std::accumulate(intExp.begin(), intExp.end(), 0.0);
    double intTheorTotal = std::accumulate(intTheor.begin(), intTheor.end(), 0.0);

    OpenSwath::normalize(intExp, intExpTotal, intExp);
    OpenSwath::normalize(intTheor, intTheorTotal, intTheor);

    manhattan = OpenSwath::manhattanDist(intExp.begin(), intExp.end(), intTheor.begin());

    //std::cout << spectrumWIso.size() << std::endl;
    DIAHelpers::addPreisotopeWeights(transitionmzs, spectrumWIsoNegPreIso, 2, -0.5);
    std::vector<double>  theorint2;
    DIAHelpers::extractSecond(spectrumWIso, theorint2);
    std::transform(theorint2.begin(), theorint2.end(), theorint2.begin(), OpenSwath::mySqrt());

    intExpTotal = OpenSwath::norm(intExp.begin(), intExp.end());
    intTheorTotal = OpenSwath::norm(theorint2.begin(), theorint2.end());

    OpenSwath::normalize(intExp, intExpTotal, intExp);
    OpenSwath::normalize(theorint2, intTheorTotal, theorint2);

    //    std::copy(intExp.begin(), intExp.end(), std::ostream_iterator<double>(std::cout, ", "));
    //    std::cout << std::endl;
    //    std::copy(theorint2.begin(), theorint2.end(), std::ostream_iterator<double>(std::cout, ", "));
    //    std::cout << std::endl;
    dotprod = OpenSwath::dotProd(intExp.begin(), intExp.end(), theorint2.begin());
  }
*/

  void DiaPrescore::updateMembers_()
  {
    dia_extract_window_ = (double) param_.getValue(
      "dia_extraction_window");
    nr_isotopes_ = (int) param_.getValue("nr_isotopes");
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
