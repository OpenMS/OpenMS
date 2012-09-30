// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Witold Wolski $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_DIAPRESCORING2_H_
#define OPENMS_ANALYSIS_OPENSWATH_DIAPRESCORING2_H_

// TODO (wolski): adhere to OpenMS coding conventions http://www-bs2.informatik.uni-tuebingen.de/services/OpenMS/OpenMS-release/html/coding__conventions.html
// TODO (wolski): comment the class and each method

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>

//#include "OpenSwath/Utils/DataFrameWriter.h"
//#include "OpenSwath/DataAccess/ITrans2Trans.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ITrans2Trans.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataFrameWriter.h"

//#include "OpenSwath/DataAccess/ISpectrumAccess.h"
//#include "OpenSwath/DataAccess/TransitionExperiment.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h"

#include <OpenMS/ANALYSIS/OPENSWATH/OpenMSHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/DIAHelpers.h>

#include "OpenMS/DATASTRUCTURES/DefaultParamHandler.h"
#include <iterator>

namespace OpenMS
{

  using namespace OpenSwath;

  struct mySqrt :
    std::unary_function<double, double>
  {
    double operator()(double x)
    {
      return sqrt(x);
    }

  };

  struct DiaPrescore2 :
    DefaultParamHandler
  {
    double dia_extract_window_;   //done
    double dia_centroided_;   //done
    double dia_byseries_intensity_min_;   //done
    double dia_byseries_ppm_diff_;   //done
    double dia_nr_isotopes_;
    double dia_nr_charges_;
    // parameters
    int nr_isotopes;
    int nr_charges;

    DiaPrescore2() :
      DefaultParamHandler("DIAPrescore2"), nr_isotopes(4), nr_charges(4)
    {
      defineDefaults();
    }

    void defineDefaults()
    {
      defaults_.setValue("dia_extraction_window", 0.05,
                         "DIA extraction window in Th.");
      defaults_.setMinFloat("dia_extraction_window", 0.0);   //done
      defaultsToParam_();
    }

    void getParams()
    {
      dia_extract_window_ = (DoubleReal) param_.getValue(
        "dia_extraction_window");
    }

public:
    void score(OpenSwath::SpectrumPtr spec,
               const std::vector<OpenSwath::LightTransition> & lt,
               double & dotprod,
               double & manhattan)
    {
      std::vector<std::pair<double, double> > res;
      getMZIntensityFromTransition(lt, res);
      std::vector<double> firstIstotope, theomasses;
      extractFirst(res, firstIstotope);
      std::vector<std::pair<double, double> > spectrum, spectrum2;
      addIsotopes2Spec(res, spectrum);
      spectrum2.resize(spectrum.size());
      std::copy(spectrum.begin(), spectrum.end(), spectrum2.begin());
      //std::cout << spectrum.size() << std::endl;
      addPreisotopeWeights(firstIstotope, spectrum, 2, 0.0);
      //extracts masses from spectrum
      extractFirst(spectrum, theomasses);
      std::vector<double>  theorint;
      extractSecond(spectrum, theorint);
      std::vector<double> intExp, mzExp;
      integrateWindows(spec, theomasses, dia_extract_window_, intExp,
                       mzExp);
      std::transform(intExp.begin(), intExp.end(), intExp.begin(), mySqrt());
      std::transform(theorint.begin(), theorint.end(), theorint.begin(), mySqrt());

      double intExptotal = std::accumulate(intExp.begin(), intExp.end(), 0.0);
      double intTheorTotal = std::accumulate(theorint.begin(), theorint.end(), 0.0);

      OpenSwath::normalize(intExp, intExptotal, intExp);
      OpenSwath::normalize(theorint, intTheorTotal, theorint);

      std::cout <<  "Exp " << intExp.size() << std::endl;
      std::cout <<  "theor" << theorint.size() << std::endl;
      //std::copy(theorint.begin(),theorint.end(),std::ostream_iterator<double>(std::cout,", "));
      //std::cout << std::endl;

      manhattan = OpenSwath::manhattanDist(intExp.begin(), intExp.end(), theorint.begin());

      //std::cout << spectrum.size() << std::endl;
      addPreisotopeWeights(firstIstotope, spectrum2, 2, -0.5);
      std::vector<double>  theorint2;
      extractSecond(spectrum, theorint2);
      std::transform(theorint2.begin(), theorint2.end(), theorint2.begin(), mySqrt());

      intExptotal = OpenSwath::norm(intExp.begin(), intExp.end());
      intTheorTotal = OpenSwath::norm(theorint2.begin(), theorint2.end());

      OpenSwath::normalize(intExp, intExptotal, intExp);
      OpenSwath::normalize(theorint2, intTheorTotal, theorint2);
      std::copy(intExp.begin(), intExp.end(), std::ostream_iterator<double>(std::cout, ", "));
      std::cout << std::endl;
      std::copy(theorint2.begin(), theorint2.end(), std::ostream_iterator<double>(std::cout, ", "));
      std::cout << std::endl;
      dotprod = OpenSwath::dotProd(intExp.begin(), intExp.end(), theorint2.begin());
    }

protected:

    void getMZIntensityFromTransition(const std::vector<OpenSwath::LightTransition> & trans, std::vector<std::pair<double, double> > & res)
    {
      for (std::size_t i = 0; i < trans.size(); ++i)
      {
        res.push_back(std::make_pair(trans[i].product_mz, trans[i].library_intensity));
      }
    }

public:
    void Prescore2(OpenSwath::SpectrumAccessPtr swath_ptr,
                   OpenSwath::LightTargetedExperiment & transition_exp_used,
                   OpenSwath::IDataFrameWriter * ivw)
    {
      getParams();
      int trans = transition_exp_used.transitions.size();
      std::cout << swath_ptr->getNrSpectra() << std::endl;

      typedef std::map<std::string, std::vector<OpenSwath::LightTransition> > Mmap;
      Mmap transmap;
      convert(transition_exp_used, transmap);
      std::cout << "nr peptides : " << transmap.size() << std::endl;

      Mmap::iterator beg = transmap.begin();
      Mmap::iterator end = transmap.end();
      std::vector<std::string> transitionsNames;

      for (; beg != end; ++beg)
      {
        transitionsNames.push_back(beg->first);
      }

      ivw->colnames(transitionsNames);
      //iterate over spectra

      for (uint32_t i = 0; i < swath_ptr->getNrSpectra(); ++i)
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
          //std::cout << "analysing transtion" << xx << beg->second.size()
          //    << " " << beg->first << std::endl;
          double score1;
          double score2;
          //OpenSwath::LightPeptide pep;
          score(spec, beg->second, score1, score2);

          score1v.push_back(score1);
          score2v.push_back(score2);
        }           //end of forloop over transitions

        std::string ispectrum = boost::lexical_cast<std::string>(i);
        std::string specRT = boost::lexical_cast<std::string>(specmeta.RT);
        ivw->store("score1_" + specRT, score1v);
        ivw->store("score2_" + specRT, score2v);
      }           //end of forloop over spectra
    }

    //iterate over all spectra and compute the DIA scores
    void Prescore(OpenSwath::SpectrumAccessPtr swath_ptr,
                  OpenSwath::LightTargetedExperiment & transition_exp_used,
                  OpenSwath::IDataFrameWriter * ivw)
    {
      getParams();
      int trans = transition_exp_used.transitions.size();
      std::cout << swath_ptr->getNrSpectra() << std::endl;

      typedef std::map<std::string, std::vector<OpenSwath::LightTransition> > Mmap;
      Mmap transmap;
      convert(transition_exp_used, transmap);
      std::cout << "nr peptides : " << transmap.size() << std::endl;

      Mmap::iterator beg = transmap.begin();
      Mmap::iterator end = transmap.end();
      std::vector<std::string> transitionsNames;

      for (; beg != end; ++beg)
      {
        transitionsNames.push_back(beg->first);
      }

      ivw->colnames(transitionsNames);
      //iterate over spectra
      for (uint32_t i = 0; i < swath_ptr->getNrSpectra(); ++i)
      {
        std::cout << "processing spec  " << i << std::endl;
        OpenSwath::SpectrumPtr spec = swath_ptr->getSpectrumById(i);
        //iterate over spectra
        size_t xx = 0;
        Mmap::iterator beg = transmap.begin();
        Mmap::iterator end = transmap.end();
        std::vector<double> score1v;
        std::vector<double> score2v;
        for (; beg != end; ++beg, ++xx)
        {
          //std::cout << "analysing transtion" << xx << beg->second.size()
          //    << " " << beg->first << std::endl;
          double score1;
          double score2;
          OpenSwath::LightPeptide pep;
          findPeptide(transition_exp_used, beg->first, pep);

          double bseries_score = 0, yseries_score = 0;
          OpenMS::AASequence aas(pep.sequence);
          for (std::vector<OpenSwath::LightModification>::const_iterator it =
                 pep.modifications.begin(); it != pep.modifications.end();
               ++it)
          {
            aas.setModification(it->location, "UniMod:" + it->unimod_id);
          }
          std::vector<double> firstIstotope, theomasses;
          std::vector<std::pair<double, double> > spectrum;

          simulateSpectrumFromAASequence(aas, firstIstotope, spectrum);
          //addPreisotopeWeights(firstIstotope, spectrum);

          extractFirst(spectrum, theomasses);
          std::vector<double> intExp, mzExp;
          integrateWindows(spec, theomasses, dia_extract_window_, intExp,
                           mzExp);


          std::vector<double>  theorint;
          extractSecond(spectrum, theorint);
          for (int k = 0; k < intExp.size(); ++k)
          {
            //std::cout << intExp[k]<< " " << theorint[k] << std::endl;
          }
          score1 = OpenSwath::dotProd(intExp.begin(), intExp.end(), theorint.begin());
          //std::cout << "score1 :" << score1 << std::endl;
          score1v.push_back(score1);


          //score2v.push_back(score2);
        }   //end of forloop over transitions

        std::string ispectrum = boost::lexical_cast<std::string>(i);
        ivw->store("score1_" + ispectrum, score1v);
        //ivw->store("score2_" + ispectrum, score2v);
      }   //end of forloop over spectra
    }

    void getNormalizedLibraryIntensities(
      const std::vector<OpenSwath::LightTransition> & transitions,
      std::vector<double> & normalizedLibraryIntensities     //normalized intensities
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

  };

}

#endif
