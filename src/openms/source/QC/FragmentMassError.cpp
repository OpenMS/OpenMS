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
// $Maintainer: Chris Bielow$
// $Authors: Patricia Scheil, Swenja Wagner$
// --------------------------------------------------------------------------

#include <OpenMS/QC/FragmentMassError.h>
#include <string>
#include <iostream>
#include <plugin-api.h>

namespace OpenMS
{
  void FragmentMassError::compute(FeatureMap& fmap, const MSExperiment& exp, const double tolerance, const bool tolerance_unit_ppm)
  {
    FMEStatistics result;
    //result.average_ppm = 0.;
    //result.variance_ppm = 0.;

    //accumulates ppm errors over all first PeptideHits
    double accumulator_ppm{};

    //counts number of ppm errors
    UInt32 counter_ppm{};

    float rt_tolerance = 0.05;

    //---------------------------------------------------------------------
    // Prepare MSExperiment
    //---------------------------------------------------------------------

    if (!exp.isSorted())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MSExperiment is not sorted by ascending RT");
    }

    // filter settings
    WindowMower window_mower_filter;
    Param filter_param = window_mower_filter.getParameters();
    filter_param.setValue("windowsize", 100.0, "The size of the sliding window along the m/z axis.");
    filter_param.setValue("peakcount", 6, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    window_mower_filter.setParameters(filter_param);

    MSExperiment exp_filtered(exp);

    for (MSSpectrum& spec : exp_filtered)
    {
      window_mower_filter.filterPeakSpectrum(spec);
    }


    auto lamCompPPM = [&exp_filtered, rt_tolerance, tolerance, tolerance_unit_ppm, &accumulator_ppm, &counter_ppm](PeptideIdentification& pep_id)
    {
      if (pep_id.getHits().empty())
      {
        LOG_WARN << "PeptideHits of PeptideIdentification with RT: " << pep_id.getRT() << " and MZ: " << pep_id.getMZ() << " is empty.";
        return;
      }

      //---------------------------------------------------------------------
      // FIND DATA FOR THEORETICAL SPECTRUM
      //---------------------------------------------------------------------

      //sequence
      AASequence seq = pep_id.getHits()[0].getSequence();

      //charge
      double mass = seq.getMonoWeight();
      double mz = pep_id.getMZ();
      double z = mass/mz;
      Int charge = round(z);

      //theoretical peak spectrum
      PeakSpectrum theo_spectrum;

      //---------------------------------------------------------------------
      // CREATE THEORETICAL SPECTRUM
      //---------------------------------------------------------------------

      //initialize a TheoreticalSpectrumGenerator
      TheoreticalSpectrumGenerator theo_gen;

      //get current parameters (default)
      Param theo_gen_settings = theo_gen.getParameters();

      //default: b- and y-ions?
      theo_gen_settings.setValue("add_a_ions", "true");
      //theo_settings.setValue("add_b_ions", "true");
      theo_gen_settings.setValue("add_c_ions", "true");
      theo_gen_settings.setValue("add_x_ions", "true");
      //theo_settings.setValue("add_y_ions", "true");
      theo_gen_settings.setValue("add_z_ions", "true");

      //set changed parameters
      theo_gen.setParameters(theo_gen_settings);

      //generate a-, b- and y-ion spectrum of peptide seq with charge
      theo_gen.getSpectrum(theo_spectrum, seq, charge, charge);

      //-----------------------------------------------------------------------
      // GET EXPERIMENTAL SPECTRUM MATCHING TO PEPTIDEIDENTIFICTION
      //-----------------------------------------------------------------------

      double rt_pep  = pep_id.getRT();

      MSExperiment::ConstIterator it = exp_filtered.RTBegin(rt_pep - rt_tolerance);
      if (it == exp_filtered.end())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The retention time of the mzML and featureXML file does not match.");
      }

      const auto& exp_spectrum = *it;

      if (exp_spectrum.getRT() - rt_pep > rt_tolerance)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "PeptideID with RT " + std::to_string(rt_pep) + " s does not have a matching MS2 Spectrum. Closest RT was " + std::to_string(exp_spectrum.getRT()) + ", which seems to far off.");
      }
      if (exp_spectrum.getMSLevel() != 2)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The matching retention time of the mzML is not a MS2 Spectrum.");
      }

      //-----------------------------------------------------------------------
      // COMPARE THEORETICAL AND EXPERIMENTAL SPECTRUM
      //-----------------------------------------------------------------------

      if (exp_spectrum.empty() || theo_spectrum.empty())
      {
        LOG_WARN << "The spectrum with " + std::to_string(exp_spectrum.getRT()) + " is empty." << "\n";
        return;
      }

      //stores ppms for one spectrum
      DoubleList ppms{};

      double inf = std::numeric_limits<double>::infinity();

      //exp_peak matching to previous theo_peak
      double current_exp = inf;

      //minimal ppm
      double ppm = inf;

      for (const Peak1D& peak : theo_spectrum)
      //for (UInt32 i = 0; i < theo_spectrum.size(); ++i)
      {
        const double theo_mz = peak.getMZ();
        Size index = exp_spectrum.findNearest(theo_mz);
        const double exp_mz = exp_spectrum[index].getMZ();

        const double mz_tolerance = tolerance_unit_ppm ?  theo_mz * tolerance * 1e-6 : tolerance;


        //found peak match
        if (std::abs(theo_mz-exp_mz) < mz_tolerance)
        {
          auto current_ppm = theo_mz - exp_mz;

          //first peak in tolerance range
          if (current_exp == inf)
          {
            ppm = current_ppm;
            current_exp = exp_mz;
          }

          //theo_peak matches to a exp_peak that is already matched
          //&& ppm is smaller than before
          if (current_exp == exp_mz && abs(current_ppm) < abs(ppm))
          {
            ppm = current_ppm;
            std::cout << "theo_mz: " << theo_mz<< std::endl;
            std::cout << "exp_mz: " << exp_mz<< std::endl;
          }

          //theo_peak matches to another exp_peaks
          if (current_exp != exp_mz)
          {
            // TODO change to ppm
            ppms.push_back(ppm);

            std::cout << "ppm " << ppm << std::endl;
            //std::cout << "Theo  " << theo_mz << "    Exp   " << exp_mz << std::endl;

            ++ counter_ppm;

            std::cout<< "c:  " << counter_ppm << std::endl;
            accumulator_ppm += ppm;
            ppm = current_ppm;
            current_exp = exp_mz;
          }

         }
      }

      //last peak doesn't have a successor so it has to be added manually
      ppms.push_back(ppm);
      accumulator_ppm += ppm;
      ++ counter_ppm;


      //-----------------------------------------------------------------------
      // WRITE PPM ERROR IN PEPTIDEHIT
      //-----------------------------------------------------------------------

      pep_id.getHits()[0].setMetaValue("ppm_errors", ppms);

    };

    auto lamVar = [&result, &counter_ppm](const PeptideIdentification& pep_id)
    {
      for (auto ppm : (pep_id.getHits()[0].getMetaValue("ppm_errors")).toDoubleList())
      {
        result.variance_ppm += (pow((ppm - result.average_ppm),2) / counter_ppm);
        std::cout << "counter_ppm: " << counter_ppm << std::endl;
      }
    };

    //computation of ppms
    QCBase::iterateFeatureMap(fmap, lamCompPPM);

    //if there are no matching peaks, the counter is zero and it is not possible to find ppms
    if (counter_ppm == 0)
    {
      results_.push_back(result);
      return;
    }

    //computes average
    result.average_ppm = accumulator_ppm/counter_ppm;

    //computes variance
    QCBase::iterateFeatureMap(fmap, lamVar);

    results_.push_back(result);

  }


  std::vector<FragmentMassError::FMEStatistics> FragmentMassError::getResults() const
  {
    return results_;
  }


  QCBase::Status FragmentMassError::requires() const
  {
    return QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT;
  }
} //namespace OpenMS