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

#include <assert.h>
#include <string>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>

namespace OpenMS
{
  void FragmentMassError::compute(FeatureMap& fmap, const MSExperiment& exp, const double tolerance, const ToleranceUnit tolerance_unit)
  {
    FMEStatistics result;

    // accumulates ppm errors over all first PeptideHits
    double accumulator_ppm{};

    // counts number of ppm errors
    UInt32 counter_ppm{};

    const float rt_tolerance = 0.05f;

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

    // computes the FragmentMassError
    auto lamCompPPM = [&exp, rt_tolerance, tolerance, tolerance_unit, &accumulator_ppm, &counter_ppm, &window_mower_filter](PeptideIdentification& pep_id)
    {
      if (pep_id.getHits().empty())
      {
        OPENMS_LOG_WARN << "PeptideHits of PeptideIdentification with RT: " << pep_id.getRT() << " and MZ: " << pep_id.getMZ() << " is empty.";
        return;
      }

      //---------------------------------------------------------------------
      // FIND DATA FOR THEORETICAL SPECTRUM
      //---------------------------------------------------------------------

      // sequence
      const AASequence& seq = pep_id.getHits()[0].getSequence();

      // charge
      Int charge = static_cast<Int>(round(seq.getMonoWeight() / pep_id.getMZ()));

      // if computed charge and the given charge in PeptideHits is not equal programm is terminated
      assert(charge == pep_id.getHits()[0].getCharge());


      //-----------------------------------------------------------------------
      // GET EXPERIMENTAL SPECTRUM MATCHING TO PEPTIDEIDENTIFICTION
      //-----------------------------------------------------------------------

      double rt_pep  = pep_id.getRT();

      MSExperiment::ConstIterator it = exp.RTBegin(rt_pep - rt_tolerance);
      if (it == exp.end())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The retention time of the mzML and featureXML file does not match.");
      }

      const auto& exp_spectrum = *it;

      if (exp_spectrum.getRT() - rt_pep > rt_tolerance)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "PeptideID with RT " + String(rt_pep) + " s does not have a matching MS2 Spectrum. Closest RT was "
        + String(exp_spectrum.getRT()) + ", which seems to far off.");
      }
      if (exp_spectrum.getMSLevel() != 2)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The matching retention time of the mzML is not a MS2 Spectrum.");
      }


      //---------------------------------------------------------------------
      // CREATE THEORETICAL SPECTRUM
      //---------------------------------------------------------------------

      // theoretical peak spectrum
      PeakSpectrum theo_spectrum;

      // initialize a TheoreticalSpectrumGenerator
      TheoreticalSpectrumGenerator theo_gen;

      // get current parameters (default)
      // default with b and y ions
      Param theo_gen_settings = theo_gen.getParameters();


      if (exp_spectrum.getPrecursors().empty() || exp_spectrum.getPrecursors()[0].getActivationMethods().empty())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No fragmentation method given.");
      }

      Precursor::ActivationMethod fm = (*exp_spectrum.getPrecursors()[0].getActivationMethods().begin());

      if (fm == Precursor::ActivationMethod::CID || fm == Precursor::ActivationMethod::HCID)
      {
        theo_gen_settings.setValue("add_b_ions", "true");
        theo_gen_settings.setValue("add_y_ions", "true");
      }

      else if (fm == Precursor::ActivationMethod::ECD || fm == Precursor::ActivationMethod::ETD)
      {
        theo_gen_settings.setValue("add_c_ions", "true");
        theo_gen_settings.setValue("add_z_ions", "true");
        theo_gen_settings.setValue("add_b_ions", "false");
        theo_gen_settings.setValue("add_y_ions", "false");
      }

      else
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Fragmentation method is not supported.");
      }

      // set changed parameters
      theo_gen.setParameters(theo_gen_settings);

      // generate a-, b- and y-ion spectrum of peptide seq with charge
      theo_gen.getSpectrum(theo_spectrum, seq, charge, charge);

      //-----------------------------------------------------------------------
      // COMPARE THEORETICAL AND EXPERIMENTAL SPECTRUM
      //-----------------------------------------------------------------------

      if (exp_spectrum.empty() || theo_spectrum.empty())
      {
        OPENMS_LOG_WARN << "The spectrum with RT: " + String(exp_spectrum.getRT()) + " is empty." << "\n";
        return;
      }

      auto exp_spectrum_filtered(exp_spectrum);
      window_mower_filter.filterPeakSpectrum(exp_spectrum_filtered);

      // stores ppms for one spectrum
      DoubleList ppms{};
      DoubleList dalton{};

      // exp_peak matching to previous theo_peak
      double current_exp = std::numeric_limits<double>::max();
      
      // max ppm
      double ppm = std::numeric_limits<double>::max();

      //max da
      double da = std::numeric_limits<double>::max();

      for (const Peak1D& peak : theo_spectrum)
      {
        const double theo_mz = peak.getMZ();
        Size index = exp_spectrum_filtered.findNearest(theo_mz);
        const double exp_mz = exp_spectrum_filtered[index].getMZ();
        
        const double mz_tolerance = (tolerance_unit==ToleranceUnit::PPM) ?  Math::ppmToMass(tolerance, theo_mz) : tolerance;

        // found peak match
        if (std::abs(theo_mz-exp_mz) < mz_tolerance)
        {
          auto current_ppm = Math::getPPM(exp_mz, theo_mz);
          auto current_da = exp_mz - theo_mz;

          // first peak in tolerance range
          if (current_exp == std::numeric_limits<double>::max())
          {
            ppm = current_ppm;
            da = current_da;
            current_exp = exp_mz;
          }

          // theo_peak matches to a exp_peak that is already matched
          // && ppm is smaller than before
          if (current_exp == exp_mz && abs(current_ppm) < abs(ppm))
          {
            ppm = current_ppm;
            da = current_da;
          }

          // theo_peak matches to another exp_peaks
          if (current_exp != exp_mz)
          {
            ppms.push_back(ppm);
            dalton.push_back(da);

            ++ counter_ppm;

            accumulator_ppm += ppm;
            ppm = current_ppm;
            da = current_da;
            current_exp = exp_mz;
          }

         }
      }

      // last peak doesn't have a successor so it has to be added manually
      ppms.push_back(ppm);
      dalton.push_back(da);
      accumulator_ppm += ppm;
      ++ counter_ppm;


      //-----------------------------------------------------------------------
      // WRITE PPM ERROR IN PEPTIDEHIT
      //-----------------------------------------------------------------------

      pep_id.getHits()[0].setMetaValue("fragment_mass_error_ppm", ppms);
      pep_id.getHits()[0].setMetaValue("fragment_mass_error_da", dalton);


    };

    auto lamVar = [&result](const PeptideIdentification& pep_id)
    {
      if (pep_id.getHits().empty())
      {
        OPENMS_LOG_WARN << "There is a Peptideidentification(RT: " << pep_id.getRT() << ", MZ: " << pep_id.getMZ() <<  ") without PeptideHits. " << "\n";
        return;
      }
      for (auto ppm : (pep_id.getHits()[0].getMetaValue("fragment_mass_error_ppm")).toDoubleList())
      {
        result.variance_ppm += pow((ppm - result.average_ppm),2);
      }
    };

    // computation of ppms
    QCBase::iterateFeatureMap(fmap, lamCompPPM);
    // if there are no matching peaks, the counter is zero and it is not possible to find ppms
    if (counter_ppm == 0)
    {
      results_.push_back(result);
      return;
    }

    // computes average
    result.average_ppm = accumulator_ppm / counter_ppm;

    // computes variance
    QCBase::iterateFeatureMap(fmap, lamVar);

    result.variance_ppm = result.variance_ppm / counter_ppm;

    results_.push_back(result);

  }

  const String& FragmentMassError::getName() const
  {
    return name_;
  }

  const std::vector<FragmentMassError::FMEStatistics>& FragmentMassError::getResults() const
  {
    return results_;
  }


  QCBase::Status FragmentMassError::requires() const
  {
    return QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT;
  }
} // namespace OpenMS



