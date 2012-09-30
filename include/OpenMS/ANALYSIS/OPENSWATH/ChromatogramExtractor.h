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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_EXTRACTCHROMATOGRAM_H
#define OPENMS_ANALYSIS_OPENSWATH_EXTRACTCHROMATOGRAM_H

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>

#ifdef _OPENMP
#include <omp.h>
#endif


namespace OpenMS
{

  /**
  @brief The ChromatogramExtractor extracts chromatograms from a mzML file.

  It will take as input a set of (TraML) transitions and will extract the
  signal of the provided map at the product ion m/z values specified by the
  transitions. The map is thus assumed to be an MS2 map from a SWATH / DIA
  experiment.

  */
  class ChromatogramExtractor :
    public ProgressLogger
  {

public:

    /// Constructor
    ChromatogramExtractor()
    {
    }

    /// Destructor
    ~ChromatogramExtractor()
    {
    }

    /// Extract chromatograms defined by the TargetedExperiment from the input map and write them to the output map
    template <typename ExperimentType>
    void extractChromatograms(const ExperimentType & input, ExperimentType & output, OpenMS::TargetedExperiment & transition_exp, double & extract_window, bool ppm,
                              TransformationDescription & trafo, double rt_extraction_window, String filter)
    {

      // invert the trafo because we want to transform nRT values to "real" RT values
      trafo.invert();

      Size input_size = input.size();
      if (input_size < 1)
      {
        return;
      }
      SpectrumSettings settings = input[0];
      int used_filter = -1;
      if (filter == "tophat")
      {
        used_filter = 1;
      }
      else if (filter == "bartlett")
      {
        used_filter = 2;
      }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                         "Filter either needs to be tophat or bartlett");
      }

      // Store the peptide retention times in an intermediate map
      PeptideRTMap.clear();
      for (Size i = 0; i < transition_exp.getPeptides().size(); i++)
      {
        const TargetedExperiment::Peptide & pep = transition_exp.getPeptides()[i];
        if (pep.rts.empty() || pep.rts[0].getCVTerms()["MS:1000896"].empty())
        {
          // we dont have retention times -> this is only a problem if we actually
          // wanted to use the RT limit feature.
          if (rt_extraction_window >= 0)
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                             "Error: Peptide " + pep.id + " does not have normalized retention times (term 1000896) which are necessary to perform an RT-limited extraction");
          }
          continue;
        }
        PeptideRTMap[pep.id] = pep.rts[0].getCVTerms()["MS:1000896"][0].getValue().toString().toDouble();
      }

      // sort the transition experiment by product mass
      // this is essential because the algorithm assumes sorted transitions!
      transition_exp.sortTransitionsByProductMZ();

      // prepare all the spectra (but leave them empty)
      std::vector<typename ExperimentType::ChromatogramType> chromatograms;
      prepare_spectra(settings, chromatograms, transition_exp);

      //go through all spectra
      startProgress(0, input_size, "Extracting chromatograms");
      for (Size scan_idx = 0; scan_idx < input_size; ++scan_idx)
      {
        setProgress(scan_idx);

        if (input[scan_idx].size() == 0)
          continue;

        Size peak_idx = 0;

        double mz;
        double integrated_intensity = 0;

        // go through all transitions / chromatograms which are sorted by
        // ProductMZ. We can use this to step through the spectrum and at the
        // same time step through the transitions. We increase the peak counter
        // until we hit the next transition and then extract the signal.
        for (Size k = 0; k < chromatograms.size(); ++k)
        {

          double current_rt = input[scan_idx].getRT();
          if (outside_extraction_window(transition_exp.getTransitions()[k], current_rt, trafo, rt_extraction_window))
          {
            continue;
          }

          typename ExperimentType::ChromatogramType::PeakType p;
          mz = transition_exp.getTransitions()[k].getProductMZ();

          if (used_filter == 1)
          {
            extract_value_tophat(input[scan_idx], mz, peak_idx, integrated_intensity, extract_window, ppm);
          }
          else if (used_filter == 2)
          {
            extract_value_bartlett(input[scan_idx], mz, peak_idx, integrated_intensity, extract_window, ppm);
          }

          p.setRT(current_rt);
          p.setIntensity(integrated_intensity);
          chromatograms[k].push_back(p);
        }
      }
      endProgress();

      // add all the chromatograms to the output
      output.setChromatograms(chromatograms);
    }

private:

    /// This populates the chromatograms vector with empty chromatograms (but sets their meta-information)
    template <class SpectrumSettings, class ExperimentType>
    void prepare_spectra(SpectrumSettings & settings, std::vector<ExperimentType> & chromatograms, OpenMS::TargetedExperiment & transition_exp)
    {

      // first prepare all the spectra (but leave them empty)
      for (Size i = 0; i < transition_exp.getTransitions().size(); i++)
      {
        const ReactionMonitoringTransition * transition = &transition_exp.getTransitions()[i];

        ExperimentType chrom;
        // Create precursor and set
        // 1) the target m/z
        // 2) the isolation window (upper/lower)
        // 3) the peptide sequence
        Precursor prec;
        prec.setMZ(transition->getPrecursorMZ());
        if (settings.getPrecursors().size() > 0)
        {
          prec.setIsolationWindowLowerOffset(settings.getPrecursors()[0].getIsolationWindowLowerOffset());
          prec.setIsolationWindowUpperOffset(settings.getPrecursors()[0].getIsolationWindowUpperOffset());
        }

        //set precursor sequence
        String pepref = transition->getPeptideRef();
        for (Size pep_idx = 0; pep_idx < transition_exp.getPeptides().size(); pep_idx++)
        {
          const OpenMS::TargetedExperiment::Peptide * pep = &transition_exp.getPeptides()[pep_idx];
          if (pep->id == pepref)
          {
            prec.setMetaValue("peptide_sequence", pep->sequence);
            break;
          }
        }
        // add precursor to spectrum
        chrom.setPrecursor(prec);

        // Create product and set its m/z
        Product prod;
        prod.setMZ(transition->getProductMZ());
        chrom.setProduct(prod);

        // Set the rest of the meta-data
        chrom.setInstrumentSettings(settings.getInstrumentSettings());
        chrom.setAcquisitionInfo(settings.getAcquisitionInfo());
        chrom.setSourceFile(settings.getSourceFile());

        // Set the id of the chromatogram, using the id of the transition (this gives directly the mapping of the two)
        chrom.setNativeID(transition->getNativeID());
        chrom.setChromatogramType(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
        chromatograms.push_back(chrom);
      }

    }

    template <typename SpectrumType>
    void extract_value_tophat(const SpectrumType & input, const double & mz, Size & peak_idx, double & integrated_intensity, const double & extract_window, const bool ppm)
    {
      // calculate extraction window
      double left, right;
      if (ppm)
      {
        left  = mz - mz * extract_window / 2.0 * 1.0e-6;
        right = mz + mz * extract_window / 2.0 * 1.0e-6;
      }
      else
      {
        left  = mz - extract_window / 2.0;
        right = mz + extract_window / 2.0;
      }

      Size walker;

      // advance the peak_idx until we hit the m/z value of the next transition
      while (peak_idx < input.size() && input[peak_idx].getMZ() < mz)
      {
        peak_idx++;
      }

      integrated_intensity = 0;

      // walk right and left and add to our intensity
      walker = peak_idx;
      // if we moved past the end of the spectrum, we need to try the last peak of the spectrum (it could still be within the window)
      if (peak_idx >= input.size())
      {
        walker = input.size() - 1;
      }

      // add the current peak if it is between right and left
      if (input[walker].getMZ() > left && input[walker].getMZ() < right)
      {
        integrated_intensity += input[walker].getIntensity();
      }

      // walk to the right until we go outside the window, then walk to the left until we are outside the window
      walker = peak_idx - 1;
      while (walker > 0 && input[walker].getMZ() > left && input[walker].getMZ() < right)
      {
        integrated_intensity += input[walker].getIntensity(); walker--;
      }
      walker = peak_idx + 1;
      while (walker<input.size() && input[walker].getMZ()> left &&  input[walker].getMZ() < right)
      {
        integrated_intensity += input[walker].getIntensity(); walker++;
      }
    }

    template <typename SpectrumType>
    void extract_value_bartlett(const SpectrumType & input, const double & mz, Size & peak_idx, double & integrated_intensity, const double & extract_window, const bool ppm)
    {
      // calculate extraction window
      double left, right, half_window_size, weight;
      if (ppm)
      {
        half_window_size = mz * extract_window / 2.0 * 1.0e-6;
        left  = mz - mz * extract_window / 2.0 * 1.0e-6;
        right = mz + mz * extract_window / 2.0 * 1.0e-6;
      }
      else
      {
        half_window_size = extract_window / 2.0;
        left  = mz - extract_window / 2.0;
        right = mz + extract_window / 2.0;
      }

      Size walker;

      // advance the peak_idx until we hit the m/z value of the next transition
      while (peak_idx < input.size() && input[peak_idx].getMZ() < mz)
      {
        peak_idx++;
      }

      integrated_intensity = 0;

      // walk right and left and add to our intensity
      walker = peak_idx;
      // if we moved past the end of the spectrum, we need to try the last peak of the spectrum (it could still be within the window)
      if (peak_idx >= input.size())
      {
        walker = input.size() - 1;
      }

      // add the current peak if it is between right and left
      if (input[walker].getMZ() > left && input[walker].getMZ() < right)
      {
        weight =  1 - fabs(input[walker].getMZ() - mz) / half_window_size;
        integrated_intensity += input[walker].getIntensity() * weight;
      }

      // walk to the right until we go outside the window, then walk to the left until we are outside the window
      walker = peak_idx - 1;
      while (walker > 0 && input[walker].getMZ() > left && input[walker].getMZ() < right)
      {
        weight =  1 - fabs(input[walker].getMZ() - mz) / half_window_size;
        integrated_intensity += input[walker].getIntensity() * weight; walker--;
      }
      walker = peak_idx + 1;
      while (walker<input.size() && input[walker].getMZ()> left &&  input[walker].getMZ() < right)
      {
        weight = 1 - fabs(input[walker].getMZ() - mz) / half_window_size;
        integrated_intensity += input[walker].getIntensity() * weight; walker++;
      }
    }

    void extract_value_tophat(const std::vector<double>::const_iterator & mz_start, std::vector<double>::const_iterator & mz_it,
                              const std::vector<double>::const_iterator & mz_end, std::vector<double>::const_iterator & int_it,
                              const double & mz, double & integrated_intensity, double & extract_window, bool ppm)
    {
      // calculate extraction window
      double left, right;
      if (ppm)
      {
        left  = mz - mz * extract_window / 2.0 * 1.0e-6;
        right = mz + mz * extract_window / 2.0 * 1.0e-6;
      }
      else
      {
        left  = mz - extract_window / 2.0;
        right = mz + extract_window / 2.0;
      }

      std::vector<double>::const_iterator mz_walker;
      std::vector<double>::const_iterator int_walker;

      // advance the mz / int iterator until we hit the m/z value of the next transition
      while (mz_it != mz_end && (*mz_it) < mz)
      {
        mz_it++; int_it++;
      }

      integrated_intensity = 0;

      // walk right and left and add to our intensity
      mz_walker  = mz_it;
      int_walker = int_it;

      // if we moved past the end of the spectrum, we need to try the last peak of the spectrum (it could still be within the window)
      if (mz_it == mz_end)
      {
        mz_walker--; int_walker--;
      }

      // add the current peak if it is between right and left
      if ((*mz_walker) > left && (*mz_walker) < right)
      {
        integrated_intensity += (*int_walker);
      }

      // walk to the right until we go outside the window, then walk to the left until we are outside the window
      mz_walker  = mz_it;
      int_walker = int_it;
      mz_walker--;
      int_walker--;
      while (mz_walker != mz_start && (*mz_walker) > left && (*mz_walker) < right)
      {
        integrated_intensity += (*int_walker); mz_walker--; int_walker--;
      }
      mz_walker  = mz_it;
      int_walker = int_it;
      mz_walker++;
      int_walker++;
      while (mz_walker != mz_end && (*mz_walker) > left && (*mz_walker) < right)
      {
        integrated_intensity += (*int_walker); mz_walker++; int_walker++;
      }
    }

    bool outside_extraction_window(const ReactionMonitoringTransition & transition, double current_rt,
                                   const TransformationDescription & trafo, double rt_extraction_window)
    {
      if (rt_extraction_window < 0)
      {
        return false;
      }

      // Get the expected retention time, apply the RT-transformation
      // (which describes the normalization) and then take the difference.
      // Note that we inverted the transformation in the beginning because
      // we want to transform from normalized to real RTs here and not the
      // other way round.
      double expected_rt = PeptideRTMap[transition.getPeptideRef()];
      double de_normalized_experimental_rt = trafo.apply(expected_rt);
      if (current_rt < de_normalized_experimental_rt - rt_extraction_window || current_rt > de_normalized_experimental_rt + rt_extraction_window)
      {
        return true;
      }
      return false;
    }

    std::map<OpenMS::String, double> PeptideRTMap;

  };

}

#endif
