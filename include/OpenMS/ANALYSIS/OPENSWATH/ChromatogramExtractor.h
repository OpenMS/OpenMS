// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#ifndef OPENMS_ANALYSIS_OPENSWATH_CHROMATOGRAMEXTRACTOR_H
#define OPENMS_ANALYSIS_OPENSWATH_CHROMATOGRAMEXTRACTOR_H

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>

namespace OpenMS
{

  /**
   * @brief The ChromatogramExtractor extracts chromatograms from a spectra file.
   *
   * It will take as input a set of transitions coordinates and will extract
   * the signal of the provided map at the product ion m/z and retention time
   * (rt) values specified by the extraction coordinates. There are two
   * interfaces, the old interface will take a full TargetedExperiment and
   * assume that one wants to extract at the m/z of the transitions present in
   * the TargetedExperiment. The new interface (see also the
   * ChromatogramExtractorAlgorithm class) only expects a set of coordinates
   * which are up to the user to fill but a convenient prepare_coordinates
   * function is provided to create the coordinates for the most common case of
   * an MS2 and MS1 extraction.
   * 
   * In the case of MS2 extraction, the map is assumed to originate from a SWATH
   * (data-independent acquisition or DIA) experiment.
   *
  */
  class OPENMS_DLLAPI ChromatogramExtractor :
    public ProgressLogger
  {

public:

    typedef ChromatogramExtractorAlgorithm::ExtractionCoordinates ExtractionCoordinates;

    /**
     * @brief Extract chromatograms defined by the TargetedExperiment from the input map and write them to the output map.
     *
     * @param mz_extraction_window Extracts a window of this size in m/z
     * dimension (e.g. a window of 50 ppm means an extraction of 25 ppm on
     * either side)
     * @param rt_extraction_window Extracts a window of this size in RT
     * dimension (e.g. a window of 600 seconds means an extraction of 300
     * seconds on either side)
     *
     * @note: it will replace chromatograms in the output map, not append to
     * them!
     * @note: whenever possible, please use the ChromatogramExtractorAlgorithm
     * implementation since it is more flexible
    */
    template <typename ExperimentT>
    void extractChromatograms(const ExperimentT& input, ExperimentT& output, 
        OpenMS::TargetedExperiment& transition_exp, double mz_extraction_window, bool ppm,
        TransformationDescription trafo, double rt_extraction_window, String filter)
    {
      // invert the trafo because we want to transform nRT values to "real" RT values
      trafo.invert();

      Size input_size = input.size();
      if (input_size < 1)
      {
        return;
      }

      int used_filter = getFilterNr(filter);
      populatePeptideRTMap_(transition_exp, rt_extraction_window);

      // sort the transition experiment by product mass
      // this is essential because the algorithm assumes sorted transitions!
      transition_exp.sortTransitionsByProductMZ();

      // prepare all the spectra (but leave them empty)
      SpectrumSettings settings = input[0];
      std::vector<typename ExperimentT::ChromatogramType> chromatograms;
      prepareSpectra_(settings, chromatograms, transition_exp);

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
          if (outsideExtractionWindow_(transition_exp.getTransitions()[k], current_rt, trafo, rt_extraction_window))
          {
            continue;
          }

          typename ExperimentT::ChromatogramType::PeakType p;
          mz = transition_exp.getTransitions()[k].getProductMZ();

          if (used_filter == 1)
          {
            extract_value_tophat(input[scan_idx], mz, peak_idx, integrated_intensity, mz_extraction_window, ppm);
          }
          else if (used_filter == 2)
          {
            extract_value_bartlett(input[scan_idx], mz, peak_idx, integrated_intensity, mz_extraction_window, ppm);
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

    /**
     * @brief Extract chromatograms at the m/z and RT defined by the ExtractionCoordinates.
     *
     * @note: whenever possible, please use this ChromatogramExtractorAlgorithm implementation
     *
    */
    void extractChromatograms(const OpenSwath::SpectrumAccessPtr input, 
        std::vector< OpenSwath::ChromatogramPtr >& output, 
        std::vector<ExtractionCoordinates> extraction_coordinates,
        double mz_extraction_window, bool ppm, String filter)
    {
      ChromatogramExtractorAlgorithm().extractChromatograms(input, output, 
          extraction_coordinates, mz_extraction_window, ppm, filter);
    }

public:

    /**
     * @brief Prepare the extraction coordinates from a TargetedExperiment 
     *
     * Will fill the coordinates vector with the appropriate extraction
     * coordinates (transitions for MS2 extraction, peptide m/z for MS1
     * extraction). The output will be sorted by m/z.
     *
     * @param output_chromatograms An empty vector which will be initialized correctly
     * @param coordinates An empty vector which will be filled with the
     *   appropriate extraction coordinates in m/z and rt and sorted by m/z (to
     *   be used as input to extractChromatograms)
     * @param transition_exp The transition experiment used as input (is constant)
     * @param rt_extraction_window Full RT extraction window (rt_end - rt_start
     *   will equal this window size). Enforces the presence of retention times
     *   if larger than zero (throws an exception), if less than zero, rt_end
     *   will be set to -1 and rt_start to 0.
     * @param ms1 Whether to extract for MS1 (peptide level) or MS2 (transition level)
     *
    */
    void prepare_coordinates(std::vector< OpenSwath::ChromatogramPtr > & output_chromatograms,
      std::vector< ExtractionCoordinates > & coordinates,
      OpenMS::TargetedExperiment & transition_exp,
      const double rt_extraction_window,
      const bool ms1) const;

    /**
     * @brief This converts the ChromatogramPtr to MSChromatogram and adds meta-information.
     *
     * It sets
     * 1) the target m/z
     * 2) the isolation window (upper/lower)
     * 3) the peptide sequence
     * 4) the fragment m/z
     * 5) the meta-data, e.g. InstrumentSettings, AcquisitionInfo, 
     *     sourceFile and DataProcessing
     * 6) the native ID from the transition
     *
     */
    template <typename TransitionExpT>
    static void return_chromatogram(std::vector< OpenSwath::ChromatogramPtr > & chromatograms,
      std::vector< ChromatogramExtractor::ExtractionCoordinates > & coordinates,
      TransitionExpT& transition_exp_used, SpectrumSettings settings,
      std::vector<OpenMS::MSChromatogram<> > & output_chromatograms, bool ms1)
    {
      typedef std::map<String, const typename TransitionExpT::Transition* > TransitionMapType;
      TransitionMapType trans_map;
      for (Size i = 0; i < transition_exp_used.getTransitions().size(); i++)
      {
        trans_map[transition_exp_used.getTransitions()[i].getNativeID()] = &transition_exp_used.getTransitions()[i];
      }

      for (Size i = 0; i < chromatograms.size(); i++)
      { 
        const OpenSwath::ChromatogramPtr & chromptr = chromatograms[i];
        const ChromatogramExtractor::ExtractionCoordinates & coord = coordinates[i];

        typename TransitionExpT::Peptide pep;
        typename TransitionExpT::Transition transition;
        OpenMS::MSChromatogram<> chrom;

        // copy data
        OpenSwathDataAccessHelper::convertToOpenMSChromatogram(chrom, chromptr);
        chrom.setNativeID(coord.id);

        // Create precursor and set
        // 1) the target m/z
        // 2) the isolation window (upper/lower)
        // 3) the peptide sequence
        Precursor prec;
        if (ms1) 
        {
          pep = transition_exp_used.getPeptideByRef(coord.id); 
          prec.setMZ(coord.mz);
          chrom.setChromatogramType(ChromatogramSettings::BASEPEAK_CHROMATOGRAM);
        }
        else 
        {
          transition = (*trans_map[coord.id]);
          pep = transition_exp_used.getPeptideByRef(transition.getPeptideRef()); 

          prec.setMZ(transition.getPrecursorMZ());
          if (settings.getPrecursors().size() > 0)
          {
            prec.setIsolationWindowLowerOffset(settings.getPrecursors()[0].getIsolationWindowLowerOffset());
            prec.setIsolationWindowUpperOffset(settings.getPrecursors()[0].getIsolationWindowUpperOffset());
          }

          // Create product and set its m/z
          Product prod;
          prod.setMZ(transition.getProductMZ());
          chrom.setProduct(prod);
          chrom.setChromatogramType(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
        }
        prec.setMetaValue("peptide_sequence", pep.sequence);
        chrom.setPrecursor(prec);

        // Set the rest of the meta-data
        chrom.setInstrumentSettings(settings.getInstrumentSettings());
        chrom.setAcquisitionInfo(settings.getAcquisitionInfo());
        chrom.setSourceFile(settings.getSourceFile());

        for (Size i = 0; i < settings.getDataProcessing().size(); ++i)
        {
          DataProcessing dp = settings.getDataProcessing()[i];
          dp.setMetaValue("performed_on_spectra", "true");
          chrom.getDataProcessing().push_back(dp);
        }
        output_chromatograms.push_back(chrom);
      }
    }

    template <typename SpectrumT>
    void extract_value_tophat(const SpectrumT& input, const double& mz, Size& peak_idx,
        double& integrated_intensity, const double& extract_window, const bool ppm)
    {
      integrated_intensity = 0;
      if (input.size() == 0)
      {
        return;
      }

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

      // walk to the left until we go outside the window, then start walking to the right until we are outside the window
      walker = peak_idx;
      if (walker > 0)
      {
        walker--;
      }
      while (walker > 0 && input[walker].getMZ() > left && input[walker].getMZ() < right)
      {
        integrated_intensity += input[walker].getIntensity(); walker--;
      }
      walker = peak_idx;
      if (walker < input.size() )
      {
        walker++;
      }
      while (walker<input.size() && input[walker].getMZ()> left &&  input[walker].getMZ() < right)
      {
        integrated_intensity += input[walker].getIntensity(); walker++;
      }
    }

    template <typename SpectrumT>
    void extract_value_bartlett(const SpectrumT& input, const double& mz, Size& peak_idx,
        double& integrated_intensity, const double& extract_window, const bool ppm)
    {
      integrated_intensity = 0;
      if (input.size() == 0)
      {
        return;
      }

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

      // walk to the left until we go outside the window, then start walking to the right until we are outside the window
      walker = peak_idx;
      if (walker > 0 )
      {
        walker--;
      }
      while (walker > 0 && input[walker].getMZ() > left && input[walker].getMZ() < right)
      {
        weight =  1 - fabs(input[walker].getMZ() - mz) / half_window_size;
        integrated_intensity += input[walker].getIntensity() * weight; walker--;
      }
      walker = peak_idx;
      if (walker < input.size() )
      {
        walker++;
      }
      while (walker<input.size() && input[walker].getMZ()> left &&  input[walker].getMZ() < right)
      {
        weight = 1 - fabs(input[walker].getMZ() - mz) / half_window_size;
        integrated_intensity += input[walker].getIntensity() * weight; walker++;
      }
    }

private:

    /**
     * @brief This populates the chromatograms vector with empty chromatograms
     * (but sets their meta-information)
     *
     * It extracts
     * 1) the target m/z
     * 2) the isolation window (upper/lower)
     * 3) the peptide sequence
     * 4) the fragment m/z
     * 5) Copy the meta-data, e.g. InstrumentSettings, AcquisitionInfo, 
     *     sourceFile and DataProcessing
     * 6) the native ID from the transition
     *
     */
    template <class SpectrumSettingsT, class ChromatogramT>
    void prepareSpectra_(SpectrumSettingsT& settings, std::vector<ChromatogramT>& chromatograms, OpenMS::TargetedExperiment& transition_exp)
    {

      // first prepare all the spectra (but leave them empty)
      for (Size i = 0; i < transition_exp.getTransitions().size(); i++)
      {
        const ReactionMonitoringTransition* transition = &transition_exp.getTransitions()[i];

        // 1) and 2) Extract precursor m/z and isolation window
        ChromatogramT chrom;
        Precursor prec;
        prec.setMZ(transition->getPrecursorMZ());
        if (settings.getPrecursors().size() > 0)
        {
          prec.setIsolationWindowLowerOffset(settings.getPrecursors()[0].getIsolationWindowLowerOffset());
          prec.setIsolationWindowUpperOffset(settings.getPrecursors()[0].getIsolationWindowUpperOffset());
        }

        // 3) set precursor peptide sequence
        String pepref = transition->getPeptideRef();
        for (Size pep_idx = 0; pep_idx < transition_exp.getPeptides().size(); pep_idx++)
        {
          const OpenMS::TargetedExperiment::Peptide* pep = &transition_exp.getPeptides()[pep_idx];
          if (pep->id == pepref)
          {
            prec.setMetaValue("peptide_sequence", pep->sequence);
            break;
          }
        }
        // add precursor to spectrum
        chrom.setPrecursor(prec);

        // 4) Create product and set its m/z
        Product prod;
        prod.setMZ(transition->getProductMZ());
        chrom.setProduct(prod);

        // 5) Set the rest of the meta-data
        chrom.setInstrumentSettings(settings.getInstrumentSettings());
        chrom.setAcquisitionInfo(settings.getAcquisitionInfo());
        chrom.setSourceFile(settings.getSourceFile());

        for (Size i = 0; i < settings.getDataProcessing().size(); ++i)
        {
          DataProcessing dp = settings.getDataProcessing()[i];
          dp.setMetaValue("performed_on_spectra", "true");
          chrom.getDataProcessing().push_back(dp);
        }

        // Set the id of the chromatogram, using the id of the transition (this gives directly the mapping of the two)
        chrom.setNativeID(transition->getNativeID());
        chrom.setChromatogramType(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
        chromatograms.push_back(chrom);
      }

    }

    bool outsideExtractionWindow_(const ReactionMonitoringTransition& transition, double current_rt,
                                   const TransformationDescription& trafo, double rt_extraction_window);

    int getFilterNr(String filter);

    void populatePeptideRTMap_(OpenMS::TargetedExperiment& transition_exp, double rt_extraction_window);

    std::map<OpenMS::String, double> PeptideRTMap_;

  };

}

#endif
