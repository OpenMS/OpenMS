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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>

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
     * @param input The input spectra from which to extract chromatograms
     * @param output The output vector in which to store the chromatograms
     * @param transition_exp The extraction coordinates (m/z, RT, ion mobility)
     * @param mz_extraction_window Extracts a window of this size in m/z
     *                             dimension (e.g. a window of 50 ppm means an
     *                             extraction of 25 ppm on either side)
     * @param ppm Whether mz windows in in ppm
     * @param trafo A transformation description for RT space
     * @param rt_extraction_window Extracts a window of this size in RT
     *                             dimension (e.g. a window of 600 seconds
     *                             means an extraction of 300 seconds on either
     *                             side)
     * @param filter Which filter to use (bartlett or tophat)
     *
     * @note: it will replace chromatograms in the output map, not append to them!
     * @note: TODO deprecate this function (use ChromatogramExtractorAlgorithm instead)
    */
    template <typename ExperimentT>
    void extractChromatograms(const ExperimentT& input,
                              ExperimentT& output,
                              OpenMS::TargetedExperiment& transition_exp,
                              double mz_extraction_window,
                              bool ppm, TransformationDescription trafo,
                              double rt_extraction_window,
                              const String& filter)
    {
      // invert the trafo because we want to transform nRT values to "real" RT values
      trafo.invert();

      Size input_size = input.size();
      if (input_size < 1)
      {
        return;
      }

      int used_filter = getFilterNr_(filter);
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
     * @param input The input spectra from which to extract chromatograms
     * @param output The output vector in which to store the chromatograms
     * (needs to be of the same length as the extraction coordinates, use
     * prepare_coordinates)
     * @param extraction_coordinates The extraction coordinates (m/z, RT, ion mobility)
     * @param mz_extraction_window Extracts a window of this size in m/z
     * dimension (e.g. a window of 50 ppm means an extraction of 25 ppm on
     * either side)
     * @param ppm Whether mz windows in in ppm
     * @param filter Which filter to use (bartlett or tophat)
     *
     * @note: whenever possible, please use this ChromatogramExtractorAlgorithm implementation
     *
    */
    void extractChromatograms(const OpenSwath::SpectrumAccessPtr input, 
                              std::vector< OpenSwath::ChromatogramPtr >& output,
                              const std::vector<ExtractionCoordinates>& extraction_coordinates,
                              double mz_extraction_window,
                              bool ppm,
                              const String& filter)
    {
      ChromatogramExtractorAlgorithm().extractChromatograms(input, output, 
          extraction_coordinates, mz_extraction_window, ppm, -1, filter);
    }

    /**
     * @brief Extract chromatograms at the m/z and RT defined by the ExtractionCoordinates.
     *
     * @param input The input spectra from which to extract chromatograms
     * @param output The output vector in which to store the chromatograms
     * (needs to be of the same length as the extraction coordinates, use
     * prepare_coordinates)
     * @param extraction_coordinates The extraction coordinates (m/z, RT, ion mobility)
     * @param mz_extraction_window Extracts a window of this size in m/z
     * dimension (e.g. a window of 50 ppm means an extraction of 25 ppm on
     * either side)
     * @param ppm Whether mz windows in in ppm
     * @param im_extraction_window Extracts a window of this size in ion mobility
     * @param filter Which filter to use (bartlett or tophat)
     *
     * @note: whenever possible, please use this ChromatogramExtractorAlgorithm implementation
     *
    */
    void extractChromatograms(const OpenSwath::SpectrumAccessPtr input,
                              std::vector< OpenSwath::ChromatogramPtr >& output,
                              const std::vector<ExtractionCoordinates>& extraction_coordinates,
                              double mz_extraction_window,
                              bool ppm,
                              double im_extraction_window,
                              const String& filter) 
    {
      ChromatogramExtractorAlgorithm().extractChromatograms(input, output, 
          extraction_coordinates, mz_extraction_window, ppm, im_extraction_window, filter);
    }

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
     * @param rt_extraction_window If non-negative, full RT extraction window,
     *   centered on the first RT value (@p rt_end - @p rt_start will equal this
     *   window size). If negative, @p rt_end will be set to -1 and @p rt_start
     *   to 0 (i.e. full RT range). If NaN, exactly two RT entries are expected
     *   - the first is used as @p rt_start and the second as @p rt_end.
     * @param ms1 Whether to extract for MS1 (peptide level) or MS2 (transition level)
     *
     * @throw Exception::IllegalArgument if RT values are expected (depending on @p rt_extraction_window) but not provided
    */
    static void prepare_coordinates(std::vector< OpenSwath::ChromatogramPtr > & output_chromatograms,
                                    std::vector< ExtractionCoordinates > & coordinates,
                                    const OpenMS::TargetedExperiment & transition_exp,
                                    const double rt_extraction_window,
                                    const bool ms1 = false,
                                    const int ms1_isotopes = 0);

    static void prepare_coordinates(std::vector< OpenSwath::ChromatogramPtr > & output_chromatograms,
                                    std::vector< ExtractionCoordinates > & coordinates,
                                    const OpenSwath::LightTargetedExperiment & transition_exp_used,
                                    const double rt_extraction_window,
                                    const bool ms1 = false,
                                    const int ms1_isotopes = 0);

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
     * 7) ion mobility extraction target and window (lower/upper)
     *
     */
    template <typename TransitionExpT>
    static void return_chromatogram(const std::vector< OpenSwath::ChromatogramPtr > & chromatograms,
                                    const std::vector< ChromatogramExtractor::ExtractionCoordinates > & coordinates,
                                    TransitionExpT& transition_exp_used,
                                    SpectrumSettings settings,
                                    std::vector<OpenMS::MSChromatogram > & output_chromatograms,
                                    bool ms1,
                                    double im_extraction_width = 0.0)
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

        // copy data
        OpenMS::MSChromatogram chrom;
        OpenSwathDataAccessHelper::convertToOpenMSChromatogram(chromptr, chrom);
        chrom.setNativeID(coord.id);

        // Create precursor and set
        // 1) the target m/z
        // 2) the isolation window (upper/lower)
        // 3) the peptide sequence
        Precursor prec;
        if (ms1) 
        {
          prec.setMZ(coord.mz);
          chrom.setChromatogramType(ChromatogramSettings::BASEPEAK_CHROMATOGRAM);

          // extract compound / peptide id from transition and store in
          // more-or-less default field
          String transition_group_id = OpenSwathHelper::computeTransitionGroupId(coord.id);
          if (!transition_group_id.empty())
          {
            int prec_charge = 0;
            String r = extract_id_(transition_exp_used, transition_group_id, prec_charge);
            prec.setCharge(prec_charge);
            prec.setMetaValue("peptide_sequence", r);
          }
        }
        else 
        {
          typename TransitionExpT::Transition transition = (*trans_map[coord.id]);

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

          // extract compound / peptide id from transition and store in
          // more-or-less default field
          if (!transition.getPeptideRef().empty())
          {
            int prec_charge = 0;
            String r = extract_id_(transition_exp_used, transition.getPeptideRef(), prec_charge);
            prec.setCharge(prec_charge);
            prec.setMetaValue("peptide_sequence", r);
          }
          else
          {
            int prec_charge = 0;
            String r = extract_id_(transition_exp_used, transition.getCompoundRef(), prec_charge);
            prec.setCharge(prec_charge);
            prec.setMetaValue("peptide_sequence", r);
          }
        }

        if (coord.ion_mobility >= 0 && im_extraction_width > 0.0)
        {
          prec.setDriftTime(coord.ion_mobility);
          prec.setDriftTimeWindowLowerOffset(im_extraction_width / 2.0);
          prec.setDriftTimeWindowUpperOffset(im_extraction_width / 2.0);
        }
        chrom.setPrecursor(prec);

        // Set the rest of the meta-data
        chrom.setInstrumentSettings(settings.getInstrumentSettings());
        chrom.setAcquisitionInfo(settings.getAcquisitionInfo());
        chrom.setSourceFile(settings.getSourceFile());

        for (Size j = 0; j < settings.getDataProcessing().size(); ++j)
        {
          settings.getDataProcessing()[j]->setMetaValue("performed_on_spectra", "true");
          chrom.getDataProcessing().push_back(settings.getDataProcessing()[j]);
        }
        output_chromatograms.push_back(chrom);
      }
    }

     /// @note: TODO deprecate this function (use ChromatogramExtractorAlgorithm instead)
    template <typename SpectrumT>
    void extract_value_tophat(const SpectrumT& input,
                              const double mz,
                              Size& peak_idx,
                              double& integrated_intensity,
                              const double extract_window,
                              const bool ppm)
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
      // if we moved past the end of the spectrum, we need to try the last peak
      // of the spectrum (it could still be within the window)
      if (peak_idx >= input.size())
      {
        walker = input.size() - 1;
      }

      // add the current peak if it is between right and left
      if (input[walker].getMZ() > left && input[walker].getMZ() < right)
      {
        integrated_intensity += input[walker].getIntensity();
      }

      // (i) Walk to the left one step and then keep walking left until we go
      // outside the window. Note for the first step to the left we have to
      // check for the walker becoming zero.
      walker = peak_idx;
      if (walker > 0)
      {
        walker--;
        // special case: walker is now zero
        if (walker == 0 && input[walker].getMZ() > left && input[walker].getMZ() < right)
        {
          integrated_intensity += input[walker].getIntensity();
        }
      }
      while (walker > 0 && input[walker].getMZ() > left && input[walker].getMZ() < right)
      {
        integrated_intensity += input[walker].getIntensity(); walker--;
      }

      // (ii) Walk to the right one step and then keep walking right until we
      // are outside the window
      walker = peak_idx;
      if (walker < input.size() )
      {
        walker++;
      }
      while (walker < input.size() && input[walker].getMZ() > left &&  input[walker].getMZ() < right)
      {
        integrated_intensity += input[walker].getIntensity(); walker++;
      }
    }

     /// @note: TODO deprecate this function (use ChromatogramExtractorAlgorithm instead)
    template <typename SpectrumT>
    void extract_value_bartlett(const SpectrumT& input,
                                const double mz,
                                Size& peak_idx,
                                double& integrated_intensity,
                                const double extract_window,
                                const bool ppm)
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
      // if we moved past the end of the spectrum, we need to try the last peak
      // of the spectrum (it could still be within the window)
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

      // (i) Walk to the left one step and then keep walking left until we go
      // outside the window. Note for the first step to the left we have to
      // check for the walker becoming zero.
      walker = peak_idx;
      if (walker > 0)
      {
        walker--;
        // special case: walker is now zero
        if (walker == 0 && input[walker].getMZ() > left && input[walker].getMZ() < right)
        {
          integrated_intensity += input[walker].getIntensity();
        }
      }
      while (walker > 0 && input[walker].getMZ() > left && input[walker].getMZ() < right)
      {
        weight =  1 - fabs(input[walker].getMZ() - mz) / half_window_size;
        integrated_intensity += input[walker].getIntensity() * weight; walker--;
      }

      // (ii) Walk to the right one step and then keep walking right until we
      // are outside the window
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
     * @brief Extracts id (peptide sequence or compound name) for a compound
     *
     * @param transition_exp The transition experiment used as input (is constant) and either of type LightTargetedExperiment or TargetedExperiment
     * @param id The identifier of the compound or peptide
     * @param prec_charge The charge state of the precursor
     *
     */
    template <typename TransitionExpT>
    static String extract_id_(TransitionExpT& transition_exp_used, const String& id, int& prec_charge);

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
    void prepareSpectra_(SpectrumSettingsT& settings,
                         std::vector<ChromatogramT>& chromatograms,
                         OpenMS::TargetedExperiment& transition_exp)
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

        // 3) set precursor peptide sequence / compound id in more-or-less default field
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
        String compref = transition->getCompoundRef();
        for (Size comp_idx = 0; comp_idx < transition_exp.getCompounds().size(); comp_idx++)
        {
          const OpenMS::TargetedExperiment::Compound* comp = &transition_exp.getCompounds()[comp_idx];
          if (comp->id == compref)
          {
            prec.setMetaValue("peptide_sequence", String(comp->id) );
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

        for (Size j = 0; j < settings.getDataProcessing().size(); ++j)
        {
          settings.getDataProcessing()[j]->setMetaValue("performed_on_spectra", "true");
          chrom.getDataProcessing().push_back(settings.getDataProcessing()[j]);
        }

        // Set the id of the chromatogram, using the id of the transition (this gives directly the mapping of the two)
        chrom.setNativeID(transition->getNativeID());
        chrom.setChromatogramType(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
        chromatograms.push_back(chrom);
      }

    }

     /// @note: TODO deprecate this function (use ChromatogramExtractorAlgorithm instead)
    bool outsideExtractionWindow_(const ReactionMonitoringTransition& transition,
                                  double current_rt,
                                  const TransformationDescription& trafo,
                                  double rt_extraction_window);

     /// @note: TODO deprecate this function (use ChromatogramExtractorAlgorithm instead)
    int getFilterNr_(const String& filter);

     /// @note: TODO deprecate this function (use ChromatogramExtractorAlgorithm instead)
    void populatePeptideRTMap_(OpenMS::TargetedExperiment& transition_exp,
                               double rt_extraction_window);

    std::map<OpenMS::String, double> PeptideRTMap_;

  };
    
  // Specialization for template (LightTargetedExperiment)
  template<>
  inline String ChromatogramExtractor::extract_id_<OpenSwath::LightTargetedExperiment>(OpenSwath::LightTargetedExperiment& transition_exp_used,
                                                                                       const String& id,
                                                                                       int & prec_charge)
  {
    const OpenSwath::LightCompound comp = transition_exp_used.getCompoundByRef(id);
    prec_charge = comp.charge;
    if (!comp.sequence.empty())
    {
      return comp.sequence;
    }
    else
    {
      return comp.compound_name;
    }
  }


  // Specialization for template (TargetedExperiment)
  template<>
  inline String ChromatogramExtractor::extract_id_<OpenMS::TargetedExperiment>(OpenMS::TargetedExperiment& transition_exp_used,
                                                                               const String& id,
                                                                               int & prec_charge)
  {
    if (transition_exp_used.hasPeptide(id))
    {
      const TargetedExperiment::Peptide p = transition_exp_used.getPeptideByRef(id);
      if (p.hasCharge()) {prec_charge = p.getChargeState();}
      return p.sequence;
    }
    else if (transition_exp_used.hasCompound(id))
    {
      const TargetedExperiment::Compound c = transition_exp_used.getCompoundByRef(id);
      if (c.hasCharge()) {prec_charge = c.getChargeState();}
      return c.id;
    }
    else
    {
      return "";
    }
  }
}

