// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
   * @brief The ChromatogramExtractor extracts chromatograms (intensity vs retention time) from mass spectrometry data.
   *
   * This class provides functionality to extract chromatographic traces from mass spectrometry data
   * based on specified coordinates (m/z, retention time, and optionally ion mobility values).
   * 
   * The extractor supports two main interfaces:
   * 1. Legacy interface: Takes a TargetedExperiment object containing transitions and extracts 
   *    chromatograms at the m/z values specified in those transitions.
   * 2. Modern interface: Takes a set of ExtractionCoordinates that specify the exact coordinates 
   *    for extraction. This provides more flexibility and control over the extraction process.
   *    The prepare_coordinates() helper function can generate these coordinates for common 
   *    MS1 and MS2 extraction scenarios.
   *
   * Key features:
   * - Supports both MS1 and MS2 level extractions
   * - Configurable extraction window sizes in m/z dimension (absolute or ppm)
   * - Multiple filter types available (Bartlett, tophat) for signal processing
   * - Handles ion mobility data when available
   * - Optimized for SWATH/DIA (Data Independent Acquisition) experiments
   * - Progress logging capabilities through ProgressLogger base class
   *
   * For MS2 extractions, the input data is expected to come from a SWATH/DIA experiment
   * where precursor ions are fragmented in wide isolation windows, allowing for extraction
   * of fragment ion chromatograms.
   *
   * @see ChromatogramExtractorAlgorithm For the underlying extraction algorithm implementation
   * @see ExtractionCoordinates For the coordinate specification format
   */

  class OPENMS_DLLAPI ChromatogramExtractor :
    public ProgressLogger
  {

public:

    typedef ChromatogramExtractorAlgorithm::ExtractionCoordinates ExtractionCoordinates;


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
     * @param ms1_isotopes Number of isotopes to include in @p coordinates when in MS1 mode
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
          if (!settings.getPrecursors().empty())
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


private:
    /**
     * @brief Extracts id (peptide sequence or compound name) for a compound
     *
     * @param transition_exp_used The transition experiment used as input (is constant) and either of type LightTargetedExperiment or TargetedExperiment
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

