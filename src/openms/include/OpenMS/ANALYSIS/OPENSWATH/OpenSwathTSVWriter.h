// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

// Interfaces
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <OpenMS/KERNEL/FeatureMap.h>

#include <fstream>

namespace OpenMS
{

  /**
    @brief Class to write out an OpenSwath TSV output (mProphet input).

    The output is organized as a set of rows where each row describes a single
    peak in a chromatogram with a start and end retention time. Per
    chromatogram multiple rows can be reported if more than one potential peak
    was found. See also OpenSwathOSWWriter for another output format.

    The class can take a FeatureMap and create a set of string from it
    suitable for output to tsv using the prepareLine() function.
    These lines can also be directly written to a file using writeLines().

    The output format written by this class is a TSV file (tab-separated plain
    text file) with the following columns:

      <table>
        <tr> <td BGCOLOR="#EBEBEB">transition_group_id</td> <td>string</td> <td> designates the transition group [e.g. peptide/molecule] to which this transition belongs</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">peptide_group_label</td> <td>string</td> <td> designates to which peptide label group (as defined in MS:1000893) the peptide belongs to<sup>2</sup></td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">run_id</td> <td>integer</td> <td> LC-MS/MS run (currently always 0)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">filename</td> <td>string</td> <td> filename of the raw LC-MS/MS file </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">RT</td> <td>float</td> <td> peak group retention time (apex)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">id</td> <td>string</td> <td> unique peak group identifier </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">Sequence</td> <td>string</td> <td> Peptide sequence </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">MC</td> <td>int</td> <td> Number of missed cleavages </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">FullPeptideName</td> <td>string</td> <td> Full peptide sequence including modifications <sup>1</sup></td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">Charge</td> <td>string</td> <td> Peptide (analyte) precursor charge </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">m/z</td> <td>string</td> <td> Peptide (analyte) precursor m/z </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">Intensity</td> <td>float</td> <td> Peptide (analyte) intensity (sum over all fragment ion intensities)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">ProteinName</td> <td>string</td> <td> Protein Identifier </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">GeneName</td> <td>string</td> <td> Gene identifier </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">decoy</td> <td>string</td> <td> Whether peak group was found in a decoy chromatogram (0 = false) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">assay_rt</td> <td>string</td> <td> The retention time at which the peptide (analyte) was expected to elute based on the retention time calibration</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">delta_rt</td> <td>string</td> <td> The difference in retention between expected retention time (assay_rt) and peak group retention time (RT) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">leftWidth</td> <td>float</td> <td> Retention time start of the peak (left width) in seconds</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">main_var_xx_swath_prelim_score</td> <td>float</td> <td> Preliminary separation for pyProphet initialization </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">norm_RT</td> <td>string</td> <td> The position of the peak group in the normalized retention time space (e.g. fx(RT) where fx describes the transformation fx)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">nr_peaks</td> <td>int</td> <td> The number of transitions used </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">peak_apices_sum</td> <td>float</td> <td> The sum of the peak apices intensities </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">potentialOutlier</td> <td>string</td> <td> Potential outlier transition </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">initialPeakQuality</td> <td>float</td> <td> Initial peak quality score (if computing peak quality was enabled) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">rightWidth</td> <td>string</td> <td> Retention time end of the peak (left width) in seconds </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">rt_score</td> <td>string</td> <td> sequence </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">sn_ratio</td> <td>float</td> <td> Signal-to-Noise ratio </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">total_xic</td> <td>float</td> <td> Total XIC </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">var_...</td> <td>float</td> <td> A variable used for the post-processing and scoring</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">aggr_prec_Peak_Area</td> <td>float-list</td> <td> MS1 %Precursor peak area (for each isotope) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">aggr_prec_Peak_Apex</td> <td>float-list</td> <td> MS1 %Precursor peak apex (for each isotope) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">aggr_prec_Annotation</td> <td>string-list</td> <td> MS1 %Precursor annotation (for each isotope) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">aggr_Peak_Area</td> <td>float-list</td> <td> Fragment ion peak area (for each transition) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">aggr_Peak_Apex</td> <td>float-list</td> <td> Fragment ion peak apex (for each transition) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">aggr_Fragment_Annotation</td> <td>string-list</td> <td> Fragment ion annotation (for each transition) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">masserror_ppm</td> <td>string</td> <td> Fragment-level mass error for each transition (see aggr_Fragment_Annotation for order)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">rt_fwhm</td> <td>string</td> <td> Fragment-level FWHM (full width half maximum) for each individual transition (see aggr_Fragment_Annotation for order)</td> </tr>
      </table>

      <p>
      Remarks:
      </p>
      <ul>
        <li>
          1. modifications should be supplied inside the sequence using UniMod
            identifiers or freetext identifiers that are understood by %OpenMS. See also @ref OpenMS::AASequence for more information. For example:
            <ul>
            <li> PEPT(Phosphorylation)IDE(UniMod:27)A ) </li>
            </ul>
        </li>
        <li>
          2. peptide label groups designate groups of peptides that are isotopically
          modified forms of the same peptide species. For example, the heavy and
          light forms of the same peptide will both be assigned the same peptide
          group label. For example:
            <ul>
            <li> PEPTIDEAK -> gets label "PEPTIDEAK_gr1"  </li>
            <li> PEPTIDEAK[+8] -> gets label "PEPTIDEAK_gr1"  </li>
            <li> PEPT(Phosphorylation)IDEAK -> gets label "PEPTIDEAK_gr2"  </li>
            <li> PEPT(Phosphorylation)IDEAK[+8] -> gets label "PEPTIDEAK_gr2"  </li>
            </ul>
        </li>
      </ul>
      </p>

   */
  class OPENMS_DLLAPI OpenSwathTSVWriter
  {
    std::ofstream ofs;
    String input_filename_;
    bool doWrite_;
    bool use_ms1_traces_;

  public:

    OpenSwathTSVWriter(const String& output_filename,
                       const String& input_filename = "inputfile",
                       bool ms1_scores = false);

    bool isActive() const;

    /**
     * @brief Initializes file by writing TSV header
     *
     */
    void writeHeader();

    /**
     * @brief Prepare a single line (feature) for output
     *
     * The result can be flushed to disk using writeLines (either line by line
     * or after collecting several lines).
     *
     * @param pep The compound (peptide/metabolite) used for extraction
     * @param transition The transition used for extraction
     * @param output The feature map containing all features (each feature will generate one entry in the output)
     * @param id The transition group identifier (peptide/metabolite id)
     *
     * @returns A string to be written using writeLines
     *
     */
    String prepareLine(const OpenSwath::LightCompound& pep,
        const OpenSwath::LightTransition * transition,
        const FeatureMap& output, const String& id) const;

    /**
     * @brief Write data to disk
     *
     * Takes a set of pre-prepared data statements from prepareLine and flushes them to disk
     *
     * @param to_output Statements generated by prepareLine to be written to a file
     *
     * @note Only call inside an OpenMP critical section
     *
     */
    void writeLines(const std::vector<String>& to_output);

  };

}


