// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger $
// --------------------------------------------------------------------------

#pragma once

// Interfaces
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <OpenMS/CONCEPT/UniqueIdGenerator.h>

#include <OpenMS/KERNEL/FeatureMap.h>

#include <fstream>

namespace OpenMS
{

  /**
    @brief Class to write out an OpenSwath OSW SQLite output (PyProphet input).

    The class can take a FeatureMap and create a set of string from it
    suitable for output to OSW using the prepareLine function. The SQL data is
    directly linked to the PQP file format described in the TransitionPQPFile class.
    See also OpenSwathTSVWriter for another output format.

    The file format has the following tables:

      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>RUN</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">ID</td> <td>INT</td> <td> Primary Key (run id)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">FILENAME</td> <td>TEXT</td> <td> Original filename associated with the run </td> </tr>
      </table>

      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>FEATURE</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">ID</td> <td>INT</td> <td> Primary Key (feature id)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">RUN_ID</td> <td>INT</td> <td> Foreign Key (RUN.ID)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">PRECURSOR_ID</td> <td>INT</td> <td> Foreign Key (TransitionPQPFile PRECURSOR.ID) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">EXP_RT</td> <td>REAL</td> <td>Experimental RT (retention time) of the feature </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">NORM_RT</td> <td>REAL</td> <td>Normalized RT (retention time) of the feature. The position of the peak group in the normalized retention time space (e.g. fx(RT) where fx describes the transformation fx)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">DELTA_RT</td> <td>REAL</td> <td>The difference in retention between expected retention time of the assay and the measured feature retention time (EXP_RT) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">LEFT_WIDTH</td> <td>REAL</td> <td>Retention time start of the peak (left width) in seconds</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">RIGHT_WIDTH</td> <td>REAL</td> <td>Retention time end of the peak (right width) in seconds</td> </tr>
      </table>

      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>FEATURE_MS1</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">FEATURE_ID</td> <td>INT</td> <td>Foreign Key (FEATURE.ID)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">AREA_INTENSITY</td> <td>REAL</td> <td>%Precursor intensity (area) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">APEX_INTENSITY</td> <td>REAL</td> <td>%Precursor intensity (apex) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">VAR_...</td> <td>REAL</td> <td>%Precursor score used in pyProphet </td> </tr>
      </table>

      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>FEATURE_MS2</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">FEATURE_ID</td> <td>INT</td> <td> Foreign Key (FEATURE.ID)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">AREA_INTENSITY</td> <td>REAL</td> <td>Summed fragment ion intensity (area) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">TOTAL_AREA_INTENSITY </td> <td>REAL</td> <td>Summed total XIC of the whole chromatogram </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">APEX_INTENSITY</td> <td>REAL</td> <td>Summed fragment ion intensity (apex) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">TOTAL_MI</td> <td>REAL</td> <td>Total mutual information (MI)  </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">VAR_...</td> <td>REAL</td> <td>Fragment ion score used in pyProphet </td> </tr>
      </table>

      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>FEATURE_PRECURSOR</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">FEATURE_ID</td> <td>INT</td> <td> Foreign Key (FEATURE.ID)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">ISOTOPE</td> <td>INT</td> <td>Isotope identifier </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">AREA_INTENSITY</td> <td>REAL</td> <td>%Precursor isotope ion intensity (area) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">APEX_INTENSITY</td> <td>REAL</td> <td>%Precursor isotope ion intensity (apex) </td> </tr>
      </table>

      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>FEATURE_TRANSITION</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">FEATURE_ID</td> <td>INT</td> <td> Foreign Key (FEATURE.ID)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">TRANSITION_ID</td> <td>INT</td> <td> Foreign Key (transition identifier)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">AREA_INTENSITY</td> <td>REAL</td> <td>Fragment ion intensity (area) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">TOTAL_AREA_INTENSITY </td> <td>REAL</td> <td>Total XIC of the whole chromatogram </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">APEX_INTENSITY</td> <td>REAL</td> <td>Fragment ion intensity (apex) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">TOTAL_MI</td> <td>REAL</td> <td>Total mutual information (MI)  </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">VAR_...</td> <td>REAL</td> <td>Fragment ion score used in pyProphet  </td> </tr>
      </table>

   */
  class OPENMS_DLLAPI OpenSwathOSWWriter
  {
    String output_filename_;
    String input_filename_;
    OpenMS::UInt64 run_id_;
    bool doWrite_;
    bool enable_uis_scoring_;

  public:

    OpenSwathOSWWriter(const String& output_filename,
                       const UInt64 run_id,
                       const String& input_filename = "inputfile",
                       bool uis_scores = false);

    bool isActive() const;

    /**
     * @brief Initializes file by generating SQLite tables
     *
     */
    void writeHeader();

    /**
     * @brief Prepare scores for SQLite insertion
     *
     * Some scores might not be defined, those are reported as NULL
     *
     * @param feature The feature being evaluated
     * @param score_name The name of the queried score
     *
     * @returns A string with the queried score
     *
     */
    String getScore(const Feature& feature, const std::string& score_name) const;

    /**
     * @brief Prepare concatenated scores for SQLite insertion
     *
     * Some scores might not be defined, those are reported as NULL
     *
     * @param feature The feature being evaluated
     * @param score_name The name of the queried score
     *
     * @returns A vector of strings with the queried scores
     *
     */
    std::vector<String> getSeparateScore(const Feature& feature, const std::string& score_name) const;

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
        const OpenSwath::LightTransition* transition,
        const FeatureMap& output, const String& id) const;

    /**
     * @brief Write data to disk
     *
     * Takes a set of pre-prepared data statements from prepareLine and flushes them to disk
     * 
     * @param to_osw_output Statements generated by prepareLine
     *
     * @note Try to call this function as little as possible (it opens a new
     * database connection each time)
     *
     * @note Only call inside an OpenMP critical section
     *
     */
    void writeLines(const std::vector<String>& to_osw_output);

  };

}

