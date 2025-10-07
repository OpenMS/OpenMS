// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Oliver Alka $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <OpenMS/METADATA/Precursor.h>

#include <string>
#include <vector>
#include <set>

namespace OpenMS
{
  /**
  @brief This class provides methods for precursor correction.

  Supported methods:
  getPrecursors: Extract precursors and associated information (mz, scan index).
  writeHist: Write output .csv for validation purposes (corrected, uncorrected).
  correctToNearestMS1Peak: Correct to the peak in closest proximity in a certain mass range.
  correctToHighestIntensityMS1Peak: Correct to the peak with the highest intensity in a certain mass range.
  correctToNearestFeature: Use feature information to re-annotate a precursor (e.g. falsely assigned to non mono-isotopic trace).
  */
class OPENMS_DLLAPI PrecursorCorrection
{
     public:

     static const std::string csv_header;

     /**
     @brief Extract precursors and associated information (precursor retention time and precursor scan index).
     @param exp: Spectra with precursors
     @param[out] precursors: vector of all precursors in @p exp (can be more than one per MSn spectrum)
     @param[out] precursors_rt: vector double of precursors retention time (same length as @p precursors)
     @param[out] precursor_scan_index: Indices into @p exp, which have a precursor
     */
     static void getPrecursors(const MSExperiment & exp,
                              std::vector<Precursor> & precursors,
                              std::vector<double> & precursors_rt,
                              std::vector<Size> & precursor_scan_index);


     /**
     @brief Writer can be used in association with correctToNearestMS1Peak or correctToHighestIntensityMS1Peak.
     A csv file with additional information (RT, uncorrectedMZ, correctedMZ, deltaMZ).

     Format:
     RT	    uncorrectedMZ	correctedMZ	deltaMZ
     100.1	509.9999	    510	         0.0001
     180.9	610.0001	    610	        -0.0001
     183.92	611.0035	    611.0033	  -0.0002
     
     @param out_csv: constant String for csv output.
     @param delta_mzs: delta m/z column values.
     @param mzs: m/z column vector (uncorrectedMZ)
     @param rts: retention time column vector 
     
     */
     static void writeHist(const String& out_csv,
                           const std::vector<double> & delta_mzs,
                           const std::vector<double> & mzs,
                           const std::vector<double> & rts);
     /**
     @brief Selection of the peak in closest proximity as corrected precursor mass in a given mass range (e.g. precursor mass +/- 0.2 Da).

     For each MS2 spectrum the corresponding MS1 spectrum is determined by using the rt information of the precursor.
     In the MS1, the peak closest to the uncorrected precursor m/z is selected and used as corrected precursor m/z.

     @param exp: MSExperiment.
     @param mz_tolerance: double tolerance used for precursor correction in mass range.
     @param ppm: bool enables usage of ppm.
     @param delta_mzs: vector double delta mass to charge.
     @param mzs: vector double mass to charge.
     @param rts: vector double retention time.
     @return set of Size with corrected precursor information.
     */
     static std::set<Size> correctToNearestMS1Peak(MSExperiment & exp,
                                                   double mz_tolerance,
                                                   bool ppm,
                                                   std::vector<double> & delta_mzs,
                                                   std::vector<double> & mzs,
                                                   std::vector<double> & rts);

     /**
     @brief Selection of the peak with the highest intensity as corrected precursor mass in a given mass range (e.g. precursor mass +/- 0.2 Da)

     For each MS2 spectrum the corresponding MS1 spectrum is determined by using the rt information of the precursor.
     In the MS1, the peak with the highest intensity in a given mass range to the uncorrected precursor m/z is selected and used as corrected precursor m/z.

     @param exp: MSExperiment.
     @param mz_tolerance: double tolerance used for precursor correction in mass range.
     @param ppm: bool enables usage of ppm.
     @param delta_mzs: vector double delta mass to charge.
     @param mzs: vector double mass to charge.
     @param rts: vector double retention time.
     @return set of Size with corrected precursor information.
     */
     static std::set<Size> correctToHighestIntensityMS1Peak(MSExperiment & exp,
                                                            double mz_tolerance,
                                                            bool ppm,
                                                            std::vector<double> & delta_mzs,
                                                            std::vector<double> & mzs,
                                                            std::vector<double> & rts);


     /**
     @brief Reassigns a precursor to the nearest feature in a given rt and mass range.
     Wrong assignment of the mono-isotopic mass for precursors are assumed:
     - if precursor_mz matches the mz of a non-monoisotopic feature mass trace
     - and in the case that believe_charge is true: if feature_charge matches the precursor_charge
     In the case of wrong mono-isotopic assignment several options for correction are available:
     keep_original will create a copy of the precursor and tandem spectrum for the new mono-isotopic mass trace and retain the original one.
     all_matching_features does this not for only the closest feature but all features in a question.

     @param features: constant FeatureMap.
     @param exp: MSExperiment.
     @param rt_tolerance_s: double retention time tolerance in seconds.
     @param mz_tolerance: double tolerance used for precursor correction in mass range.
     @param ppm: bool enables usage of ppm.
     @param believe_charge: bool only add features that match the precursor charge.
     @param keep_original: bool this will create a copy of the precursor and tandem spectrum for the new mono-isotopic trace and retain the original one.
     @param all_matching_features: bool correction is performed for all features in question not only the closest one.
     @param max_trace: integer maximal number of traces used.
     @param debug_level: integer debug level.
     @return set of Size with corrected precursor information.


     */
     static std::set<Size> correctToNearestFeature(const FeatureMap& features,
                                                   MSExperiment & exp,
                                                   double rt_tolerance_s = 0.0,
                                                   double mz_tolerance = 0.0,
                                                   bool ppm = true,
                                                   bool believe_charge = false,
                                                   bool keep_original = false,
                                                   bool all_matching_features = false,
                                                   int max_trace = 2,
                                                   int debug_level = 0);

  protected:

      /**
      @brief Check if precursor is located in the bounding box of a features convex hull.
      Here the bounding box of the feature is extended by the retention time tolerance and
      afterwards the precursor location is validated.

      @param feature: constant Feature.
      @param rt: constant double retention time.
      @param pc_mz: constant double precursor mass to charge.
      @param rt_tolerance: constant double retention time tolerance in seconds.
      @return static boolean to check if the precursor is located in the bounding box of a features convex hull.
      */
      static bool overlaps_(const Feature& feature,
                            const double rt,
                            const double pc_mz,
                            const double rt_tolerance);

      /**
      @brief Check precursor and feature compatibility
      If the precursor mz is in one of the masstraces the feature is compatible.
      Dependent on 13C mass difference and charge.

      @param feature: constant Feature.
      @param pc_mz: double precursor mass to charge.
      @param mz_tolerance: double mass to charge tolerance.
      @param max_trace_number: Size maximum number of mass traces.
      @param debug_level: integer debug level.
      @return static boolean if the precursor mass to charge is in one of the features masstraces.
      */
      static bool compatible_(const Feature& feature,
                              double pc_mz,
                              double mz_tolerance,
                              Size max_trace_number = 2,
                              int debug_level = 0);

};
}
