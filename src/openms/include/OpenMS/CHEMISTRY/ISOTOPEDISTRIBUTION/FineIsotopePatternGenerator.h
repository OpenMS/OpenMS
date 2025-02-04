// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Rost $
// $Authors: Hannes Rost $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopePatternGenerator.h>

namespace OpenMS
{

  /**
    * @ingroup Chemistry
    *
    * @brief Isotope pattern generator for fine isotope distributions.
    *
    * This algorithm implements IsotopePatternGenerator and generates
    * theoretical pattern distributions for empirical formulas with high
    * resolution (while the CoarseIsotopePatternGenerator will generate
    * low-resolution patterns). The output is a list of pairs containing
    * isotope probabilities paired with the accurate m/z for the analyte
    * isotopic composition.
    *
    * For example, for a C100H202 molecule (at 0.01 threshold), you will get:
    *
    * @code
    *     m/z 1403.5806564438 : INT 0.333207070827484
    *     m/z 1404.5840114438 : INT 0.360387712717056
    *     m/z 1404.5869331919 : INT 0.00774129061028361
    *     m/z 1405.5873664438 : INT 0.19294385612011
    *     m/z 1405.5902881919 : INT 0.00837276969105005
    *     m/z 1406.5907214438 : INT 0.0681697279214859
    *     m/z 1406.5936431919 : INT 0.00448260130360723
    *     m/z 1407.5940764438 : INT 0.0178796537220478
    *     m/z 1407.5969981919 : INT 0.00158376491162926
    *     ...
    * @endcode
    *
    * For comparison, the CoarseIsotopePatternGenerator will produce the
    * following result for a C100H202 molecule:
    *
    * @code
    *     m/z 1403.58 INT: 0.333489
    *     m/z 1404.58 INT: 0.36844
    *     m/z 1405.59 INT: 0.201576
    *     m/z 1406.59 INT: 0.0728113
    *     m/z 1407.59 INT: 0.0195325
    *     ...
    * @endcode
    *
    * From the above example, we can see that the CoarseIsotopePatternGenerator
    * will generate a single peak at nominal mass 1404 which sums up the
    * probability of both the 13C and the 2H (deuterium) peak, while the
    * FineIsotopePatternGenerator will generate two peaks at 1404.5840 (for
    * 13C) and at 1404.5869 (for 2H). The probabilities of 36.0% and 0.77% add
    * up to 36.8% which is the same as the sum reported by the
    * CoarseIsotopePatternGenerator for the nominal mass at 1404. Note that for
    * the peak at 1405 the FineIsotopePatternGenerator only reports two out of
    * the three probabilities due to the chosen probability cutoffs.
    *
    * One important value to set is the threshold with tells the algorithm when
    * to stop calculating isotopic peaks to calculate. The default stop
    * condition is to stop when only a small portion (such as 0.01) of the
    * total probability is unexplained and the reported values cover most of
    * the probability (e.g. 0.99).
    *
    * Another way to stop the search is when any new peak would be less than
    * 0.01 in height (absolute) or when it would be less than 0.01 of the
    * highest isotopic peak (relative). This is how the stop_condition
    * parameter is interpreted when use_total_prob is set to false.
    *
    * @note Computation of fine isotope patterns can be slow for large
    *       molecules, if you don't need fine isotope distributions consider using
    *       CoarseIsotopePatternGenerator.
    *
    * @note Consider using IsoSpec directly or the OpenMS IsoSpecWrapper /
    *       IsoSpecGeneratorWrapper classes defined in IsoSpecWrapper.h for
    *       increased performance since this class will sort the resuly by m/z
    *       while the wrapper will not; sorting substantially decreases
    *       performance.
    *
    * The computation is based on the IsoSpec algorithm, please cite
    *
    * @code
    * Łącki MK, Startek M, Valkenborg D, Gambin A.
    * IsoSpec: Hyperfast Fine Structure Calculator.
    * Anal Chem. 2017 Mar 21;89(6):3272-3277. doi: 10.1021/acs.analchem.6b01459.
    * @endcode
    *
    * See also method run()
    **/
  class OPENMS_DLLAPI FineIsotopePatternGenerator
    : public IsotopePatternGenerator
  {

 public:

    /**
      * @brief Default Constructor
      *
      **/
    FineIsotopePatternGenerator() = default;

    /**
      * @brief Constructor
      *
      * @param stop_condition The total probability (if use_total_prob == true) or
      *        threshold (if use_total_prob is false) (see class docu)
      *
      * @param use_total_prob Whether the stop_condition should be interpreted as a
      *        probability threshold (only configurations with intensity above this
      *        threshold will be returned) or as a total probability that the distribution
      *        should cover.
      *
      * @param absolute Whether threshold is absolute or relative (ignored if use_total_prob is true, see class docu)
      *
      **/
    FineIsotopePatternGenerator(double stop_condition, bool use_total_prob = true, bool absolute = false) :
      stop_condition_(stop_condition),
      absolute_(absolute),
      use_total_prob_(use_total_prob)
    {}

    /**
      * @brief Creates an isotope distribution from an empirical sum formula
      *
      * Iterates through all elements, convolves them according to the number
      * of atoms from that element and sums up the result.
      *
      * @note The constructed isotope distribution is sorted by m/z which slows
      * down processing, consider using IsoSpec (IsoSpecWrapper /
      * IsoSpecGeneratorWrapper) directly for increased performance.
      *
      **/
    IsotopeDistribution run(const EmpiricalFormula&) const override;

    /// Set probability stop condition (lower values generate fewer results)
    void setThreshold(double stop_condition)
    {
      stop_condition_ = stop_condition;
    }

    /// Get probability stop condition (lower values generate fewer results)
    double getThreshold() const
    {
      return stop_condition_;
    }

    /// Set whether threshold is absolute or relative probability (ignored if use_total_prob is true, see class docu)
    void setAbsolute(bool absolute)
    {
      absolute_ = absolute;
    }

    /// Returns whether threshold is absolute or relative probability (ignored if use_total_prob is true, see class docu)
    bool getAbsolute() const
    {
      return absolute_;
    }

    /// Set whether total probability should be computed (see FineIsotopePatternGenerator() )
    void setTotalProbability(bool total)
    {
      use_total_prob_ = total;
    }

    /// Returns whether total probability should be computed (see FineIsotopePatternGenerator() )
    bool getTotalProbability() const
    {
      return use_total_prob_;
    }

 protected:
    double stop_condition_ = 0.01;
    bool absolute_ = false;
    bool use_total_prob_ = true;

  };

} // namespace OpenMS

