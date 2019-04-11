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
    * This algorithm generates theoretical pattern distributions for empirical
    * formulas with high resolution (while the CoarseIsotopePatternGenerator
    * will generate low-resolution patterns). The output is a list of pairs
    * containing isotope probabilities paired with the accurate m/z for the
    * analyte isotopic composition.
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
    *       increased performance.
    *
    * The computation is based on the IsoSpec algorithm
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
    IsotopeDistribution run(const EmpiricalFormula&) const;

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

    /// Set whether total probability should be computed
    void setTotalProbability(bool total)
    {
      use_total_prob_ = total;
    }

    /// Returns whether total probability should be computed
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

