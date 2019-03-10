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
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>


namespace OpenMS
{

  /**
    * @ingroup Chemistry
    * @brief Isotope pattern generator for fine isotope distributions.
    *
    * NOTE: This interface is DEPRECATED, and kept for backward compatibility
    * only. Please use the classes defined in IsoSpecWrapper.h directly instead.
    *
    * This algorithm generates theoretical pattern distributions for empirical
    * formulas with high resolution. The output is a list of pairs containing
    * isotope probabilities paired with the accurate m/z for the analyte
    * isotopic composition.
    *
    * For example, for a C100 molecule, you will get:
    *
    * @code
    *     1200            : 0.338014274835587
    *     1201.0033548352 : 0.368628561496735
    *     1202.0067096704 : 0.198997721076012
    *     1203.0100645056 : 0.0708935707807541
    * @endcode
    *
    * One important value to set is the threshold with tells the algorithm when
    * to stop calculating isotopic peaks to calculate. Here, a threshold of
    * 0.01 would mean that the algorithm either stops calculating when any new
    * peak would be less than 0.01 in height (absolute) or when it would be
    * less than 0.01 of the highest isotopic peak (relative). This is how the
    * stop_condition parameter is interpreted when use_total_prob is set to false.
    *
    * Another way to parametrise the search is to cover a certain portion
    * of total probability, such as 0.99. This is what most users will want to use;
    * see the documentation on IsoSpecTotalProbWrapper and IsoSpecThresholdWrapper
    * for the details (also, you should use those classes directly for increased
    * performance). This is how the stop_condition parameter is interpreted when
    * use_total_prob is set to true.
    *
    * @note Although use_total_prob is set to false by default for backward
    *       compatibility, most users should set it to true and parametrise by
    *       total probability needed.
    *
    * @note absolute parameter is ignored when use_total_prob is set to true
    *
    * @note When use_total_prob is set to false, common sensible values of
    *       stop_condition will be close to 0.0 (like 0.01), when it's set
    *       to true - close to 1.0, like 0.99.
    *
    * @note Computation of fine isotope patterns can be slow for large
    * molecules, if you don't need fine isotope distributions consider using
    * CoarseIsotopePatternGenerator.
    * @note Consider using IsoSpec directly for increased performance.
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
      * @param absolute Whether threshold is absolute or relative (see class docu)
      * @param use_total_prob Whether the stop_condition should be interpreted as a
      *        probability threshold (only configurations with intensity above this
      *        threshold will be returned) or as a total probability that the distribution
      *        should cover.
      *
      **/
    FineIsotopePatternGenerator(double stop_condition, bool absolute = false, bool use_total_prob = false) :
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
      * down processing, consider using IsoSpec directly for increased
      * performance.
      *
      **/
    IsotopeDistribution run(const EmpiricalFormula&) const;

    /// Set probability stop condition
    void setThreshold(double stop_condition)
    {
      stop_condition_ = stop_condition;
    }

    /// Get probability stop condition
    double getThreshold()
    {
      return stop_condition_;
    }

    /// Set whether threshold is absolute or relative probability
    void setAbsolute(bool absolute)
    {
      absolute_ = absolute;
    }

    /// Returns whether threshold is absolute or relative probability
    bool getAbsolute()
    {
      return absolute_;
    }

 protected:
    double stop_condition_ = 0.01;
    bool absolute_ = false;
    bool use_total_prob_ = false;

  };

} // namespace OpenMS

