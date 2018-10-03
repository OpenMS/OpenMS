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
    * This algorithm generates theoretical pattern distributions for empirical
    * formulas with high resolution. The output is a list of pairs containing
    * isotope probabilities paired with the accurate m/z for the analyte
    * isotopic composition.
    *
    * For example, for a C100 molecule, you will get:
    *
    * 1200            : 0.338014274835587
    * 1201.0033548352 : 0.368628561496735
    * 1202.0067096704 : 0.198997721076012
    * 1203.0100645056 : 0.0708935707807541
    *
    * One important value to set is the threshold with tells the algorithm how
    * many isotopic peaks to calculate, e.g. a threshold of 0.01 would only
    * calculate peaks that contribute at least 1% to the total isotopic
    * abundance.
    *
    * See also method run()
    **/

  class OPENMS_DLLAPI FineIsotopePatternGenerator 
    : public IsotopePatternGenerator
  {

 public:
    FineIsotopePatternGenerator() = default;
    FineIsotopePatternGenerator(double threshold) : threshold_(threshold) {}

    /**
      * @brief Creates an isotope distribution from an empirical sum formula
      *
      * Iterates through all elements, convolves them according to the number
      * of atoms from that element and sums up the result.
      *
      **/
    IsotopeDistribution run(const EmpiricalFormula&) const;

    void setThreshold(double threshold)
    {
      threshold_ = threshold;
    }

    double getThreshold()
    {
      return threshold_;
    }

 protected:
    double threshold_ = 0.01;

  };

} // namespace OpenMS


