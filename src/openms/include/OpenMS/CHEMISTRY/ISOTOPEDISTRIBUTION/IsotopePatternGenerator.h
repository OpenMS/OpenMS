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
// $Maintainer: Timo Sachsenberg $
// $Authors: Nikos Patikos $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_ISOTOPEPATTERNGENERATOR_H
#define OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_ISOTOPEPATTERNGENERATOR_H

#include <OpenMS/config.h>

namespace OpenMS
{
  class EmpiricalFormula;
  class IsotopeDistribution;

  /** 
      @brief Provides an interface for different isotope pattern generator methods.
      
      The IsotopePatternGenerator interface  allows the developer integrate various 
      isotope pattern generator methods in the OpenMS code. It provides a run() method 
      that generates but  does not hold any generated isotope distribution data in 
      the class. Instead it returns an IsotopeDistribution to the caller.

   */
  class OPENMS_DLLAPI IsotopePatternGenerator
  {
 public:
    IsotopePatternGenerator();
    IsotopePatternGenerator(double probability_cutoff);
    
    /** 
        @brief interface that is being used by the Isotope Pattern Generator methods.
        
        Method that calculates the isotope distribution for the given formula.

     */
    virtual IsotopeDistribution run(const EmpiricalFormula&) const = 0;
    virtual ~IsotopePatternGenerator();

 protected:
    double min_prob_;
  };
}

#endif
