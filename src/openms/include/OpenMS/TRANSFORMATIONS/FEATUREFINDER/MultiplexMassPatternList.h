// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXMASSPATTERNLIST_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXMASSPATTERNLIST_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexMassPattern.h>

#include <vector>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
  /**
   * @brief complete list of all possible mass shifts due to isotopic labelling
   * 
   * Isotopic labelling results in the shift of peptide masses. For example
   * in a Lys8/Arg10 SILAC labelled sample, some peptides (the ones with one
   * Arg in their sequence) will show a relative mass shift between light and
   * heavy partners of 10 Da. This class constructs the complete list of all
   * possible mass shifts that arise from isotopic labelling.
   */
  class OPENMS_DLLAPI MultiplexMassPatternList
  {
    public:

    /**
     * @brief constructor
     * 
     * @param labels    isotopic labels
     * @param missed_cleavages    maximum number of missed cleavages due to incomplete digestion
     * @param knock_out    Do we expect some peptides in the multiplets to be absent?
     * For example du to knock-outs in one of the samples.
     */
    MultiplexMassPatternList(String labels, int missed_cleavages, bool knock_out);
    
    /**
     * @brief returns mass shift at position i
     */
    String getLabels() const;
    
    private:
   
    /**
     * @brief isotopic labels
     */
    String labels_;

    /**
     * @brief maximum number of missed cleavages
     */
    int missed_cleavages_;

    /**
     * @brief Do we expect some peptides to be absent?
     */
    bool knock_out_;

    /**
     * @brief list of all possible mass shift patterns
     */
    std::vector<MultiplexMassPattern> mass_pattern_list_;
      
 };
  
}

#endif /* MULTIPLEXMASSPATTERNLIST_H */
