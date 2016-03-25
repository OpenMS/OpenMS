// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXDELTAMASSESGENERATOR_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXDELTAMASSESGENERATOR_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>

#include <vector>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
  /**
   * @brief generates complete list of all possible mass shifts due to isotopic labelling
   * 
   * Isotopic labelling results in the shift of peptide masses. For example
   * in a Lys8/Arg10 SILAC labelled sample, some peptides (the ones with one
   * Arg in their sequence) will show a relative mass shift between light and
   * heavy partners of 10 Da. This class constructs the complete list of all
   * possible mass shifts that arise from isotopic labelling.
   */
  class OPENMS_DLLAPI MultiplexDeltaMassesGenerator
  {
    public:

    /**
     * @brief constructor
     * 
     * @param labels    isotopic labels
     * @param missed_cleavages    maximum number of missed cleavages due to incomplete digestion
     * For example due to knock-outs in one of the samples.
     */
    MultiplexDeltaMassesGenerator(String labels, int missed_cleavages, std::map<String,double> label_mass_shift);
        
    /**
     * @brief generate all mass shifts that can occur due to the absence of one or multiple peptides
     * (e.g. for a triplet experiment generate the doublets and singlets that might be present)
     */
    void generateKnockoutDeltaMasses();

    /**
     * @brief write the list of labels for each sample
     */
    void printLabelsList() const;
    
    /**
     * @brief write the list of all mass patterns
     */
    void printMassPatternList() const;
    
    /**
     * @brief returns the list of mass shift patterns
     */
    std::vector<MultiplexDeltaMasses> getMassPatternList() const;
    
    private:
   
    /**
     * @brief isotopic labels
     */
    String labels_;
    
    /**
     * @brief list of samples with their corresponding labels
     */
    std::vector<std::vector<String> > samples_labels_;
    
    /**
     * @brief maximum number of missed cleavages
     */
    int missed_cleavages_;

    /**
     * @brief list of all possible mass shift patterns
     */
    std::vector<MultiplexDeltaMasses> mass_pattern_list_;
      
    /**
     * @brief mapping from single label to mass shift
     * e.g. "Arg10" -> 10.0082686
     */
    std::map<String,double> label_mass_shift_;
    
 };
  
}

#endif /* MULTIPLEXDELTAMASSESGENERATOR_H */
