// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>

#include <vector>
#include <algorithm>
#include <iosfwd>


namespace OpenMS
{
  /**
   * @brief generates complete list of all possible mass shifts due to isotopic labelling
   * 
   * Isotopic labelling results in the shift of peptide masses.
   * 
   * For example in a Lys8/Arg10 SILAC labelled sample, some peptides (the ones with one
   * Arg in their sequence) will show a relative mass shift between light and heavy
   * partners of 10 Da. This class constructs the complete list of all possible mass
   * shifts that arise from isotopic labelling.
   */
  class OPENMS_DLLAPI MultiplexDeltaMassesGenerator :
    public DefaultParamHandler
  {
    public:

    /**
     * @brief complete label information
     */
    struct OPENMS_DLLAPI Label
    {
      String short_name;
      String long_name;
      String description;
      double delta_mass;
      
      Label(String sn, String ln, String d, double dm);
    };
   
    /**
     * @brief constructor
     */
    MultiplexDeltaMassesGenerator();
    
    /**
     * @brief constructor
     * 
     * @param labels    string describing the labels used in each sample. [...] specifies the labels for a single sample. For example
     * For example, [][Lys8,Arg10] describes a standard SILAC experiment. In the "light" sample, none of the amino acids are labelled [].
     * In the "heavy" sample, lysines and arginines are isotopically labelled [Lys8,Arg10].
     * @param missed_cleavages    maximum number of missed cleavages due to incomplete digestion
     * @param label_mass_shift    name of labels (e.g. Lys8) and their corresponding mass shifts (e.g. 8.0141988132)
     */
    MultiplexDeltaMassesGenerator(String labels, int missed_cleavages, std::map<String,double> label_mass_shift);
     
    /**
     * @brief generate all mass shifts that can occur due to the absence of one or multiple peptides
     * (e.g. for a triplet experiment generate the doublets and singlets that might be present)
     */
    void generateKnockoutDeltaMasses();

    /**
     * @brief write the list of labels for each of the sample
     * 
     * For example in a standard SILAC experiment, sample 1 (light) is unlabelled and sample 2 (heavy) contains Lys8 and Arg 10 labels.
     * sample 1:    no_label    
     * sample 2:    Lys8    Arg10
     */
    void printSamplesLabelsList(std::ostream &stream) const;
    
    /**
     * @brief write the list of all mass patterns
     * 
     * For example in a standard SILAC experiment allowing for one missed cleavage, five mass shift patterns are possible.
     * mass shift 1:    0 (no_label)    8.0142 (Lys8)    
     * mass shift 2:    0 (no_label)    10.0083 (Arg10)    
     * mass shift 3:    0 (no_label)    16.0284 (Lys8,Lys8)    
     * mass shift 4:    0 (no_label)    18.0225 (Arg10,Lys8)    
     * mass shift 5:    0 (no_label)    20.0165 (Arg10,Arg10)   
     * 
     * @param stream    output stream 
     */
    void printDeltaMassesList(std::ostream &stream) const;
    
    /**
     * @brief returns the list of mass shift patterns
     * 
     * @param stream    output stream 
     */
    std::vector<MultiplexDeltaMasses> getDeltaMassesList();
    
    /**
     * @brief returns the list of mass shift patterns
     */
    const std::vector<MultiplexDeltaMasses>& getDeltaMassesList() const;
    
    /**
     * @brief returns the list of samples with their corresponding labels
     * 
     * For example in a standard SILAC experiment:
     * sample 1:    no_label    
     * sample 2:    Lys8    Arg10
     */
    std::vector<std::vector<String> > getSamplesLabelsList();
    
    /**
     * @brief returns the list of samples with their corresponding labels
     * 
     * For example in a standard SILAC experiment:
     * sample 1:    no_label    
     * sample 2:    Lys8    Arg10
    */
    const std::vector<std::vector<String> >& getSamplesLabelsList() const;
    
    /**
     * @brief returns the short label string
     * 
     * @param label    long label, UniMod name as it appears in peptide sequences, e.g. "Label:13C(6)15N(4)"
     */
    String getLabelShort(String label);
    
    /**
     * @brief returns the long label string
     * 
     * @param label    short label, as it appears in the "labels" parameter, e.g. "Arg10"
     */
    String getLabelLong(String label);
    
    /**
     * @brief extract the label set from the sequence
     *
     * @param sequence    amino acid sequence
     * 
     * For example, the sequence VLSEEEIDDNFK(Label:13C(6)15N(2))AQR(Label:13C(6)15N(4))
     * contains a set of two labels, Lys8 and Arg10.
     */
    MultiplexDeltaMasses::LabelSet extractLabelSet(AASequence sequence);
    
    private:
   
    /**
     * @brief isotopic labels
     */
    String labels_;
    
    /**
     * @brief flat list of all occurring isotopic labels
     */
    std::vector<String> labels_list_;
    
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
    std::vector<MultiplexDeltaMasses> delta_masses_list_;
      
    /**
     * @brief master list of all labels
     */
    std::vector<Label> label_master_list_;
    
    /**
     * @brief mapping from single label to delta mass
     * e.g. "Arg10" -> 10.0082686
     */
    std::map<String, double> label_delta_mass_;
    
    /**
     * @brief mapping from a short label (as in the user params) to a long label (as in PSI-MS name)
     * e.g. "Arg10" -> "Label:13C(6)15N(4)"
     */
    std::map<String, String> label_short_long_;
    
    /**
     * @brief mapping from a long label (as in PSI-MS name) to a short label (as in the user params)
     * e.g. "Label:13C(6)15N(4)" -> "Arg10"
     */
    std::map<String, String> label_long_short_;
    
    /**
     * @brief fill label master list
     */
    void fillLabelMasterList_();
 };
  
}

