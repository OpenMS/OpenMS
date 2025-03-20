// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
#include <OpenMS/FEATUREFINDER/MultiplexDeltaMasses.h>

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
    String getLabelShort(const String& label);
    
    /**
     * @brief returns the long label string
     * 
     * @param label    short label, as it appears in the "labels" parameter, e.g. "Arg10"
     */
    String getLabelLong(const String& label);
    
    /**
     * @brief extract the label set from the sequence
     *
     * @param sequence    amino acid sequence
     * 
     * For example, the sequence VLSEEEIDDNFK(Label:13C(6)15N(2))AQR(Label:13C(6)15N(4))
     * contains a set of two labels, Lys8 and Arg10.
     */
    MultiplexDeltaMasses::LabelSet extractLabelSet(const AASequence& sequence);
    
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

