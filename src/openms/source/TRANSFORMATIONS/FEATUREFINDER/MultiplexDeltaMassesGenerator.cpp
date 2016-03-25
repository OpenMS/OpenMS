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

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

namespace OpenMS
{

  MultiplexDeltaMassesGenerator::MultiplexDeltaMassesGenerator(String labels, int missed_cleavages, std::map<String,double> label_mass_shift) :
    labels_(labels), samples_labels_(), missed_cleavages_(missed_cleavages), label_mass_shift_(label_mass_shift)
  {
    // split the labels_ string
    String temp_labels(labels_);
    std::vector<String> temp_samples;
    
    boost::replace_all(temp_labels, "[]", "no_label");
    boost::replace_all(temp_labels, "()", "no_label");
    boost::replace_all(temp_labels, "{}", "no_label");
    boost::split(temp_samples, temp_labels, boost::is_any_of("[](){}")); // any bracket allowed to separate samples
    
    for (unsigned i = 0; i < temp_samples.size(); ++i)
    {
      if (!temp_samples[i].empty())
      {
        if (temp_samples[i]=="no_label")
        {
          vector<String> temp_labels;
          temp_labels.push_back("no_label");
          samples_labels_.push_back(temp_labels);
        }
        else
        {
          vector<String> temp_labels;
          boost::split(temp_labels, temp_samples[i], boost::is_any_of(",;: ")); // various separators allowed to separate labels
          samples_labels_.push_back(temp_labels);
        }
      }
    }
    
    if (samples_labels_.empty())
    {
      vector<String> temp_labels;
      temp_labels.push_back("no_label");
      samples_labels_.push_back(temp_labels);
    }
    
    
    // What kind of labelling do we have?
    // SILAC, Dimethyl, ICPL or no labelling ??

    bool labelling_SILAC = ((labels_.find("Arg") != std::string::npos) || (labels_.find("Lys") != std::string::npos));
    bool labelling_Dimethyl = (labels_.find("Dimethyl") != std::string::npos);
    bool labelling_ICPL = (labels_.find("ICPL") != std::string::npos);
    bool labelling_none = labels_.empty() || (labels_ == "[]") || (labels_ == "()") || (labels_ == "{}");

    bool SILAC = (labelling_SILAC && !labelling_Dimethyl && !labelling_ICPL && !labelling_none);
    bool Dimethyl = (!labelling_SILAC && labelling_Dimethyl && !labelling_ICPL && !labelling_none);
    bool ICPL = (!labelling_SILAC && !labelling_Dimethyl && labelling_ICPL && !labelling_none);
    bool none = (!labelling_SILAC && !labelling_Dimethyl && !labelling_ICPL && labelling_none);

    if (!(SILAC || Dimethyl || ICPL || none))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Unknown labelling. Neither SILAC, Dimethyl nor ICPL.");
    }
    
    // check if the labels are included in advanced section "labels"
    String all_labels = "Arg6 Arg10 Lys4 Lys6 Lys8 Dimethyl0 Dimethyl4 Dimethyl6 Dimethyl8 ICPL0 ICPL4 ICPL6 ICPL10 no_label";
    for (unsigned i = 0; i < samples_labels_.size(); i++)
    {
      for (unsigned j = 0; j < samples_labels_[i].size(); ++j)
      {
        if (all_labels.find(samples_labels_[i][j]) == std::string::npos)
        {
          std::stringstream stream;
          stream << "The label " << samples_labels_[i][j] << " is unknown.";
          throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, stream.str());
        }
      }
    }

    // generate mass pattern list
    if (SILAC)
    {
      // SILAC
      // We assume the first sample to be unlabelled. Even if the "[]" for the first sample in the label string has not been specified.
      
      for (unsigned ArgPerPeptide = 0; ArgPerPeptide <= (unsigned) missed_cleavages_ + 1; ArgPerPeptide++)
      {
        for (unsigned LysPerPeptide = 0; LysPerPeptide <= (unsigned) missed_cleavages_ + 1; LysPerPeptide++)
        {
          if (ArgPerPeptide + LysPerPeptide <= (unsigned) missed_cleavages_ + 1)
          {
            MultiplexDeltaMasses delta_masses_temp;    // single mass shift pattern
            delta_masses_temp.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0, "no_label"));
            for (unsigned i = 0; i < samples_labels_.size(); i++)
            {
              double mass_shift = 0;
              MultiplexDeltaMasses::LabelSet label_set;
              // Considering the case of an amino acid (e.g. LysPerPeptide != 0) for which no label is present (e.g. Lys4There && Lys6There && Lys8There == false) makes no sense. Therefore each amino acid will have to give its "Go Ahead" before the shift is calculated.
              bool goAhead_Lys = false;
              bool goAhead_Arg = false;

              for (unsigned j = 0; j < samples_labels_[i].size(); ++j)
              {
                bool Arg6There = (samples_labels_[i][j].find("Arg6") != std::string::npos); // Is Arg6 in the SILAC label?
                bool Arg10There = (samples_labels_[i][j].find("Arg10") != std::string::npos);
                bool Lys4There = (samples_labels_[i][j].find("Lys4") != std::string::npos);
                bool Lys6There = (samples_labels_[i][j].find("Lys6") != std::string::npos);
                bool Lys8There = (samples_labels_[i][j].find("Lys8") != std::string::npos);

                // construct label set
                for (unsigned k = 1; k < Arg6There * (ArgPerPeptide + 1); ++k)
                {
                  label_set.insert("Arg6");
                }
                for (unsigned k = 1; k < Arg10There * (ArgPerPeptide + 1); ++k)
                {
                  label_set.insert("Arg10");
                }
                for (unsigned k = 1; k < Lys4There * (LysPerPeptide + 1); ++k)
                {
                  label_set.insert("Lys4");
                }
                for (unsigned k = 1; k < Lys6There * (LysPerPeptide + 1); ++k)
                {
                  label_set.insert("Lys6");
                }
                for (unsigned k = 1; k < Lys8There * (LysPerPeptide + 1); ++k)
                {
                  label_set.insert("Lys8");
                }

                mass_shift = mass_shift + ArgPerPeptide * (Arg6There * label_mass_shift_["Arg6"] + Arg10There * label_mass_shift_["Arg10"]) + LysPerPeptide * (Lys4There * label_mass_shift_["Lys4"] + Lys6There * label_mass_shift_["Lys6"] + Lys8There * label_mass_shift_["Lys8"]);

                goAhead_Arg = goAhead_Arg || !(ArgPerPeptide != 0 && !Arg6There && !Arg10There);
                goAhead_Lys = goAhead_Lys || !(LysPerPeptide != 0 && !Lys4There && !Lys6There && !Lys8There);
              }

              if (goAhead_Arg && goAhead_Lys && (mass_shift != 0))
              {
                delta_masses_temp.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(mass_shift, label_set));
              }
            }

            if (delta_masses_temp.getDeltaMasses().size() > 1)
            {
              mass_pattern_list_.push_back(delta_masses_temp);
            }
          }
        }
      }

    }
    else if (Dimethyl || ICPL)
    {
      // Dimethyl or ICPL
      // We assume each sample to be labelled only once.

      for (unsigned mc = 0; mc <= (unsigned) missed_cleavages_; ++mc)
      {
        MultiplexDeltaMasses delta_masses_temp;    // single mass shift pattern
        for (unsigned i = 0; i < samples_labels_.size(); i++)
        {
          double mass_shift = (mc + 1) * (label_mass_shift_[samples_labels_[i][0]] - label_mass_shift_[samples_labels_[0][0]]);
          MultiplexDeltaMasses::LabelSet label_set;
          // construct label set
          for (unsigned k = 1; k < (mc + 2); ++k)
          {
            label_set.insert(samples_labels_[i][0]);
          }

         delta_masses_temp.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(mass_shift, label_set));
        }
        mass_pattern_list_.push_back(delta_masses_temp);
      }

    }
    else
    {
      // none (singlet detection)
      MultiplexDeltaMasses delta_masses_temp;
      delta_masses_temp.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0, "no_label"));
      mass_pattern_list_.push_back(delta_masses_temp);
    }
    
    // sort mass patterns
    // (from small mass shifts to larger ones, i.e. few miscleavages = simple explanation first)
    std::sort(mass_pattern_list_.begin(), mass_pattern_list_.end());

  }

  void MultiplexDeltaMassesGenerator::generateKnockoutDeltaMasses()
  {
    if (mass_pattern_list_.empty())
    {
      // Even in the case of a singlet search, there should be one mass shift (zero mass shift) in the list.
      throw OpenMS::Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, 0);
    }
    
    unsigned n = mass_pattern_list_[0].getDeltaMasses().size();    // n=1 for singlets, n=2 for doublets, n=3 for triplets, n=4 for quadruplets
    unsigned m = mass_pattern_list_.size();    // number of mass shift patterns before extension of the list
    if (n == 1)
    {
      throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Knock-outs for singlet detection not relevant.");
    }
    else if (n == 2)
    {
      // add singlets
      MultiplexDeltaMasses dm;
      dm.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0,"any_label_set"));    // There are two singlets with different label sets. But only a single singlet with "any_label_set" is added.
      mass_pattern_list_.push_back(dm);
    }
    else if (n == 3)
    {
      for (unsigned i = 0; i < m; ++i)
      {
        // add doublets
        MultiplexDeltaMasses doublet1;
        doublet1.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[0]);
        doublet1.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[1]);
        mass_pattern_list_.push_back(doublet1);

        MultiplexDeltaMasses doublet2;
        doublet2.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[1]);
        doublet2.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[2]);
        mass_pattern_list_.push_back(doublet2);

        MultiplexDeltaMasses doublet3;
        doublet3.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[0]);
        doublet3.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[2]);
        mass_pattern_list_.push_back(doublet3);
      }
      
      // add singlets
      MultiplexDeltaMasses dm;
      dm.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0, "any_label_set"));    // There are three singlets with different label sets. But only a single singlet with "any_label_set" is added.
      mass_pattern_list_.push_back(dm);
    }
    else if (n == 4)
    {
      for (unsigned i = 0; i < m; ++i)
      {
        // add triplets
        MultiplexDeltaMasses triplet1;
        triplet1.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[1]);
        triplet1.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[2]);
        triplet1.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[3]);
        mass_pattern_list_.push_back(triplet1);
        
        MultiplexDeltaMasses triplet2;
        triplet2.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[0]);
        triplet2.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[2]);
        triplet2.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[3]);
        mass_pattern_list_.push_back(triplet2);
        
        // Knockout combination previously forgotten. Will be un-commented in final FFM/MultiplexResolver version.
        /*MultiplexDeltaMasses triplet3;
        triplet3.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[0]);
        triplet3.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[1]);
        triplet3.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[3]);
        mass_pattern_list_.push_back(triplet3);*/
        
        MultiplexDeltaMasses triplet4;
        triplet4.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[0]);
        triplet4.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[1]);
        triplet4.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[2]);
        mass_pattern_list_.push_back(triplet4);


        // add doublets
        MultiplexDeltaMasses doublet1;
        doublet1.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[0]);
        doublet1.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[1]);
        mass_pattern_list_.push_back(doublet1);

        MultiplexDeltaMasses doublet2;
        doublet2.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[0]);
        doublet2.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[2]);
        mass_pattern_list_.push_back(doublet2);

        MultiplexDeltaMasses doublet3;
        doublet3.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[0]);
        doublet3.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[3]);
        mass_pattern_list_.push_back(doublet3);

        MultiplexDeltaMasses doublet4;
        doublet4.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[1]);
        doublet4.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[2]);
        mass_pattern_list_.push_back(doublet4);

        MultiplexDeltaMasses doublet5;
        doublet5.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[1]);
        doublet5.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[3]);
        mass_pattern_list_.push_back(doublet5);

        MultiplexDeltaMasses doublet6;
        doublet6.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[2]);
        doublet6.getDeltaMasses().push_back(mass_pattern_list_[i].getDeltaMasses()[3]);
        mass_pattern_list_.push_back(doublet6);
      }

      // add singlets
      MultiplexDeltaMasses dm;
      dm.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0,"any_label_set"));
      mass_pattern_list_.push_back(dm);
    }
    else if (n > 4)
    {
      throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Knock-outs for multiplex experiments with more than 4 samples not supported.");
    }
    
    // sort mass patterns
    // (from small mass shifts to larger ones, i.e. few miscleavages = simple explanation first)
    std::sort(mass_pattern_list_.begin(), mass_pattern_list_.end());

  }
  
  void MultiplexDeltaMassesGenerator::printLabelsList() const
  {
    cout << "\n";
    for (unsigned i = 0; i < samples_labels_.size(); ++i)
    {
      cout << "sample " << (i + 1) << ":    ";
      for (unsigned j = 0; j < samples_labels_[i].size(); ++j)
      {
        cout << samples_labels_[i][j] << "    ";
      }
      cout << "\n";
    }
  }
  
  void MultiplexDeltaMassesGenerator::printMassPatternList() const
  {
    cout << "\n";
    for (unsigned i = 0; i < mass_pattern_list_.size(); ++i)
    {
      std::cout << "mass shift " << (i + 1) << ":    ";
      for (unsigned j = 0; j < mass_pattern_list_[i].getDeltaMasses().size(); ++j)
      {
        double mass_shift = mass_pattern_list_[i].getDeltaMasses()[j].delta_mass;
        MultiplexDeltaMasses::LabelSet label_set = mass_pattern_list_[i].getDeltaMasses()[j].label_set;
        
        std::cout << mass_shift << " (";
        for (std::multiset<String>::iterator it = label_set.begin(); it != label_set.end(); ++it)
        {
          if (it != label_set.begin())
          {
            std::cout << ",";
          }
          std::cout << *it;
        }
        std::cout << ")    ";
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }
  
  std::vector<MultiplexDeltaMasses> MultiplexDeltaMassesGenerator::getMassPatternList() const
  {
    return mass_pattern_list_;
  }

}
