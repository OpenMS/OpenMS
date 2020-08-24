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

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <iostream>
#include <ostream>

using namespace std;

namespace OpenMS
{
  MultiplexDeltaMassesGenerator::Label::Label(String sn, String ln, String d, double dm) :
    short_name(sn),
    long_name(ln),
    description(d),
    delta_mass(dm)
  {
  }

  MultiplexDeltaMassesGenerator::MultiplexDeltaMassesGenerator() :
    DefaultParamHandler("labels"),
    labels_(),
    labels_list_(),
    samples_labels_(),
    missed_cleavages_(),
    label_delta_mass_()
  {
    // fill label master list
    fillLabelMasterList_();

    // set user parameters
    for (std::vector<MultiplexDeltaMassesGenerator::Label>::const_iterator it = label_master_list_.begin(); it != label_master_list_.end(); ++it)
    {
      defaults_.setValue(it->short_name, it->delta_mass, it->description);
      defaults_.setMinFloat(it->short_name, 0);    
    }
    defaultsToParam_();    
  }
  
  MultiplexDeltaMassesGenerator::MultiplexDeltaMassesGenerator(String labels, int missed_cleavages, std::map<String,double> label_delta_mass) :
    DefaultParamHandler("labels"),
    labels_(labels),
    labels_list_(),
    samples_labels_(),
    missed_cleavages_(missed_cleavages),
    label_delta_mass_(label_delta_mass)
  {
    // fill label master list
    fillLabelMasterList_();
    
    // generate short/long label mappings
    for (std::vector<MultiplexDeltaMassesGenerator::Label>::const_iterator it = label_master_list_.begin(); it != label_master_list_.end(); ++it)
    {
      label_short_long_.insert(make_pair(it->short_name, it->long_name));
      label_long_short_.insert(make_pair(it->long_name, it->short_name));
    }
    
    // split the labels_ string
    String temp_labels_string(labels_);
    std::vector<String> temp_samples;
    
    boost::replace_all(temp_labels_string, "[]", "no_label");
    boost::replace_all(temp_labels_string, "()", "no_label");
    boost::replace_all(temp_labels_string, "{}", "no_label");
    boost::split(temp_samples, temp_labels_string, boost::is_any_of("[](){}")); // any bracket allowed to separate samples
    
    for (String::size_type i = 0; i < temp_samples.size(); ++i)
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
    // SILAC, Leu, Dimethyl, ICPL, numeric labelling or no labelling ??

    bool labelling_SILAC = ((labels_.find("Arg") != std::string::npos) || (labels_.find("Lys") != std::string::npos));
    bool labelling_Leu = (labels_.find("Leu") != std::string::npos);
    bool labelling_Dimethyl = (labels_.find("Dimethyl") != std::string::npos);
    bool labelling_ICPL = (labels_.find("ICPL") != std::string::npos);
    // Check whether each label string represents a double. If yes, use these doubles as mass shifts.
    bool labelling_numeric = true;
    for (size_t i = 0; i < samples_labels_.size(); i++)
    {
      for (size_t j = 0; j < samples_labels_[i].size(); j++)
      {
        try
        {
          double mass_shift = std::stod(samples_labels_[i][j]);
          
          // For numeric mass shifts, long and short label names as well as the numerical mass shift are trivial.
          // For example, long label name ("3.1415"), short label name ("3.1415") and numerical mass shift (3.1415).
          label_delta_mass_.insert(make_pair(samples_labels_[i][j], mass_shift));
          label_short_long_.insert(make_pair(samples_labels_[i][j], samples_labels_[i][j]));
          label_long_short_.insert(make_pair(samples_labels_[i][j], samples_labels_[i][j]));
        }
        catch(...)
        {
          labelling_numeric = false;
        }
      }
    }
    bool labelling_none = labels_.empty() || (labels_ == "[]") || (labels_ == "()") || (labels_ == "{}");

    bool SILAC = (labelling_SILAC && !labelling_Leu && !labelling_Dimethyl && !labelling_ICPL && !labelling_numeric && !labelling_none);
    bool Leu = (!labelling_SILAC && labelling_Leu && !labelling_Dimethyl && !labelling_ICPL && !labelling_numeric && !labelling_none);
    bool Dimethyl = (!labelling_SILAC && !labelling_Leu && labelling_Dimethyl && !labelling_ICPL && !labelling_numeric && !labelling_none);
    bool ICPL = (!labelling_SILAC && !labelling_Leu && !labelling_Dimethyl && labelling_ICPL && !labelling_numeric && !labelling_none);
    bool numeric = (!labelling_SILAC && !labelling_Leu && !labelling_Dimethyl && !labelling_ICPL && labelling_numeric && !labelling_none);
    bool none = (!labelling_SILAC && !labelling_Leu && !labelling_Dimethyl && !labelling_ICPL && !labelling_numeric && labelling_none);

    if (!(SILAC || Leu || Dimethyl || ICPL || numeric || none))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown labelling. Neither SILAC, Leu, Dimethyl nor ICPL.");
    }
    
    // Check if the labels are included in advanced section "labels"
    // unless the labelling is numeric.
    if (!numeric)
    {
      String all_labels = "Arg6 Arg10 Lys4 Lys6 Lys8 Leu3 Dimethyl0 Dimethyl4 Dimethyl6 Dimethyl8 ICPL0 ICPL4 ICPL6 ICPL10 no_label";
      for (std::vector<std::vector<String> >::size_type i = 0; i < samples_labels_.size(); i++)
      {
        for (std::vector<String>::size_type j = 0; j < samples_labels_[i].size(); ++j)
        {
          if (all_labels.find(samples_labels_[i][j]) == std::string::npos)
          {
            std::stringstream stream;
            stream << "The label " << samples_labels_[i][j] << " is unknown.";
            throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, stream.str());
          }
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

                mass_shift = mass_shift + ArgPerPeptide * (Arg6There * label_delta_mass_["Arg6"] + Arg10There * label_delta_mass_["Arg10"]) + LysPerPeptide * (Lys4There * label_delta_mass_["Lys4"] + Lys6There * label_delta_mass_["Lys6"] + Lys8There * label_delta_mass_["Lys8"]);

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
              delta_masses_list_.push_back(delta_masses_temp);
            }
          }
        }
      }

    }
    else if (Leu)
    {
      // Leu
      // We assume each sample to be labelled only once. Hence, we only consider samples_labels_[...][0] below.
      // Unlike in classical SILAC where two labels with two different specificities are available, in Leu labelling
      // there is only one specificity Leu. Hence, [Lys8,Arg10] but [Leu3].
      
      for (unsigned mc = 0; mc <= (unsigned) missed_cleavages_; ++mc)
      {
        MultiplexDeltaMasses delta_masses_temp;    // single mass shift pattern
        
        delta_masses_temp.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0, "no_label"));
        
        double mass_shift = (mc + 1) * (label_delta_mass_[samples_labels_[1][0]] - label_delta_mass_[samples_labels_[0][0]]);
        MultiplexDeltaMasses::LabelSet label_set;
        // construct label set
        for (unsigned k = 1; k < (mc + 2); ++k)
        {
          label_set.insert(samples_labels_[1][0]);
        }
        delta_masses_temp.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(mass_shift, label_set));
        
        delta_masses_list_.push_back(delta_masses_temp);
      }
      
    }
    else if (Dimethyl || ICPL)
    {
      // Dimethyl or ICPL
      // We assume each sample to be labelled only once. Hence, we only consider samples_labels_[...][0] below.
      // Unlike in classical SILAC where two labels with two different specificities are available, in Dimethyl labelling
      // there is only one specificity (Lys and N-term). Two labels [Lys8,Arg10] in a SILAC sample are fine, since
      // both have different specificities. But two labels [Dimethyl4,Dimethyl8] make no sense, since both have the same
      // specificity. With only one specificity in Dimethyl labelling available, each sample can have only one label.

      for (unsigned mc = 0; mc <= (unsigned) missed_cleavages_; ++mc)
      {
        MultiplexDeltaMasses delta_masses_temp;    // single mass shift pattern
        for (unsigned i = 0; i < samples_labels_.size(); i++)
        {
          double mass_shift = (mc + 1) * (label_delta_mass_[samples_labels_[i][0]] - label_delta_mass_[samples_labels_[0][0]]);
          MultiplexDeltaMasses::LabelSet label_set;
          // construct label set
          for (unsigned k = 1; k < (mc + 2); ++k)
          {
            label_set.insert(samples_labels_[i][0]);
          }

         delta_masses_temp.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(mass_shift, label_set));
        }
        delta_masses_list_.push_back(delta_masses_temp);
      }

    }
    else if (numeric)
    {
      for (unsigned mc = 0; mc <= (unsigned) missed_cleavages_; ++mc)
      {
        MultiplexDeltaMasses delta_masses_temp;    // single mass shift pattern
        for (unsigned i = 0; i < samples_labels_.size(); i++)
        {
          double mass_shift = (mc + 1) * (label_delta_mass_[samples_labels_[i][0]] - label_delta_mass_[samples_labels_[0][0]]);
          MultiplexDeltaMasses::LabelSet label_set;
          for (unsigned k = 1; k < (mc + 2); ++k)
          {
            label_set.insert(samples_labels_[i][0]);
          }
          
          delta_masses_temp.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(mass_shift, label_set));
        }
        delta_masses_list_.push_back(delta_masses_temp);
      }    
    }
    else
    {
      // none (singlet detection)
      MultiplexDeltaMasses delta_masses_temp;
      delta_masses_temp.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0, "no_label"));
      delta_masses_list_.push_back(delta_masses_temp);
    }
    
    // sort mass patterns
    // (from small mass shifts to larger ones, i.e. few miscleavages = simple explanation first)
    std::sort(delta_masses_list_.begin(), delta_masses_list_.end());

    // generate flat list of all occurring isotopic labels
    for (unsigned i = 0; i < samples_labels_.size(); ++i)
    {
      for (unsigned j = 0; j < samples_labels_[i].size(); ++j)
      {
        if (samples_labels_[i][j] != "no_label")
        {
          labels_list_.push_back(samples_labels_[i][j]);
        }
      }
    }

  }

  void MultiplexDeltaMassesGenerator::generateKnockoutDeltaMasses()
  {
    if (delta_masses_list_.empty())
    {
      // Even in the case of a singlet search, there should be one mass shift (zero mass shift) in the list.
      throw OpenMS::Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 0);
    }
    
    unsigned n = delta_masses_list_[0].getDeltaMasses().size();    // n=1 for singlets, n=2 for doublets, n=3 for triplets, n=4 for quadruplets
    unsigned m = delta_masses_list_.size();    // number of mass shift patterns before extension of the list
    if (n == 1)
    {
      throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Knock-outs for singlet detection not relevant.");
    }
    else if (n == 2)
    {
      // add singlets
      MultiplexDeltaMasses dm;
      dm.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0,"any_label_set"));    // There are two singlets with different label sets. But only a single singlet with "any_label_set" is added.
      delta_masses_list_.push_back(dm);
    }
    else if (n == 3)
    {
      for (unsigned i = 0; i < m; ++i)
      {
        // add doublets
        MultiplexDeltaMasses doublet1;
        doublet1.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[0]);
        doublet1.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[1]);
        delta_masses_list_.push_back(doublet1);

        MultiplexDeltaMasses doublet2;
        doublet2.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[1]);
        doublet2.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[2]);
        delta_masses_list_.push_back(doublet2);

        MultiplexDeltaMasses doublet3;
        doublet3.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[0]);
        doublet3.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[2]);
        delta_masses_list_.push_back(doublet3);
      }
      
      // add singlets
      MultiplexDeltaMasses dm;
      dm.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0, "any_label_set"));    // There are three singlets with different label sets. But only a single singlet with "any_label_set" is added.
      delta_masses_list_.push_back(dm);
    }
    else if (n == 4)
    {
      for (unsigned i = 0; i < m; ++i)
      {
        // add triplets
        MultiplexDeltaMasses triplet1;
        triplet1.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[1]);
        triplet1.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[2]);
        triplet1.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[3]);
        delta_masses_list_.push_back(triplet1);
        
        MultiplexDeltaMasses triplet2;
        triplet2.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[0]);
        triplet2.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[2]);
        triplet2.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[3]);
        delta_masses_list_.push_back(triplet2);
        
        // Knockout combination previously forgotten. Will be un-commented in final FFM/MultiplexResolver version.
        /*MultiplexDeltaMasses triplet3;
        triplet3.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[0]);
        triplet3.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[1]);
        triplet3.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[3]);
        delta_masses_list_.push_back(triplet3);*/
        
        MultiplexDeltaMasses triplet4;
        triplet4.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[0]);
        triplet4.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[1]);
        triplet4.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[2]);
        delta_masses_list_.push_back(triplet4);


        // add doublets
        MultiplexDeltaMasses doublet1;
        doublet1.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[0]);
        doublet1.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[1]);
        delta_masses_list_.push_back(doublet1);

        MultiplexDeltaMasses doublet2;
        doublet2.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[0]);
        doublet2.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[2]);
        delta_masses_list_.push_back(doublet2);

        MultiplexDeltaMasses doublet3;
        doublet3.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[0]);
        doublet3.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[3]);
        delta_masses_list_.push_back(doublet3);

        MultiplexDeltaMasses doublet4;
        doublet4.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[1]);
        doublet4.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[2]);
        delta_masses_list_.push_back(doublet4);

        MultiplexDeltaMasses doublet5;
        doublet5.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[1]);
        doublet5.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[3]);
        delta_masses_list_.push_back(doublet5);

        MultiplexDeltaMasses doublet6;
        doublet6.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[2]);
        doublet6.getDeltaMasses().push_back(delta_masses_list_[i].getDeltaMasses()[3]);
        delta_masses_list_.push_back(doublet6);
      }

      // add singlets
      MultiplexDeltaMasses dm;
      dm.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0,"any_label_set"));
      delta_masses_list_.push_back(dm);
    }
    else if (n > 4)
    {
      throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Knock-outs for multiplex experiments with more than 4 samples not supported.");
    }
    
    // sort mass patterns
    // (from small mass shifts to larger ones, i.e. few miscleavages = simple explanation first)
    std::sort(delta_masses_list_.begin(), delta_masses_list_.end());
    
  }
  
  void MultiplexDeltaMassesGenerator::printSamplesLabelsList(std::ostream &stream) const
  {
    stream << "\n";
    for (unsigned i = 0; i < samples_labels_.size(); ++i)
    {
      stream << "sample " << (i + 1) << ":    ";
      for (unsigned j = 0; j < samples_labels_[i].size(); ++j)
      {
        stream << samples_labels_[i][j] << "    ";
      }
      stream << "\n";
    }
  }
  
  void MultiplexDeltaMassesGenerator::printDeltaMassesList(std::ostream &stream) const
  {
    stream << "\n";
    for (unsigned i = 0; i < delta_masses_list_.size(); ++i)
    {
      stream << "mass shift " << (i + 1) << ":    ";
      for (unsigned j = 0; j < delta_masses_list_[i].getDeltaMasses().size(); ++j)
      {
        double mass_shift = delta_masses_list_[i].getDeltaMasses()[j].delta_mass;
        MultiplexDeltaMasses::LabelSet label_set = delta_masses_list_[i].getDeltaMasses()[j].label_set;
        
        stream << mass_shift << " (";
        for (std::multiset<String>::iterator it = label_set.begin(); it != label_set.end(); ++it)
        {
          if (it != label_set.begin())
          {
            stream << ",";
          }
          stream << *it;
        }
        stream << ")    ";
      }
      stream << "\n";
    }
    stream << "\n";
  }
  
  std::vector<MultiplexDeltaMasses> MultiplexDeltaMassesGenerator::getDeltaMassesList()
  {
    return delta_masses_list_;
  }

  const std::vector<MultiplexDeltaMasses>& MultiplexDeltaMassesGenerator::getDeltaMassesList() const
  {
    return delta_masses_list_;
  }
  
  std::vector<std::vector<String> > MultiplexDeltaMassesGenerator::getSamplesLabelsList()
  {
    return samples_labels_;
  }

  const std::vector<std::vector<String> >& MultiplexDeltaMassesGenerator::getSamplesLabelsList() const
  {
    return samples_labels_;
  }

  String MultiplexDeltaMassesGenerator::getLabelShort(String label)
  {
    return label_long_short_[label];
  }

  String MultiplexDeltaMassesGenerator::getLabelLong(String label)
  {
    return label_short_long_[label];
  }
  
  MultiplexDeltaMasses::LabelSet MultiplexDeltaMassesGenerator::extractLabelSet(AASequence sequence)
  {
    String s(sequence.toString());

    MultiplexDeltaMasses::LabelSet label_set;
    // loop over all labels that might occur
    for (std::vector<String>::size_type i = 0; i < labels_list_.size(); ++i)
    {    
      String label("(" + getLabelLong(labels_list_[i]) + ")");
      String::size_type length_label = label.size();
      
      // check if label occurs in peptide sequence
      if (s.hasSubstring(label))
      {
        String::size_type length_before = s.size();
        s.substitute(label, "");
        String::size_type length_after = s.size();
        String::size_type multiple = (length_before - length_after)/length_label;
        
        // add as many labels to the set as occured in the peptide sequence
        for (String::size_type j = 0; j < multiple; ++j)
        {
          label_set.insert(labels_list_[i]);
        }
      }
    }

    // add no_label if nothing was found
    // (either way if no_label is in the theoretical list or not)
    if (label_set.empty())
    {
      label_set.insert("no_label");
    }    
    
    return label_set;
  }

  void MultiplexDeltaMassesGenerator::fillLabelMasterList_()
  {
    label_master_list_.push_back(MultiplexDeltaMassesGenerator::Label("Arg6", "Label:13C(6)", "Label:13C(6)  |  C(-6) 13C(6)  |  unimod #188", 6.0201290268));
    label_master_list_.push_back(MultiplexDeltaMassesGenerator::Label("Arg10", "Label:13C(6)15N(4)", "Label:13C(6)15N(4)  |  C(-6) 13C(6) N(-4) 15N(4)  |  unimod #267", 10.008268600));
    label_master_list_.push_back(MultiplexDeltaMassesGenerator::Label("Lys4", "Label:2H(4)", "Label:2H(4)  |  H(-4) 2H(4)  |  unimod #481", 4.0251069836));
    label_master_list_.push_back(MultiplexDeltaMassesGenerator::Label("Lys6", "Label:13C(6)", "Label:13C(6)  |  C(-6) 13C(6)  |  unimod #188", 6.0201290268));
    label_master_list_.push_back(MultiplexDeltaMassesGenerator::Label("Lys8", "Label:13C(6)15N(2)", "Label:13C(6)15N(2)  |  C(-6) 13C(6) N(-2) 15N(2)  |  unimod #259", 8.0141988132));
    label_master_list_.push_back(MultiplexDeltaMassesGenerator::Label("Leu3", "Label:2H(3)", "Label:2H(3)  |  H(-3) 2H(3)  |  unimod #262", 3.018830));
    label_master_list_.push_back(MultiplexDeltaMassesGenerator::Label("Dimethyl0", "Dimethyl", "Dimethyl  |  H(4) C(2)  |  unimod #36", 28.031300));
    label_master_list_.push_back(MultiplexDeltaMassesGenerator::Label("Dimethyl4", "Dimethyl:2H(4)", "Dimethyl:2H(4)  |  2H(4) C(2)  |  unimod #199", 32.056407));
    label_master_list_.push_back(MultiplexDeltaMassesGenerator::Label("Dimethyl6", "Dimethyl:2H(4)13C(2)", "Dimethyl:2H(4)13C(2)  |  2H(4) 13C(2)  |  unimod #510", 34.063117));
    label_master_list_.push_back(MultiplexDeltaMassesGenerator::Label("Dimethyl8", "Dimethyl:2H(6)13C(2)", "Dimethyl:2H(6)13C(2)  |  H(-2) 2H(6) 13C(2)  |  unimod #330", 36.075670));
    label_master_list_.push_back(MultiplexDeltaMassesGenerator::Label("ICPL0", "ICPL", "ICPL  |  H(3) C(6) N O  |  unimod #365", 105.021464));
    label_master_list_.push_back(MultiplexDeltaMassesGenerator::Label("ICPL4", "ICPL:2H(4)", "ICPL:2H(4)  |  H(-1) 2H(4) C(6) N O  |  unimod #687", 109.046571));
    label_master_list_.push_back(MultiplexDeltaMassesGenerator::Label("ICPL6", "ICPL:13C(6)", "ICPL:13C(6)  |  H(3) 13C(6) N O  |  unimod #364", 111.041593));
    label_master_list_.push_back(MultiplexDeltaMassesGenerator::Label("ICPL10", "ICPL:13C(6)2H(4)", "ICPL:13C(6)2H(4)  |  H(-1) 2H(4) 13C(6) N O  |  unimod #866", 115.066700));
  }
}
