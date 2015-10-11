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

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexMassPatternList.h>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

namespace OpenMS
{

  MultiplexMassPatternList::MultiplexMassPatternList(String labels, int missed_cleavages, std::map<String,double> label_mass_shift) :
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

    // generate mass shift list
    std::vector<std::vector<double> > list;
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
            std::vector<double> temp;
            temp.push_back(0);
            for (unsigned i = 0; i < samples_labels_.size(); i++)
            {
              double mass_shift = 0;
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

                mass_shift = mass_shift + ArgPerPeptide * (Arg6There * label_mass_shift_["Arg6"] + Arg10There * label_mass_shift_["Arg10"]) + LysPerPeptide * (Lys4There * label_mass_shift_["Lys4"] + Lys6There * label_mass_shift_["Lys6"] + Lys8There * label_mass_shift_["Lys8"]);

                goAhead_Arg = goAhead_Arg || !(ArgPerPeptide != 0 && !Arg6There && !Arg10There);
                goAhead_Lys = goAhead_Lys || !(LysPerPeptide != 0 && !Lys4There && !Lys6There && !Lys8There);
              }

              if (goAhead_Arg && goAhead_Lys && (mass_shift != 0))
              {
                temp.push_back(mass_shift);
              }
            }

            if (temp.size() > 1)
            {
              list.push_back(temp);
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
        std::vector<double> temp;
        for (unsigned i = 0; i < samples_labels_.size(); i++)
        {
          temp.push_back((mc + 1) * (label_mass_shift_[samples_labels_[i][0]] - label_mass_shift_[samples_labels_[0][0]]));
        }
        list.push_back(temp);
      }

    }
    else
    {
      // none (singlet detection)
      std::vector<double> temp;
      temp.push_back(0);
      list.push_back(temp);
    }

    // sort mass patterns
    // (from small mass shifts to larger ones, i.e. few miscleavages = simple explanation first)
    std::sort(list.begin(), list.end());

    for (unsigned i = 0; i < list.size(); ++i)
    {
      mass_pattern_list_.push_back(MultiplexDeltaMasses(list[i]));
    }

  }

  void MultiplexMassPatternList::generateKnockoutMassShifts()
  {
    if (mass_pattern_list_.empty())
    {
      // Even in the case of a singlet search, there should be one mass shift (zero mass shift) in the list.
      throw OpenMS::Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, 0);
    }
    
    unsigned n = mass_pattern_list_[0].getMassShiftCount();    // n=1 for singlets, n=2 for doublets, n=3 for triplets, n=4 for quadruplets
    unsigned m = mass_pattern_list_.size();    // number of mass shift patterns before extension of the list
    if (n == 1)
    {
      throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Knock-outs for singlet detection not relevant.");
    }
    else if (n == 2)
    {
      // add singlets
      std::vector<double> singlet(1, 0);
      mass_pattern_list_.push_back(MultiplexDeltaMasses(singlet));
    }
    else if (n == 3)
    {
      for (unsigned i = 0; i < m; ++i)
      {
        // add doublets
        std::vector<double> doublet1(1, 0);
        doublet1.push_back(mass_pattern_list_[i].getMassShiftAt(1) - mass_pattern_list_[i].getMassShiftAt(0));
        mass_pattern_list_.push_back(MultiplexDeltaMasses(doublet1));

        std::vector<double> doublet2(1, 0);
        doublet2.push_back(mass_pattern_list_[i].getMassShiftAt(2) - mass_pattern_list_[i].getMassShiftAt(1));
        mass_pattern_list_.push_back(MultiplexDeltaMasses(doublet2));

        std::vector<double> doublet3(1, 0);
        doublet3.push_back(mass_pattern_list_[i].getMassShiftAt(2) - mass_pattern_list_[i].getMassShiftAt(0));
        mass_pattern_list_.push_back(MultiplexDeltaMasses(doublet3));
      }
      
      // add singlets
      std::vector<double> singlet(1, 0);
      mass_pattern_list_.push_back(MultiplexDeltaMasses(singlet));
    }
    else if (n == 4)
    {
      for (unsigned i = 0; i < m; ++i)
      {
        // add triplets
        std::vector<double> triplet1(1, 0);
        triplet1.push_back(mass_pattern_list_[i].getMassShiftAt(2) - mass_pattern_list_[i].getMassShiftAt(1));
        triplet1.push_back(mass_pattern_list_[i].getMassShiftAt(3) - mass_pattern_list_[i].getMassShiftAt(1));
        mass_pattern_list_.push_back(MultiplexDeltaMasses(triplet1));

        std::vector<double> triplet2(1, 0);
        triplet2.push_back(mass_pattern_list_[i].getMassShiftAt(2) - mass_pattern_list_[i].getMassShiftAt(0));
        triplet2.push_back(mass_pattern_list_[i].getMassShiftAt(3) - mass_pattern_list_[i].getMassShiftAt(0));
        mass_pattern_list_.push_back(MultiplexDeltaMasses(triplet2));

        std::vector<double> triplet3(1, 0);
        triplet3.push_back(mass_pattern_list_[i].getMassShiftAt(1) - mass_pattern_list_[i].getMassShiftAt(0));
        triplet3.push_back(mass_pattern_list_[i].getMassShiftAt(2) - mass_pattern_list_[i].getMassShiftAt(0));
        mass_pattern_list_.push_back(MultiplexDeltaMasses(triplet3));


        // add doublets
        std::vector<double> doublet1(1, 0);
        doublet1.push_back(mass_pattern_list_[i].getMassShiftAt(1) - mass_pattern_list_[i].getMassShiftAt(0));
        mass_pattern_list_.push_back(MultiplexDeltaMasses(doublet1));

        std::vector<double> doublet2(1, 0);
        doublet2.push_back(mass_pattern_list_[i].getMassShiftAt(2) - mass_pattern_list_[i].getMassShiftAt(0));
        mass_pattern_list_.push_back(MultiplexDeltaMasses(doublet2));

        std::vector<double> doublet3(1, 0);
        doublet3.push_back(mass_pattern_list_[i].getMassShiftAt(3) - mass_pattern_list_[i].getMassShiftAt(0));
        mass_pattern_list_.push_back(MultiplexDeltaMasses(doublet3));

        std::vector<double> doublet4(1, 0);
        doublet4.push_back(mass_pattern_list_[i].getMassShiftAt(2) - mass_pattern_list_[i].getMassShiftAt(1));
        mass_pattern_list_.push_back(MultiplexDeltaMasses(doublet4));

        std::vector<double> doublet5(1, 0);
        doublet5.push_back(mass_pattern_list_[i].getMassShiftAt(3) - mass_pattern_list_[i].getMassShiftAt(1));
        mass_pattern_list_.push_back(MultiplexDeltaMasses(doublet5));

        std::vector<double> doublet6(1, 0);
        doublet6.push_back(mass_pattern_list_[i].getMassShiftAt(3) - mass_pattern_list_[i].getMassShiftAt(2));
        mass_pattern_list_.push_back(MultiplexDeltaMasses(doublet6));
      }

      // add singlets
      std::vector<double> singlet(1, 0);
      mass_pattern_list_.push_back(MultiplexDeltaMasses(singlet));
    }
    else if (n > 4)
    {
      throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Knock-outs for multiplex experiments with more than 4 samples not supported.");
    }
  }
  
  void MultiplexMassPatternList::printLabelsList() const
  {
    cout << "\n";
    for (unsigned i = 0; i < samples_labels_.size(); ++i)
    {
      cout << "sample " << (i + 1) << ":   ";
      for (unsigned j = 0; j < samples_labels_[i].size(); ++j)
      {
        cout << samples_labels_[i][j] << " ";
      }
      cout << "\n";
    }
  }
  
  void MultiplexMassPatternList::printMassPatternList() const
  {
    cout << "\n";
    for (unsigned i = 0; i < mass_pattern_list_.size(); ++i)
    {
      std::cout << "mass shift " << (i + 1) << ":    ";
      for (unsigned j = 0; j < mass_pattern_list_[i].getMassShiftCount(); ++j)
      {
        std::cout << mass_pattern_list_[i].getMassShiftAt(j) << "  ";
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }
  
  std::vector<MultiplexDeltaMasses> MultiplexMassPatternList::getMassPatternList() const
  {
    return mass_pattern_list_;
  }

}
