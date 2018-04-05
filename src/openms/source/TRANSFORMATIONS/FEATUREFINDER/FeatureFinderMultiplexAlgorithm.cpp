// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringCentroided.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>

#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <vector>
#include <numeric>
#include <fstream>
#include <algorithm>

using namespace std;

namespace OpenMS
{
  FeatureFinderMultiplexAlgorithm::FeatureFinderMultiplexAlgorithm(MSExperiment& exp, bool centroided) :
    DefaultParamHandler("FeatureFinderMultiplexAlgorithm"), centroided_(centroided)
  {
    // parameter section: algorithm
    defaults_.setValue("algorithm:labels", "[][Lys8,Arg10]", "Labels used for labelling the samples. If the sample is unlabelled (i.e. you want to detect only single peptide features) please leave this parameter empty. [...] specifies the labels for a single sample. For example\n\n[][Lys8,Arg10]        ... SILAC\n[][Lys4,Arg6][Lys8,Arg10]        ... triple-SILAC\n[Dimethyl0][Dimethyl6]        ... Dimethyl\n[Dimethyl0][Dimethyl4][Dimethyl8]        ... triple Dimethyl\n[ICPL0][ICPL4][ICPL6][ICPL10]        ... ICPL");
    defaults_.setValue("algorithm:charge", "1:4", "Range of charge states in the sample, i.e. min charge : max charge.");
    defaults_.setValue("algorithm:isotopes_per_peptide", "3:6", "Range of isotopes per peptide in the sample. For example 3:6, if isotopic peptide patterns in the sample consist of either three, four, five or six isotopic peaks. ", ListUtils::create<String>("advanced"));
    defaults_.setValue("algorithm:rt_typical", 40.0, "Typical retention time [s] over which a characteristic peptide elutes. (This is not an upper bound. Peptides that elute for longer will be reported.)");
    defaults_.setMinFloat("algorithm:rt_typical", 0.0);
    defaults_.setValue("algorithm:rt_band", 10.0, "RT band which is taken into considerations when filtering.TODO docu");
    defaults_.setMinFloat("algorithm:rt_band", 0.0);
    defaults_.setValue("algorithm:rt_min", 2.0, "Lower bound for the retention time [s]. (Any peptides seen for a shorter time period are not reported.)");
    defaults_.setMinFloat("algorithm:rt_min", 0.0);
    defaults_.setValue("algorithm:mz_tolerance", 6.0, "m/z tolerance for search of peak patterns.");
    defaults_.setMinFloat("algorithm:mz_tolerance", 0.0);
    defaults_.setValue("algorithm:mz_unit", "ppm", "Unit of the 'mz_tolerance' parameter.");
    defaults_.setValidStrings("algorithm:mz_unit", ListUtils::create<String>("Da,ppm"));
    defaults_.setValue("algorithm:intensity_cutoff", 1000.0, "Lower bound for the intensity of isotopic peaks.");
    defaults_.setMinFloat("algorithm:intensity_cutoff", 0.0);
    defaults_.setValue("algorithm:peptide_similarity", 0.5, "Two peptides in a multiplet are expected to have the same isotopic pattern. This parameter is a lower bound on their similarity.");
    defaults_.setMinFloat("algorithm:peptide_similarity", -1.0);
    defaults_.setMaxFloat("algorithm:peptide_similarity", 1.0);
    defaults_.setValue("algorithm:averagine_similarity", 0.4, "The isotopic pattern of a peptide should resemble the averagine model at this m/z position. This parameter is a lower bound on similarity between measured isotopic pattern and the averagine model.");
    defaults_.setMinFloat("algorithm:averagine_similarity", -1.0);
    defaults_.setMaxFloat("algorithm:averagine_similarity", 1.0);
    defaults_.setValue("averagine_similarity_scaling", 0.75, "Let x denote this scaling factor, and p the averagine similarity parameter. For the detection of single peptides, the averagine parameter p is replaced by p' = p + x(1-p), i.e. x = 0 -> p' = p and x = 1 -> p' = 1. (For knock_out = true, peptide doublets and singlets are detected simulataneously. For singlets, the peptide similarity filter is irreleavant. In order to compensate for this 'missing filter', the averagine parameter p is replaced by the more restrictive p' when searching for singlets.)", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("algorithm:averagine_similarity_scaling", 0.0);
    defaults_.setMaxFloat("algorithm:averagine_similarity_scaling", 1.0);
    defaults_.setValue("algorithm:missed_cleavages", 0, "Maximum number of missed cleavages due to incomplete digestion. (Only relevant if enzymatic cutting site coincides with labelling site. For example, Arg/Lys in the case of trypsin digestion and SILAC labelling.)");
    defaults_.setMinInt("algorithm:missed_cleavages", 0);
    defaults_.setValue("algorithm:knock_out", "false", "Is it likely that knock-outs are present? (Supported for doublex, triplex and quadruplex experiments only.)", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("algorithm:knock_out", ListUtils::create<String>("true,false"));
    defaults_.setValue("algorithm:averagine_type","peptide","The type of averagine to use, currently RNA, DNA or peptide", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("algorithm:averagine_type", ListUtils::create<String>("peptide,RNA,DNA"));
    
    // parameter section: labels
    MultiplexDeltaMassesGenerator generator;
    Param p = generator.getParameters();
    for (Param::ParamIterator it = p.begin(); it != p.end(); ++it)
    {
      defaults_.setValue(("labels:" + it->name), it->value, it->description, ListUtils::create<String>("advanced"));
      defaults_.setMinFloat(it->name, 0.0);
      
      label_mass_shift_.insert(make_pair(it->name, it->value));
    }
    
    // parameter section: algorithm, get selected charge range
    String charge_string = param_.getValue("algorithm:charge");
    charge_min_ = charge_string.prefix(':').toInt();
    charge_max_ = charge_string.suffix(':').toInt();
    if (charge_min_ > charge_max_)
    {
      swap(charge_min_, charge_max_);
    }
    
    // parameter section: algorithm, get isotopes per peptide range
    String isotopes_per_peptide_string = param_.getValue("algorithm:isotopes_per_peptide");
    isotopes_per_peptide_min_ = isotopes_per_peptide_string.prefix(':').toInt();
    isotopes_per_peptide_max_ = isotopes_per_peptide_string.suffix(':').toInt();
    if (isotopes_per_peptide_min_ > isotopes_per_peptide_max_)
    {
      swap(isotopes_per_peptide_min_, isotopes_per_peptide_max_);
    }

    // check for empty experimental data
    if (exp.getSpectra().empty())
    {
      throw OpenMS::Exception::FileEmpty(__FILE__, __LINE__, __FUNCTION__, "Error: No MS1 spectra in input file.");
    }
    
    // update m/z and RT ranges
    exp.updateRanges();
    
    // sort according to RT and MZ
    exp.sortSpectra();

    // store experiment in member varaibles
    if (centroided_)
    {
      exp.swap(exp_centroid_);
      // exp_profile_ will never be used.
    }
    else
    {
      exp.swap(exp_profile_);
      // exp_centroid_ will be constructed later on.
    }
  }
  
  /**
   * @brief order of charge states
   *
   * 2+ 3+ 4+ 1+ 5+ 6+ ...
   *
   * Order charge states by the likelihood of their occurrence, i.e. we search for the most likely charge states first.
   */
  static size_t order_charge(int charge)
  {
    if ((1 < charge) && (charge < 5))
    {
      return (charge - 1);
    }
    else if (charge == 1)
    {
      return 4;
    }
    else
    {
      return charge;
    }
  }

  /**
   * @brief comparator of peak patterns
   *
   * The comperator determines in which order the peak patterns are searched for.
   * First we check the number of mass shifts (triplets before doublets before singlets).
   * Then we check the first mass shift (for example 6 Da before 12 Da i.e. misscleavage).
   * Finally we check for charges (2+ before 1+, most likely first).
   *
   * @param pattern1    first peak pattern
   * @param pattern2    second peak pattern
   *
   * @return true if pattern1 should be searched before pattern2
   */
  static bool less_pattern(const MultiplexIsotopicPeakPattern& pattern1, const MultiplexIsotopicPeakPattern& pattern2)
  {
    if (pattern1.getMassShiftCount() == pattern2.getMassShiftCount())
    {
      // The first mass shift is by definition always zero.
      if ((pattern1.getMassShiftCount() > 1) && (pattern2.getMassShiftCount() > 1))
      {
        if (pattern1.getMassShiftAt(1) == pattern2.getMassShiftAt(1))
        {
          // 2+ before 3+ before 4+ before 1+ before 5+ before 6+ etc.
          return order_charge(pattern1.getCharge()) < order_charge(pattern2.getCharge());
        }
        else
        {
          return pattern1.getMassShiftAt(1) < pattern2.getMassShiftAt(1);
        }
      }
      else
      {
        // 2+ before 3+ before 4+ before 1+ before 5+ before 6+ etc.
        return order_charge(pattern1.getCharge()) < order_charge(pattern2.getCharge());
      }
    }
    else
    {
      // triplets before doublets before singlets
      return pattern1.getMassShiftCount() > pattern2.getMassShiftCount();
    }
  }

  std::vector<MultiplexIsotopicPeakPattern> FeatureFinderMultiplexAlgorithm::generatePeakPatterns_(int charge_min, int charge_max, int peaks_per_peptide_max, std::vector<MultiplexDeltaMasses> mass_pattern_list)
  {
    std::vector<MultiplexIsotopicPeakPattern> list;
    
    // iterate over all charge states
    for (int c = charge_max; c >= charge_min; --c)
    {
      // iterate over all mass shifts
      for (unsigned i = 0; i < mass_pattern_list.size(); ++i)
      {
        MultiplexIsotopicPeakPattern pattern(c, peaks_per_peptide_max, mass_pattern_list[i], i);
        list.push_back(pattern);
      }
    }
    
    sort(list.begin(), list.end(), less_pattern);
    
    // debug output
    /*for (int i = 0; i < list.size(); ++i)
     {
     std::cout << "charge = " << list[i].getCharge() << "+    shift = " << list[i].getMassShiftAt(1) << " Da\n";
     }*/
    
    return list;
  }

  void FeatureFinderMultiplexAlgorithm::run()
  {
    /**
     * pick peaks (if input data are in profile mode)
     */
    std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_s; // peak boundaries for spectra
    std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_c; // peak boundaries for chromatograms
    
    if (!centroided_)
    {
      PeakPickerHiRes picker;
      Param param = picker.getParameters();
      picker.setLogType(getLogType());
      param.setValue("ms_levels", ListUtils::create<Int>("1"));
      param.setValue("signal_to_noise", 0.0); // signal-to-noise estimation switched off
      picker.setParameters(param);
      
      picker.pickExperiment(exp_profile_, exp_centroid_, boundaries_exp_s, boundaries_exp_c);
    }

    /**
     * filter for peak patterns
     */
    MultiplexDeltaMassesGenerator generator = MultiplexDeltaMassesGenerator(param_.getValue("algorithm:labels"), param_.getValue("algorithm:missed_cleavages"), label_mass_shift_);
    if (param_.getValue("algorithm:knock_out") == "true")
    {
      generator.generateKnockoutDeltaMasses();
    }
    generator.printSamplesLabelsList();
    generator.printDeltaMassesList();

    std::vector<MultiplexDeltaMasses> masses = generator.getDeltaMassesList();
    std::vector<MultiplexIsotopicPeakPattern> patterns = generatePeakPatterns_(param_.getValue("algorithm:charge_min"), param_.getValue("algorithm:charge_max"), param_.getValue("algorithm:isotopes_per_peptide_max"), masses);
    
    std::vector<MultiplexFilteredMSExperiment> filter_results;
    if (centroided_)
    {
      // centroided data
      MultiplexFilteringCentroided filtering(exp_centroid_, patterns, isotopes_per_peptide_min_, isotopes_per_peptide_max_, param_.getValue("algorithm:intensity_cutoff"), param_.getValue("algorithm:rt_band"), param_.getValue("algorithm:mz_tolerance"), (param_.getValue("algorithm:mz_unit") == "ppm"), param_.getValue("algorithm:peptide_similarity"), param_.getValue("algorithm:averagine_similarity"), param_.getValue("algorithm:averagine_similarity_scaling"), param_.getValue("algorithm:averagine_type"));
      filtering.setLogType(getLogType());
      filter_results = filtering.filter();
    }

  }
}
