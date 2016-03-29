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
// $Maintainer: Erhan Kenar $
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/DATAREDUCTION/FeatureFindingMetabo.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/SYSTEM/File.h>

#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <fstream>

#include <boost/dynamic_bitset.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

// #define FFM_DEBUG

namespace OpenMS
{
  FeatureHypothesis::FeatureHypothesis() :
    iso_pattern_(),
    feat_score_(),
    charge_()
  {

  }

  FeatureHypothesis::~FeatureHypothesis()
  {

  }

  FeatureHypothesis::FeatureHypothesis(const FeatureHypothesis& fh) :
    iso_pattern_(fh.iso_pattern_),
    feat_score_(fh.feat_score_),
    charge_(fh.charge_)
  {

  }

  FeatureHypothesis& FeatureHypothesis::operator=(const FeatureHypothesis& rhs)
  {
    if (this == &rhs)
      return *this;

    iso_pattern_ = rhs.iso_pattern_;
    feat_score_ = rhs.feat_score_;
    charge_ = rhs.charge_;

    return *this;
  }

  void FeatureHypothesis::addMassTrace(const MassTrace& mt_ptr)
  {
    iso_pattern_.push_back(&mt_ptr);
  }

  double FeatureHypothesis::getMonoisotopicFeatureIntensity(bool smoothed = false) const
  {
    if (iso_pattern_.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
          "FeatureHypothesis is empty, no traces contained!", String(iso_pattern_.size()));
    }
    return iso_pattern_[0]->getIntensity(smoothed);
  }

  double FeatureHypothesis::getSummedFeatureIntensity(bool smoothed = false) const
  {
    double int_sum(0.0);
    for (Size i = 0; i < iso_pattern_.size(); ++i)
    {
      int_sum += iso_pattern_[i]->getIntensity(smoothed);
    }
    return int_sum;
  }

  Size FeatureHypothesis::getNumFeatPoints() const
  {
    Size num_points(0);

    for (Size mt_idx = 0; mt_idx < iso_pattern_.size(); ++mt_idx)
    {
      num_points += iso_pattern_[mt_idx]->getSize();
    }

    return num_points;
  }

  std::vector<ConvexHull2D> FeatureHypothesis::getConvexHulls() const
  {
    std::vector<ConvexHull2D> tmp_hulls;
    for (Size mt_idx = 0; mt_idx < iso_pattern_.size(); ++mt_idx)
    {
      ConvexHull2D::PointArrayType hull_points(iso_pattern_[mt_idx]->getSize());

      Size i = 0;
      for (MassTrace::const_iterator l_it = iso_pattern_[mt_idx]->begin(); l_it != iso_pattern_[mt_idx]->end(); ++l_it)
      {
        hull_points[i][0] = (*l_it).getRT();
        hull_points[i][1] = (*l_it).getMZ();
        ++i;
      }

      ConvexHull2D hull;
      hull.addPoints(hull_points);
      tmp_hulls.push_back(hull);
    }
    return tmp_hulls;
  }

  FeatureFindingMetabo::FeatureFindingMetabo() :
    DefaultParamHandler("FeatureFindingMetabo"), ProgressLogger()
  {
    // defaults_.setValue( "name" , 1 , "description" );
    defaults_.setValue("quant_method", String(MassTrace::names_of_quantmethod[0]), "Method of quantification for mass traces. For LC data 'area' is recommended, 'median' for direct injection data.");
    defaults_.setValidStrings("quant_method", std::vector<String>(MassTrace::names_of_quantmethod, MassTrace::names_of_quantmethod +(int)MassTrace::SIZE_OF_MT_QUANTMETHOD));
    defaults_.setValue("local_rt_range", 10.0, "RT range where to look for coeluting mass traces", ListUtils::create<String>("advanced")); // 5.0
    defaults_.setValue("local_mz_range", 6.5, "MZ range where to look for isotopic mass traces", ListUtils::create<String>("advanced")); // 6.5
    defaults_.setValue("charge_lower_bound", 1, "Lowest charge state to consider"); // 1
    defaults_.setValue("charge_upper_bound", 3, "Highest charge state to consider"); // 3
    //defaults_.setValue("mass_error_ppm", 20.0, "Allowed mass error deviation in ppm");  // 20.0
    defaults_.setValue("chrom_fwhm", 5.0, "Expected chromatographic peak width (in seconds)."); // 5.0
    defaults_.setValue("report_summed_ints", "false", "Set to true for a feature intensity summed up over all traces rather than using monoisotopic trace intensity alone.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("report_summed_ints", ListUtils::create<String>("false,true"));
    defaults_.setValue("enable_RT_filtering", "true", "Require sufficient overlap in RT while assembling mass traces. Disable for direct injection data..");
    defaults_.setValidStrings("enable_RT_filtering", ListUtils::create<String>("false,true"));
    defaults_.setValue("disable_isotope_filtering", "false", "Disable isotope filtering.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("disable_isotope_filtering", ListUtils::create<String>("false,true"));
    defaults_.setValue("isotope_model", "metabolites", "Change type of isotope model.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("isotope_model", ListUtils::create<String>("metabolites,peptides"));

    defaults_.setValue("isotope_noisemodel", "5%RMS", "SVM isotope models were trained with either 2% or 5% RMS error. Select the appropriate noise model according to the quality of measurement or MS device.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("isotope_noisemodel", ListUtils::create<String>("5%RMS,2%RMS"));

    defaults_.setValue("use_smoothed_intensities", "true", "Use LOWESS intensities instead of raw intensities.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("use_smoothed_intensities", ListUtils::create<String>("false,true"));


    defaultsToParam_();

    this->setLogType(CMD);
  }

  FeatureFindingMetabo::~FeatureFindingMetabo()
  {

  }

  void FeatureFindingMetabo::updateMembers_()
  {
    // delta_ = (Size)param_.getValue( "delta" );

    local_rt_range_ = (double)param_.getValue("local_rt_range");
    local_mz_range_ = (double)param_.getValue("local_mz_range");
    // mass_error_ppm_ = (double)param_.getValue("mass_error_ppm");
    chrom_fwhm_ = (double)param_.getValue("chrom_fwhm");

    charge_lower_bound_ = (Size)param_.getValue("charge_lower_bound");
    charge_upper_bound_ = (Size)param_.getValue("charge_upper_bound");

    report_summed_ints_ = param_.getValue("report_summed_ints").toBool();
    enable_RT_filtering_ = param_.getValue("enable_RT_filtering").toBool();
    disable_isotope_filtering_ = param_.getValue("disable_isotope_filtering").toBool();
    isotope_model_ = param_.getValue("isotope_model");
    metabo_iso_noisemodel_ = (String)param_.getValue("isotope_noisemodel");
    use_smoothed_intensities_ = param_.getValue("use_smoothed_intensities").toBool();
  }

  double FeatureFindingMetabo::computeAveragineSimScore_(const std::vector<double>& hypo_ints, const double& mol_weight)
  {
    IsotopeDistribution isodist(hypo_ints.size());
    isodist.estimateFromPeptideWeight(mol_weight);
    // isodist.renormalize();

    std::vector<std::pair<Size, double> > averagine_dist = isodist.getContainer();

    // std::vector<double> hypo_ints = feat_hypo.getAllIntensities();

    double max_int(0.0), theo_max_int(0.0);
    for (Size i = 0; i < hypo_ints.size(); ++i)
    {
      if (hypo_ints[i] > max_int)
      {
        max_int = hypo_ints[i];
      }

      if (averagine_dist[i].second > theo_max_int)
      {
        theo_max_int = averagine_dist[i].second;
      }
    }

    // compute normalized intensities
    std::vector<double> averagine_ratios, hypo_isos;
    for (Size i = 0; i < hypo_ints.size(); ++i)
    {
      averagine_ratios.push_back(averagine_dist[i].second / theo_max_int);
      hypo_isos.push_back(hypo_ints[i] / max_int);
    }

    double iso_score = computeCosineSim_(averagine_ratios, hypo_isos);
    return iso_score;
  }

  bool FeatureFindingMetabo::isLegalIsotopePattern_(const FeatureHypothesis& feat_hypo) const
  {
    if (feat_hypo.getSize() == 1)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
          "Cannot compute isotope pattern on a single mass trace!", String(feat_hypo.getSize()));
    }

    std::vector<double> all_ints = feat_hypo.getAllIntensities(use_smoothed_intensities_);
    double mono_int(all_ints[0]);

    svm_node* nodes;

    nodes = new svm_node[7];


    nodes[0].index = 1;
    nodes[0].value = (feat_hypo.getCentroidMZ() - svm_feat_centers_[0]) / svm_feat_scales_[0];

    Size i = 2;

    Size feat_size(feat_hypo.getSize());

    if (feat_size > 6)
    {
      feat_size = 6;
    }

    for (; i - 1 < feat_size; ++i)
    {
      nodes[i - 1].index = static_cast<Int>(i);

      double ratio((all_ints[i - 1] / mono_int));

      // std::cout << i << " " << ratio << " " << std::flush;

      if (ratio > 1.0)
      {
        delete[] nodes;
        return false;
      }

      double tmp_val((ratio - svm_feat_centers_[i - 1]) / svm_feat_scales_[i - 1]);
      nodes[i - 1].value = tmp_val;
    }


    for (; i < 7; ++i)
    {
      nodes[i - 1].index = static_cast<Int>(i);
      nodes[i - 1].value = (-svm_feat_centers_[i - 1]) / svm_feat_scales_[i - 1];
    }

    nodes[6].index = -1;
    nodes[6].value = 0;

    double predict = svm_predict(isotope_filt_svm_, nodes);

    delete[] nodes;

    return (predict == 2.0) ? true : false;
  }

  bool FeatureFindingMetabo::isLegalIsotopePattern2_(const FeatureHypothesis& feat_hypo) const
  {
    if (feat_hypo.getSize() == 1)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
          "Cannot compute isotope pattern on a single mass trace!", String(feat_hypo.getSize()));
    }

    std::vector<double> all_ints = feat_hypo.getAllIntensities(use_smoothed_intensities_);
    double mono_int(all_ints[0]); // monoisotopic intensity

    // Limit the feature size to 4, since the model was only trained on
    // monoisotopic + 3 isotopic traces
    const Size FEAT_NUM(4);
    Size feat_size(feat_hypo.getSize());
    if (feat_size > 4)
    {
      feat_size = 4;
    }

    svm_node* nodes;
    nodes = new svm_node[FEAT_NUM + 1];

    double charge(feat_hypo.getCharge());
    double act_mass(feat_hypo.getCentroidMZ() * charge);

    // isotope model currently restricted to formulas up to 1000 Da
    if (act_mass > 1000.0)
    {
      act_mass = 1000.0;
    }

    nodes[0].index = 1;
    nodes[0].value = (act_mass - svm_feat_centers_[0]) / svm_feat_scales_[0];

    // Iterate, start with first isotopic trace (skip monoisotopic)
    Size i = 2;
    for (; i - 1 < feat_size; ++i)
    {
      nodes[i - 1].index = static_cast<Int>(i);

      // compute ratio of trace to monoisotopic intensity
      double ratio((all_ints[i - 1] / mono_int));

      // std::cout << i << " " << ratio << " " << std::endl;
      // if (ratio > 1.0)
      // {
      //     delete[] nodes;
      //     return false;
      // }

      double tmp_val((ratio - svm_feat_centers_[i - 1]) / svm_feat_scales_[i - 1]);
      nodes[i - 1].value = tmp_val;
    }

    for (; i < FEAT_NUM + 1; ++i)
    {
      nodes[i - 1].index = static_cast<Int>(i);
      nodes[i - 1].value = (-svm_feat_centers_[i - 1]) / svm_feat_scales_[i - 1];
    }

    nodes[FEAT_NUM].index = -1;
    nodes[FEAT_NUM].value = 0;

    // debug output
    //    std::cout << "isocheck for " << feat_hypo.getLabel() << " " << feat_hypo.getSize() << std::endl;
    //    for (Size i = 0; i < FEAT_NUM + 1; ++i)
    //    {
    //        std::cout << "idx: " << nodes[i].index << " val: " << nodes[i].value << std::endl;
    //    }

    // Use SVM model to predict the category in which the current trace group
    // belongs ...
    double predict = svm_predict(isotope_filt_svm_, nodes);

    // std::cout << "predict: " << predict << std::endl;
    delete[] nodes;

    return (predict == 2.0) ? true : false;
  }

  void FeatureFindingMetabo::loadIsotopeModel_(const String& model_name)
  {
    String search_name("CHEMISTRY/" + model_name);

    std::string model_filename = File::find(search_name + ".svm");
    std::string scale_filename = File::find(search_name + ".scale");

    isotope_filt_svm_ = svm_load_model(model_filename.c_str());
    if (isotope_filt_svm_ == NULL)
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__,
          "Loading " + model_filename + " failed", model_filename);
    }

    std::ifstream ifs(scale_filename.c_str());

    std::string line;
    std::stringstream str_buf;
    std::istream_iterator<double> eol;

    svm_feat_centers_.clear();
    svm_feat_scales_.clear();

    while (getline(ifs, line))
    {
      str_buf.clear();
      str_buf << line;
      std::istream_iterator<double> istr_it(str_buf);

      while (istr_it != eol)
      {
        svm_feat_centers_.push_back(*istr_it);
        ++istr_it;
        svm_feat_scales_.push_back(*istr_it);
        ++istr_it;
      }
    }

    if (svm_feat_centers_.size() != svm_feat_scales_.size())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
          "Numbers of centers and scales from file " + scale_filename + " are different!",
          String(svm_feat_centers_.size()) + " and " + String(svm_feat_scales_.size()));
    }
  }

  double FeatureFindingMetabo::scoreMZ_(const MassTrace& tr1, const MassTrace& tr2, Size iso_pos, Size charge) const
  {
    double mu((1.000857 * (double)iso_pos + 0.001091) / (double)charge);
    double sd((0.0016633 * (double)iso_pos - 0.0004751) / (double)charge);

    double mz1(tr1.getCentroidMZ());
    double mz2(tr2.getCentroidMZ());

    // double centered_mz(std::fabs(mz2 - mz1) - mu);
    double diff_mz(std::fabs(mz2 - mz1));

    double mt_sigma1(tr1.getCentroidSD());
    double mt_sigma2(tr2.getCentroidSD());
    // double mt_variances1(mt_sigma1*mt_sigma1 + mt_sigma2*mt_sigma2);
    double mt_variances(std::exp(2 * std::log(mt_sigma1)) + std::exp(2 * std::log(mt_sigma2)));
    // std::cout << "mt1: " << mt_sigma1 << " mt2: " << mt_sigma2 << " mt_variances: " << mt_variances << " old " << mt_variances1 <<  std::endl;

    // double score_sigma_old(std::sqrt(sd*sd + mt_variances));

    double score_sigma(std::sqrt(std::exp(2 * std::log(sd)) + mt_variances));

    // std::cout << std::setprecision(15) << "old " << score_sigma_old << " new " << score_sigma << std::endl;

    double sigma_mult(3.0);

    double mz_score(0.0);


    if ((diff_mz < mu + sigma_mult * score_sigma) && (diff_mz > mu - sigma_mult * score_sigma))
    {
      double tmp_exponent((diff_mz - mu) / score_sigma);
      mz_score = std::exp(-0.5 * tmp_exponent * tmp_exponent);

    }

    // std::cout << tr1.getLabel() << "_" << tr2.getLabel() << " diffmz: " << diff_mz << " charge " << charge << " isopos: " << iso_pos << " score: " << mz_score << std::endl ;

    return mz_score;
  }

  double FeatureFindingMetabo::scoreMZ2_(const MassTrace& tr1, const MassTrace& tr2, Size iso_pos, Size charge) const
  {
    double mu((1.003355 * (double)iso_pos) / (double)charge); 
    double sd(0.01 / (double)charge);
    // only difference to above: model parameters

    double mz1(tr1.getCentroidMZ());
    double mz2(tr2.getCentroidMZ());

    // double centered_mz(std::fabs(mz2 - mz1) - mu);
    double diff_mz(std::fabs(mz2 - mz1));

    double mt_sigma1(tr1.getCentroidSD());
    double mt_sigma2(tr2.getCentroidSD());
    // double mt_variances1(mt_sigma1*mt_sigma1 + mt_sigma2*mt_sigma2);
    double mt_variances(std::exp(2 * std::log(mt_sigma1)) + std::exp(2 * std::log(mt_sigma2)));
    // std::cout << "mt1: " << mt_sigma1 << " mt2: " << mt_sigma2 << " mt_variances: " << mt_variances << " old " << mt_variances1 <<  std::endl;

    // double score_sigma_old(std::sqrt(sd*sd + mt_variances));

    double score_sigma(std::sqrt(std::exp(2 * std::log(sd)) + mt_variances));

    // std::cout << std::setprecision(15) << "old " << score_sigma_old << " new " << score_sigma << std::endl;

    double sigma_mult(3.0);

    double mz_score(0.0);


    if ((diff_mz < mu + sigma_mult * score_sigma) && (diff_mz > mu - sigma_mult * score_sigma))
    {
      double tmp_exponent((diff_mz - mu) / score_sigma);
      mz_score = std::exp(-0.5 * tmp_exponent * tmp_exponent);

    }

    // std::cout << tr1.getLabel() << "_" << tr2.getLabel() << " diffmz: " << diff_mz << " charge " << charge << " isopos: " << iso_pos << " score: " << mz_score << std::endl ;

    return mz_score;
  }

    /// Not used any more ???  -> seems to be old model with special treatment of isotopic position 1
    /// TODO: remove

//double FeatureFindingMetabo::scoreMZ_(const MassTrace& tr1, const MassTrace& tr2, Size iso_pos, Size charge)
//{
//    double sigma_mult(3.0);
//    double center((1.001156 * (double)iso_pos + 0.001349) / (double)charge);

//    double mz1(tr1.getCentroidMZ());
//    double mz2(tr2.getCentroidMZ());

//    double centered_mz(std::fabs(mz2 - mz1) - center);

//    // setup gaussian mixture model for valid isotopic m/z distances
//    //double iso_mz_diff(1.002245);

//    // std::cout << "iso_pos: " << iso_pos << " charge: " << charge << " mz_diff" << centered_mz << std::endl;

//    double mu1(0.000981383);
//    double sigma1(0.001657985);
//    double mu2(-0.00559452);
//    double sigma2(0.00118085);

//    double mt_sigma1(tr1.getCentroidSD());
//    double mt_sigma2(tr2.getCentroidSD());
//    // double mt_variances1(mt_sigma1*mt_sigma1 + mt_sigma2*mt_sigma2);
//    double mt_variances(std::exp(2 * std::log(mt_sigma1)) + std::exp(2 * std::log(mt_sigma2)));
//    // std::cout << "mt1: " << mt_sigma1 << " mt2: " << mt_sigma2 << " mt_variances: " << mt_variances << " old " << mt_variances1 <<  std::endl;


//    double score_sigma1(std::sqrt(sigma1 * sigma1 + mt_variances));
//    double score_sigma2(std::sqrt(sigma2 * sigma2 + mt_variances));

//    // std::cout << "score_sigma1: " << score_sigma1 << std::endl;
//    // std::cout << "score_sigma2: " << score_sigma2 << std::endl;


//    double mz_score(0.0);


//    if (iso_pos == 1)
//    {
//        if ((centered_mz < mu1 + sigma_mult * score_sigma1) && (centered_mz > mu1 - sigma_mult * score_sigma1))
//        {
//            double tmp_exponent1((centered_mz - mu1) / score_sigma1);
//            mz_score = std::exp(-0.5 * tmp_exponent1 * tmp_exponent1);

//        }
//    }
//    else
//    {
//        if ((centered_mz < mu1 + sigma_mult * score_sigma1) && (centered_mz > mu2 - sigma_mult * score_sigma2))
//        {
//            double tmp_exponent1((centered_mz - mu1) / score_sigma1);
//            double tmp_exponent2((centered_mz - mu2) / score_sigma2);

//            double mz_score1(std::exp(-0.5 * tmp_exponent1 * tmp_exponent1));
//            double mz_score2(std::exp(-0.5 * tmp_exponent2 * tmp_exponent2));

//            mz_score = (mz_score1 > mz_score2) ? mz_score1 : mz_score2;

//        }
//    }
//    // std::cout<< tr1.getLabel() << "_" << tr2.getLabel() << " mz score: " << mz_score << std::endl;

//    // double diff_mz(mz2-mz1);
//    // std::cout << tr1.getLabel() << "_" << tr2.getLabel() << " diffmz: " << diff_mz << " charge " << charge << " isopos: " << iso_pos << " score: " << mz_score << std::endl ;

//    return mz_score;
//}

  double FeatureFindingMetabo::scoreRT_(const MassTrace& tr1, const MassTrace& tr2) const
  {
    // return success if this filter is disabled
    if (!enable_RT_filtering_) return 1.0;

    // continue to check overlap and cosine similarity
    // ...
    std::map<double, std::vector<double> > coinciding_rts;

    std::pair<Size, Size> tr1_fwhm_idx(tr1.getFWHMborders());
    std::pair<Size, Size> tr2_fwhm_idx(tr2.getFWHMborders());

    //    std::cout << tr1_fwhm_idx.first << " " << tr1_fwhm_idx.second << std::endl;
    //    std::cout << tr2_fwhm_idx.first << " " << tr2_fwhm_idx.second << std::endl;

    //    Size tr1_fwhm_size(tr1_fwhm_idx.second - tr1_fwhm_idx.first);
    //    Size tr2_fwhm_size(tr2_fwhm_idx.second - tr2_fwhm_idx.first);

    //    double max_length = (tr1_fwhm_size > tr2_fwhm_size) ? tr1_fwhm_size : tr2_fwhm_size;

    double tr1_length(tr1.getFWHM());
    double tr2_length(tr2.getFWHM());
    double max_length = (tr1_length > tr2_length) ? tr1_length : tr2_length;

    // std::cout << "tr1 " << tr1_length << " tr2 " << tr2_length << std::endl;

    // Extract peak shape between FWHM borders for both peaks
    for (Size i = tr1_fwhm_idx.first; i <= tr1_fwhm_idx.second; ++i)
    {
      coinciding_rts[tr1[i].getRT()].push_back(tr1[i].getIntensity());
    }
    for (Size i = tr2_fwhm_idx.first; i <= tr2_fwhm_idx.second; ++i)
    {
      coinciding_rts[tr2[i].getRT()].push_back(tr2[i].getIntensity());
    }

    // Look at peaks at the same RT 
    // TODO: this only works if both traces are sampled with equal rate at the same RT
    std::vector<double> x, y, overlap_rts;
    for (std::map<double, std::vector<double> >::const_iterator m_it = coinciding_rts.begin(); m_it != coinciding_rts.end(); ++m_it)
    {
      if (m_it->second.size() == 2)
      {
        x.push_back(m_it->second[0]);
        y.push_back(m_it->second[1]);
        overlap_rts.push_back(m_it->first);
      }
    }

    //    if (x.size() < std::floor(0.8*max_length))
    //        {
    //            return 0.0;
    //        }
    // double rt_range(0.0)
    // if (coinciding_rts.size() > 0)
    // {
    //     rt_range = std::fabs(coinciding_rts.rbegin()->first - coinciding_rts.begin()->first);
    // }


    double overlap(0.0);
    if (overlap_rts.size() > 0)
    {
      double start_rt(*(overlap_rts.begin())), end_rt(*(overlap_rts.rbegin()));
      overlap = std::fabs(end_rt - start_rt);
    }

    double proportion(overlap / max_length);
    if (proportion < 0.7)
    {
      return 0.0;
    }
    return computeCosineSim_(x, y);
  }

  double FeatureFindingMetabo::computeCosineSim_(const std::vector<double>& x, const std::vector<double>& y) const
  {
    if (x.size() != y.size())
    {
      return 0.0;
    }

    double mixed_sum(0.0);
    double x_squared_sum(0.0);
    double y_squared_sum(0.0);

    for (Size i = 0; i < x.size(); ++i)
    {
      mixed_sum += x[i] * y[i];
      x_squared_sum += x[i] * x[i];
      y_squared_sum += y[i] * y[i];
    }

    double denom(std::sqrt(x_squared_sum) * std::sqrt(y_squared_sum));
    return (denom > 0.0) ? mixed_sum / denom : 0.0;
  }

  double FeatureFindingMetabo::computeOLSCoeff_(const std::vector<double>& x, const std::vector<double>& y) const
  {
    if (x.size() != y.size())
    {
      return 0.0;
    }

    double mixed_sum(0.0);
    double x_squared_sum(0.0);

    for (Size i = 0; i < x.size(); ++i)
    {
      mixed_sum += x[i] * y[i];
      x_squared_sum += x[i] * x[i];
    }

    return (x_squared_sum > 0.0) ? mixed_sum / x_squared_sum : 0.0;
  }

  void FeatureFindingMetabo::findLocalFeatures_(const std::vector<const MassTrace*>& candidates, std::vector<FeatureHypothesis>& output_hypos)
  {
    // single Mass trace hypothesis
    FeatureHypothesis tmp_hypo;
    tmp_hypo.addMassTrace(*candidates[0]);
    tmp_hypo.setScore((candidates[0]->getIntensity(use_smoothed_intensities_)) / total_intensity_);

#ifdef _OPENMP
#pragma omp critical (OPENMS_FFMetabo_output_hypos)
#endif
    {
      // pushing back to shared vector needs to be synchronized
      output_hypos.push_back(tmp_hypo);
    }

    for (Size charge = charge_lower_bound_; charge <= charge_upper_bound_; ++charge)
    {
      FeatureHypothesis fh_tmp;
      fh_tmp.addMassTrace(*candidates[0]);
      fh_tmp.setScore((candidates[0]->getIntensity(use_smoothed_intensities_)) / total_intensity_);

      // double mono_iso_rt(candidates[0]->getCentroidRT());
      // double mono_iso_mz(candidates[0]->getCentroidMZ());
      // double mono_iso_int(candidates[0]->computePeakArea());

      Size last_iso_idx(0);
      Size iso_pos_max(std::floor(charge * local_mz_range_));
      for (Size iso_pos = 1; iso_pos <= iso_pos_max; ++iso_pos)
      {

        // Find mass trace that best agrees with current hypothesis of charge
        // and isotopic position
        double best_so_far(0.0);
        Size best_idx(0);
        for (Size mt_idx = last_iso_idx + 1; mt_idx < candidates.size(); ++mt_idx)
        {
          // double tmp_iso_rt(candidates[mt_idx]->getCentroidRT());
          // double tmp_iso_mz(candidates[mt_idx]->getCentroidMZ());
          // double tmp_iso_int(candidates[mt_idx]->computePeakArea());

#ifdef FFM_DEBUG
          std::cout << "scoring " << candidates[0]->getLabel() << " " << candidates[0]->getCentroidMZ() << 
            " with " << candidates[mt_idx]->getLabel() << " " << candidates[mt_idx]->getCentroidMZ() << std::endl;
#endif

          // Score current mass trace candidates against hypothesis
          double rt_score(scoreRT_(*candidates[0], *candidates[mt_idx]));
          double mz_score(scoreMZ_(*candidates[0], *candidates[mt_idx], iso_pos, charge));

          // disable intensity scoring for now...
          double int_score(1.0);
          // double int_score((candidates[0]->getIntensity(use_smoothed_intensities_))/total_weight + (candidates[mt_idx]->getIntensity(use_smoothed_intensities_))/total_weight);

          if (isotope_model_ == "peptides")
          {
            std::vector<double> tmp_ints(fh_tmp.getAllIntensities());
            tmp_ints.push_back(candidates[mt_idx]->getIntensity(use_smoothed_intensities_));
            int_score = computeAveragineSimScore_(tmp_ints, candidates[mt_idx]->getCentroidMZ() * charge);
          }

#ifdef FFM_DEBUG
          std::cout << fh_tmp.getLabel() << "_" << candidates[mt_idx]->getLabel() << 
            "\t" << "ch: " << charge << " isopos: " << iso_pos << " rt: " << 
            rt_score << "mz: " << mz_score << "int: " << int_score << std::endl;
#endif

          double total_pair_score(0.0);
          if (rt_score > 0.0 && mz_score > 0.0 && int_score > 0.0)
          {
            total_pair_score = std::exp(std::log(rt_score) + log(mz_score) + log(int_score));
          }
          if (total_pair_score > best_so_far)
          {
            best_so_far = total_pair_score;
            best_idx = mt_idx;
          }
        } // end mt_idx

        // Store mass trace that best agrees with current hypothesis of charge
        // and isotopic position
        if (best_so_far > 0.0)
        {
          fh_tmp.addMassTrace(*candidates[best_idx]);
          double weighted_score(((candidates[best_idx]->getIntensity(use_smoothed_intensities_)) * best_so_far) / total_intensity_);

          fh_tmp.setScore(fh_tmp.getScore() + weighted_score);
          fh_tmp.setCharge(charge);
          last_iso_idx = best_idx;

#ifdef _OPENMP
#pragma omp critical (OPENMS_FFMetabo_output_hypos)
#endif
          {
            // pushing back to shared vector needs to be synchronized
            output_hypos.push_back(fh_tmp);
          }
        }
        else
        {
          break;
        }
      } // end for iso_pos

#ifdef FFM_DEBUG
      std::cout << "best found for ch " << charge << ":" << fh_tmp.getLabel() << " score: " << fh_tmp.getScore() << std::endl;
#endif
    } // end for charge
  } // end of findLocalFeatures_(...)

  void FeatureFindingMetabo::run(std::vector<MassTrace>& input_mtraces, FeatureMap& output_featmap)
  {
    if (input_mtraces.empty()) 
    {
      return;
    }

    // mass traces must be sorted by their centroid MZ
    std::sort(input_mtraces.begin(), input_mtraces.end(), CmpMassTraceByMZ());

    this->startProgress(0, input_mtraces.size(), "assembling mass traces to features");

    // *********************************************************** //
    // Step 1 configure quantification method
    // *********************************************************** //
    MassTrace::MT_QUANTMETHOD method = MassTrace::getQuantMethod((String)param_.getValue("quant_method"));
    for (std::vector<MassTrace>::iterator it = input_mtraces.begin();
      it != input_mtraces.end();
      ++it)
    {
      it->setQuantMethod(method);
    }

    // *********************************************************** //
    // Step 2 initialize SVM model for isotope ratio filtering
    // *********************************************************** //
    if (metabo_iso_noisemodel_ == "2%RMS")
    {
      LOG_INFO << "Loading metabolite isotope model with 2% RMS error" << std::endl;
      loadIsotopeModel_("MetaboliteIsoModelNoised2");
    }
    else
    {
      LOG_INFO << "Loading metabolite isotope model with 5% RMS error" << std::endl;
      loadIsotopeModel_("MetaboliteIsoModelNoised5");
    }

    total_intensity_ = 0.0;
    for (Size i = 0; i < input_mtraces.size(); ++i)
    {
      total_intensity_ += input_mtraces[i].getIntensity(use_smoothed_intensities_);
    }

    // *********************************************************** //
    // Step 3 Iterate through all mass traces to find likely matches 
    // and generate isotopic / charge hypotheses
    // *********************************************************** //
    std::vector<FeatureHypothesis> feat_hypos;
    Size progress(0);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize i = 0; i < (SignedSize)input_mtraces.size(); ++i)
    {
      IF_MASTERTHREAD this->setProgress(progress);
#ifdef _OPENMP
#pragma omp atomic
#endif
      ++progress;

      std::vector<const MassTrace*> local_traces;
      double ref_trace_mz(input_mtraces[i].getCentroidMZ());
      double ref_trace_rt(input_mtraces[i].getCentroidRT());

      local_traces.push_back(&input_mtraces[i]);

      for (Size ext_idx = i + 1; ext_idx < input_mtraces.size(); ++ext_idx)
      {
        // traces are sorted by m/z, so we can break when we leave the allowed window
        double diff_mz = std::fabs(input_mtraces[ext_idx].getCentroidMZ() - ref_trace_mz);
        if (diff_mz > local_mz_range_) break;

        double diff_rt = std::fabs(input_mtraces[ext_idx].getCentroidRT() - ref_trace_rt);
        if (diff_rt <= local_rt_range_)
        {
          // std::cout << " accepted!" << std::endl;
          local_traces.push_back(&input_mtraces[ext_idx]);
        }
      }
      findLocalFeatures_(local_traces, feat_hypos);
    }
    this->endProgress();

    // sort feature candidates by their score
    std::sort(feat_hypos.begin(), feat_hypos.end(), CmpHypothesesByScore());

#ifdef FFM_DEBUG
    std::cout << "size of hypotheses: " << feat_hypos.size() << std::endl;
    // output all hypotheses:
    for (Size hypo_idx = 0; hypo_idx < feat_hypos.size(); ++ hypo_idx)
    {
      bool legal(false);
      if (feat_hypos[hypo_idx].getSize() > 1)
      {
        legal = isLegalIsotopePattern_(feat_hypos[hypo_idx]);
      }
      std::cout << feat_hypos[hypo_idx].getLabel() << " ch: " << feat_hypos[hypo_idx].getCharge() << 
        " score: " << feat_hypos[hypo_idx].getScore() << " legal: " << legal << std::endl;
    }
#endif

    // *********************************************************** //
    // Step 4 Iterate through all hypotheses, starting with the highest 
    // scoring one. Accept them if they do not contain traces that have 
    // already been used by a higher scoring hypothesis.
    // *********************************************************** //
    std::map<String, bool> trace_excl_map;
    for (Size hypo_idx = 0; hypo_idx < feat_hypos.size(); ++hypo_idx)
    {
      // std::cout << "score now: " <<  feat_hypos[hypo_idx].getScore() << std::endl;
      std::vector<String> labels(feat_hypos[hypo_idx].getLabels());
      bool trace_coll = false;   // trace collision?
      for (Size lab_idx = 0; lab_idx < labels.size(); ++lab_idx)
      {
        if (trace_excl_map.find(labels[lab_idx]) != trace_excl_map.end())
        {
          trace_coll = true;
          break;
        }
      }

#ifdef FFM_DEBUG
      if (feat_hypos[hypo_idx].getSize() > 1)
      {
        std::cout << "check for collision: " << trace_coll << " " << 
          feat_hypos[hypo_idx].getLabel() << " " << isLegalIsotopePattern_(feat_hypos[hypo_idx]) << 
          " " << feat_hypos[hypo_idx].getScore() << std::endl;
      }
#endif

      // Skip hypotheses that contain a mass trace that has already been used
      if (trace_coll) 
      {
        continue;
      }

      // Check whether the trace  passes the intensity filter (metabolites
      // only). This is based on a pre-trained SVM model of isotopic
      // intensities.
      bool pass_isotope_filter = true;
      if (feat_hypos[hypo_idx].getSize() > 1)
      {
        if (!disable_isotope_filtering_)
        {
          if (isotope_model_ == "metabolites")
          {
            pass_isotope_filter = isLegalIsotopePattern2_(feat_hypos[hypo_idx]);
          }
          else if (isotope_model_ == "peptides")
          {
            pass_isotope_filter = true;
          }
        }
        // std::cout << "\nlegal iso? " << feat_hypos[hypo_idx].getLabel() << " score: " << feat_hypos[hypo_idx].getScore() << " " << result << std::endl;
      }

      if (!pass_isotope_filter) 
      {
        continue;
      }

      //
      // Now accept hypothesis
      //

      Feature f;
      f.setRT(feat_hypos[hypo_idx].getCentroidRT());
      f.setMZ(feat_hypos[hypo_idx].getCentroidMZ());

      if (report_summed_ints_)
      {
        f.setIntensity(feat_hypos[hypo_idx].getSummedFeatureIntensity(use_smoothed_intensities_));
      }
      else
      {
        f.setIntensity(feat_hypos[hypo_idx].getMonoisotopicFeatureIntensity(use_smoothed_intensities_));
      }

      f.setWidth(feat_hypos[hypo_idx].getFWHM());
      f.setCharge(feat_hypos[hypo_idx].getCharge());
      f.setMetaValue(3, feat_hypos[hypo_idx].getLabel());

      // store isotope intensities
      std::vector<double> all_ints(feat_hypos[hypo_idx].getAllIntensities(use_smoothed_intensities_));

      f.setMetaValue("num_of_masstraces", all_ints.size());
      for (Size int_idx = 0; int_idx < all_ints.size(); ++int_idx)
      {
        f.setMetaValue("masstrace_intensity_" + String(int_idx), all_ints[int_idx]);
      }

      // TODO add flag to skip this step, blows up in memory and increases size of featureXML ...
      f.setConvexHulls(feat_hypos[hypo_idx].getConvexHulls());
      f.setOverallQuality(feat_hypos[hypo_idx].getScore());

      output_featmap.push_back(f);

      // add used traces to exclusion map
      for (Size lab_idx = 0; lab_idx < labels.size(); ++lab_idx)
      {
        trace_excl_map[labels[lab_idx]] = true;
      }
    }
  } // end of FeatureFindingMetabo::run

}
