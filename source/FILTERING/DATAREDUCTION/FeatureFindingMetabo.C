// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Erhan Kenar $
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/DATAREDUCTION/FeatureFindingMetabo.h>

#include <OpenMS/SYSTEM/File.h>

#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <fstream>


#include <boost/dynamic_bitset.hpp>

namespace OpenMS
{
FeatureHypothesis::FeatureHypothesis()
    : iso_pattern_(),
      feat_score_(),
      charge_()
{

}

FeatureHypothesis::~FeatureHypothesis()
{

}


FeatureHypothesis::FeatureHypothesis(const FeatureHypothesis& fh)
    : iso_pattern_(fh.iso_pattern_),
      feat_score_(fh.feat_score_),
      charge_(fh.charge_)
{

}


FeatureHypothesis& FeatureHypothesis::operator=(const FeatureHypothesis& rhs)
{
    if (this==&rhs) return *this;
    iso_pattern_ = rhs.iso_pattern_;
    feat_score_ = rhs.feat_score_;
    charge_ = rhs.charge_;

    return *this;
}


void FeatureHypothesis::addMassTrace(MassTrace& mt_ptr)
{
    iso_pattern_.push_back(&mt_ptr);

    return ;
}


DoubleReal FeatureHypothesis::getMonoisotopicFeatureIntensity(bool smoothed = false)
{
    if (iso_pattern_.empty())
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "FeatureHypothesis is empty, no traces contained!", String(iso_pattern_.size()));
    }

    DoubleReal result;

    if (smoothed)
    {
        result = iso_pattern_[0]->computeSmoothedPeakArea();
    }
    else
    {
        result = iso_pattern_[0]->computePeakArea();
    }

    return result;
}


DoubleReal FeatureHypothesis::getSummedFeatureIntensity(bool smoothed = false)
{
    DoubleReal int_sum(0.0);

    for (Size i = 0; i < iso_pattern_.size(); ++i)
    {
        if (smoothed)
        {
            int_sum += iso_pattern_[i]->computeSmoothedPeakArea();
        }
        else
        {
            int_sum += iso_pattern_[i]->computePeakArea();
        }
    }

    return int_sum;
}


Size FeatureHypothesis::getNumFeatPoints() const
{
    Size num_points(0);

    for (Size mt_idx = 0; mt_idx < iso_pattern_.size(); ++mt_idx) {
        num_points += iso_pattern_[mt_idx]->getSize();
    }

    return num_points;
}


std::vector<ConvexHull2D> FeatureHypothesis::getConvexHulls() const
{
    std::vector<ConvexHull2D> tmp_hulls;

    for (Size mt_idx = 0; mt_idx < iso_pattern_.size(); ++mt_idx) {

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


FeatureFindingMetabo::FeatureFindingMetabo()
    : DefaultParamHandler("FeatureFindingMetabo"), ProgressLogger()
{
    // defaults_.setValue( "name" , 1 , "descript" );
    defaults_.setValue("local_rt_range" , 5.0 , "RT range where to look for coeluting mass traces", StringList::create("advanced")); // 5.0
    defaults_.setValue("local_mz_range" , 6.5 , "MZ range where to look for isotopic mass traces", StringList::create("advanced"));    // 6.5
    defaults_.setValue("charge_lower_bound" , 1 , "Lowest charge state to consider");    // 1
    defaults_.setValue("charge_upper_bound" , 5 , "Highest charge state to consider");    // 5
    //defaults_.setValue("mass_error_ppm", 20.0, "Allowed mass error deviation in ppm");  // 20.0
    //defaults_.setValue("chrom_fwhm" , 3.0 , "Expected chromatographic peak width (in seconds).");    // 3.0
    defaults_.setValue("report_summed_ints", "false", "Set to true for a feature intensity summed up over all traces rather than using monoisotopic trace intensity alone.", StringList::create("advanced"));
    defaults_.setValidStrings("report_summed_ints", StringList::create(("false,true")));
    defaults_.setValue("disable_isotope_filtering", "false", "Disable metabolite isotope filtering.", StringList::create("advanced"));
    defaults_.setValidStrings("disable_isotope_filtering", StringList::create(("false,true")));
    defaults_.setValue("use_smoothed_intensities", "false", "Use LOWESS intensities instead of raw intensities.", StringList::create("advanced"));
    defaults_.setValidStrings("use_smoothed_intensities", StringList::create(("false,true")));


    defaultsToParam_();

    this->setLogType(CMD);
}


FeatureFindingMetabo::~FeatureFindingMetabo()
{

}



void FeatureFindingMetabo::updateMembers_()
{
    // delta_ = (Size)param_.getValue( "delta" );

    local_rt_range_ = (DoubleReal)param_.getValue("local_rt_range");
    local_mz_range_ = (DoubleReal)param_.getValue("local_mz_range");
    // mass_error_ppm_ = (DoubleReal)param_.getValue("mass_error_ppm");
    // chrom_fwhm_ = (DoubleReal)param_.getValue("chrom_fwhm");

    charge_lower_bound_ = (Size)param_.getValue("charge_lower_bound");
    charge_upper_bound_ = (Size)param_.getValue("charge_upper_bound");

    report_summed_ints_ = param_.getValue("report_summed_ints").toBool();
    disable_isotope_filtering_ = param_.getValue("disable_isotope_filtering").toBool();
    use_smoothed_intensities_ = param_.getValue("use_smoothed_intensities").toBool();
}


bool FeatureFindingMetabo::isLegalIsotopePattern_(FeatureHypothesis& feat_hypo)
{
    if (feat_hypo.getSize() == 1)
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Cannot compute isotope pattern on a single mass trace!", String(feat_hypo.getSize()));
    }

    std::vector<DoubleReal> all_ints = feat_hypo.getAllIntensities(use_smoothed_intensities_);

    DoubleReal mono_int(all_ints[0]);

    svm_node* nodes;

    nodes = new svm_node[7];


    nodes[0].index = 1;
    nodes[0].value = (feat_hypo.getCentroidMZ() - svm_feat_centers[0])/svm_feat_scales[0];

    Size i = 2;

    Size feat_size(feat_hypo.getSize());

    if (feat_size > 6)
    {
        feat_size = 6;
    }

    for ( ; i - 1 < feat_size; ++i)
    {
        nodes[i - 1].index = i;

        DoubleReal ratio((all_ints[i - 1] / mono_int));

        // std::cout << i << " " << ratio << " " << std::flush;

        if (ratio > 1.0)
        {
            return false;
        }

        DoubleReal tmp_val((ratio - svm_feat_centers[i - 1])/svm_feat_scales[i - 1]);
        nodes[i - 1].value = tmp_val;
    }


    for ( ; i < 7; ++i)
    {
        nodes[i - 1].index = i;
        nodes[i - 1].value = (-svm_feat_centers[i - 1])/svm_feat_scales[i - 1];
    }

    nodes[6].index = -1;
    nodes[6].value = 0;

    DoubleReal predict = svm_predict(isotope_filt_svm, nodes);

    delete(nodes);

    return (predict == 2.0)? true : false;
}

void FeatureFindingMetabo::loadIsotopeModel_()
{
    std::string model_filename = File::find("CHEMISTRY/MetaboliteIsoModel.svm");
    std::string scale_filename = File::find("CHEMISTRY/MetaboliteIsoModel.scale");

    isotope_filt_svm = svm_load_model(model_filename.c_str());

    std::ifstream ifs(scale_filename.c_str());

    std::string line;
    std::stringstream str_buf;
    std::istream_iterator<DoubleReal> eol;

    svm_feat_centers.clear();
    svm_feat_scales.clear();

    while (getline(ifs, line))
    {
        str_buf.clear();
        str_buf << line;
        std::istream_iterator<DoubleReal> istr_it(str_buf);

        while(istr_it != eol)
        {
            svm_feat_centers.push_back(*istr_it);
            ++istr_it;
            svm_feat_scales.push_back(*istr_it);
            ++istr_it;
        }
    }

    if (svm_feat_centers.size() != svm_feat_scales.size())
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Numbers of centers and scales from file " + scale_filename + " are different!", String(svm_feat_centers.size()) + " and " + String(svm_feat_scales.size()));
    }

    return ;
}



//DoubleReal FeatureFindingMetabo::scoreMZ_(const MassTrace& tr1, const MassTrace& tr2, Size iso_pos, Size charge)
//{
//    // DoubleReal mu(std::pow(1.0029316*iso_pos, -0.0002107)/charge);

//    // DoubleReal mu((1.003355*iso_pos)/charge);
//    // DoubleReal mu((DoubleReal)iso_pos/(DoubleReal)charge);

//    // DoubleReal mu((1.001527*(DoubleReal)iso_pos)/(DoubleReal)charge);
//    //DoubleReal mu((1.003355*(DoubleReal)iso_pos)/(DoubleReal)charge);

//    // std::cout << "checking mu " << mu << " " << iso_pos << " "  << charge << std::endl;



//    // DoubleReal err_ppm1((mz1/1000000)*mass_error_ppm_);
//    // DoubleReal err_ppm2((mz2/1000000)*mass_error_ppm_);

//    // DoubleReal sigma = 0.02017855;

//    DoubleReal sigma_mult(3.0);
//    DoubleReal C13_diff(1.003355);

//    DoubleReal mz1(tr1.getCentroidMZ());
//    DoubleReal mz2(tr2.getCentroidMZ());

//    DoubleReal centered_mz(std::fabs(mz2 - mz1) - (C13_diff*(DoubleReal)iso_pos)/(DoubleReal)charge);

//    // setup gaussian mixture model for valid isotopic m/z distances
//    //DoubleReal iso_mz_diff(1.002245);

//    // std::cout << "iso_pos: " << iso_pos << " charge: " << charge << " mz_diff" << centered_mz << std::endl;

//    DoubleReal mu1(-0.001210044);
//    DoubleReal sigma1(0.0009562774);
//    DoubleReal mu2(-0.008932877);
//    DoubleReal sigma2(0.0021871884);

//    DoubleReal mt_sigma1(tr1.getCentroidSD());
//    DoubleReal mt_sigma2(tr2.getCentroidSD());
//    DoubleReal mt_variances(mt_sigma1*mt_sigma1 + mt_sigma2*mt_sigma2);

//    //std::cout << "mt_variances: " << mt_variances << std::endl;


//    DoubleReal score_sigma1(std::sqrt(sigma1*sigma1 + mt_variances));
//    DoubleReal score_sigma2(std::sqrt(sigma2*sigma2 + mt_variances));

//    // std::cout << "score_sigma1: " << score_sigma1 << std::endl;
//    // std::cout << "score_sigma2: " << score_sigma2 << std::endl;


//    DoubleReal mz_score(0.0);


//    if ((centered_mz < mu1 + sigma_mult*score_sigma1) && (centered_mz > mu2 - sigma_mult*score_sigma2))
//    {
//        DoubleReal tmp_exponent1((centered_mz - mu1)/score_sigma1);
//        DoubleReal tmp_exponent2((centered_mz - mu2)/score_sigma2);

//        DoubleReal mz_score1(std::exp(-0.5*tmp_exponent1*tmp_exponent1));
//        DoubleReal mz_score2(std::exp(-0.5*tmp_exponent2*tmp_exponent2));

//        mz_score = (mz_score1 > mz_score2)? mz_score1 : mz_score2;

//    }




//    // std::cout << tr1.getLabel() << "_" << tr2.getLabel() <<  "mass: " << mz1 << " diffppm: " << " diffmz1: " << diff_mz-mu1 << " diffmz2: " << diff_mz-mu2 << " 3sigma1: " << sigma_mult*score_sigma1 << " 3sigma2: " << sigma_mult*score_sigma2 << " score: " << mz_score << std::endl ;


//    return mz_score;
//}

//DoubleReal FeatureFindingMetabo::scoreRT_(const MassTrace& tr1, const MassTrace& tr2)
//{
//    DoubleReal rt1(tr1.getCentroidRT());
//    DoubleReal rt2(tr2.getCentroidRT());

//    std::vector<DoubleReal> x, y;
//    std::map<DoubleReal, std::vector<DoubleReal> > overlap;

//    for (MassTrace::const_iterator mt_it = tr1.begin(); mt_it != tr1.end(); mt_it++)
//    {
//        overlap[mt_it->getRT()].push_back(mt_it->getIntensity());
//    }

//    for (MassTrace::const_iterator mt_it = tr2.begin(); mt_it != tr2.end(); mt_it++)
//    {
//        overlap[mt_it->getRT()].push_back(mt_it->getIntensity());
//    }

//    for (std::map<DoubleReal, std::vector<DoubleReal> >::const_iterator m_it = overlap.begin(); m_it != overlap.end(); ++m_it)
//    {
//        if (m_it->second.size() == 2)
//        {
//            x.push_back(m_it->second[0]);
//            y.push_back(m_it->second[1]);
//        }
//    }

//    DoubleReal diff_rt(std::fabs(rt2 - rt1));

//    //    DoubleReal sigma1(tr1.estimateFWHM(true)/2.3548);
//    //    DoubleReal sigma2(tr2.estimateFWHM(true)/2.3548);

//    DoubleReal sigma1(chrom_fwhm_/2.3548);
//    DoubleReal sigma2(chrom_fwhm_/2.3548);
//    DoubleReal sigma(std::sqrt(sigma1*sigma1 + sigma2*sigma2));

//    // std::cout << "==> RT: " << diff_rt << " " << std::exp(-0.5*((diff_rt)/sigma)*((diff_rt)/sigma)) << " olsCoeff " << computeCosineSim(x, y) << std::endl;

//    //    if (diff_rt > sigma)
//    //    {
//    //        return 0.0;
//    //    }

//    //    return std::exp(-0.5*((diff_rt)/sigma)*((diff_rt)/sigma));
//    return computeCosineSim(x, y);
//}


DoubleReal FeatureFindingMetabo::scoreMZ_(const MassTrace& tr1, const MassTrace& tr2, Size iso_pos, Size charge) {

    DoubleReal mz1(tr1.getCentroidMZ());
    DoubleReal mz2(tr2.getCentroidMZ());

    DoubleReal diff_mz(std::fabs(mz2 - mz1));
    DoubleReal avg_diff(1.001881);

    DoubleReal center((avg_diff*(DoubleReal)iso_pos)/(DoubleReal)charge);

    DoubleReal diff_sigma(0.004631371);

    DoubleReal sigma1(tr1.getCentroidSD());
    DoubleReal sigma2(tr2.getCentroidSD());

    // DoubleReal sigma = std::sqrt(sigma1*sigma1 + sigma2*sigma2);
    DoubleReal sigma(std::sqrt(diff_sigma*diff_sigma + sigma1*sigma1 + sigma2*sigma2));


    // std::cout << "mass: " << mz1 << " diffmz: " << diff_mz-center << " sigma: " << sigma << std::endl ;


    DoubleReal mz_score(0.0);

    if (std::fabs(diff_mz - center) < 3*sigma)
    {
        mz_score = std::exp(-0.5*((diff_mz - center)/sigma)*((diff_mz - center)/sigma));
    }

    return mz_score;
}

DoubleReal FeatureFindingMetabo::scoreRT_(const MassTrace& tr1, const MassTrace& tr2)
{
    DoubleReal rt1(tr1.getCentroidRT());
    DoubleReal rt2(tr2.getCentroidRT());


    DoubleReal diff_rt(std::fabs(rt2 - rt1));

//    DoubleReal sigma1(chrom_fwhm_/2.3548);
//    DoubleReal sigma(std::sqrt(2*sigma1*sigma1));

    DoubleReal sigma(1.0);

    if (diff_rt > 3*sigma)
    {
        return 0.0;
    }

    return std::exp(-0.5*((diff_rt)/sigma)*((diff_rt)/sigma));
}


DoubleReal FeatureFindingMetabo::computeCosineSim(const std::vector<DoubleReal>& x, const std::vector<DoubleReal>& y)
{
    if (x.size() != y.size())
    {
        return 0.0;
    }

    DoubleReal mixed_sum(0.0);
    DoubleReal x_squared_sum(0.0);
    DoubleReal y_squared_sum(0.0);


    for (Size i = 0; i < x.size(); ++i)
    {
        mixed_sum += x[i]*y[i];
        x_squared_sum += x[i]*x[i];
        y_squared_sum += y[i]*y[i];
    }

    DoubleReal denom(std::sqrt(x_squared_sum)*std::sqrt(y_squared_sum));

    return (denom > 0.0) ? mixed_sum/denom : 0.0;
}

DoubleReal FeatureFindingMetabo::computeOLSCoeff(const std::vector<DoubleReal>& x, const std::vector<DoubleReal>& y)
{
    if (x.size() != y.size())
    {
        return 0.0;
    }

    DoubleReal mixed_sum(0.0);
    DoubleReal x_squared_sum(0.0);

    for (Size i = 0; i < x.size(); ++i)
    {
        mixed_sum += x[i]*y[i];
        x_squared_sum += x[i]*x[i];
    }

    return (x_squared_sum > 0.0) ? mixed_sum/x_squared_sum : 0.0;
}

//DoubleReal FeatureFindingMetabo::scoreTraceSim_(MassTrace a, MassTrace b)
//{
//    std::map<DoubleReal, std::vector<DoubleReal> > intersect;

//    for (MassTrace::const_iterator c_it = a.begin(); c_it != a.end(); ++c_it)
//    {
//        intersect[c_it->getRT()].push_back(c_it->getIntensity());
//    }

//    for (MassTrace::const_iterator c_it = b.begin(); c_it != b.end(); ++c_it)
//    {
//        intersect[c_it->getRT()].push_back(c_it->getIntensity());
//    }

//    std::map<DoubleReal, std::vector<DoubleReal> >::const_iterator m_it = intersect.begin();

//    std::vector<DoubleReal> x, y;

//    for ( ; m_it != intersect.end(); ++m_it)
//    {
//        if (m_it->second.size() == 2)
//        {
//            x.push_back(m_it->second[0]);
//            y.push_back(m_it->second[1]);
//        }
//    }

//    if ( x.empty() || y.empty() )
//    {
//        return 0.0;
//    }

//    DoubleReal x_mean(0.0), y_mean(0.0);

//    x_mean = accumulate(x.begin(), x.end(), x_mean)/x.size();
//    y_mean = accumulate(y.begin(), y.end(), y_mean)/y.size();

//    DoubleReal counter(0.0), denom_x(0.0), denom_y(0.0);

//    for (Size i = 0; i < x.size(); ++i)
//    {
//        counter += (x[i] - x_mean)*(y[i] - y_mean);
//    }

//    for (Size i = 0; i < x.size(); ++i)
//    {
//        denom_x += (x[i] - x_mean)*(x[i] - x_mean);
//        denom_y += (y[i] - y_mean)*(y[i] - y_mean);
//    }

//    DoubleReal sim_score(counter/sqrt(denom_x*denom_y));

//    if (sim_score > 0.0)
//    {
//        return sim_score;
//    }

//    return 0.0;
//}



//DoubleReal FeatureFindingMetabo::scoreIntRatio_(DoubleReal int1, DoubleReal int2, Size iso_pos)
//{
//    DoubleReal int_ratio(0.0);

//    if (int2 > 0.0)
//    {
//        int_ratio = int2/int1;
//    }

//    DoubleReal mu(0.0), sigma(1.0);

//    switch (iso_pos)
//    {
//    case 1: mu = 0.4102466; sigma = 0.128907; break;
//    case 2: mu = 0.1034883; sigma = 0.04742052; break;
//    case 3: mu = 0.01910963; sigma = 0.01197569; break;
//    case 4: mu = 0.00286942; sigma  = 0.002266673; break;
//    default: mu = 0.0; sigma = 0.0003450974; break;
//    }

//    if (std::fabs(int_ratio - mu) > 2*sigma)
//    {
//        return 0.0;
//    }


//    DoubleReal int_score(std::exp(-0.5*((int_ratio - mu)/sigma)*((int_ratio - mu)/sigma)));

//    return int_score;
//}


void FeatureFindingMetabo::findLocalFeatures_(std::vector<MassTrace*>& candidates, std::vector<FeatureHypothesis>& output_hypos)
{
    FeatureHypothesis tmp_hypo;
    tmp_hypo.addMassTrace(*candidates[0]);
    tmp_hypo.setScore(0.0);

    output_hypos.push_back(tmp_hypo);

    for (Size charge = charge_lower_bound_; charge <= charge_upper_bound_; ++charge)
    {
        FeatureHypothesis fh_tmp;
        fh_tmp.setScore(0.0);

        fh_tmp.addMassTrace(*candidates[0]);

        //        DoubleReal mono_iso_rt(candidates[0]->getCentroidRT());
        //        DoubleReal mono_iso_mz(candidates[0]->getCentroidMZ());
        // DoubleReal mono_iso_int(candidates[0]->computePeakArea());

        Size last_iso_idx(0);

        // Size iso_pos_max(std::floor(charge * local_mz_range_));
        Size iso_pos_max(6);

        // std::cout << "isoposmax: " << iso_pos_max << std::endl;

        for (Size iso_pos = 1; iso_pos <= iso_pos_max; ++iso_pos) {
            DoubleReal best_so_far(0.0);
            Size best_idx(0);

            for (Size mt_idx = last_iso_idx + 1; mt_idx < candidates.size(); ++mt_idx)
            {
                // DoubleReal tmp_iso_rt(candidates[mt_idx]->getCentroidRT());
                // DoubleReal tmp_iso_mz(candidates[mt_idx]->getCentroidMZ());
                // DoubleReal tmp_iso_int(candidates[mt_idx]->computePeakArea());


                DoubleReal rt_score(scoreRT_(*candidates[0], *candidates[mt_idx]));

                DoubleReal mz_score(scoreMZ_(*candidates[0], *candidates[mt_idx], iso_pos, charge));

                // disable intensity scoring for now...
                DoubleReal int_score(1.0);

                DoubleReal total_pair_score(0.0);

                if (rt_score > 0.7 && mz_score > 0.0 && int_score > 0.0)
                {
                    total_pair_score = std::exp(std::log(rt_score) + log(mz_score) + log(int_score));
                }

                if (total_pair_score > best_so_far)
                {
                    best_so_far = total_pair_score;
                    best_idx = mt_idx;

                }
            } // end mt_idx

            if (best_so_far > 0.0)
            {
                fh_tmp.addMassTrace(*candidates[best_idx]);
                fh_tmp.setScore(fh_tmp.getScore() + best_so_far);
                fh_tmp.setCharge(charge);
                output_hypos.push_back(fh_tmp);
                last_iso_idx = best_idx;
            }
            else
            {
                break;
            }


        } // end for iso_pos

    } // end for charge

    return ;

} // end of findLocalFeatures_(...)


void FeatureFindingMetabo::run(std::vector<MassTrace>& input_mtraces, FeatureMap<>& output_featmap)
{
    // mass traces must be sorted by their centroid MZ
    std::sort(input_mtraces.begin(), input_mtraces.end(), CmpMassTraceByMZ());

    std::vector<FeatureHypothesis> feat_hypos;

    this->startProgress(0, input_mtraces.size(), "assembling mass traces to features");

    // initialize SVM model for isotope ratio filtering
    loadIsotopeModel_();

    if (input_mtraces.size() > 0) {
        for (Size i = 0; i < input_mtraces.size(); ++i)
        {
            this->setProgress(i);
            std::vector<MassTrace*> local_traces;

            DoubleReal ref_trace_mz(input_mtraces[i].getCentroidMZ());
            DoubleReal ref_trace_rt(input_mtraces[i].getCentroidRT());

            local_traces.push_back(&input_mtraces[i]);

            DoubleReal diff_mz(0.0), diff_rt(0.0);
            Size ext_idx(i + 1);

            while (diff_mz <= local_mz_range_ && ext_idx < input_mtraces.size())
            {
                // update diff_mz and diff_rt
                diff_mz = std::fabs(input_mtraces[ext_idx].getCentroidMZ() - ref_trace_mz);
                diff_rt = std::fabs(input_mtraces[ext_idx].getCentroidRT() - ref_trace_rt);

                if (diff_mz <= local_mz_range_ && diff_rt <= local_rt_range_)
                {
                    local_traces.push_back(&input_mtraces[ext_idx]);
                }

                ++ext_idx;
            }

            findLocalFeatures_(local_traces, feat_hypos);
        }
        this->endProgress();

        // sort feature candidates by their score
        std::sort(feat_hypos.begin(), feat_hypos.end(), CmpHypothesesByScore());

        std::map<String, bool> trace_excl_map;

        std::cout << "size of hypotheses: " << feat_hypos.size() << std::endl;

        for (Size hypo_idx = 0; hypo_idx < feat_hypos.size(); ++hypo_idx)
        {

            // std::cout << "score now: " <<  feat_hypos[hypo_idx].getScore() << std::endl;
            std::vector<String> labels(feat_hypos[hypo_idx].getLabels());

            bool trace_coll = false;

            for (Size lab_idx = 0; lab_idx < labels.size(); ++lab_idx)
            {
                if (trace_excl_map[labels[lab_idx]])
                {
                    trace_coll = true;
                }
            }

            //            if (feat_hypos[hypo_idx].getSize() > 1)
            //            {
            //                std::cout << "check for collision: " << trace_coll << " " << feat_hypos[hypo_idx].getLabel() << " " << isLegalIsotopePattern_(feat_hypos[hypo_idx]) << " " << feat_hypos[hypo_idx].getScore() << std::endl;
            //            }

            if (!trace_coll)
            {
                bool result = true;

                if (feat_hypos[hypo_idx].getSize() > 1)
                {
                    DoubleReal mono_int(feat_hypos[hypo_idx].getAllIntensities()[0]);


                    result = disable_isotope_filtering_ ? true : isLegalIsotopePattern_(feat_hypos[hypo_idx]);

                    // std::cout << "\nlegal iso? " << feat_hypos[hypo_idx].getLabel() << " score: " << feat_hypos[hypo_idx].getScore() << " " << result << std::endl << std::endl;
                }

                if (result) {
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

                    f.setWidth(feat_hypos[hypo_idx].getFWHM(true));
                    f.setCharge(feat_hypos[hypo_idx].getCharge());
                    f.setMetaValue(3,feat_hypos[hypo_idx].getLabel());
                    f.setConvexHulls(feat_hypos[hypo_idx].getConvexHulls());
                    f.setOverallQuality(feat_hypos[hypo_idx].getScore());

                    output_featmap.push_back(f);

                    for (Size lab_idx = 0; lab_idx < labels.size(); ++lab_idx)
                    {
                        trace_excl_map[labels[lab_idx]] = true;
                    }
                }

            }
        }
    }

    return ;
} // end of FeatureFindingMetabo::run



}
