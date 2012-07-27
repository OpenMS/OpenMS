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


#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <sstream>

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


DoubleReal FeatureHypothesis::getMonoisotopicFeatureIntensity()
{
    if (iso_pattern_.empty())
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "FeatureHypothesis is empty, no traces contained!", String(iso_pattern_.size()));
    }

    return iso_pattern_[0]->computePeakArea();

}


DoubleReal FeatureHypothesis::getSummedFeatureIntensity()
{
    DoubleReal int_sum(0.0);

    for (Size i = 0; i < iso_pattern_.size(); ++i)
    {
        int_sum += iso_pattern_[i]->computePeakArea();
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
    defaults_.setValue("mass_error_ppm", 20.0, "Allowed mass error deviation in ppm");  // 20.0
    defaults_.setValue("chrom_sigma" , 3.0 , "Expected shift in RT (in seconds) between the apeces of mass traces that putatively belong to the same feature.");    // 3.0
    defaults_.setValue("report_summed_ints", "false", "Set to true for a feature intensity summed up over all traces rather than using monoisotopic trace intensity alone.", StringList::create("advanced"));
    defaults_.setValidStrings("report_summed_ints", StringList::create(("false,true")));

    defaultsToParam_();

    this->setLogType(CMD);
}


FeatureFindingMetabo::~FeatureFindingMetabo()
{

}


void FeatureFindingMetabo::run(std::vector<MassTrace>& input_mtraces, FeatureMap<>& output_featmap)
{
    // mass traces must be sorted by their centroid MZ
    std::sort(input_mtraces.begin(), input_mtraces.end(), CmpMassTraceByMZ());

    std::vector<FeatureHypothesis> feat_hypos;

    this->startProgress(0, input_mtraces.size(), "assembling mass traces to features");

    if (input_mtraces.size() > 1) {
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

        for (Size hypo_idx = 0; hypo_idx < feat_hypos.size(); ++hypo_idx)
        {
            std::vector<String> labels(feat_hypos[hypo_idx].getLabels());

            bool trace_coll = false;

            for (Size lab_idx = 0; lab_idx < labels.size(); ++lab_idx)
            {
                if (trace_excl_map[labels[lab_idx]])
                {
                    trace_coll = true;
                }
            }

            if (!trace_coll)
            {
                Feature f;
                f.setRT(feat_hypos[hypo_idx].getCentroidRT());
                f.setMZ(feat_hypos[hypo_idx].getCentroidMZ());

                if (report_summed_ints_)
                {
                    f.setIntensity(feat_hypos[hypo_idx].getSummedFeatureIntensity());
                }
                else
                {
                    f.setIntensity(feat_hypos[hypo_idx].getMonoisotopicFeatureIntensity());
                }

                f.setWidth(feat_hypos[hypo_idx].getFWHM());
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

    return ;
} // end of FeatureFindingMetabo::run






void FeatureFindingMetabo::updateMembers_()
{
    // delta_ = (Size)param_.getValue( "delta" );

    local_rt_range_ = (DoubleReal)param_.getValue("local_rt_range");
    local_mz_range_ = (DoubleReal)param_.getValue("local_mz_range");
    mass_error_ppm_ = (DoubleReal)param_.getValue("mass_error_ppm");
    chrom_sigma_ = (DoubleReal)param_.getValue("chrom_sigma");

    charge_lower_bound_ = (Size)param_.getValue("charge_lower_bound");
    charge_upper_bound_ = (Size)param_.getValue("charge_upper_bound");

    report_summed_ints_ = param_.getValue("report_summed_ints").toBool();
}


DoubleReal FeatureFindingMetabo::scoreMZ_(DoubleReal mz1, DoubleReal mz2, Size iso_pos, Size charge) {

    DoubleReal diff_mz(std::fabs(mz2 - mz1)/iso_pos);

    DoubleReal mu(std::pow(1.0029316*iso_pos, -0.0002107)/charge);
    DoubleReal err_ppm((mz2/1000000)*mass_error_ppm_);
    DoubleReal sigma((4*err_ppm)/2.3548);

    if (diff_mz - mu > 2*sigma)
    {
        return 0.0;
    }

    DoubleReal mz_score(std::exp(-0.5*((diff_mz - mu)/sigma)*((diff_mz - mu)/sigma)));

    return mz_score;
}

DoubleReal FeatureFindingMetabo::scoreRT_(DoubleReal rt1, DoubleReal rt2) {

    DoubleReal diff_rt(std::fabs(rt2 - rt1));

    DoubleReal sigma(chrom_sigma_/2.3548);

    if (diff_rt > 2*sigma)
    {
        return 0.0;
    }

    return std::exp(-0.5*((diff_rt)/sigma)*((diff_rt)/sigma));
}

DoubleReal FeatureFindingMetabo::scoreTraceSim_(MassTrace a, MassTrace b)
{
    std::map<DoubleReal, std::vector<DoubleReal> > intersect;

    for (MassTrace::const_iterator c_it = a.begin(); c_it != a.end(); ++c_it)
    {
        intersect[c_it->getRT()].push_back(c_it->getIntensity());
    }

    for (MassTrace::const_iterator c_it = b.begin(); c_it != b.end(); ++c_it)
    {
        intersect[c_it->getRT()].push_back(c_it->getIntensity());
    }

    std::map<DoubleReal, std::vector<DoubleReal> >::const_iterator m_it = intersect.begin();

    std::vector<DoubleReal> x, y;

    for ( ; m_it != intersect.end(); ++m_it)
    {
        if (m_it->second.size() == 2)
        {
            x.push_back(m_it->second[0]);
            y.push_back(m_it->second[1]);
        }
    }

    if ( x.empty() || y.empty() )
    {
        return 0.0;
    }

    DoubleReal x_mean(0.0), y_mean(0.0);

    x_mean = accumulate(x.begin(), x.end(), x_mean)/x.size();
    y_mean = accumulate(y.begin(), y.end(), y_mean)/y.size();

    DoubleReal counter(0.0), denom_x(0.0), denom_y(0.0);

    for (Size i = 0; i < x.size(); ++i)
    {
        counter += (x[i] - x_mean)*(y[i] - y_mean);
    }

    for (Size i = 0; i < x.size(); ++i)
    {
        denom_x += (x[i] - x_mean)*(x[i] - x_mean);
        denom_y += (y[i] - y_mean)*(y[i] - y_mean);
    }

    DoubleReal sim_score(counter/sqrt(denom_x*denom_y));

    if (sim_score > 0.0)
    {
        return sim_score;
    }

    return 0.0;
}



DoubleReal FeatureFindingMetabo::scoreIntRatio_(DoubleReal int1, DoubleReal int2, Size iso_pos)
{
    DoubleReal int_ratio(0.0);

    if (int2 > 0.0)
    {
        int_ratio = int2/int1;
    }

    DoubleReal mu(0.0), sigma(1.0);

    switch (iso_pos)
    {
    case 1: mu = 0.4102466; sigma = 0.128907; break;
    case 2: mu = 0.1034883; sigma = 0.04742052; break;
    case 3: mu = 0.01910963; sigma = 0.01197569; break;
    case 4: mu = 0.00286942; sigma  = 0.002266673; break;
    default: mu = 0.0; sigma = 0.0003450974; break;
    }

    if (std::fabs(int_ratio - mu) > 2*sigma)
    {
        return 0.0;
    }


    DoubleReal int_score(std::exp(-0.5*((int_ratio - mu)/sigma)*((int_ratio - mu)/sigma)));

    return int_score;
}


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

        DoubleReal mono_iso_rt(candidates[0]->getCentroidRT());
        DoubleReal mono_iso_mz(candidates[0]->getCentroidMZ());
        // DoubleReal mono_iso_int(candidates[0]->computePeakArea());

        Size last_iso_idx(0);

        Size iso_pos_max(std::floor(charge_upper_bound_ * local_mz_range_));

        for (Size iso_pos = 1; iso_pos <= iso_pos_max; ++iso_pos) {
            DoubleReal best_so_far(0.0);
            Size best_idx(0);

            for (Size mt_idx = last_iso_idx + 1; mt_idx < candidates.size(); ++mt_idx)
            {
                DoubleReal tmp_iso_rt(candidates[mt_idx]->getCentroidRT());
                DoubleReal tmp_iso_mz(candidates[mt_idx]->getCentroidMZ());
                // DoubleReal tmp_iso_int(candidates[mt_idx]->computePeakArea());


                DoubleReal rt_score(scoreRT_(mono_iso_rt, tmp_iso_rt));

                DoubleReal mz_score(scoreMZ_(mono_iso_mz, tmp_iso_mz, iso_pos, charge));

                // disable intensity scoring for now...
                DoubleReal int_score(1.0);

                DoubleReal total_pair_score(0.0);

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



}
