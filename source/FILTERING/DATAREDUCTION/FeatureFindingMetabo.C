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
#include <sstream>

#include <boost/dynamic_bitset.hpp>

namespace OpenMS
{
    FeatureHypothesis::FeatureHypothesis()
        : iso_pattern_(),
        feat_score_()
    {

    }

    FeatureHypothesis::~FeatureHypothesis()
    {

    }


    FeatureHypothesis::FeatureHypothesis(const FeatureHypothesis& fh)
        : iso_pattern_(fh.iso_pattern_),
        feat_score_(fh.feat_score_)
    {

    }

    FeatureHypothesis& FeatureHypothesis::operator=(const FeatureHypothesis& rhs)
                                                   {
        if (this==&rhs) return *this;
        iso_pattern_ = rhs.iso_pattern_;
        feat_score_ = rhs.feat_score_;

        return *this;
    }



    void FeatureHypothesis::addMassTrace(MassTrace& mt_ptr)
    {
        iso_pattern_.push_back(&mt_ptr);

        return ;
    }

    DoubleReal FeatureHypothesis::computeFeatureIntensity()
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
        defaults_.setValue( "local_rt_range" , 3.0 , "RT range where to look for coeluting mass traces");
        defaults_.setValue( "local_mz_range" , 5.0 , "MZ range where to look for isotopic mass traces");
        defaults_.setValue( "mass_error_ppm", 20.0, "Allowed mass error deviation in ppm");

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

        if (input_mtraces.size() > 1) {
            for (Size i = 0; i < input_mtraces.size(); ++i)
            {
                // std::cout << input_mtraces[i].getCentroidMZ() << " " << input_mtraces[i].getSmoothedMaxRT() << std::endl;
                std::vector<MassTrace*> local_traces;

                DoubleReal ref_trace_mz(input_mtraces[i].getCentroidMZ());
                DoubleReal ref_trace_rt(input_mtraces[i].getSmoothedMaxRT());

                local_traces.push_back(&input_mtraces[i]);

                DoubleReal diff_mz(0.0), diff_rt(0.0);
                Size ext_idx(i + 1);

                while (diff_mz <= local_mz_range_ && ext_idx < input_mtraces.size())
                {
                    // update diff_mz and diff_rt
                    diff_mz = std::fabs(input_mtraces[ext_idx].getCentroidMZ() - ref_trace_mz);
                    diff_rt = std::fabs(input_mtraces[ext_idx].getSmoothedMaxRT() - ref_trace_rt);

                    if (diff_mz <= local_mz_range_ && diff_rt <= local_rt_range_)
                    {
                        local_traces.push_back(&input_mtraces[ext_idx]);
                    }

                    ++ext_idx;
                }


                // std::cout << "look at " << local_traces.size() << std::endl;

                findLocalFeatures_(local_traces, feat_hypos);
            }


            // sort feature candidates by their score
            std::sort(feat_hypos.begin(), feat_hypos.end(), CmpHypothesesByScore());

            std::map<String, bool> trace_excl_map;

            for (Size hypo_idx = 0; hypo_idx < feat_hypos.size(); ++hypo_idx)
            {
                std::cout << feat_hypos[hypo_idx].getLabel() << " " << feat_hypos[hypo_idx].getScore() << std::endl;

                std::vector<String> labels(feat_hypos[hypo_idx].getLabels());

                bool trace_coll = false;

                for (Size lab_idx = 0; lab_idx < labels.size(); ++lab_idx)
                {
                    if (trace_excl_map[labels[lab_idx]])
                    {
                        trace_coll = true;
                    }
                }

                if (!trace_coll) {
                    Feature f;
                    f.setRT(feat_hypos[hypo_idx].getCentroidRT());
                    f.setMZ(feat_hypos[hypo_idx].getCentroidMZ());
                    f.setIntensity(feat_hypos[hypo_idx].computeFeatureIntensity());
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

    }


    DoubleReal FeatureFindingMetabo::scoreMZ_(DoubleReal mz1, DoubleReal mz2, Size iso_pos, Size charge) {

        DoubleReal diff_mz(std::fabs(mz2 - mz1)/iso_pos);

        // use a piecewise linear function for first isotope trace
#if 0
        if (iso_pos == 1)
        {
            if (diff_mz < 1.000716)
            {
                return 0.0;
            }
            else if (diff_mz < 1.001746)
            {
                return 87.53355*diff_mz - 87.59622;
            }
            else if (diff_mz < 1.003294)
            {
                return 587.741*diff_mz - 588.6771;
            }
            else if (diff_mz < 1.003608)
            {
                return -3092.83*diff_mz + 3104.016;
            }
            else
            {
                return 0.0;
            }


        }
#endif

        DoubleReal mu(std::pow(1.0029316*iso_pos, -0.0002107)/charge);
        DoubleReal err_ppm((mz2/1000000)*mass_error_ppm_);
        DoubleReal sigma((2*err_ppm)/2.3548);



       // DoubleReal mu = (1.003355/charge)*iso_pos;

      //  DoubleReal sigma(0.01);

        DoubleReal mz_score(std::exp(-0.5*((diff_mz - mu)/sigma)*((diff_mz - mu)/sigma)));

        return mz_score;
    }

    DoubleReal FeatureFindingMetabo::scoreRT_(DoubleReal rt1, DoubleReal rt2) {

        DoubleReal diff_rt(std::fabs(rt2 - rt1));

        // DoubleReal mu_rt = (1.003355/charge)*iso_pos;

        DoubleReal sigma(2.0);

        DoubleReal rt_score(std::exp(-0.5*((diff_rt)/sigma)*((diff_rt)/sigma)));

        return rt_score;
    }


    void FeatureFindingMetabo::findLocalFeatures_(std::vector<MassTrace*>& candidates, std::vector<FeatureHypothesis>& output_hypos)
    {
        // check for singleton traces
//        if (candidates.size() == 1)
//        {
            FeatureHypothesis tmp_hypo;
            tmp_hypo.addMassTrace(*candidates[0]);
            tmp_hypo.setScore(0.0);

            output_hypos.push_back(tmp_hypo);

//            return ;
//        }


        bool singleton_trace = true;

        for (Size charge = 1; charge < 4; ++charge) {
            FeatureHypothesis fh_tmp;
            fh_tmp.setScore(0.0);

            fh_tmp.addMassTrace(*candidates[0]);

            DoubleReal mono_iso_rt(candidates[0]->getSmoothedMaxRT());
            DoubleReal mono_iso_mz(candidates[0]->getCentroidMZ());
            DoubleReal mono_iso_int(candidates[0]->computePeakArea());

            std::cout << candidates[0]->getLabel() << " with ";

            Size last_iso_idx(0);

            for (Size iso_pos = 1; iso_pos < 5; ++iso_pos) {
                DoubleReal best_so_far(0.0);
                Size best_idx(0);

                for (Size mt_idx = last_iso_idx + 1; mt_idx < candidates.size(); ++mt_idx)
                {
                    DoubleReal tmp_iso_rt(candidates[mt_idx]->getSmoothedMaxRT());
                    DoubleReal tmp_iso_mz(candidates[mt_idx]->getCentroidMZ());
                    DoubleReal tmp_iso_int(candidates[mt_idx]->computePeakArea());

                    DoubleReal rt_score(scoreRT_(mono_iso_rt, tmp_iso_rt));
                    DoubleReal mz_score(scoreMZ_(mono_iso_mz, tmp_iso_mz, iso_pos, charge));

                    DoubleReal total_pair_score(rt_score * mz_score);

                    if (total_pair_score > best_so_far)
                    {
                        best_so_far = total_pair_score;
                        best_idx = mt_idx;
                        //                        fh_tmp.addMassTrace(*candidates[mt_idx]);
                        //                        fh_tmp.setScore(fh_tmp.getScore() + total_pair_score);
                    }



                    std::cout << candidates[mt_idx]->getLabel() << " score: " << mz_score << " " << rt_score << std::endl;

                } // end mt_idx

                if (best_so_far > 0.01)
                {
                    fh_tmp.addMassTrace(*candidates[best_idx]);
                    fh_tmp.setScore(fh_tmp.getScore() + best_so_far);

                    output_hypos.push_back(fh_tmp);

                    last_iso_idx = best_idx;
                    singleton_trace = false;
                }
                else
                {
                    break;
                }


            } // end for iso_pos

//            if (fh_tmp.getSize() > 1)
//            {
//                output_hypos.push_back(fh_tmp);
//            }
            // std::cout << fh_tmp.getLabel() << " " << fh_tmp.getScore() << std::endl;
        } // end for charge


//        if (singleton_trace)
//        {
//            FeatureHypothesis fh_tmp;
//            fh_tmp.setScore(0.0);

//            fh_tmp.addMassTrace(*candidates[0]);

//            output_hypos.push_back(fh_tmp);
//        }


        return ;

    } // end of findLocalFeatures_(...)



}
