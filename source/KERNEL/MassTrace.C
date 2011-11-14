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

#include <OpenMS/KERNEL/MassTrace.h>

#include <boost/dynamic_bitset.hpp>

namespace OpenMS
{
    MassTrace::MassTrace()
        : trace_peaks_(),
        centroid_mz_(),
        centroid_rt_(),
        label_(),
        smoothed_intensities_(),
        rough_fwhm_points_(),
        prev_counter_(),
        prev_denom_()
    {
    }

    MassTrace::~MassTrace()
    {
    }

    MassTrace::MassTrace(const MassTrace& mt)
        : trace_peaks_(mt.trace_peaks_),
        centroid_mz_(mt.centroid_mz_),
        centroid_rt_(mt.centroid_rt_),
        label_(mt.label_),
        smoothed_intensities_(mt.smoothed_intensities_),
        rough_fwhm_points_(mt.rough_fwhm_points_),
        prev_counter_(mt.prev_counter_),
        prev_denom_(mt.prev_denom_)
    {
    }

    MassTrace& MassTrace::operator= (const MassTrace& rhs)
                                    {
        if (this==&rhs) return *this;

        trace_peaks_ = rhs.trace_peaks_;
        centroid_mz_ = rhs.centroid_mz_;
        centroid_rt_ = rhs.centroid_rt_;
        label_ = rhs.label_;
        smoothed_intensities_ = rhs.smoothed_intensities_;
        rough_fwhm_points_ = rhs.rough_fwhm_points_;
        prev_counter_ = rhs.prev_counter_;
        prev_denom_ = rhs.prev_denom_;

        return *this;
    }

    void MassTrace::prependPeak(PeakType p)
    {
        trace_peaks_.push_front(p);
        updateMedianMZ_();
        // updateIterativeWeightedMeanMZ_(p);
        return ;
    }

    void MassTrace::appendPeak(PeakType p)
    {
        trace_peaks_.push_back(p);
        updateMedianMZ_();
        // updateIterativeWeightedMeanMZ_(p);
        return ;
    }

    DoubleReal MassTrace::computeWeightedMeanMZ()
    {
        DoubleReal weighted_sum(0.0);
        DoubleReal total_weight(0.0);

        for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
        {
            DoubleReal w_i = (*l_it).getIntensity();
            total_weight += w_i;
            weighted_sum += w_i * (*l_it).getMZ();
        }

        if (total_weight == 0.0)
            return 0.0;

        return weighted_sum/total_weight;
    }

    DoubleReal MassTrace::computePeakArea() {
        DoubleReal peak_area(0.0);

        if (trace_peaks_.size() == 0)
            return peak_area;

        for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
        {
            peak_area += (*l_it).getIntensity();
        }

        return peak_area;
    }


    Size MassTrace::findSmoothedMaxIdx()
    {
        if (trace_peaks_.size() != smoothed_intensities_.size())
        {
            throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace was not smoothed before! Aborting...", String(smoothed_intensities_.size()));
        }


        DoubleReal max_int(smoothed_intensities_[0]);
        DoubleReal max_idx(0);

        for (Size i = 0; i < smoothed_intensities_.size(); ++i)
        {
            if (smoothed_intensities_[i] > max_int)
            {
                max_int = smoothed_intensities_[i];
                max_idx = i;
            }
        }

        return max_idx;
    }

    DoubleReal MassTrace::findMaxPeakRT()
    {
        DoubleReal max_rt((*(trace_peaks_.begin())).getRT());

        DoubleReal max_int(0.0);

        for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
        {
            if ((*l_it).getIntensity() > max_int) {
                max_int = (*l_it).getIntensity();
                max_rt = (*l_it).getRT();
            }
        }

        return max_rt;
    }


    DoubleReal MassTrace::estimateFWHM()
    {
        Size max_idx(this->findSmoothedMaxIdx());
        DoubleReal half_max_int(smoothed_intensities_[max_idx]/2.0);

        Size left_border(max_idx), right_border(max_idx);

        while (left_border > 0)
        {
            --left_border;

            if (smoothed_intensities_[left_border] <= half_max_int)
            {
                break ;
            }
        }

        while (right_border+1 < smoothed_intensities_.size())
        {
            ++right_border;

            if (smoothed_intensities_[right_border] <= half_max_int)
            {
                break ;
            }
        }
        MassTrace::const_iterator l_it = trace_peaks_.begin();
        MassTrace::const_iterator r_it = trace_peaks_.begin();
        std::advance(l_it, left_border);

        if (right_border < smoothed_intensities_.size())
        {
            std::advance(r_it, right_border);
        }
        else
        {
            std::advance(r_it, right_border - 1);
        }



        return std::fabs(r_it->getRT() - l_it->getRT());
    }





    void MassTrace::findLocalExtrema(const Size& num_neighboring_peaks, std::vector<Size>& chrom_maxes, std::vector<Size>& chrom_mins)
    {
        Size mt_length(smoothed_intensities_.size());

        if (mt_length != trace_peaks_.size())
        {
            throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace was not smoothed before! Aborting...", String(smoothed_intensities_.size()));
        }

        // Extract RTs from the chromatogram and store them into into vectors for index access

        //  Store indices along with smoothed_ints to keep track of the peak order
        std::multimap<DoubleReal, Size> intensity_indices;
        boost::dynamic_bitset<> used_idx(mt_length);

        for (Size i = 0; i < mt_length; ++i)
        {
            intensity_indices.insert(std::make_pair(smoothed_intensities_[i], i));
        }


        for(std::multimap<DoubleReal, Size>::const_iterator c_it = intensity_indices.begin(); c_it != intensity_indices.end(); ++c_it)
        {
            DoubleReal ref_int = c_it->first;
            Size ref_idx = c_it->second;

            if (!(used_idx[ref_idx]))
            {
                bool real_max = true;

                Size start_idx(0);

                if(ref_idx > num_neighboring_peaks)
                {
                    start_idx = ref_idx - num_neighboring_peaks;
                }

                Size end_idx = ref_idx + num_neighboring_peaks;

                if(end_idx > mt_length)
                {
                    end_idx = mt_length;
                }

                for(Size j = start_idx; j < end_idx; ++j)
                {
                    if(used_idx[j])
                    {
                        real_max = false;
                        break;
                    }

                    if(j == ref_idx)
                    {
                        continue;
                    }

                    if(smoothed_intensities_[j] > ref_int)
                    {
                        real_max = false;
                    }
                }

                if(real_max)
                {
                    chrom_maxes.push_back(ref_idx);

                    for(Size j = start_idx; j < end_idx; ++j)
                    {
                        used_idx[j] = true;
                    }
                }

            }
        }

        std::sort(chrom_maxes.begin(), chrom_maxes.end());


        for (Size i = 0; i < chrom_maxes.size() - 1; ++i)
        {
#if 0
            // linear search for minimum
            DoubleReal max_int((smoothed_intensities_[chrom_maxes[i]] > smoothed_intensities_[chrom_maxes[i + 1]]) ? smoothed_intensities_[chrom_maxes[i]] : smoothed_intensities_[chrom_maxes[i + 1]]);
            Size min_idx(chrom_maxes[i] + 1);

            Size linear_steps(0);

            for (Size j = chrom_maxes[i] + 1; j < chrom_maxes[i + 1]; ++j)
            {
                if (smoothed_intensities_[j] < max_int)
                {
                    max_int = smoothed_intensities_[j];
                    min_idx = j;
                }
                ++linear_steps;
            }
#endif

            // binary search
            Size left_bound(chrom_maxes[i] + 1);
            Size right_bound(chrom_maxes[i + 1] - 1);

            // Size binary_steps(0);

            while ((left_bound + 1) < right_bound)
            {
                //                std::cout << left_bound << "___" << right_bound << std::endl;
                DoubleReal mid_dist((right_bound - left_bound)/2.0);

                Size mid_element_idx(left_bound + std::floor(mid_dist));

                DoubleReal mid_element_int = smoothed_intensities_[mid_element_idx];

                if (mid_element_int <= smoothed_intensities_[mid_element_idx + 1])
                {
                    right_bound = mid_element_idx;
                }
                else // or to the right...
                {
                    left_bound = mid_element_idx;
                }

                // ++binary_steps;
            }

            Size min_rt((smoothed_intensities_[left_bound] < smoothed_intensities_[right_bound]) ? left_bound : right_bound);
            chrom_mins.push_back(min_rt);
        }


        return ;
    }

    ConvexHull2D MassTrace::getConvexhull() const
    {
        ConvexHull2D::PointArrayType hull_points(trace_peaks_.size());

        Size i = 0;
        for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
        {
            hull_points[i][0] = (*l_it).getRT();
            hull_points[i][1] = (*l_it).getMZ();
            ++i;
        }

        ConvexHull2D hull;
        hull.addPoints(hull_points);

        return hull;
    }

    void MassTrace::updateWeightedMeanRT_()
    {
        DoubleReal trace_area(this->computePeakArea());

        DoubleReal wmean_rt(0.0);

        for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
        {
            wmean_rt += ((*l_it).getIntensity() * (*l_it).getRT())/trace_area;
        }

        centroid_rt_ = wmean_rt;
    }


    void MassTrace::updateMedianRT_()
    {
        if (trace_peaks_.size() == 0)
        {
            return ;
        }

        if (trace_peaks_.size() == 1)
        {
            centroid_rt_ = (*(trace_peaks_.begin())).getRT();

            return ;
        }

        // copy mz values to temp vec
        std::vector<DoubleReal> temp_rt;

        for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
        {
            temp_rt.push_back((*l_it).getRT());
        }

        std::sort(temp_rt.begin(), temp_rt.end());

        Size temp_mz_size = temp_rt.size();

        if ((temp_mz_size % 2) == 0)
        {
            centroid_rt_ = (temp_rt[std::floor(temp_mz_size/2.0) - 1] +  temp_rt[std::floor(temp_mz_size/2.0)])/2;
        }
        else
        {
            centroid_rt_ = temp_rt[std::floor(temp_mz_size/2.0)];
        }


        return ;
    }


    void MassTrace::updateMedianMZ_()
    {
        if (trace_peaks_.size() == 0)
        {
            return ;
        }

        if (trace_peaks_.size() == 1)
        {
            centroid_mz_ = (*(trace_peaks_.begin())).getMZ();

            return ;
        }

        // copy mz values to temp vec
        std::vector<DoubleReal> temp_mz;

        for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
        {
            temp_mz.push_back((*l_it).getMZ());
        }

        std::sort(temp_mz.begin(), temp_mz.end());

        Size temp_mz_size = temp_mz.size();

        if ((temp_mz_size % 2) == 0)
        {
            centroid_mz_ = (temp_mz[std::floor(temp_mz_size/2.0) - 1] +  temp_mz[std::floor(temp_mz_size/2.0)])/2;
        }
        else
        {
            centroid_mz_ = temp_mz[std::floor(temp_mz_size/2.0)];
        }

        return ;
    }


    void MassTrace::updateMeanMZ_() {
        Size trace_size = trace_peaks_.size();

        if (trace_size == 0)
            return ;

        DoubleReal sum_mz(0.0);

        for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
        {
            sum_mz += (*l_it).getMZ();
        }

        centroid_mz_ = sum_mz/trace_size;

        return ;
    }

    void MassTrace::updateIterativeWeightedMeanMZ_(const PeakType& added_peak) {
        Size trace_size = trace_peaks_.size();

        if (trace_size == 1) {
            centroid_mz_ = (*(trace_peaks_.begin())).getMZ();
            prev_counter_ = (*(trace_peaks_.begin())).getIntensity() * (*(trace_peaks_.begin())).getMZ();
            prev_denom_ = (*(trace_peaks_.begin())).getIntensity();

            return ;
        }

        DoubleReal new_weight = added_peak.getIntensity();
        DoubleReal new_mz = added_peak.getMZ();

        DoubleReal counter_tmp = (1 + (new_weight*new_mz)/prev_counter_);
        DoubleReal denom_tmp = (1 + (new_weight)/prev_denom_);
        centroid_mz_ = centroid_mz_ * (counter_tmp/denom_tmp);
        prev_counter_ = prev_counter_ * counter_tmp;
        prev_denom_ = prev_denom_ * denom_tmp;

        return ;
    }



} // end of MassTrace.C
