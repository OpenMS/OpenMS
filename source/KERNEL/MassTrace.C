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
      fwhm_num_scans_()
{
}

MassTrace::MassTrace(const std::list<PeakType>& tmp_lst)
{
    trace_peaks_.clear();

    for (std::list<PeakType>::const_iterator l_it = tmp_lst.begin(); l_it != tmp_lst.end(); ++l_it)
    {
        trace_peaks_.push_back((*l_it));
    }
}

MassTrace::MassTrace(const std::vector<PeakType>& tmp_vec)
{
    trace_peaks_ = tmp_vec;
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
      fwhm_num_scans_(mt.fwhm_num_scans_)
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
    fwhm_num_scans_ = rhs.fwhm_num_scans_;

    return *this;
}


DoubleReal MassTrace::computePeakArea() {
    DoubleReal peak_area(0.0);

    if (trace_peaks_.empty())
        return peak_area;

    for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
    {
        peak_area += (*l_it).getIntensity();
    }

    return peak_area;
}


Size MassTrace::findMaxByIntPeak(bool use_smoothed_ints = false)
{
    if (use_smoothed_ints && smoothed_intensities_.empty())
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace was not smoothed before! Aborting...", String(smoothed_intensities_.size()));
    }

    if (trace_peaks_.empty())
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace appears to be empty! Aborting...", String(trace_peaks_.size()));
    }

    DoubleReal max_int;
    Size max_idx(0);

    if (use_smoothed_ints)
    {
        max_int = smoothed_intensities_[0];
    }
    else
    {
        max_int = trace_peaks_.begin()->getIntensity();
    }

    for (Size i = 0; i < trace_peaks_.size(); ++i)
    {
        DoubleReal act_int = use_smoothed_ints ? smoothed_intensities_[i] : trace_peaks_[i].getIntensity();

        if (act_int > max_int)
        {
            max_int = act_int;
            max_idx = i;
        }

    }

    return max_idx;
}

DoubleReal MassTrace::estimateFWHM(bool use_smoothed_ints = false)
{
    Size max_idx(this->findMaxByIntPeak(use_smoothed_ints ? true : false));

    std::vector<DoubleReal> tmp_ints;

    if (use_smoothed_ints)
    {
        tmp_ints = smoothed_intensities_;
    }
    else
    {
        for (Size vec_idx = 0; vec_idx < trace_peaks_.size(); ++vec_idx)
        {
            tmp_ints.push_back(trace_peaks_[vec_idx].getIntensity());
        }
    }

    DoubleReal half_max_int(tmp_ints[max_idx]/2.0);

    Size left_border(max_idx), right_border(max_idx);

    while (left_border > 0 && tmp_ints[left_border] >= half_max_int)
    {
        --left_border;
    }

    while (right_border + 1 < tmp_ints.size() && tmp_ints[right_border] >= half_max_int)
    {
        ++right_border;
    }

    // side effect: record number of peaks/scans that span the fwhm of the mass trace; useful for smoothing techniques (window size)

    fwhm_num_scans_ = right_border - left_border + 1;

    return std::fabs(trace_peaks_[right_border].getRT() - trace_peaks_[left_border].getRT());
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

        // bisection
        Size left_bound(chrom_maxes[i] + 1);
        Size right_bound(chrom_maxes[i + 1] - 1);


        while ((left_bound + 1) < right_bound)
        {
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

void MassTrace::updateWeightedMeanRT()
{
    if (trace_peaks_.empty())
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace is empty... centroid RT undefined!", String(trace_peaks_.size()));
    }

    DoubleReal trace_area(this->computePeakArea());

    if (trace_area < std::numeric_limits<DoubleReal>::epsilon())
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Peak area equals to zero... impossible to compute weights!", String(trace_peaks_.size()));
    }

    DoubleReal wmean_rt(0.0);

    for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
    {
        wmean_rt += ((*l_it).getIntensity() * (*l_it).getRT())/trace_area;
    }

    centroid_rt_ = wmean_rt;
}


void MassTrace::updateMedianRT()
{
    if (trace_peaks_.empty())
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace is empty... centroid RT undefined!", String(trace_peaks_.size()));
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


void MassTrace::updateMedianMZ()
{
    if (trace_peaks_.empty())
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace is empty... centroid MZ undefined!", String(trace_peaks_.size()));
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


void MassTrace::updateMeanMZ()
{
    if (trace_peaks_.empty())
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace is empty... centroid MZ undefined!", String(trace_peaks_.size()));
    }

    Size trace_size = trace_peaks_.size();

    DoubleReal sum_mz(0.0);

    for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
    {
        sum_mz += (*l_it).getMZ();
    }

    centroid_mz_ = sum_mz/trace_size;

    return ;
}

void MassTrace::updateWeightedMeanMZ()
{
    if (trace_peaks_.empty())
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "MassTrace is empty... centroid MZ undefined!", String(trace_peaks_.size()));
    }

    DoubleReal weighted_sum(0.0);
    DoubleReal total_weight(0.0);

    for (MassTrace::const_iterator l_it = trace_peaks_.begin(); l_it != trace_peaks_.end(); ++l_it)
    {
        DoubleReal w_i = (*l_it).getIntensity();
        total_weight += w_i;
        weighted_sum += w_i * (*l_it).getMZ();
    }

    if (total_weight < std::numeric_limits<DoubleReal>::epsilon())
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "All weights were equal to zero! Empty trace? Aborting...", String(total_weight));
    }

    centroid_mz_ = weighted_sum/total_weight;
}





} // end of MassTrace.C
