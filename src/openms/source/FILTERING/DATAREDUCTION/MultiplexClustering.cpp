// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/FILTERING/DATAREDUCTION/PeakPattern.h>
#include <OpenMS/FILTERING/DATAREDUCTION/FilterResult.h>
#include <OpenMS/FILTERING/DATAREDUCTION/FilterResultRaw.h>
#include <OpenMS/FILTERING/DATAREDUCTION/FilterResultPeak.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplinePackage.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MultiplexFiltering.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MultiplexClustering.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>
#include <OpenMS/COMPARISON/CLUSTERING/Cluster.h>
#include <OpenMS/COMPARISON/CLUSTERING/LocalClustering.h>

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

namespace OpenMS
{

	MultiplexClustering::MultiplexClustering(MSExperiment<Peak1D> exp_profile, MSExperiment<Peak1D> exp_picked, std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries, double rt_typical, double rt_minimum, bool debug)
    : rt_typical_(rt_typical), rt_minimum_(rt_minimum), debug_(debug)
	{
        if (exp_picked.size() != boundaries.size())
        {
            throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Centroided data and the corresponding list of peak boundaries do not contain same number of spectra.");
        }
        
        // ranges of the experiment
        double mz_min = exp_profile.getMin().getX();
        double mz_max = exp_profile.getMax().getX();
        double rt_min = exp_profile.getMin().getY();
        double rt_max = exp_profile.getMax().getY();
        
        // generate hash grid spacing
        PeakWidthEstimator estimator(exp_picked, boundaries, 40);
        for (double mz = mz_min; mz < mz_max; mz = mz + estimator.getPeakWidth(mz) / 10)
        {
            // We assume that the jitter of the peak centres are less than 1/10 of the peak width.
            // The factor 1/10 ensures that two neighbouring peaks at the same RT cannot be in the same cluster. 
            grid_spacing_mz_.push_back(mz);
        }
        grid_spacing_mz_.push_back(mz_max);
        
        for (double rt = rt_min; rt < rt_max; rt = rt + rt_typical)
        {
            grid_spacing_rt_.push_back(rt);
        }
        grid_spacing_rt_.push_back(rt_max);
        
        // determine RT scaling
        std::vector<double> mz;
        MSExperiment<Peak1D>::Iterator it_rt;
        for (it_rt = exp_picked.begin(); it_rt < exp_picked.end(); ++it_rt)
        {
            MSSpectrum<Peak1D>::Iterator it_mz;
            for (it_mz = it_rt->begin(); it_mz < it_rt->end(); ++it_mz)
            {
                mz.push_back(it_mz->getMZ());
             }
        }
        std::sort(mz.begin(), mz.end());
        // RT scaling = peak width at the median of the m/z distribuation / RT threshold
        rt_scaling_ = estimator.getPeakWidth(mz[(int) mz.size() / 2]) / rt_typical_;
        
	}
    
    std::vector<std::map<int,Cluster> > MultiplexClustering::cluster(std::vector<FilterResult> filter_results)
    {
        std::vector<std::map<int,Cluster> > cluster_results;
        
        // loop over patterns i.e. cluster each of the corresponding filter results
        for (unsigned i = 0; i < filter_results.size(); ++i)
        {
            LocalClustering clustering(filter_results[i].getMz(), filter_results[i].getRt(), grid_spacing_mz_, grid_spacing_rt_, rt_scaling_);
            clustering.cluster();
            //clustering.extendClustersY();
            //clustering.removeSmallClustersY(rt_minimum_);
            cluster_results.push_back(clustering.getResults());
            
            // debug output
            vector<DebugPoint> debug_clustered;
            if (debug_)
            {
                std::map<int,Cluster> cluster_result = clustering.getResults();
                FilterResult filter_result = filter_results[i];
                
                int cluster_id = 0;
                for(std::map<int,Cluster>::iterator it = cluster_result.begin(); it != cluster_result.end(); ++it) {
                    std::vector<int> points = (it->second).getPoints();
                    for (std::vector<int>::iterator it2 = points.begin(); it2 != points.end(); ++it2)
                    {
                        DebugPoint data_point;
                        data_point.rt = filter_result.getRt(*it2);
                        data_point.mz = filter_result.getMz(*it2);
                        data_point.cluster = cluster_id;
                        debug_clustered.push_back(data_point);
                    }
                    ++cluster_id;
                }
            }
            std::cout << "    debug size = " << debug_clustered.size() << "\n";
            writeDebug(debug_clustered, i);

        }
        
        return cluster_results;
    }
    
    MultiplexClustering::PeakWidthEstimator::PeakWidthEstimator(MSExperiment<Peak1D> exp_picked, std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries, int quantiles)
    {
        if (exp_picked.size() != boundaries.size())
        {
            throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Centroided data and the corresponding list of peak boundaries do not contain same number of spectra.");
        }

        std::vector<double> mz;
        std::vector<double> peak_width;

        MSExperiment<Peak1D>::Iterator it_rt;
        vector<vector<PeakPickerHiRes::PeakBoundary> >::const_iterator it_rt_boundaries;
        for (it_rt = exp_picked.begin(), it_rt_boundaries = boundaries.begin();
            it_rt < exp_picked.end() && it_rt_boundaries < boundaries.end();
            ++it_rt, ++it_rt_boundaries)
        {
            MSSpectrum<Peak1D>::Iterator it_mz;
            vector<PeakPickerHiRes::PeakBoundary>::const_iterator it_mz_boundary;
            for (it_mz = it_rt->begin(), it_mz_boundary = it_rt_boundaries->begin();
                 it_mz < it_rt->end(), it_mz_boundary < it_rt_boundaries->end();
                 ++it_mz, ++it_mz_boundary)
            {
                mz.push_back(it_mz->getMZ());
                peak_width.push_back((*it_mz_boundary).mz_max - (*it_mz_boundary).mz_min);
            }
        }
        std::sort(mz.begin(), mz.end());
        std::sort(peak_width.begin(), peak_width.end());
        
        std::vector<double> mz_quantiles;
        std::vector<double> peak_width_quantiles;
        for (int i = 1; i < quantiles; ++i)
        {
            mz_quantiles.push_back(mz[(int) mz.size() * i / quantiles]);
            peak_width_quantiles.push_back(peak_width[(int) peak_width.size() * i / quantiles]);
        }
        
        mz_min_ = mz_quantiles.front();
        mz_max_ = mz_quantiles.back();
        
        spline_ = new CubicSpline2d(mz_quantiles, peak_width_quantiles);
    }
    
    double MultiplexClustering::PeakWidthEstimator::getPeakWidth(double mz) const
    {
        if (mz < mz_min_)
        {
            return (*spline_).eval(mz_min_);
        }
        else if (mz > mz_max_)
        {
            return (*spline_).eval(mz_max_);
        }
        else
        {
            return (*spline_).eval(mz);
        }
    }
    
    void MultiplexClustering::writeDebug(vector<DebugPoint> points, int pattern) const
    {
        // fill consensus map
        ConsensusMap map;
        for (std::vector<DebugPoint>::const_iterator it = points.begin(); it != points.end(); ++it)
        {
            ConsensusFeature consensus;
            consensus.setRT((*it).rt);
            consensus.setMZ((*it).mz);
            consensus.setIntensity((*it).cluster);
            consensus.setCharge(1);    // dummy
            std::cout << "colour = " << getColour((*it).cluster) << "\n";
            consensus.setMetaValue("color", getColour((*it).cluster));
            consensus.setMetaValue("Cluster ID", (*it).cluster);
            
            FeatureHandle feature;
            feature.setRT((*it).rt);
            feature.setMZ((*it).mz);
            feature.setUniqueId((*it).cluster);
            consensus.insert(feature);
            
            map.getFileDescriptions()[0].size++;
            map.push_back(consensus);
        }
        
        ConsensusMap::FileDescription & desc = map.getFileDescriptions()[0];
        desc.filename = "debug";
        desc.label = "Cluster";
        
        map.sortByPosition();
        map.applyMemberFunction(&UniqueIdInterface::setUniqueId);
        map.setExperimentType("multiplex");

        // write consensus file
        ConsensusXMLFile file;
        String file_name = "debug_clustered_";
        file_name = file_name + pattern + ".consensusXML";
        file.store(file_name, map);
    }
    
    String MultiplexClustering::getColour(int c) const
    {
        // 15 HTML colors
        static const String colours[] =
        {
          "#00FFFF", "#000000", "#0000FF", "#FF00FF", "#008000",
          "#808080", "#00FF00", "#800000", "#000080", "#808000",
          "#800080", "#FF0000", "#C0C0C0", "#008080", "#FFFF00",
        };
        
        return colours[c % 15];
    }
    
}
