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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_OPENSWATH_TEST_H
#define OPENMS_OPENSWATH_TEST_H

#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>

#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

namespace OpenSWATH_Test
{

  using namespace OpenMS;

  // Above are the definitions using real chromatograms.
  // Below are the definitions using spectra objects to store the data.
  //
  // This is necessary as long as peak picking cannot be done on chromatograms
  // natively.  the classes needed at the moment are SavitzkyGolayFilter,
  // GaussFilter, PeakPickerHiRes -- SignalToNoiseEstimatorMedian seems to work
  // already.
#if 0
  typedef MSChromatogram<ChromatogramPeak> RichPeakChromatogram;
  typedef MRMTransitionGroup<MSChromatogram, ChromatogramPeak> MRMTransitionGroupType;
#else
  typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram;
  typedef OpenSwath::LightTransition TransitionType;
  //typedef ReactionMonitoringTransition TransitionType;
  typedef MRMTransitionGroup<RichPeakChromatogram, TransitionType> MRMTransitionGroupType;
#endif

  MRMFeature createMockFeature()
  {

    MRMFeature feature;
    feature.setRT(3120);
    feature.setIntensity(static_cast<float>(973.122));

    {
      Feature f;
      static const double arr1[] = 
      {
        3103.13,
        3106.56,
        3109.98,
        3113.41,
        3116.84,
        3120.26,
        3123.69,
        3127.11,
        3130.54,
        3133.97,
        3137.4 
      };
      std::vector<double> mz (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
      static const double arr2[] = 
      {
        5.97544,
        4.27492,
        3.33018,
        4.08597,
        5.50307,
        5.24327,
        8.40812,
        2.8342 ,
        6.94379,
        7.69957,
        4.08597
      };
      std::vector<double> intensity (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

      ConvexHull2D::PointArrayType hull_points;
      for (Size i = 0; i < mz.size(); i++)
      {
        DPosition<2> p(mz[i],intensity[i]);
        hull_points.push_back(p);
      }
      ConvexHull2D hull;
      hull.setHullPoints(hull_points);
      f.getConvexHulls().push_back(hull);
      f.setIntensity(static_cast<float>(58.38450));
      feature.addFeature(f, "tr3");
    }

    {
      Feature f;
      static const double arr1[] = 
      {
        3103.13,     
        3106.56,
        3109.98,
        3113.41,
        3116.84,
        3120.26,
        3123.69,
        3127.11,
        3130.54,
        3133.97,
        3137.4 
      };
      std::vector<double> mz (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
      static const double arr2[] =
      {
        15.8951,     
        41.5446,
        76.0746,
        109.069,
        111.904,
        169.792,
        121.044,
        63.0137,
        44.615 ,
        21.4927,
       7.93576 
      };
      std::vector<double> intensity (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

      ConvexHull2D::PointArrayType hull_points;
      for (Size i = 0; i < mz.size(); i++)
      {
        DPosition<2> p(mz[i],intensity[i]);
        hull_points.push_back(p);
      }
      ConvexHull2D hull;
      hull.setHullPoints(hull_points);
      f.setIntensity(782.38073);
      f.getConvexHulls().push_back(hull);
      feature.addFeature(f, "tr1");
    }

    {  
      Feature f;
      static const double arr1[] = 
      {
        3103.13,
        3106.56,
        3109.98,
        3113.41,
        3116.84,
        3120.26,
        3123.69,
        3127.11,
        3130.54,
        3133.97,
        3137.4 
      };
      std::vector<double> mz (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
      static const double arr2[] = 
      {
        5.73925,
        6.7076 ,
        2.85782,
        5.0307 ,
        8.95135,
        14.4544,
        20.9731,
        24.3033,
        20.6897,
        13.7459,
        8.90411
      };
      std::vector<double> intensity (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

      ConvexHull2D::PointArrayType hull_points;
      for (Size i = 0; i < mz.size(); i++)
      {
        DPosition<2> p(mz[i],intensity[i]);
        hull_points.push_back(p);
      }
      ConvexHull2D hull;
      hull.setHullPoints(hull_points);
      f.setIntensity(static_cast<float>(58.38450));
      f.getConvexHulls().push_back(hull);
      feature.addFeature(f, "tr5");
    }

    return feature;
  }

  MRMTransitionGroupType createMockTransitionGroup()
  {

    MRMTransitionGroupType transition_group;
    {
      String native_id = "tr3";
      RichPeakChromatogram chrom;
      chrom.setNativeID(native_id);
      transition_group.addChromatogram(chrom, native_id );

      TransitionType tr;
      tr.library_intensity = 10000;
      tr.product_mz = 618.31;
      tr.product_charge = 1;
      tr.transition_name = native_id;
      transition_group.addTransition(tr, native_id );
    }
    {
      String native_id = "tr1";
      RichPeakChromatogram chrom;
      chrom.setNativeID(native_id);
      transition_group.addChromatogram(chrom, native_id );

      TransitionType tr;
      tr.library_intensity = 1;
      tr.product_mz = 628.435;
      tr.product_charge = 1;
      tr.transition_name = native_id;
      transition_group.addTransition(tr, native_id );
    }
    {
      String native_id = "tr5";
      RichPeakChromatogram chrom;
      chrom.setNativeID(native_id);
      transition_group.addChromatogram(chrom, native_id );

      TransitionType tr;
      tr.library_intensity = 2000;
      tr.product_mz = 628.435;
      tr.product_charge = 1;
      tr.transition_name = native_id;
      transition_group.addTransition(tr, native_id );
    }

    return transition_group;
  }

}

#endif
