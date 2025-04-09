// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>

#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

namespace OpenSWATH_Test
{

  using namespace OpenMS;

  typedef OpenSwath::LightTransition TransitionType;
  //typedef ReactionMonitoringTransition TransitionType;
  typedef MRMTransitionGroup<MSChromatogram, TransitionType> MRMTransitionGroupType;

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
      MSChromatogram chrom;
      chrom.setNativeID(native_id);
      transition_group.addChromatogram(chrom, native_id );

      TransitionType tr;
      tr.library_intensity = 10000;
      tr.product_mz = 618.31;
      tr.fragment_charge = 1;
      tr.transition_name = native_id;
      transition_group.addTransition(tr, native_id );
    }
    {
      String native_id = "tr1";
      MSChromatogram chrom;
      chrom.setNativeID(native_id);
      transition_group.addChromatogram(chrom, native_id );

      TransitionType tr;
      tr.library_intensity = 1;
      tr.product_mz = 628.435;
      tr.fragment_charge = 1;
      tr.transition_name = native_id;
      transition_group.addTransition(tr, native_id );
    }
    {
      String native_id = "tr5";
      MSChromatogram chrom;
      chrom.setNativeID(native_id);
      transition_group.addChromatogram(chrom, native_id );

      TransitionType tr;
      tr.library_intensity = 2000;
      tr.product_mz = 628.435;
      tr.fragment_charge = 1;
      tr.transition_name = native_id;
      transition_group.addTransition(tr, native_id );
    }

    return transition_group;
  }

}

