// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: ClusterFactory.C,v 1.3 2006/03/29 12:30:29 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer:  $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/CLUSTERING/ClusterFactory.h>

//Filter Functions
#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ComplementFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/GoodDiffFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/IntensityBalanceFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/IntensityDistBins.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NeutralLossDiffFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeDiffFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/KellerQuality.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ParentFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PeakDensityFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PeakDiffBins.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PeakPosBins.h>
#include <OpenMS/FILTERING/TRANSFORMERS/TradSeqQuality.h>

//Mower hopefully improve the Quality of the Spectra
#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <OpenMS/FILTERING/TRANSFORMERS/MarkerMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NeutralLossMarker.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ComplementMarker.h>
#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeMarker.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SqrtMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ParentPeakMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Scaler.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/BernNorm.h>

//Similarity Functors
#include <OpenMS/COMPARISON/SPECTRA/CompareFunctor.h>
//work on BinnedReps
#include <OpenMS/COMPARISON/SPECTRA/BinnedRepSpectrumContrastAngle.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedRepSharedPeakCount.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedRepSumAgreeingIntensities.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedRepMutualInformation.h>
//work on 1D Spectra
#include <OpenMS/COMPARISON/SPECTRA/SpectrumCheapDPCorr.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumPrecursorComparator.h>
#include <OpenMS/COMPARISON/SPECTRA/SequestCompareFunctor.h>

//Clustering Algorithms
#include <OpenMS/COMPARISON/CLUSTERING/ClusterFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/LinkageCluster.h>

//Analysis Functors
#include <OpenMS/COMPARISON/CLUSTERING/AnalysisFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterCompareFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/DistanceAnalyzer.h>

#include <OpenMS/CONCEPT/Exception.h>

using namespace std;

namespace OpenMS
{
 
  ClusterFactory* ClusterFactory::instancep_ = 0;
  
  //singleton
  ClusterFactory* ClusterFactory::instance()
  {
    if (!instancep_) 
    {
      instancep_ = new ClusterFactory();
      instancep_->init();
    }
    return instancep_;
  }
  
  ClusterFactory::ClusterFactory()
    :inventory_()
  {
  }

  ClusterFactory::~ClusterFactory()
  {
  }

  FactoryProduct* ClusterFactory::create(String name) const
  {
    if (inventory_.find(name) != inventory_.end())
    {
      return (*(inventory_.find(name)->second))();
    }
    else throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"unregistered FactoryProduct",name.c_str());
  }

  void ClusterFactory::registerfp(String name,FactoryProduct*(*fp)() )
  {
    if ( inventory_.find(name) != inventory_.end() )
    {
      throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"already registered",name);
    }
    inventory_.insert(make_pair(name,fp));
  }

  // all FactoryProducts should be registered here
  void ClusterFactory::init()
  {
    // Filters
	
		// TODO register childs automatically
		//registerfp(FilterFunctor::getName(), &FilterFunctor::create);
	
		//@TODO registering of the names
    registerfp(ComplementFilter::getName(), &ComplementFilter::create);
    registerfp(GoodDiffFilter::getName(),&GoodDiffFilter::create);
    registerfp(IntensityBalanceFilter::getName(),&IntensityBalanceFilter::create);
    registerfp(IntensityDistBins::getName(),&IntensityDistBins::create);
    registerfp(NeutralLossDiffFilter::getName(),&NeutralLossDiffFilter::create);
    registerfp(IsotopeDiffFilter::getName(),&IsotopeDiffFilter::create);
    registerfp(KellerQuality::getName(),&KellerQuality::create);
    registerfp(ParentFilter::getName(),&ParentFilter::create);
    registerfp(TICFilter::getName(),&TICFilter::create);
    registerfp(PeakDensityFilter::getName(),&PeakDensityFilter::create);
    registerfp(PeakDiffBins::getName(),&PeakDiffBins::create);
    registerfp(PeakPosBins::getName(),&PeakPosBins::create);
    registerfp(TradSeqQuality::getName(),&TradSeqQuality::create);
    
    // PeakMarker
		//registerfp(MarkerFunctor::getName(), &MarkerFunctor::create);
    registerfp(NeutralLossMarker::getName(),&NeutralLossMarker::create);
    registerfp(IsotopeMarker::getName(),&IsotopeMarker::create);
    registerfp(ComplementMarker::getName(),&ComplementMarker::create);
    
    
    // Preprocessing
		//registerfp(PreprocessingFunctor::getName(), &PreprocessingFunctor::create);
    registerfp(ThresholdMower::getName(),&ThresholdMower::create);
    registerfp(WindowMower::getName(),&WindowMower::create);
    registerfp(Scaler::getName(),&Scaler::create);
    registerfp(NLargest::getName(),&NLargest::create);
    registerfp(BernNorm::getName(),&BernNorm::create);
    registerfp(MarkerMower::getName(),&MarkerMower::create);
    registerfp(SqrtMower::getName(),&SqrtMower::create);
    registerfp(Normalizer::getName(),&Normalizer::create);
    registerfp(ParentPeakMower::getName(),&ParentPeakMower::create);
    
    // Comparison functions
		//registerfp(CompareFunctor::getName(), &CompareFunctor::create);
    registerfp(BinnedRepSpectrumContrastAngle::getName(),&BinnedRepSpectrumContrastAngle::create);
    registerfp(BinnedRepMutualInformation::getName(),&BinnedRepMutualInformation::create);
    registerfp(BinnedRepSumAgreeingIntensities::getName(),&BinnedRepSumAgreeingIntensities::create);
    registerfp(BinnedRepSharedPeakCount::getName(),&BinnedRepSharedPeakCount::create);
    registerfp(SpectrumCheapDPCorr::getName(),&SpectrumCheapDPCorr::create);
    registerfp(SpectrumPrecursorComparator::getName(),&SpectrumPrecursorComparator::create);
    registerfp(SequestCompareFunctor::getName(),&SequestCompareFunctor::create);
    
    // Clustering
		//registerfp(ClusterFunctor::getName(), &ClusterFunctor::create);
    registerfp(LinkageCluster::getName(),&LinkageCluster::create);
    
    // Analysis
		//registerfp(AnalysisFunctor::getName(), &AnalysisFunctor::create);
    registerfp(ClusterCompareFunctor::getName(),&ClusterCompareFunctor::create);
    registerfp(DistanceAnalyzer::getName(),&DistanceAnalyzer::create);
  }

  /** 
  code needs to be changed for new types <br>
  */
  vector<String> ClusterFactory::catalogue(String type) const
  {
	
    vector<String> result;
    for ( map<String, FactoryProduct*(*)()>::const_iterator cmit = inventory_.begin(); cmit != inventory_.end(); ++cmit )
    {
      FactoryProduct* tmp = create(cmit->first);
      if ( type == "CompareFunctor") 
      {
        tmp = dynamic_cast<CompareFunctor*>(tmp);
      }
      else if ( type == "ClusterFunctor" )
      {
        tmp = dynamic_cast<ClusterFunctor*>(tmp);
      }
      else if ( type == "FilterFunctor" )
      {
        tmp = dynamic_cast<FilterFunctor*>(tmp);
      }
      else if ( type == "PreprocessingFunctor" )
      {
        tmp = dynamic_cast<PreprocessingFunctor*>(tmp);
      }
      else if ( type == "AnalysisFunctor" )
      {
        tmp = dynamic_cast<AnalysisFunctor*>(tmp);
      }
      else if ( type == "PeakMarker" )
      {
        tmp = dynamic_cast<PeakMarker*>(tmp);
      }
      else if ( type == "FactoryProduct")
      {
        tmp = dynamic_cast<FactoryProduct*>(tmp);
      }
      else
      {
        throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"unknown type, check ClusterFactory::catalogue",type.c_str());
      }
      if ( tmp ) 
      {
        result.push_back(cmit->first);
      }
    }
    return result;
  }

  // copying of FactoryProducts is discouraged <br>
  // (see Gotcha #76 in "C++ Gotchas" by Stephen C. Dewhurst ), <br>
  // use this instead <br>
  FactoryProduct* ClusterFactory::duplicate(const FactoryProduct* tmplate ) const
  {
    FactoryProduct* copy = create(tmplate->getName());
    if ( dynamic_cast<CompareFunctor*>(copy) )
    {
      CompareFunctor* cfp = dynamic_cast<CompareFunctor*>(copy);
      const CompareFunctor* cftmp = dynamic_cast<const CompareFunctor*>(tmplate);
      *cfp = *cftmp;
    }
    else if ( dynamic_cast<FilterFunctor*>(copy) )
    {
      FilterFunctor* ffp = dynamic_cast<FilterFunctor*>(copy);
      const FilterFunctor* fftmp = dynamic_cast<const FilterFunctor*>(tmplate);
      *ffp = *fftmp;
    }
    else if ( dynamic_cast<PreprocessingFunctor*>(copy) )
    {
      PreprocessingFunctor* mfp = dynamic_cast<PreprocessingFunctor*>(copy);
      const PreprocessingFunctor* mftmp = dynamic_cast<const PreprocessingFunctor*>(tmplate);
      *mfp = *mftmp;
    }
    else if ( dynamic_cast<AnalysisFunctor*>(copy) )
    {
      AnalysisFunctor* afp = dynamic_cast<AnalysisFunctor*>(copy);
      const AnalysisFunctor* aftmp = dynamic_cast<const AnalysisFunctor*>(tmplate);
      *afp = *aftmp;
    }
    else if ( dynamic_cast<PeakMarker*>(copy) )
    {
      PeakMarker* pmp = dynamic_cast<PeakMarker*>(copy);
      const PeakMarker* pmtmp = dynamic_cast<const PeakMarker*>(tmplate);
      *pmp = *pmtmp;
    }
    return copy;
  }
}
