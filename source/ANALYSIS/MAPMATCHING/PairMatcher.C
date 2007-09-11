// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/PairMatcher.h>

#include <iomanip>

using namespace std;

namespace OpenMS
{
		const double PairMatcher::sqrt2_ = sqrt(2);


		PairMatcher::PairMatcher(FeatureMapType& features)
			: FactoryProduct(PairMatcher::getProductName()), 
				features_(features), 
				pairs_(), 
				best_pairs_()
		{
			defaults_.setValue("rt_stdev_low", 0.22, "standard deviation below optimal retention time distance", false);
			defaults_.setValue("rt_stdev_high", 0.65, "standard deviation above optimal retention time distance", false);
			defaults_.setValue("mz_stdev", 0.025, "standard deviation from optimal m/z distance\n", false);
			defaults_.setValue("mz_pair_dist", 4.0, "optimal pair distance in m/z [Th] for features with charge +1 (adapted to +2, +3, .. by division through charge)", false);
			defaults_.setValue("rt_pair_dist", 0.3, "optimal pair distance in RT [sec]", false);
			
			defaultsToParam_();
		}

		PairMatcher::~PairMatcher()
		{	
		}

    PairMatcher& PairMatcher::operator = (const PairMatcher& source)
    {
			FactoryProduct::operator = (source);
			features_ = source.features_;
			return *this;
		}

		PairMatcher::PairMatcher(const PairMatcher& source)
		: FactoryProduct(source), 
			features_(source.features_),
			pairs_(source.pairs_), 
			best_pairs_()
		{	
		}

		const PairMatcher::PairVectorType& PairMatcher::run()
		{
			//RT settings
			double rt_pair_dist = param_.getValue("rt_pair_dist");
			double rt_stdev_low = param_.getValue("rt_stdev_low");
			double rt_stdev_high = param_.getValue("rt_stdev_high");
			//MZ settings
			double mz_stdev = param_.getValue("mz_stdev");
			double mz_pair_dist = param_.getValue("mz_pair_dist");

			//cout << "MZ Window: " << mz_pair_dist << " +/- " << mz_stdev << endl;						
			//cout << "RT Window: " << rt_pair_dist << " + " << rt_stdev_high << " - " << rt_stdev_low << endl;
			
			pairs_.clear();

			// sort features by RT (and MZ) to speed up searching afterwards
			features_.sortByPosition();
			
			// set id for each feature
			int id = -1;
			for (FeatureMapType::Iterator it = features_.begin(); it != features_.end(); ++it)
			{
				it->setMetaValue(11,++id);
				//cout <<"Feature " << id << ": " << it->getRT() << " / " << it->getMZ() << endl;
			}

			// check each feature
			for (FeatureMapType::const_iterator it=features_.begin(); it!=features_.end(); ++it)
			{
				//cout << "*****************************************************************" << endl;
				//cout << "Testing feature: " << it->getRT() << " / " << it->getMZ() << endl;
				//cout << "RT range: " << it->getRT()+rt_pair_dist - 2.0*rt_stdev_low << " - " << it->getRT()+rt_pair_dist + 2.0*rt_stdev_high << endl;
				//cout << "MZ range: " << it->getMZ()+mz_pair_dist/it->getCharge()-2.0*mz_stdev << " - " << it->getMZ()+mz_pair_dist/it->getCharge()+2.0*mz_stdev << endl;
				//cout << "*****************************************************************" << endl;
				FeatureMapType::const_iterator range = lower_bound(features_.begin(),features_.end(),it->getRT()+rt_pair_dist - 2.0*rt_stdev_low, Feature::NthPositionLess<0>());
				while (range!=features_.end() && range->getRT() <= it->getRT()+rt_pair_dist + 2.0*rt_stdev_high)
				{
					//cout << "Checking: " << range->getRT() << " / " << range->getMZ() << endl;
					if (range->getCharge() == it->getCharge()
						&& range->getMZ() >=it->getMZ()+mz_pair_dist/it->getCharge()-2.0*mz_stdev 
						&& range->getMZ() <=it->getMZ()+mz_pair_dist/it->getCharge()+2.0*mz_stdev)
					{
						
						DoubleReal score =  PValue_(range->getMZ() - it->getMZ(), mz_pair_dist/it->getCharge(), mz_stdev, mz_stdev)
															* PValue_(range->getRT() - it->getRT(), rt_pair_dist, rt_stdev_low, rt_stdev_high)
															* range->getOverallQuality()
															* it->getOverallQuality();
						//cout << "HIT: " << score << endl;
						pairs_.push_back(PairType( *it, *range, score));
					}
					++range;
				}
			}
			return pairs_;
		}

		const PairMatcher::PairVectorType& PairMatcher::getBestPairs()
		{
			best_pairs_.clear();
			
			typedef std::list< PairType* > Feature2PairList;
			typedef vector<Feature2PairList> ListVector;
			ListVector feature2pair(features_.size());

			std::sort(pairs_.begin(), pairs_.end(), PairMatcher::Comparator());

			for (PairVectorType::iterator it=pairs_.begin(); it!=pairs_.end(); ++it)
			{
				it->getFirst().setMetaValue(12,0);
				int id1 = it->getFirst().getMetaValue(11);
				int id2 = it->getSecond().getMetaValue(11);

				feature2pair[id1].push_back( &(*it) );
				feature2pair[id2].push_back( &(*it) );
			}

			for (PairVectorType::iterator pair=pairs_.begin(); pair!=pairs_.end(); ++pair)
			{
				// Pair still in set
				if (static_cast<int>(pair->getFirst().getMetaValue(12))==0)
				{
					int id1 = pair->getFirst().getMetaValue(11);
					int id2 = pair->getSecond().getMetaValue(11);
					// 'Remove' (by setting the flag) all additional pairs the features belongs to
					for (Feature2PairList::const_iterator it=feature2pair[id1].begin();
							it!=feature2pair[id1].end(); ++it)
						(*it)->getFirst().setMetaValue(12,1);

					for (Feature2PairList::const_iterator it=feature2pair[id2].begin();
							it!=feature2pair[id2].end(); ++it)
						(*it)->getFirst().setMetaValue(12,1);

					// Add pair into vector of best pairs
					best_pairs_.push_back(*pair);
				}
			}
			return best_pairs_;
		}

		void PairMatcher::printInfo(std::ostream& out, const PairVectorType& pairs)
		{
			out << "Found the following " << pairs.size() << " pairs:\n"
					<< "Quality\tFirst[RT]\tFirst[MZ]\tFirst[Int]\tFirst[Corr]"
					<< "\tSecond[RT]\tSecond[MZ]\tSecond[Int]\tSecond[Corr]"
					<< "\tRatio\tCharge\tDiff[RT]\tDiff[MZ]\n";
			for (UInt i=0; i<pairs.size(); ++i)
			{
				DPosition<2> diff = pairs[i].getSecond().getPosition()-pairs[i].getFirst().getPosition();
				out << setiosflags(ios::fixed) << setprecision(2)
						<< pairs[i].getQuality() << "\t" << pairs[i].getFirst().getRT() << "\t" 
						<< pairs[i].getFirst().getMZ() << "\t" << pairs[i].getFirst().getIntensity() << "\t" 
						<< pairs[i].getFirst().getOverallQuality() << "\t" << pairs[i].getSecond().getRT() << "\t"
						<< pairs[i].getSecond().getMZ() << "\t" << pairs[i].getSecond().getIntensity() << "\t"
						<< pairs[i].getSecond().getOverallQuality() << "\t" << pairs[i].getFirst().getIntensity()/pairs[i].getSecond().getIntensity() << "\t"
						<< pairs[i].getFirst().getCharge() << "\t" << diff[RawDataPoint2D::RT] << "\t"
						<< diff[RawDataPoint2D::MZ] << endl;
			}
		}
	}
