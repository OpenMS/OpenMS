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
			defaults_.setValue("rt_stdev_low",0.22);
			defaults_.setValue("rt_stdev_high",0.65);
			defaults_.setValue("mz_stdev",0.025);
			defaults_.setValue("mz_pair_dist",4.0);
			defaults_.setValue("rt_pair_dist",0.3);
			
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
			features_
			(source.features_),
			pairs_(source.pairs_), 
			best_pairs_()
		{	
		}

		const PairMatcher::PairVectorType& PairMatcher::run()
		{
			double rt_stdev_low = param_.getValue("rt_stdev_low");
			double rt_stdev_high = param_.getValue("rt_stdev_high");
			double mz_stdev = param_.getValue("mz_stdev");
			double mz_pair_dist = param_.getValue("mz_pair_dist");
			double rt_pair_dist = param_.getValue("rt_pair_dist");

			double rt_min = rt_pair_dist - 2.0*rt_stdev_low;
			double rt_max = rt_pair_dist + 2.0*rt_stdev_high;
			double mz_diff = 2.0*mz_stdev;

			pairs_.clear();

			// calculate area of map
			features_.updateRanges();

			// fill tree
			QuadTreeType tree(DRange<2>(features_.getMin()[0], features_.getMax()[0], features_.getMin()[1], features_.getMax()[1]));
			for (Size i=0; i<features_.size(); ++i)
			{
				try
				{
					tree.insert(features_[i].getPosition(), &features_[i] );
				}
				catch(Exception::IllegalTreeOperation e)
				{
					cout << "Warning: Multiple identical feature positions in given feature map!" << endl;
				}
			}

			// set id for each feature
			int id = -1;
			for (FeatureMapType::Iterator it = features_.begin(); it != features_.end(); ++it)
			{
				it->setMetaValue(ID,++id);
			}

			// check each feature
			DRange<2> local;
			for (FeatureMapType::const_iterator it=features_.begin(); it!=features_.end(); ++it)
			{
				//cout << "Testing feature " << it->getPosition()[0] << " " << it->getPosition()[1] << endl;
				// set up local area to search for feature partner
				int charge = it->getCharge();
				double mz_opt = mz_pair_dist/charge;
				local.setMinX(it->getPosition()[0]-rt_max);
				local.setMaxX(it->getPosition()[0]-rt_min);
				local.setMinY(it->getPosition()[1]+mz_opt-mz_diff);
				local.setMaxY(it->getPosition()[1]+mz_opt+mz_diff);
				
				//cout << local << endl;
				
				for (QuadTreeType::Iterator check=tree.begin(local); check!=tree.end(); ++check)
				{
					//cout << "  Testing point" << endl;
					if ( check->second->getCharge() == charge)
					{
						//cout << "    charge ok" << endl;
						// calculate score
						double diff[2];
						diff[MZ] = fabs( it->getPosition()[MZ] - check->second->getPosition()[MZ] );
						diff[RT] = it->getPosition()[RT] - check->second->getPosition()[RT];

						double score =  PValue_(diff[MZ], mz_opt, mz_stdev, mz_stdev)
													* PValue_(diff[RT],rt_pair_dist, rt_stdev_low, rt_stdev_high)
													* check->second->getOverallQuality()
													* it->getOverallQuality();

						pairs_.push_back(PairType( *it, *check->second, score));
					}
				}
			}

			return pairs_;
		}

		const PairMatcher::PairVectorType& PairMatcher::getBestPairs()
		{
			typedef std::list< PairType* > Feature2PairList;
			typedef vector<Feature2PairList> ListVector;
			ListVector feature2pair(features_.size());

			std::sort(pairs_.begin(), pairs_.end(), PairMatcher::Comparator());

			for (PairVectorType::iterator it=pairs_.begin(); it!=pairs_.end(); ++it)
			{
				it->getFirst().setMetaValue(LOW_QUALITY,0);
				int id1 = it->getFirst().getMetaValue(ID);
				int id2 = it->getSecond().getMetaValue(ID);

				feature2pair[id1].push_back( &(*it) );
				feature2pair[id2].push_back( &(*it) );
			}

			for (PairVectorType::iterator pair=pairs_.begin(); pair!=pairs_.end(); ++pair)
			{
				// Pair still in set
				if (static_cast<int>(pair->getFirst().getMetaValue(LOW_QUALITY))==0)
				{
					int id1 = pair->getFirst().getMetaValue(ID);
					int id2 = pair->getSecond().getMetaValue(ID);
					// 'Remove' (by setting the flag) all additional pairs the features belongs to
					for (Feature2PairList::const_iterator it=feature2pair[id1].begin();
							it!=feature2pair[id1].end(); ++it)
						(*it)->getFirst().setMetaValue(LOW_QUALITY,1);

					for (Feature2PairList::const_iterator it=feature2pair[id2].begin();
							it!=feature2pair[id2].end(); ++it)
						(*it)->getFirst().setMetaValue(LOW_QUALITY,1);

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
			for (Size i=0; i<pairs.size(); ++i)
			{
				DPosition<2> diff = pairs[i].getFirst().getPosition()-pairs[i].getSecond().getPosition();
				out << setiosflags(ios::fixed) << setprecision(2)
						<< pairs[i].getQuality() << "\t" << pairs[i].getFirst().getPosition()[0] << "\t" 
						<< pairs[i].getFirst().getPosition()[1] << "\t" << pairs[i].getFirst().getIntensity() << "\t" 
						<< pairs[i].getFirst().getOverallQuality() << "\t" << pairs[i].getSecond().getPosition()[0] << "\t"
						<< pairs[i].getSecond().getPosition()[1] << "\t" << pairs[i].getSecond().getIntensity() << "\t"
						<< pairs[i].getSecond().getOverallQuality() << "\t" << pairs[i].getFirst().getIntensity()/pairs[i].getSecond().getIntensity() << "\t"
						<< pairs[i].getFirst().getCharge() << "\t" << diff[RT] << "\t"
						<< diff[MZ] << endl;
			}
		}
	}
