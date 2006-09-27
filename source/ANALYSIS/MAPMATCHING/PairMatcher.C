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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/PairMatcher.h>
#include <qstring.h>

using namespace std;

namespace OpenMS
{

		const double PairMatcher::sqrt2 = sqrt(2);

		PairMatcher::PairMatcher(FeatureMapType& features)
		: FactoryProduct(), features_(features), pairs_(), best_pairs_()
		{
			name_ = PairMatcher::getName();
			defaults_.setValue("rt_stdev_low",0.22);
			defaults_.setValue("rt_stdev_high",0.65);
			defaults_.setValue("mz_stdev",0.025);
			defaults_.setValue("mz_pair_dist",4.0);
			defaults_.setValue("rt_pair_dist",0.3);
			param_ = defaults_;
		}

		PairMatcher::~PairMatcher(){}

    PairMatcher& PairMatcher::operator = (const PairMatcher& source)
    {
			FactoryProduct::operator = (source);
			features_ = source.features_;
			return *this;
		}

		PairMatcher::PairMatcher(const PairMatcher& source)
		: FactoryProduct(source), features_(source.features_),
			pairs_(source.pairs_), best_pairs_()
		{}

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
			float	min_x = numeric_limits<float>::max();
			float max_x = -1 * numeric_limits<float>::max();
			float	min_y = numeric_limits<float>::max();
			float max_y = -1 * numeric_limits<float>::max();
			for (FeatureMapType::ConstIterator it = features_.begin(); it != features_.end(); ++it)
			{
				if (it->getPosition()[0] < min_x) min_x = it->getPosition()[0];
				if (it->getPosition()[0] > max_x) max_x = it->getPosition()[0];
				if (it->getPosition()[1] < min_y) min_y = it->getPosition()[1];
				if (it->getPosition()[1] > max_y) max_y = it->getPosition()[1];
			}

			// fill tree
			QuadTreeType tree(DRange<2>(min_x, max_x, min_y, max_y));
			for (Size i=0; i<features_.size(); ++i)
			{
				try{
					tree.insert(features_[i].getPosition(), &features_[i] );
				}catch(Exception::IllegalTreeOperation e)
				{
					cout << "Warning: Multiple identical feature positions in given feature map!" << endl;
				}
			}

			// clear convex hulls and meta value, set id for each feature
			int id = -1;
			for (FeatureMapType::Iterator it = features_.begin(); it != features_.end(); ++it)
			{
				it->getConvexHulls().clear();
				it->setMetaValue(3,std::string(""));
				it->setMetaValue(ID,++id);
			}

			// check each feature
			DRange<2> local;
			for (FeatureMapType::const_iterator it=features_.begin(); it!=features_.end(); ++it)
			{
				// set up local area to search for feature partner
				int charge = it->getCharge();
				double mz_opt = mz_pair_dist/charge;
				local.setMinX(it->getPosition()[0]-rt_max);
				local.setMaxX(it->getPosition()[0]-rt_min);
				local.setMinY(it->getPosition()[1]+mz_opt-mz_diff);
				local.setMaxY(it->getPosition()[1]+mz_opt+mz_diff);

				for (QuadTreeType::Iterator check=tree.begin(local); check!=tree.end(); ++check)
				{
					if ( check->second->getCharge() == charge)
					{
						// calculate score
						double diff[2];
						diff[MZ] = fabs( it->getPosition()[MZ] - check->second->getPosition()[MZ] );
						diff[RT] = it->getPosition()[RT] - check->second->getPosition()[RT];

						double score =  p_value_(diff[MZ], mz_opt, mz_stdev, mz_stdev)
													* p_value_(diff[RT],rt_pair_dist, rt_stdev_low, rt_stdev_high)
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

		void PairMatcher::fillFeatureMap(FeatureMapType& map, const PairVectorType& pairs)
		{
			map.clear();
			// put pairs in feature map, make sure every feature is unique
			// to avoid problems with quad tree in TOPPView
			typedef std::map< DPosition<2> , DFeature<2> > UniqueFeatureMap;
			UniqueFeatureMap tmp;

			for (Size i=0; i<pairs.size(); ++i)
			{
				tmp.insert( make_pair( pairs[i].getFirst().getPosition() , pairs[i].getFirst()) );
				tmp.insert( make_pair( pairs[i].getSecond().getPosition(), pairs[i].getSecond()) );

				UniqueFeatureMap::iterator first = tmp.find( pairs[i].getFirst().getPosition() );

				DFeature<2>::ConvexHullType hull;
				const double fak = 0.3;
				hull.push_back( pairs[i].getFirst().getPosition() - DPosition<2>(fak,0) );
				hull.push_back( pairs[i].getFirst().getPosition() - DPosition<2>(0,fak) );
				hull.push_back( pairs[i].getFirst().getPosition() + DPosition<2>(fak,0) );
				hull.push_back( pairs[i].getFirst().getPosition() + DPosition<2>(0,fak) );
				hull.push_back( pairs[i].getSecond().getPosition() + DPosition<2>(fak,0) );
				hull.push_back( pairs[i].getSecond().getPosition() + DPosition<2>(0,fak) );
				hull.push_back( pairs[i].getSecond().getPosition() - DPosition<2>(fak,0) );
				hull.push_back( pairs[i].getSecond().getPosition() - DPosition<2>(0,fak) );
				first->second.getConvexHulls().push_back( hull );

				DPosition<2> diff
					= pairs[i].getFirst().getPosition()-pairs[i].getSecond().getPosition();
				double ratio
					= pairs[i].getFirst().getIntensity()/pairs[i].getSecond().getIntensity();
				QString meta = QString("%1Quality: %2, Intens.Ratio: %3, Dist.: RT %4, MZ %5; ")
											.arg(string(first->second.getMetaValue(3)))
											.arg(pairs[i].getQuality(),4,'f',2)
											.arg(ratio,4,'f',2)
											.arg(fabs(diff[RT]),4,'f',2)
											.arg(fabs(diff[MZ]),4,'f',2);
				first->second.setMetaValue(3,meta.ascii() );
			}

			for (UniqueFeatureMap::const_iterator it=tmp.begin(); it!=tmp.end(); ++it)
				map.push_back(it->second);
		}

		void PairMatcher::printInfo(std::ostream& out, const PairVectorType& pairs)
		{
			out << "Found the following " << pairs.size() << " pairs:\n"
					<< "Quality\tFirst[RT]\tFirst[MZ]\tFirst[Int]\tFirst[Corr]"
					<< "\tSecond[RT]\tSecond[MZ]\tSecond[Int]\tSecond[Corr]"
					<< "\tRatio\tCharge\tDiff[RT]\tDiff[MZ]\n";
			for (Size i=0; i<pairs.size(); ++i)
			{
				DPosition<2> diff
					= pairs[i].getFirst().getPosition()-pairs[i].getSecond().getPosition();
				out <<	QString("%1\t%2\t%3\t%4\t%5\t%6\t%7\t%8\t%9")
										.arg(pairs[i].getQuality(),7,'f',2)
										.arg(pairs[i].getFirst().getPosition()[0],7,'f',2)
										.arg(pairs[i].getFirst().getPosition()[1],7,'f',2)
										.arg(pairs[i].getFirst().getIntensity(),7,'f',2)
										.arg(pairs[i].getFirst().getOverallQuality(),7,'f',2)
										.arg(pairs[i].getSecond().getPosition()[0],7,'f',2)
										.arg(pairs[i].getSecond().getPosition()[1],7,'f',2)
										.arg(pairs[i].getSecond().getIntensity(),7,'f',2)
										.arg(pairs[i].getSecond().getOverallQuality(),7,'f',2)
						<< QString("\t%1\t%2\t%3\t%4\n")
										.arg(pairs[i].getFirst().getIntensity()/
												 pairs[i].getSecond().getIntensity(),7,'f',2)
										.arg(pairs[i].getFirst().getCharge())
										.arg(diff[RT],4,'f',2)
										.arg(diff[MZ],4,'f',2);
			}
		}

}
