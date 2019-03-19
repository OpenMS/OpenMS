//#include <OpenMS/QC/QualityControl.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/QC/ChargeDistribution.h>
#include <algorithm>

using namespace OpenMS;
using namespace std;

ChargeDistribution::ChargeDistribution(FeatureMap& features) :features_(features)
{
}

vector<pair<uint64_t, uint64_t>> ChargeDistribution:: calculate (FeatureMap features_) const
{
		vector<uint32_t> all_charges;
		for (uint64_t i=0; i<=features_.size(); i++)
		{ 
			all_charges.push_back(features_[i].getCharge);
		}
		uint32_t * max = max_element(all_charges.begin(), all_charges.end());
		if (max != all_charges.end())
			int maximum = *max;
		else
			int maximum = 0;

		const vector<pair <uint64_t, uint64_t>> ChargeDistribution::getCharges() const
		{
			return charges;
		}


};

	

