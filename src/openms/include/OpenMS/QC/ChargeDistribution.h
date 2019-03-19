#pragma once
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CONCEPT/Types.h>

class ChargeDistribution
{
public:
	ChargeDistribution(OpenMS::FeatureMap& features);
	const std::vector<std::pair <UInt32, Size>> getCharges() const;
	std::vector<std::pair <UInt32, Size>> calculate (OpenMS::FeatureMap& features) const;

private:
	
	const OpenMS::FeatureMap& features_;
	std::vector<std::pair <UInt32, Size>> charges{};

};