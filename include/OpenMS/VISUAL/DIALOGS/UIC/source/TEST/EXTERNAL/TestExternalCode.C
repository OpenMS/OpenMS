#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

using namespace OpenMS;

int main(int argc, char *argv[])
{

	FeatureMap<> fm;
	Feature feature;
	fm.push_back(feature);
	
	std::cout << "All good and well!\n";
	
	return 0;
}
