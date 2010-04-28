#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include "ExampleLibraryFile.h"

using namespace OpenMS;
using namespace OpenMSExternal;

int main(int argc, char *argv[])
{

	FeatureMap<> fm;
	Feature feature;
	fm.push_back(feature);
	std::string s = ExampleLibraryFile::printSomething();
	std::cout << "From external lib: " << s << "\n";
	std::cout << "All good and well!\n";
	
	return 0;
}
