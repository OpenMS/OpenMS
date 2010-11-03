#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletTransform.h>

namespace OpenMS
{
bool myCudaComparator (const cudaHelp& a, const cudaHelp& b)
{
	return (a.intens < b.intens);
}
}
