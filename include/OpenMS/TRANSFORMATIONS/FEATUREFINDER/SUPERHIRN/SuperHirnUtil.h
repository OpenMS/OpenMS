/*
 * SuperHirnUtil.h
 *
 *  Created on: Oct 10, 2011
 *      Author: pkunszt
 */

#ifndef SUPERHIRNUTIL_H_
#define SUPERHIRNUTIL_H_

#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{

	class OPENMS_DLLAPI SuperHirnUtil
	{
		public:

			/*
			 * @brief Compare two masses at the PPM value and decide if they fall into the m/z tolerance window
			 *
			 */
			static bool compareMassValuesAtPPMLevel(double massA, double massB, double ppmTolerance);

			/*
			 * @brief Get mass error at PPM value - basically ppmTolerance * m/z /10^6
			 */
			static double getMassErrorAtPPMLevel(double mz, double ppmTolerance);
	};
//////////////////////////////////////////////////////
// compare two masses at the PPM value and decide
// if they fall into the m/z tolerance window
	inline bool SuperHirnUtil::compareMassValuesAtPPMLevel(double mzA, double mzB, double ppmTol)
	{

		// take the average mass, parts per million:
		double avMassPPM = (mzA + mzB) / 2000000.0;
		double ppmDeltaTol = avMassPPM * ppmTol;

		double deltaMass = fabs(mzA - mzB);
		if (deltaMass > ppmDeltaTol)
		{
			return false;
		}

		return true;
	}

//////////////////////////////////////////////////////
// get the mass error at the PPM value
	inline double SuperHirnUtil::getMassErrorAtPPMLevel(double mz, double ppmTol)
	{
		return mz * ppmTol / 1000000.00;
	}

}
#endif /* SUPERHIRNUTIL_H_ */
