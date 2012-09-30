/*
 * DIAHelpers.C
 *
 *  Created on: Aug 24, 2012
 *      Author: witold
 */

#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DIAHelpers.h"
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace OpenSwath
{

	void normalize(
			const std::vector<double> & intensities,
			double normalizer,
			std::vector<double> & normalized_intensities )
	{	//compute total intensities
		//normalize intensities
		double res = 0;
		if (normalizer >0) {
			std::transform(intensities.begin(), intensities.end(),
					normalized_intensities.begin(), std::bind2nd(std::divides<double>(), normalizer));
		}
	}

	void integrateWindows(const OpenSwath::SpectrumPtr spectrum,
			const std::vector<double> & windowsCenter, double width,
			std::vector<double> & integratedWindowsIntensity,
			std::vector<double> & integratedWindowsMZ,
			bool remZero)
	{
		std::vector<double>::const_iterator beg = windowsCenter.begin();
		std::vector<double>::const_iterator end = windowsCenter.end();

		double mz, intensity;
		for (; beg != end; ++beg)
		{
			double left = *beg - width / 2.0;
			double right = *beg + width / 2.0;
			if(integrateWindow(spectrum, left, right, mz, intensity, false)){
				integratedWindowsIntensity.push_back(intensity);
				integratedWindowsMZ.push_back(mz);
			}else if(!remZero){
				integratedWindowsIntensity.push_back(0.);
				integratedWindowsMZ.push_back(*beg);
			}else{}
		}
	}


	/* integrate all masses in window */
	bool integrateWindow(const OpenSwath::SpectrumPtr spectrum, double mz_start,
			double mz_end, double & mz, double & intensity, bool centroided)
	{
		//check precondtion
		if (std::adjacent_find(spectrum->getMZArray()->data.begin(),
				spectrum->getMZArray()->data.end(), std::greater<double>())
				!= spectrum->getMZArray()->data.end()) {
			throw std::runtime_error("Precondition MZ vector needs to be sorted!");
		}

		intensity = 0;
		if (!centroided) {
			// get the weighted average for noncentroided data.
			// TODO this is not optimal if there are two peaks in this window (e.g. if
			// the window is too large)
			mz = 0;
			intensity = 0;

			std::vector<double>::const_iterator mz_arr_end =
					spectrum->getMZArray()->data.end();
			std::vector<double>::const_iterator int_it =
					spectrum->getIntensityArray()->data.begin();

			// this assumes that the spectra are sorted!
			std::vector<double>::const_iterator mz_it = std::lower_bound(
					spectrum->getMZArray()->data.begin(),
					spectrum->getMZArray()->data.end(), mz_start);
			std::vector<double>::const_iterator mz_it_end = std::lower_bound(mz_it,
					mz_arr_end, mz_end);

			// also advance intensity iterator now
			int iterator_pos =
					std::distance(
							(std::vector<double>::const_iterator) spectrum->getMZArray()->data.begin(),
							mz_it);
			std::advance(int_it, iterator_pos);

			for (; mz_it != mz_it_end; ++mz_it, ++int_it) {
				intensity += (*int_it);
				mz += (*int_it) * (*mz_it);
			}

			if (intensity > 0.) {
				mz /= intensity;
				return true;
			} else {
				mz = -1;
				intensity = 0;
				return false;
			}

		} else {
			// not implemented
			throw "Not implemented";
		}
	}

}

