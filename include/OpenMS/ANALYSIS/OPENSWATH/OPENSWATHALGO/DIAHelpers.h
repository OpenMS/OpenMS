/*
 * DIAHelpers.h
 *
 *  Created on: Aug 24, 2012
 *      Author: witold
 */

#ifndef DIAHELPERS_H_
#define DIAHELPERS_H_

#include <cmath>
#include <vector>
#include <numeric>
#include <boost/bind.hpp>
#include <complex>
#include <algorithm>
#include <cmath>

#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h"
namespace OpenSwath
{

	void normalize(const std::vector<double> & intensities,
			double normalization_factor,
			std::vector<double> & normalized_intensities);

	template<typename T>
	double norm(T beg, T end)
	{
		double res = 0.0;
		for (; beg != end; ++beg) {
			double tmp = *beg;
			res += tmp * tmp;
		}
		return sqrt(res);
	}

	//integrate Window
	bool integrateWindow(const OpenSwath::SpectrumPtr spectrum, double mz_start,
			double mz_end, double & mz, double & intensity, bool centroided = false);
	//integrate Window
	void integrateWindows(const OpenSwath::SpectrumPtr spectrum,
			const std::vector<double> & windowsCenter, double width,
			std::vector<double> & integratedWindowsIntensity,
			std::vector<double> & integratedWindowsMZ, bool remZero = false);

	template<typename Texp, typename Ttheo>
	double dotProd(Texp intExpBeg, Texp intExpEnd, Ttheo intTheo)
	{
		std::vector<double> res(std::distance(intExpBeg, intExpEnd));
		std::transform(intExpBeg, intExpEnd, intTheo, res.begin(),
				std::multiplies<double>());
		double sum = std::accumulate(res.begin(), res.end(), 0.);
		return sum;
	}

	inline double dotprodScoring(std::vector<double> intExp,
			std::vector<double> theorint)
	{

		for(int i = 0 ; i < intExp.size(); ++i){

			intExp[i] = sqrt(intExp[i]);
			theorint[i] = sqrt(theorint[i]);
			//std::transform(intExp.begin(), intExp.end(), intExp.begin(), sqrt);
			//std::transform(theorint.begin(), theorint.end(), theorint.begin(), sqrt);
		}

		double intExptotal = norm(intExp.begin(), intExp.end());
		double intTheorTotal = norm(theorint.begin(), theorint.end());
		OpenSwath::normalize(intExp, intExptotal, intExp);
		OpenSwath::normalize(theorint, intTheorTotal, theorint);
		double score2 = OpenSwath::dotProd(intExp.begin(), intExp.end(),
				theorint.begin());
		return score2;
	}

	template<typename Texp, typename Ttheo>
	double manhattanDist(Texp itExpBeg, Texp itExpEnd, Ttheo itTheo)
	{
		double sum = 0.0;
		for (std::size_t i = 0; itExpBeg < itExpEnd; ++itExpBeg, ++itTheo, ++i) {
			double x = *itExpBeg - *itTheo;
			x = fabs(x);
			sum += x;
		}
		return sum;
	}

	inline double manhattanScoring(std::vector<double> intExp,
			std::vector<double> theorint)
	{

		for(int i = 0 ; i < intExp.size(); ++i){

			intExp[i] = sqrt(intExp[i]);
			theorint[i] = sqrt(theorint[i]);
			//std::transform(intExp.begin(), intExp.end(), intExp.begin(), sqrt);
			//std::transform(theorint.begin(), theorint.end(), theorint.begin(), sqrt);
		}

		double intExptotal = std::accumulate(intExp.begin(), intExp.end(), 0.0);
		double intTheorTotal = std::accumulate(theorint.begin(), theorint.end(),
				0.0);
		OpenSwath::normalize(intExp, intExptotal, intExp);
		OpenSwath::normalize(theorint, intTheorTotal, theorint);
		double score2 = OpenSwath::manhattanDist(intExp.begin(), intExp.end(),
				theorint.begin());
		return score2;
	}
}

#endif /* DIAHELPERS_H_ */
