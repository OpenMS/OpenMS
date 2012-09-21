// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Florian Zeller $
// $Authors: Markus Mueller $
// --------------------------------------------------------------------------
//
/*
 *  CentroidPeak.h
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *  
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */

#ifndef _CENTROIDPEAK_H_
#define _CENTROIDPEAK_H_

#include <OpenMS/CONCEPT/Types.h>

#include <ostream>
#include <cmath>
#include <vector>

namespace OpenMS
{

	class OPENMS_DLLAPI CentroidPeak
	{
		public:

//		static int sfCentroidWindowWidth;

			CentroidPeak();
			CentroidPeak(double, double);
			CentroidPeak(double, double, double);
			CentroidPeak(const CentroidPeak&);
			// Copy constructor
			CentroidPeak(const CentroidPeak*);

			CentroidPeak& operator=(const CentroidPeak&);

			bool operator<(const CentroidPeak&);

			virtual ~CentroidPeak();

			// getters and setters
			double getMass();
			double getIntensity();
			int getIsotopIdx();
			double getSignalToNoise();
			double getFittedIntensity();
			double getOrgIntensity();
			std::string getExtraPeakInfo();
			double getRetentionTime();

			void setMass(double pMass);
			void setIntensity(double pIntensity);
			void setIsotopIdx(double pIsotopIdx);
			void setSignalToNoise(double in);
			void setFittedIntensity(double pFittedIntensity);
			void setOrgIntensity(double pOrgIntensity);
			void setExtraPeakInfo(std::string in);
			void setRetentionTime(double in);

			void show_info();
			void subtractIntensity(double);

		protected:

			int isotopIdx_;
			double mass_;
			double intensity_;
			double fittedIntensity_;
			double orgIntensity_;
			double tr_;
			double signalToNoise_;
			std::string extraPeakInfo_;
	};

	// Class for deconvoluted isotopic patterns
	class DeconvPeak: public CentroidPeak
	{
		public:

			DeconvPeak();
			DeconvPeak(double, double, int, int, double, double);
			DeconvPeak(const DeconvPeak&);
			DeconvPeak(const DeconvPeak*);

			DeconvPeak& operator=(const DeconvPeak&);

			virtual ~DeconvPeak();

			// shows the info of the peak:
			void show_info();

			// getters and setters
			int getCharge();
			int getNrIsotopes();
			double getC13MassError();
			double getScore();
			std::vector<CentroidPeak> getIsotopicPeaks();

			void setCharge(int pCharge);
			void setNrIsotopes(int pNrIsotopes);
			void setC13MassError(double pC13MassError);
			void setScore(double pScore);
			void setIsotopicPeaks(std::vector<CentroidPeak> pIsotopicPeaks);

		protected:

			int charge_;
			int nrIsotopes_;
			double c13MassError_;
			double score_;
			std::vector<CentroidPeak> isotopicPeaks_;
	};

	// stream operators
	std::ostream& operator<<(std::ostream&, CentroidPeak&);
	std::ostream& operator<<(std::ostream&, DeconvPeak&);

	// inline implementation of getters and setters

	inline double CentroidPeak::getMass()
	{
		return mass_;
	}

	inline double CentroidPeak::getIntensity()
	{
		return intensity_;
	}

	inline int CentroidPeak::getIsotopIdx()
	{
		return isotopIdx_;
	}

	inline double CentroidPeak::getSignalToNoise()
	{
		return signalToNoise_;
	}

	inline double CentroidPeak::getFittedIntensity()
	{
		return fittedIntensity_;
	}

	inline double CentroidPeak::getOrgIntensity()
	{
		return orgIntensity_;
	}

	inline std::string CentroidPeak::getExtraPeakInfo()
	{
		return extraPeakInfo_;
	}

	inline double CentroidPeak::getRetentionTime()
	{
		return tr_;
	}

	inline void CentroidPeak::setMass(double pMass)
	{
		mass_ = pMass;
	}

	inline void CentroidPeak::setIntensity(double pIntensity)
	{
		intensity_ = pIntensity;
	}

	inline void CentroidPeak::setIsotopIdx(double pIsotopIdx)
	{
		isotopIdx_ = (int) pIsotopIdx;
	}

	inline void CentroidPeak::setSignalToNoise(double in)
	{
		signalToNoise_ = in;
	}

	inline void CentroidPeak::setFittedIntensity(double pFittedIntensity)
	{
		fittedIntensity_ = pFittedIntensity;
	}

	inline void CentroidPeak::setOrgIntensity(double pOrgIntensity)
	{
		orgIntensity_ = pOrgIntensity;
	}

	inline void CentroidPeak::setExtraPeakInfo(std::string in)
	{
		extraPeakInfo_ = in;
	}

	inline void CentroidPeak::setRetentionTime(double in)
	{
		tr_ = in;
	}

	//

	inline int DeconvPeak::getCharge()
	{
		return charge_;
	}

	inline int DeconvPeak::getNrIsotopes()
	{
		return nrIsotopes_;
	}

	inline double DeconvPeak::getC13MassError()
	{
		return c13MassError_;
	}

	inline double DeconvPeak::getScore()
	{
		return score_;
	}

	inline std::vector<CentroidPeak> DeconvPeak::getIsotopicPeaks()
	{
		return isotopicPeaks_;
	}

	inline void DeconvPeak::setCharge(int pCharge)
	{
		charge_ = pCharge;
	}

	inline void DeconvPeak::setC13MassError(double pC13MassError)
	{
		c13MassError_ = pC13MassError;
	}

	inline void DeconvPeak::setNrIsotopes(int pNrIsotopes)
	{
		nrIsotopes_ = pNrIsotopes;
	}

	inline void DeconvPeak::setScore(double pScore)
	{
		score_ = pScore;
	}

	inline void DeconvPeak::setIsotopicPeaks(std::vector<CentroidPeak> pIsotopicPeaks)
	{
		isotopicPeaks_ = pIsotopicPeaks;
	}

} // ns

#endif
