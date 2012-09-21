// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Florian Zeller $
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////
//
//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
//
//  written by Markus Mueller, markus.mueller@imsb.biol.ethz.ch
//  and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch
//  October 2005
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
// 

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <stdio.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MSPeak.h>

namespace OpenMS
{

	using namespace std;

////////////////////////////////////////////////
// constructor for the object MSPeak:
	MSPeak::MSPeak()
	{
		MZ = 0;
		INTENSITY = 0;
		SCAN = 0;
		TR = 0;
		CHRG = 0;
		NRISOTOPES = 0;
		SCORE = 0;
		precursorMZ = 0;
		SignalToNoise = 1;
		precursorMass = false;
		childScan = -1;

	}

////////////////////////////////////////////////
// constructor for the object MSPeak:
	MSPeak::MSPeak(int IN_scan, double IN_mass, float IN_intens)
	{
		MZ = IN_mass;
		INTENSITY = IN_intens;
		SCAN = IN_scan;
		TR = 0;
		NRISOTOPES = 0;
		SCORE = 0;
		CHRG = 0;
		precursorMZ = 0;
		SignalToNoise = 1;
		precursorMass = false;
		childScan = -1;

	}

////////////////////////////////////////////////
// constructor for the object MSPeak:
	MSPeak::MSPeak(int IN_scan, double IN_mass, float IN_intens, unsigned int IN_CHRG, unsigned int IN_NRISOTOPES,
			float IN_SCORE, vector<CentroidPeak> IN_ISOPEAKS)
	{
		MZ = IN_mass;
		INTENSITY = IN_intens;
		SCAN = IN_scan;
		TR = 0;
		precursorMZ = 0;
		SignalToNoise = 1;

		CHRG = IN_CHRG;
		NRISOTOPES = IN_NRISOTOPES;
		SCORE = IN_SCORE;
		ISOPEAKS = IN_ISOPEAKS;
		precursorMass = false;
		childScan = -1;

	}

//////////////////////////////////////////////////
// class desctructor of MSPeak
	MSPeak::~MSPeak()
	{
		MZ = 0;
		INTENSITY = 0;
		SCAN = 0;
		TR = 0;
		CHRG = 0;
		NRISOTOPES = 0;
		SCORE = 0;
		precursorMZ = 0;
		SignalToNoise = 0;
		precursorMass = false;
		childScan = -1;

	}

//////////////////////////////////////////////////
// copy constructor:
	MSPeak::MSPeak(const MSPeak& tmp)
	{
		MZ = tmp.MZ;
		SignalToNoise = tmp.SignalToNoise;
		INTENSITY = tmp.INTENSITY;
		SCAN = tmp.SCAN;
		TR = tmp.TR;
		CHRG = tmp.CHRG;
		SCORE = tmp.SCORE;
		NRISOTOPES = tmp.NRISOTOPES;
		ISOPEAKS = tmp.ISOPEAKS;
		precursorMZ = tmp.precursorMZ;
		precursorMass = tmp.precursorMass;
		childScan = tmp.childScan;
		extraMSPeakInfo = tmp.extraMSPeakInfo;

	}

//////////////////////////////////////////////////
// copy constructor:
	MSPeak::MSPeak(const MSPeak* tmp)
	{
		MZ = tmp->MZ;
		SignalToNoise = tmp->SignalToNoise;
		INTENSITY = tmp->INTENSITY;
		SCAN = tmp->SCAN;
		TR = tmp->TR;
		CHRG = tmp->CHRG;
		SCORE = tmp->SCORE;
		NRISOTOPES = tmp->NRISOTOPES;
		ISOPEAKS = tmp->ISOPEAKS;
		precursorMZ = tmp->precursorMZ;
		precursorMass = tmp->precursorMass;
		childScan = tmp->childScan;
		extraMSPeakInfo = tmp->extraMSPeakInfo;

	}

//////////////////////////////////////////////////
// copy constructor:
	MSPeak& MSPeak::operator=(const MSPeak& tmp)
	{
		MZ = tmp.MZ;
		SignalToNoise = tmp.SignalToNoise;
		INTENSITY = tmp.INTENSITY;
		SCAN = tmp.SCAN;
		TR = tmp.TR;
		CHRG = tmp.CHRG;
		SCORE = tmp.SCORE;
		NRISOTOPES = tmp.NRISOTOPES;
		ISOPEAKS = tmp.ISOPEAKS;
		precursorMZ = tmp.precursorMZ;
		precursorMass = tmp.precursorMass;
		childScan = tmp.childScan;
		extraMSPeakInfo = tmp.extraMSPeakInfo;

		return *this;
	}

//////////////////////////////////////////////////
// show info:
	void MSPeak::show_info()
	{
		printf("mz=%0.4f,int=%0.1f,scan=%d,tr=%0.2f,+%d", MZ, INTENSITY, SCAN, TR, CHRG); // precursorMZ

		if (precursorMZ > 1.0)
		{
			printf(",preMZ=%0.4f\n", precursorMZ);
		}
		else
		{
			printf("\n");
		}

		if (!getExtraPeakInfo().empty())
		{
			cout << getExtraPeakInfo() << endl;
		}

		// print the isotope pattern:

		if (!ISOPEAKS.empty())
		{
			printf("\t");
			vector<CentroidPeak>::iterator I = ISOPEAKS.begin();
			while (I != ISOPEAKS.end())
			{
				printf("%0.4f(%0.0f[%0.0f]) ", (*I).getMass(), (*I).getFittedIntensity(), (*I).getOrgIntensity());
				I++;
			}
			printf("\n");
		}

	}

// copied from simple_math
	double simple_math_getMassErrorAtPPMLevel(double mz, double PPM_TOLERANCE)
	{
		double ppmValue = mz / 1000000.00;
		return ppmValue * PPM_TOLERANCE;
	}

// copied from simple_math
	bool simple_math_compareMassValuesAtPPMLevel2(double mzA, double mzB, double PPM_TOLERANCE)
	{

		// take the average mass:
		double avMass = (mzA + mzB) / 2.0;

		// define the parts per million:
		double ppmValue = avMass / 1000000.00;
		double ppmDeltaTol = ppmValue * PPM_TOLERANCE;

		double deltaMass = fabs(mzA - mzB);
		if (deltaMass > ppmDeltaTol)
		{
			return false;
		}

		return true;
	}

//////////////////////////////////////////////////
// check if the input mass matches one of the isotopic masses
	bool MSPeak::checkIsotopeBelongingAndAdjustMass(double mass, double mzTolerance)
	{

		// check if the mass is really smaller then the mono isotope,
		// then it cant bellong to this peak:
		double deltaSmall = get_MZ() - mass;
		deltaSmall -= simple_math_getMassErrorAtPPMLevel(mass, mzTolerance);
		if (deltaSmall > 0)
		{
			return false;
		}

		// ok if not, then check the isotopes:
		// highest isotope to consider:
		int max = 2;
		// check now also the isotopes
		if (!ISOPEAKS.empty())
		{
			vector<CentroidPeak>::iterator I = ISOPEAKS.begin();
			int i = 1;
			while (I != ISOPEAKS.end())
			{
				if (simple_math_compareMassValuesAtPPMLevel2(mass, (*I).getMass(), mzTolerance))
				{
					return true;
				}
				i++;
				I++;
				if (i > max)
				{
					break;
				}
			}
		}
		else
		{
			if (simple_math_compareMassValuesAtPPMLevel2(mass, MZ, mzTolerance))
			{
				return true;
			}
		}

		return false;
	}

//////////////////////////////////////////////////
// store the MS/MS scan number and activate this peak as precursor peak:
	void MSPeak::activateAsPrecursorPeak(int in)
	{

		precursorMass = true;
		childScan = in;

	}

}
