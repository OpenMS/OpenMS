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
// $Maintainer: Peter Kunszt $
// $Authors: Peter Kunszt $
// --------------------------------------------------------------------------
// 

#ifndef OPENMS_FEATUREFINDER_SUPERHIRN_SUPERHIRNPARAMETERS_H
#define OPENMS_FEATUREFINDER_SUPERHIRN_SUPERHIRNPARAMETERS_H

namespace OpenMS
{

	/**
	 @brief SuperHirn parameters singleton class containing all static configuration variables
	 */

	class SuperHirnParameters
	{
		public:
			static SuperHirnParameters* instance();

			~SuperHirnParameters();

			double getBackgroundIntensityBinsTR();
			double getBackgroundIntensityBinsMZ();
			double getBackgroundIntensityBinsIntens();
			int getBackgroundIntensityBinsMinBinCount();

			double getMinTR();
			double getMaxTR();
			double getMinFeatureMZ();
			double getMaxFeatureMZ();
			int getMinFeatureChrg();
			int getMaxFeatureChrg();

			float getIntensityThreshold();
			double getToleranceMZ();
			double getMaxInterScanRetentionTimeDistance();
			int getMinNbClusterMembers();
			std::map<int, float>* getScanTRIndex();
			//    bool getMonoIsoDebugging();
			//    double getDebugMonoIsoMassMin();
			//    double getDebugMonoIsoMassMax();
			//    double getMS1IntensityApexPercentilCutoff();
			double getMS1TRResolution();
			bool centroidDataModus();
			int getCentroidWindowWidth();

			double getDetectableIsotopeFactor();
			/*
			 * @brief Maximal deviation between expected and measured isotopic intensities
			 */
			double getIntensityCV();

			/*
			 * @brief  Mass tolerance in ppm between isotopes
			 */
			double getMassTolPpm();

			/*
			 * @brief Mass tolerance in Da between isotopes - total mass to = mass*fMassTolPpm/1000000 + fMassTolDa
			 */
			double getMassTolDa();

			/*
			 * @brief Peak below this values are not considered as monoisotopic peaks
			 */
			double getMinIntensity();

			/*
			 * @brief Intensities below this value are considered as 0
			 */
			double getIntensityFloor();

			/*
			 * @brief M/Z tolerance in Parts per Million
			 */
			double getMzTolPpm();

			/*
			 * @brief TR tolerance 
			 */
			double getTrTol();

			double getPeptideProbabilityThreshold();
			bool storeAllLowProbabilityMS2Scans();

			bool createFeatureElutionProfiles();
			bool ms1FeatureClustering();

			double getMs1PeakAreaTrResolution();
			double getInitialTrTolerance();
			double getMs1FeatureMergingTrTolerance();
			double getPercentageIntensityElutionBorderVariation();
			double getPpmToleranceForMZClustering();

		double getLowIntensityMSSignalThreshold();

		bool isInitIsotopeDist();
		void setInitIsotopeDist();

			friend class FeatureFinderAlgorithmSHCtrl;

		private:
			/*
			 * @brief The constructors are private as only the class itself
			 * can construct an instance. Use the static instance() method to access the parameters.
			 */
			SuperHirnParameters();
			SuperHirnParameters(const SuperHirnParameters &);
			SuperHirnParameters & operator =(const SuperHirnParameters &);

			static SuperHirnParameters *instance_;  // the singleton instance
			static bool haveInstance_;

			double backgroundIntensityBinsTR_;
			double backgroundIntensityBinsMZ_;
			double backgroundIntensityBinsIntens_;
			int backgroundIntensityBinsMinBinCount_;

			double minTR_;
			double maxTR_;
			double minFeatureMZ_;
			double maxFeatureMZ_;
			int minFeatureChrg_;
			int maxFeatureChrg_;

			float intensityThreshold_; // minimal intensity level:  NEVER USED
			double toleranceMZ_;   // m/z tolerance value: NEVER CONFIGURED
			double maxInterScanRetentionTimeDistance_;   // max_distance from next elution peak member in min.
			int minNbClusterMembers_;   // define minimal number of members in LC elution peaks cluster

			std::map<int, float> scanTRIndex_;

			//    bool monoIsoDebugging_;   // to track detected monoistopic mass for debugging
			//    double debugMonoIsoMassMin_;
			//    double debugMonoIsoMassMax_;
			//    double ms1IntensityApexPercentilCutoff_;
			double ms1TRResolution_;
			bool centroidDataModus_;   // if data are in centroid form or not
			int centroidWindowWidth_;

			//    int	reportMonoPeaks_; // 1 if info about monoisotopic peaks should be written to mono_peaks.txt
			//    std::string debugDirectory_; // Directory where peak detection debug files are written
			//    int	reportScanNumber_; // if sfReportMonoPeaks is set to 1, details about this spectrum will be written to debug files

			//    int ms1BaseInterScanDistance_;
			//    int ms2BaseInterScanDistance;
			//    bool ms2PeakProcessing_;

			//    std::vector<double> peakExtractionScanLevels_;
			//    std::vector<double> fragmentMassScanLevels_;

			double detectableIsotopeFactor_;
			double intensityCV_;

			double massTolPpm_;  // mass tolerance in ppm between isotopes
			double massTolDa_; // mass tolerance in Da between isotopes - total mass to = mass*fMassTolPpm/1000000 + fMassTolDa
			double minIntensity_; // peak below this values are not considered as monoisotopic peaks
			double intensityFloor_; // intensities below this value are considered as 0


		  double mzTolPpm_;		// tolerance in m/z and TR:
			double trTol_;
			double peptideProbabilityThreshold_;
			bool storeAllLowProbabilityMS2Scans_;

			bool createFeatureElutionProfiles_;
		  bool ms1FeatureClustering_;

		  double ms1PeakAreaTrResolution_;
		  double initialTrTolerance_;
		  double ms1FeatureMergingTrTolerance_;
		  double percentageIntensityElutionBorderVariation_;
		  double ppmToleranceForMZClustering_;

			double lowIntensityMSSignalThreshold_;
			bool initIsotopeDist_;

	};

//------------------------- inline methods --------------------------

	inline SuperHirnParameters::SuperHirnParameters(const SuperHirnParameters &)
	{
	}

	inline SuperHirnParameters & SuperHirnParameters::operator =(const SuperHirnParameters &)
	{
		return *this;
	}

	inline SuperHirnParameters* SuperHirnParameters::instance()
	{
		if (haveInstance_)
		{
			return instance_;
		}
		instance_ = new SuperHirnParameters();
		haveInstance_ = true;
		return instance_;
	}

	inline double SuperHirnParameters::getBackgroundIntensityBinsTR()
	{
		return backgroundIntensityBinsTR_;
	}

	inline double SuperHirnParameters::getBackgroundIntensityBinsMZ()
	{
		return backgroundIntensityBinsMZ_;
	}

	inline double SuperHirnParameters::getBackgroundIntensityBinsIntens()
	{
		return backgroundIntensityBinsIntens_;
	}

	inline int SuperHirnParameters::getBackgroundIntensityBinsMinBinCount()
	{
		return backgroundIntensityBinsMinBinCount_;
	}

	inline double SuperHirnParameters::getMinTR()
	{
		return minTR_;
	}

	inline double SuperHirnParameters::getMaxTR()
	{
		return maxTR_;
	}

	inline double SuperHirnParameters::getMinFeatureMZ()
	{
		return minFeatureMZ_;
	}

	inline double SuperHirnParameters::getMaxFeatureMZ()
	{
		return maxFeatureMZ_;
	}

	inline int SuperHirnParameters::getMinFeatureChrg()
	{
		return minFeatureChrg_;
	}

	inline int SuperHirnParameters::getMaxFeatureChrg()
	{
		return maxFeatureChrg_;
	}

	inline float SuperHirnParameters::getIntensityThreshold()
	{
		return intensityThreshold_;
	}

	inline double SuperHirnParameters::getToleranceMZ()
	{
		return toleranceMZ_;
	}

	inline double SuperHirnParameters::getMaxInterScanRetentionTimeDistance()
	{
		return maxInterScanRetentionTimeDistance_;
	}

	inline int SuperHirnParameters::getMinNbClusterMembers()
	{
		return minNbClusterMembers_;
	}

	inline std::map<int, float>* SuperHirnParameters::getScanTRIndex()
	{
		return &scanTRIndex_;
	}

	/*
	 inline bool SuperHirnParameters::getMonoIsoDebugging()
	 {
	 return monoIsoDebugging_;
	 }

	 inline double SuperHirnParameters::getDebugMonoIsoMassMin()
	 {
	 return debugMonoIsoMassMin_;
	 }

	 inline double SuperHirnParameters::getDebugMonoIsoMassMax()
	 {
	 return debugMonoIsoMassMax_;
	 }

	 inline double SuperHirnParameters::getMS1IntensityApexPercentilCutoff()
	 {
	 return ms1IntensityApexPercentilCutoff_;
	 }
	 */

	inline double SuperHirnParameters::getMS1TRResolution()
	{
		return ms1TRResolution_;
	}

	inline bool SuperHirnParameters::centroidDataModus()
	{
		return centroidDataModus_;
	}

	inline int SuperHirnParameters::getCentroidWindowWidth()
	{
		return centroidWindowWidth_;
	}

	inline double SuperHirnParameters::getDetectableIsotopeFactor()
	{
		return detectableIsotopeFactor_;
	}

	inline double SuperHirnParameters::getIntensityCV()
	{
		return intensityCV_;
	}

	inline double SuperHirnParameters::getMassTolPpm()
	{
		return massTolPpm_;
	}

	inline double SuperHirnParameters::getMassTolDa()
	{
		return massTolDa_;
	}

	inline double SuperHirnParameters::getMinIntensity()
	{
		return minIntensity_;
	}

	inline double SuperHirnParameters::getIntensityFloor()
	{
		return intensityFloor_;
	}

	inline double SuperHirnParameters::getMzTolPpm()
	{
		return mzTolPpm_;
	}

	inline double SuperHirnParameters::getTrTol()
	{
		return trTol_;
	}

	inline double SuperHirnParameters::getPeptideProbabilityThreshold()
	{
		return peptideProbabilityThreshold_;
	}

	inline bool SuperHirnParameters::storeAllLowProbabilityMS2Scans()
	{
		return storeAllLowProbabilityMS2Scans_;
	}


	inline bool SuperHirnParameters::createFeatureElutionProfiles()
	{
		return createFeatureElutionProfiles_;
	}

	inline bool SuperHirnParameters::ms1FeatureClustering()
	{
		return ms1FeatureClustering_;
	}

	inline bool SuperHirnParameters::isInitIsotopeDist()
	{
		return initIsotopeDist_;
	}

	inline void SuperHirnParameters::setInitIsotopeDist()
	{
		initIsotopeDist_ = true;
	}

	inline double SuperHirnParameters::getMs1PeakAreaTrResolution()
	{
		return ms1PeakAreaTrResolution_;
	}

	inline double SuperHirnParameters::getInitialTrTolerance()
	{
		return initialTrTolerance_;
	}

	inline double SuperHirnParameters::getMs1FeatureMergingTrTolerance()
	{
		return ms1FeatureMergingTrTolerance_;
	}

	inline double SuperHirnParameters::getPercentageIntensityElutionBorderVariation()
	{
		return percentageIntensityElutionBorderVariation_;
	}

	inline double SuperHirnParameters::getPpmToleranceForMZClustering()
	{
		return ppmToleranceForMZClustering_;
	}

	inline double SuperHirnParameters::getLowIntensityMSSignalThreshold() {
		return lowIntensityMSSignalThreshold_;
	}

}

#endif
