// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_LCMSSIM_H
#define OPENMS_SIMULATION_LCMSSIM_H

// STL includes
#include <map>
#include <vector>
#include <cmath>
#include <iostream>
#include <ctime>

// GSL includes (random number generation)
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// OpenMS includes
#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/FORMAT/DTA2DFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ProductModel.h>

#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/SIMULATION/LCMSSample.h>
#include <OpenMS/SIMULATION/IsotopeModelGeneral.h>
#include <OpenMS/SIMULATION/MixtureModel.h>
#include <OpenMS/SIMULATION/ElutionModel.h>

namespace OpenMS{

	/**
		@brief Simulation of a LC/MS experiment

	*/
	class OPENMS_DLLAPI LCMSSim
		: public DefaultParamHandler, public ProgressLogger
	{

	public:

		/** @name Typedefs
		*/
		//@{
		/// An LC-MS data point
		typedef Peak1D PointType;

		/// Coordinate (m/z or rt) of an LC-MS data point
		typedef PointType::CoordinateType CoordinateType;

		/// Intensity of an LC-MS data point
		typedef PointType::IntensityType IntensityType;

		/// Charge of a peptide
		typedef Feature::ChargeType ChargeType;

		/// The peptides
		typedef LCMSSample::PeptideSequences PeptideSequences;

		/// Retention time table (format: rt , pointer to sequence)
		typedef std::multimap< CoordinateType, PeptideSequences::const_iterator > RTTable;

		/// Data structure for LC-MS spectra
		typedef MSExperiment< PointType > LCMSmap;

		/// Possible ionization methods
		typedef enum
		{
			SIMPLE,	   // Only 1+ charges (debugging)
			ESI       		// Electro Spray Ionisation
		} IonizationType;
		//@}

		/// A posttranslational modification
		struct PTM
		{
			/// (Simplified) name
			String name_;

			/// Formula
			EmpiricalFormula formula_;

			/// Relative abundance (in %)
			double abundance_;

			/// mass shift (positive = true, negative = false)
			bool shift_;
		};

		/// Dictionary of allowed posttranslational modifications (one residue can have several modifications)
		typedef std::multimap<String,PTM> ModTable;

		/// PTM iterator
		typedef ModTable::iterator ModTableIterator;

		/** @name Constructors and Destructors
			*/
		//@{
		/// Default constructor
		LCMSSim();

		/// Copy constructor
		LCMSSim(const LCMSSim& source);

		/// Destructor
		virtual	~LCMSSim();
		//@}

		LCMSSim& operator = (const LCMSSim& source);

		/** @name Accessors
			*/
		//@{
		/// Add a peptide sample to the pool
		void setSample(LCMSSample& sample);

		/// Set file name of SVM model for RT prediction
		void setRTModelFile(OpenMS::String RTModelFile);

		//@}

		/// Run simulation
		void run();

		/// Export spectrum data as mzXML
		void exportFeatureMap(const String& filename);

		/// Export spectrum data as mzDATA
		void exportMzData(const String& filename);

		/// Returns the set of valid LC conditions
    std::vector<String> getValidColumnConditions();

	private:

		/// Synchronize members with param class
		void updateMembers_();

		/// Create table for peptide retention times
		void predictRT_(RTTable& rtTable);

		/// samples data points from an averagine model
		void samplePeptideModel_(const ProductModel<2>& pm,
																						const CoordinateType mz_start,  const CoordinateType mz_end,
																						CoordinateType rt_start,  CoordinateType rt_end);

		/// Counts the number of strongly basic residues (Arg, Lys, His) in an amino acid sequence
		unsigned int countBasicResidues_(const AASequence& seq) const;

		/// Adds common contaminants to the map
		void addContaminants_();

		/// Remove duplicate data points (e.g. with distance <= mz_sampling) in all changed scans and sum their intensities
		void removeDuplicatePoints_();

		/// Merge all duplicate data points (e.g. with distance <= mz_sampling) and sum their intensities
		UInt removeAllDuplicatePoints_();

		/// Read allowed post-translational modifications
		void readFromModFile_();

		/// Read set of contaminants present
		void readFromContaminationFile_(std::vector<EmpiricalFormula> & vef);

		/// Sample from the set of allowed post-translational modifications
		double sampleModifications_(AASequence& aas, EmpiricalFormula& ef);

		/// Inserts a 2D peptide signal
		void insertPeptideIon_(const EmpiricalFormula& ef, const CoordinateType rt, const ChargeType c, const double ab);

		/// Add shot noise to map
		void addShotNoise_();

		/// Choose an elution profile at random
		void chooseElutionProfile_(ProductModel<2>& pm, const CoordinateType rt, const double scale);

		/// Adds a baseline to each spectrum
		void addBaseline_();

		/// Check if feature given by @p pos would overlap with already existing features
		bool checkForOverlaps_(DPosition<2> pos);

		/// The sample (e.g. collection of digested peptides)
		LCMSSample sample_;

		/// Name of the svm model file
		OpenMS::String RTModelFile_;

		/// Random number generator
		gsl_rng* rand_gen_;

		/// Data structure storing the LC-MS map
		LCMSmap exp_;

		/// Length of gradient (in seconds)
		CoordinateType gradientTime_;

		/// Sampling steps in rt (e.g. time distance between consecutive scans)
		CoordinateType rt_sampling_;

		/// Full width at half maximum of simulated peaks
		CoordinateType peak_std_;

		/// Mass accuracy in ppm
 		CoordinateType msAccuracy_;
		/// Bin size
		CoordinateType msBinSize_;

		/// Mean of peak m/z error
		CoordinateType mzMeanError_;
		/// Standard deviation of peak m/z error
		CoordinateType mzStdDevError_;

		/// Mean of peak intensity error
		IntensityType intMeanError_;
		/// Standard deviation of peak intensity error
		IntensityType intStdDevError_;

		/// Maximum m/z detected by mass analyser
		CoordinateType maxMapMZ_;
		/// Minimum m/z detected by mass analyser
		CoordinateType minMapMZ_;

		/// Mean intensity scaling
		CoordinateType mean_scaling_;

		/// Number of peptide ions
		UInt ion_count_;

		/// Allowed posttranslational modifications
		ModTable allowedMods_;

		/// The current peptide or metabolite feature
		Feature current_feature_;

		/// List of simulated features
		FeatureMap< > features_;

		/// List of simulated contaminants
		FeatureMap< > contaminants_;

		/// LC conditions (noise parameter for EMG)
		DoubleReal distortion_;

		/// Upped bound of EMG symmetry ( > 0 tailed peak, < 0 fronted peak )
		DoubleReal symmetry_up_;

		/// Lower bound of EMG symmetry ( > 0 tailed peak, < 0 fronted peak )
		DoubleReal symmetry_down_;

		/// Remembers which scans were changed after the last call to removeDuplicatePoints_()
		std::vector<bool> changed_scans_;

		/// Do we allow overlapping peptide signals?
		UInt allow_overlaps_;

	};

} // namespace OpenMS

#endif // OPENMS_SIMULATION_LCMSSIM_H
