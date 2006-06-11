// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: XMLSchemes.h,v 1.13 2006/05/30 15:46:38 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Jens Joachim $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDERS_XMLSCHEMES_H
#define OPENMS_FORMAT_HANDERS_XMLSCHEMES_H

#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
  /** @brief the XML schema namespace

		 	This namespace contains the XML schemas.
		 	Schemata are sorted descending with respect to their version

			A schema consists of an array of Strings.
			<br>
			The number of elements has to
			be consistent with the number of Maps in the Handler class.
			<br>
			The first element has to contain the String the schema is
			recognized by e.g. the name of the xsd-file or the version number.
			<br>
			All other Strings in the array are used to fill maps in the corresponding
			Handler and contain strings separated by semicolons.
			These strings correspond to values of enumerations defined in the METADATA classes.
			Empty strings are uses if no the schema does not define one (e.g. for NULL).
  */
  namespace Schemes
  {
		/// Number of available MzXML schemata
		const UnsignedInt MzXML_num = 2;

    ///Schemata for MzXML
    const String MzXML[][11] = {
			///---------------- MzXML Version 2.1 --------------------------
			{
				// name of scheme
				"mzXML_idx_2.1.xsd",

				// Polarity
      	"any;+;-",

				//Ionization
      	";ESI;EI;CI;FAB;TSP;MALDI;FD;FI;PD;SI;TI;API;ISI;CID;CAD;HN;APCI;APPI;ICP",

				// Detector
      	";EMT;Daly;;Faraday Cup;;;;Channeltron",

      	// Analyzer
      	";Quadrupole;Quadrupole Ion Trap;;;TOF;Magnetic Sector;FT-ICR;",

      	// ScanMode
      	";zoom;Full",

  	    // Attribute
      	";polarity;scanType;centroided;deisotoped;chargeDeconvoluted;"
      	"retentionTime;ionisationEnergy;collisionEnergy;cidGasPressure;"
      	"startMz;endMz;lowMz;highMz;basePeakMz;basePeakIntensity;"
      	"totIonCurrent;peaksCount;num;msLevel;scanCount;fileName;fileType;"
				"version;name;type;completionTime;precursorIntensity;precursorCharge;"
				"first;last;email;phone;URI;value;category;precision;byteOrder;"
				"pairOrder;xsi:schemaLocation;spotIntegration;intensityCutoff",

	      // Tag
      	";msRun;index;offset;sha1;parentFile;msInstrument;dataProcessing;"
      	"separation;spotting;scan;scanOrigin;precursorMz;maldi;peaks;"
      	"nameValue;comment;software;indexOffset;operator;msManufacturer;"
	      "msModel;msIonisation;msMassAnalyzer;msDetector;msResolution;"
      	"mzXML;processingOperation;separationTechnique",

				// ResolutionMethod
      	";FWHM;TenPercentValley;Baseline",

				// Peak Processing
				";1;0",

 				// Precision of Base64 encoding
				";32;64"
	  	},

	    ///---------------- MzXML Version 1.01 --------------------------
  		{
				// name of scheme file
				"MsXML.xsd",

				// Polarity
      	"any;+;-",

				//Ionization
      	";ESI;EI;CI;FAB;TSP;MALDI;FD;FI;PD;SI;TI;API;ISI;CID;CAD;HN;APCI;APPI;ICP",

				// Detector
      	";EMT;Daly;;Faraday Cup;;;;Channeltron",

      	// Analyzer
      	";Quadrupole;Ion Trap;;;TOF;Magnetic Sector;FT-ICR;",

      	// ScanMode
      	";zoom;full",

	      // Attribute
      	";polarity;scanType;centroided;deisotoped;chargeDeconvoluted;"
	      "retentionTime;ionisationEnergy;collisionEnergy;cidGasPressure;"
      	"startMz;endMz;lowMz;highMz;basePeakMz;basePeakIntensity;"
      	"totIonCurrent;peaksCount;num;msLevel;scanCount;fileName;fileType;"
				"version;name;type;completionTime;precursorIntensity;precursorCharge;"
				"first;last;email;phone;URI;value;category;precision;byteOrder;"
				"pairOrder;xsi:schemaLocation;spotIntegration;intensityCutoff",

      	// Tag
	      ";msRun;index;offset;sha1;parentFile;instrument;dataProcessing;"
      	"separation;spotting;scan;scanOrigin;precursorMz;maldi;peaks;"
      	"nameValue;comment;software;indexOffset;operator;manufacturer;"
      	"model;ionisation;msType;detector;msResolution;"
      	"mzXML;processingOperation;separationTechnique",

				// ResolutionMethod
      	";FWHM;TenPercentValley;Baseline",

				// Peak Processing
				";1;0",

 				// Precision of Base64 encoding
				";32;64"
			}
		};


		/// Number of available MzData schemata
		const UnsignedInt MzData_num = 1;

    ///Schemata for MzData
    const String MzData[][9] = {
			///---------------- MzData Version 1.05 --------------------------
			{
				// name of scheme
				"1.05",

 				// Precision of Base64 encoding
				";32;64",

				// endian-type of Base64 encoding (bigEndian, littleEndian)
				";little;big",

				// EnergyUnits
				";eV;Percent",

				// ScanMode
				";SelectedIonDetection;MassScan",

				// Polarity
				";Positive;Negative",

				// ActivationMethod
				";CID;PSD;PD;SID",

				// Ontology
				";ScanMode;Polarity;TimeInMinutes;TimeInSeconds;MassToChargeRatio;ChargeState;"
				"Intensity;IntensityUnits;Method;CollisionEnergy;EnergyUnits",

				// Tags
				";mzData;description;spectrumList;spectrum;spectrumDesc;spectrumSettings;"
				"acqSpecification;acquisition;spectrumInstrument;precursorList;ionSelection;"
				"activation;precursor;supDataDesc;supDesc;supSourceFile;data;"
				"intenArrayBinary;mzArrayBinary;cvParam;userParam;acqInstrument;acqSettings;"
				"acqDesc;cvLookup;supDataArrayBinary;supDataArray;arrayName;comments;"
				"nameOfFile;pathToFile;fileType",
	  	}
		};



		/// Number of available ExperimentalSettings schemata for MzData
		const UnsignedInt MzDataExpSett_num = 1;

    ///Schemata for the ExperimentalSettings of MzData
    const String MzDataExpSett[][18] = {
			///---------------- MzData Version 1.05 --------------------------
			{
				// name of scheme
				"1.05",

				// IonizationMode
				";PositiveIonMode;NegativeIonMode",

				// ResolutionMethod
				";FWHM;TenPercentValley;Baseline",

				// ResolutionType
				";Constant;Proportional",

				// ScanFunction
				";SelectedIonDetection;MassScan",

				// ScanDirection
				";Up;Down",

				// ScanLaw
				";Exponential;Linear;Quadratic",

				// SampleState
				";Solid;Liquid;Gas;Solution;Emulsion;Suspension",

				// PeakProcessing
				";CentroidMassSpectrum;ContinuumMassSpectrum",

				// ReflectronState
				";On;Off;None",

				// AcquisitionMode
				";PulseCounting;ADC;TDC;TransientRecorder",

				// IonizationType
				";ESI;EI;CI;FAB;TSP;LD;FD;FI;PD;SI;TI;API;ISI;CID;CAD;HN;APCI;APPI;ICP",

				// InletType
				";Direct;Batch;Chromatography;ParticleBeam;MembraneSeparator;OpenSplit;"
				"JetSeparator;Septum;Reservoir;MovingBelt;MovingWire;FlowInjectionAnalysis;"
				"ElectrosprayInlet;ThermosprayInlet;Infusion;ContinuousFlowFastAtomBombardment;"
				"InductivelyCoupledPlasma",

				// TandemScanningMethod
				";ProductIonScan;PrecursorIonScan;ConstantNeutralLoss;SingleReactionMonitoring;"
				"MultipleReactionMonitoring;SingleIonMonitoring;MultipleIonMonitoring",

				// Type
				";EM;Photomultiplier;FocalPlaneArray;FaradayCup;ConversionDynodeElectronMultiplier;"
				"ConversionDynodePhotomultiplier;Multi-Collector;ChannelElectronMultiplier",

				// AnalyzerType
				";Quadrupole;PaulIonTrap;RadialEjectionLinearIonTrap;AxialEjectionLinearIonTrap;TOF;"
				"Sector;FourierTransform;IonStorage",

				// Ontology
				";SampleNumber;SampleName;SampleState;SampleMass;SampleVolume;SampleConcentration;"
				"InletType;IonizationType;IonizationMode;AnalyzerType;Resolution;ResolutionMethod;"
				"ResolutionType;Accuracy;ScanRate;ScanTime;ScanFunction;ScanDirection;ScanLaw;"
		 		"TandemScanningMethod;ReflectronState;TOFTotalPathLength;IsolationWidth;"
				"FinalMSExponent;MagneticFieldStrength;DetectorType;DetectorAcquisitionMode;"
				"DetectorResolution;ADCSamplingFrequency;Vendor;Model;Customization;"
		 		"Deisotoping;ChargeDeconvolution;PeakProcessing",

				// Tags
				";description;admin;sampleName;sampleDescription;sourceFile;nameOfFile;"
				"pathToFile;fileType;contact;name;institution;contactInfo;instrument;instrumentName;"
				"source;detector;analyzerList;analyzer;additional;dataProcessing;software;version;"
				"processingMethod;comments;cvParam;userParam"
	  	}
		};

		/// Number of available DFeatureMap schemata
		const UnsignedInt DFeatureMap_num = 1;

    ///Schemata for the DFeatureMap
    const String DFeatureMap[][2] = {
			///---------------- DFeatureMap 1.0 --------------------------
			{
				// name of scheme
				"1.0",

				// Tags
	 			";featureList;feature;position;intensity;quality;acquisition;overallquality;"
				"charge;model;param;convexhull;hullpoint;hposition;meta;description;featureMap"
	  	}
		};

		/// Number of available DFeaturePairs schemata
		const UnsignedInt DFeaturePairs_num = 1;

    ///Schemata for the DFeatureMap
    const String DFeaturePairs[][2] = {
			///---------------- DFeaturePairs 1.0 --------------------------
			{
				// name of scheme
				"1.0",

				// Tags
				";pairlist;pair;pairquality;first;second;feature;position;intensity;quality;"
				"overallquality;charge;model;param;convexhull;hullpoint;hposition"
	  	}
		};

  }
}

#endif // OPENMS_FORMAT_HANDERS_XMLSCHEMES_H

