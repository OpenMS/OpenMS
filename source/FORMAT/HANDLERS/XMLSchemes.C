// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/XMLSchemes.h>

namespace OpenMS
{
	namespace Schemes
	{
		const UnsignedInt MzXML_num = 2;

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
				"pairOrder;xsi:schemaLocation;spotIntegration;intensityCutoff;"
				"startTime;endTime;fileSha1;parentFileID;precursorScanNum;"
				"windowWideness;plateID;spotID;laserShootCount;laserFrequency;"
				"laserIntesity",

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
				"pairOrder;xsi:schemaLocation;spotIntegration;intensityCutoff;"
				"startTime;endTime;fileSha1;parentFileID;precursorScanNum;"
				"windowWideness;plateID;spotID;laserShootCount;laserFrequency;"
				"laserIntesity",

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

		const UnsignedInt MzData_num = 1;

		const String MzData[][10] = {
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
				
				// Attributes
				";name;value;id;count;spectrumType;methodOfCombination;acqNumber;msLevel;"
				"mzRangeStart;mzRangeStop;supDataArrayRef;precision;endian;length;version",
			}
		};

		const UnsignedInt MzDataExpSett_num = 1;

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

		const UnsignedInt DFeatureMap_num = 1;

		const String DFeatureMap[][3] = {
			///---------------- DFeatureMap 1.0 --------------------------
			{
				// name of scheme
				"1.0",

				// Tags
				";featureList;feature;position;intensity;quality;acquisition;overallquality;"
				"charge;model;param;convexhull;hullpoint;hposition;meta;description;featureMap",
				
				// Attributes
				";dim;name;value"
			}
		};

		const UnsignedInt DFeaturePairs_num = 1;

		const String DFeaturePairs[][3] = {
			///---------------- DFeaturePairs 1.0 --------------------------
			{
				// name of scheme
				"1.0",

				// Tags
				";featurePairs;pair;pairquality;first;second;feature;position;intensity;quality;"
				"overallquality;charge;model;param;convexhull;hullpoint;hposition",
				
				// Attributes
				";dim;name;value"
			}
		};

		const UnsignedInt ConsensusXML_num = 1;

		const String ConsensusXML[][3] =
			{
				///---------------- consensusXML 1.0 --------------------------
				{
					// name of scheme
					"1.0",

					// Tags
					";consensusXML;mapList;mapType;map;alignment;alignmentMethod;"
					"matchingAlgorithm;consensusAlgorithm;alignmentNewickTree;transformationList;"
					"transformation;cell;range;parameters;consensusElementList;consensusElement;"
					"centroid;groupedElementList;element;mappinglist;rtMapping;mzMapping;param",
					
					// Attributes
					";count;name;id;rt;mz;it;rtMin;rtMax;mzMin;mzMax;itMin;itMax;map"

				}
			};

	}
}
