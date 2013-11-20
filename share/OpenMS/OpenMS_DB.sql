SET FOREIGN_KEY_CHECKS=0;

-- --------------------------------------------------------

--
-- Table structure for table `ADMIN_Version`
--

DROP TABLE IF EXISTS ADMIN_Version;

CREATE TABLE ADMIN_Version (
  version varchar(50) NOT NULL default ''
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `META_MSExperiment`
--

DROP TABLE IF EXISTS META_MSExperiment;

CREATE TABLE META_MSExperiment (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_MetaInfo bigint(20) unsigned default NULL,
  `Date` date default NULL,
  Description text collate latin1_general_ci NOT NULL,
  FractionIdentifier text collate latin1_general_ci NOT NULL,
  PRIMARY KEY  (id),
  KEY fid_MetaInfo (fid_MetaInfo)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- --------------------------------------------------------

--
-- Table structure for table `DATA_Spectrum`
--

DROP TABLE IF EXISTS DATA_Spectrum;

CREATE TABLE DATA_Spectrum (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_MSExperiment bigint(20) unsigned NOT NULL default '0',
  fid_File bigint(20) unsigned default NULL,
  fid_MetaInfo bigint(20) unsigned default NULL,
  Description text character set latin1 collate latin1_general_ci NOT NULL,
  NativeID text character set latin1 collate latin1_general_ci NOT NULL,
  MSLevel smallint(5) unsigned NOT NULL default '0',
  MassType enum('monoisotopic','average') NOT NULL default 'monoisotopic',
  RetentionTime double unsigned default NULL,
  TotalIonCurrent int(11) unsigned default NULL,
  `Type` enum('Unknown','PeakData','RawData') NOT NULL default 'Unknown',
  PRIMARY KEY  (id),
  KEY fid_MSExperiment (fid_MSExperiment),
  KEY fid_MetaInfo (fid_MetaInfo),
  KEY fid_File (fid_File)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `DATA_Peak`
--

DROP TABLE IF EXISTS DATA_Peak;

CREATE TABLE DATA_Peak (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_Spectrum bigint(20) unsigned NOT NULL default '0',
  fid_MetaInfo bigint(20) unsigned default NULL,
  mz float unsigned NOT NULL default '0',
  Intensity float unsigned default NULL,
  PRIMARY KEY  (id),
  KEY peak_peak_list (fid_Spectrum),
  KEY fid_MetaInfo (fid_MetaInfo)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `META_MetaInfoDescription`
--

DROP TABLE IF EXISTS META_MetaInfoDescription;

CREATE TABLE META_MetaInfoDescription (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_Spectrum bigint(20) unsigned NOT NULL default '0',
  fid_MetaInfo bigint(20) unsigned default NULL,
  Name varchar(30) collate latin1_general_ci NOT NULL default '',
  PRIMARY KEY  (id),
  KEY fid_Spectrum (fid_Spectrum),
  KEY fid_MetaInfo (fid_MetaInfo)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- --------------------------------------------------------

--
-- Table structure for table `DATA_PeakMetaData`
--

DROP TABLE IF EXISTS DATA_PeakMetaData;

CREATE TABLE DATA_PeakMetaData (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_Peak bigint(20) unsigned NOT NULL default '0',
  fid_MetaInfoDescription bigint(20) unsigned NOT NULL default '0',
  Value float,
  PRIMARY KEY (id),
  KEY fid_Peak (fid_Peak),
  KEY fid_MetaInfoDescription (fid_MetaInfoDescription)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `DATA_Precursor`
--

DROP TABLE IF EXISTS DATA_Precursor;

CREATE TABLE DATA_Precursor (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_Spectrum bigint(20) unsigned NOT NULL default '0',
  fid_MetaInfo bigint(20) unsigned default NULL,
  WindowMz float unsigned NOT NULL default '0',
  WindowLow float unsigned default NULL,
  WindowUp float unsigned default NULL,
  Intensity float unsigned default NULL,
  Charge smallint(6) default NULL,
  ActivationEnergy float default NULL,
  PRIMARY KEY  (id),
  KEY peak_list (fid_Spectrum),
  KEY fid_MetaInfo (fid_MetaInfo)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `DATA_PrecursorPCS`
--

DROP TABLE IF EXISTS DATA_PrecursorPCS;

CREATE TABLE DATA_PrecursorPCS (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_Precursor bigint(20) unsigned NOT NULL default '0',
  PossibleChargeStates smallint(6) default NULL,
  PRIMARY KEY  (id),
  KEY fid_Precursor (fid_Precursor)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `DATA_PrecursorAM`
--

DROP TABLE IF EXISTS DATA_PrecursorAM;

CREATE TABLE DATA_PrecursorAM (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_Precursor bigint(20) unsigned NOT NULL default '0',
  ActivationMethods enum('UNKNOWN','CID','PSD','PD','SID','BIRD','ECD','IMD','SORI','HCID','LCID','PHD','ETD','PQD') NOT NULL default 'UNKNOWN',
  PRIMARY KEY  (id),
  KEY fid_Precursor (fid_Precursor)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `DATA_Products`
--

DROP TABLE IF EXISTS DATA_Products;

CREATE TABLE DATA_Products (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_Spectrum bigint(20) unsigned NOT NULL default '0',
  fid_MetaInfo bigint(20) unsigned default NULL,
  WindowMz float unsigned NOT NULL default '0',
  WindowLow float unsigned default NULL,
  WindowUp float unsigned default NULL,
  PRIMARY KEY  (id),
  KEY fid_Spectrum (fid_Spectrum),
  KEY fid_MetaInfo (fid_MetaInfo)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


-- --------------------------------------------------------

--
-- Table structure for table `ID_PeptideIdentification`
--

DROP TABLE IF EXISTS ID_PeptideIdentification;

CREATE TABLE ID_PeptideIdentification (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_Spectrum bigint(20) unsigned default NULL,
  fid_File bigint(20) unsigned default NULL,
  fid_MetaInfo bigint(20) unsigned default NULL,
  SignificanceThreshold float NOT NULL default '0',
  ScoreType varchar(100) NOT NULL default '',
  HigherScoreBetter tinyint(1) NOT NULL default '0',
  PRIMARY KEY  (id),
  KEY fid_Spectrum (fid_Spectrum),
  KEY fid_MetaInfo (fid_MetaInfo),
  KEY fid_File (fid_File)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `ID_PeptideHit`
--

DROP TABLE IF EXISTS ID_PeptideHit;

CREATE TABLE ID_PeptideHit (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_Identification bigint(20) unsigned NOT NULL default '0',
  fid_MetaInfo bigint(20) unsigned default NULL,
  Score float NOT NULL default '0',
  Sequence varchar(100) NOT NULL default '',
  charge tinyint(4) default NULL,
  AABefore char(1) default NULL,
  AAAfter char(1) default NULL,
  PRIMARY KEY  (id),
  KEY db_search (fid_Identification),
  KEY Score (Score),
  KEY fid_MetaInfo (fid_MetaInfo)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `META_InstrumentSettings`
--

DROP TABLE IF EXISTS META_InstrumentSettings;

CREATE TABLE META_InstrumentSettings (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_Spectrum bigint(20) unsigned NOT NULL default '0',
  fid_MetaInfo bigint(20) unsigned default NULL,
  Polarity enum('UNKNOWN','POSITIVE','NEGATIVE') default 'UNKNOWN',
  ScanMode enum('UNKNOWN','FULL','ZOOM','SELECTEDIONMONITORING','SELECTEDREACTIONMONITORING','CONSECUTIVEREACTIONMONITORING','CONSTANTNEUTRALGAIN','CONSTANTNEUTRALLOSS','PRECURSOR','PHOTODIODEARRAYDETECTOR','ENHANCEDMULTIPLYCHARGED','TIMEDELAYEDFRAGMENTATION') default 'UNKNOWN',
  ZoomScan tinyint NOT NULL default '0',
  PRIMARY KEY  (id),
  KEY fid_MetaInfo (fid_MetaInfo),
  KEY fid_Spectrum (fid_Spectrum)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `META_ScanWindows`
--

DROP TABLE IF EXISTS META_ScanWindows;

CREATE TABLE META_ScanWindows (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_Spectrum bigint(20) unsigned NOT NULL default '0',
  fid_MetaInfo bigint(20) unsigned default NULL,
  MZRangeBegin float default NULL,
  MzRangeEnd float default NULL,
  PRIMARY KEY  (id),
  KEY fid_MetaInfo (fid_MetaInfo),
  KEY fid_Spectrum (fid_Spectrum)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `META_SpectrumQuality`
--

DROP TABLE IF EXISTS META_SpectrumQuality;

CREATE TABLE META_SpectrumQuality (
  fid_Spectrum bigint(20) unsigned NOT NULL default '0',
  Score float NOT NULL default '0',
  `Type` enum('Misc') NOT NULL default 'Misc',
  PRIMARY KEY  (fid_Spectrum)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `META_AcquisitionInfo`
--

DROP TABLE IF EXISTS META_AcquisitionInfo;

CREATE TABLE META_AcquisitionInfo (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_Spectrum bigint(20) unsigned NOT NULL default '0',
  MethodOfCombination varchar(30) collate latin1_general_ci NOT NULL default '',
  PRIMARY KEY  (id),
  KEY fid_Spectrum (fid_Spectrum)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- --------------------------------------------------------

--
-- Table structure for table `META_Acquisition`
--

DROP TABLE IF EXISTS META_Acquisition;

CREATE TABLE META_Acquisition (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_AcquisitionInfo bigint(20) unsigned NOT NULL default '0',
  fid_MetaInfo bigint(20) unsigned default NULL,
  Number int(11) NOT NULL default '0',
  PRIMARY KEY  (id),
  KEY fid_AcquisitionInfo (fid_AcquisitionInfo),
  KEY fid_MetaInfo (fid_MetaInfo)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- --------------------------------------------------------

--
-- Table structure for table `META_DataProcessing`
--

DROP TABLE IF EXISTS META_DataProcessing;

CREATE TABLE META_DataProcessing (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_Spectrum bigint(20) unsigned NOT NULL default '0',
  fid_MetaInfo bigint(20) unsigned default NULL,
  CompletionTime varchar(50) NOT NULL default '0',
  PRIMARY KEY  (id),
  KEY fid_MetaInfo (fid_MetaInfo),
  KEY fid_Spectrum (fid_Spectrum)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
-- --------------------------------------------------------

--
-- Table structure for table `META_ProcessingActions`
--

DROP TABLE IF EXISTS META_ProcessingActions;

CREATE TABLE META_ProcessingActions (
  id bigint(20) unsigned NOT NULL auto_increment,
  ProcessingActionType enum('CHARGE_DECONVOLUTION','DEISOTOPING','SMOOTHING','CHARGE_CALCULATION','PRECURSOR_RECALCULATION','BASELINE_REDUCTION','PEAK_PICKING','PEAK_PICKING_SUM','PEAK_PICKING_MAX','FEATURE_FINDING','ALIGNMENT','LOW_INTENSITY_REMOVAL','HIGH_INTENSITY_REMOVAL','CONVERSION_MZDATA','CONVERSION_MZML','CONVERSION_MZXML','CONVERSION_DTA','CONVERSION_FEATUREXML','CONVERSION_CONSENSUSXML') NOT NULL default 'CONVERSION_MZML',
  fid_DataProcessing bigint(20) unsigned default NULL,
  PRIMARY KEY  (id),
  KEY fid_DataProcessing (fid_DataProcessing)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- --------------------------------------------------------

--
-- Table structure for table `ID_ProteinIdentification`
--

DROP TABLE IF EXISTS ID_ProteinIdentification;

CREATE TABLE ID_ProteinIdentification (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_MSExperiment bigint(20) unsigned default NULL,
  fid_File bigint(20) unsigned default NULL,
  fid_MetaInfo bigint(20) unsigned default NULL,
  SearchEngine varchar(40) NOT NULL default '',
  SearchEngineVersion varchar(10) NOT NULL default '',
  `Date` date default NULL,
  ScoreType varchar(100) NOT NULL default '',
  HigherScoreBetter tinyint(1) NOT NULL default '0',
  SignificanceThreshold float NOT NULL default '0',
  PRIMARY KEY  (id),
  KEY fid_MSExperiment (fid_MSExperiment),
  KEY fid_MetaInfo (fid_MetaInfo),
  KEY fid_File (fid_File)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `ID_ProteinHit`
--

DROP TABLE IF EXISTS ID_ProteinHit;

CREATE TABLE ID_ProteinHit (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_ProteinIdentification bigint(20) unsigned NOT NULL default '0',
  fid_MetaInfo bigint(20) unsigned default NULL,
  Score float NOT NULL default '0',
  Accession varchar(20) NOT NULL default '0',
  Sequence varchar(100) NOT NULL default '',
  Rank tinyint unsigned default NULL,
  PRIMARY KEY  (id),
  KEY db_search (fid_ProteinIdentification),
  KEY Score (Score),
  KEY fid_MetaInfo (fid_MetaInfo)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `ID_SearchParameters`
--

DROP TABLE IF EXISTS ID_SearchParameters;

CREATE TABLE ID_SearchParameters (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_ProteinIdentification bigint(20) unsigned NOT NULL default '0',
  fid_MetaInfo bigint(20) unsigned default NULL,
  DB varchar(40) NOT NULL default '',
  DBVersion varchar(10) NOT NULL default '',
  Taxonomy varchar(50) NOT NULL default '',
  Charges varchar(20) NOT NULL default '',
  MassType enum('average','monoisotopic') NOT NULL default 'monoisotopic',
  Enzyme enum('trypsin','pepsin_a','protease_k','chymotrypsin','no_enzyme','unknown_enzyme') NOT NULL default 'unknown_enzyme',
  MissedCleavages tinyint(4) NOT NULL default '0',
  PeakMassTolerance float NOT NULL default '0',
  PrecursorTolerance float NOT NULL default '0',
  PRIMARY KEY  (id),
  KEY fid_ProteinIdentification (fid_ProteinIdentification),
  KEY fid_MetaInfo (fid_MetaInfo)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `ID_VariableModifications`
--

DROP TABLE IF EXISTS ID_VariableModifications;

CREATE TABLE ID_VariableModifications (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_SearchParameters bigint(20) unsigned NOT NULL default '0',
  name varchar(50) NOT NULL default '',
  PRIMARY KEY  (id),
  KEY fid_SearchParameters (fid_SearchParameters)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `ID_FixedModifications`
--

DROP TABLE IF EXISTS ID_FixedModifications;

CREATE TABLE ID_FixedModifications (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_SearchParameters bigint(20) unsigned NOT NULL default '0',
  name varchar(50) NOT NULL default '',
  PRIMARY KEY  (id),
  KEY fid_SearchParameters (fid_SearchParameters)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `META_ContactPerson`
--

DROP TABLE IF EXISTS META_ContactPerson;

CREATE TABLE META_ContactPerson (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_MSExperiment bigint(20) unsigned NOT NULL default '0',
  fid_MetaInfo bigint(20) unsigned default NULL,
  PreName varchar(50) collate latin1_general_ci NOT NULL default '',
  LastName varchar(40) collate latin1_general_ci NOT NULL default '',
  Affiliation varchar(50) collate latin1_general_ci NOT NULL default '',
  Email varchar(40) collate latin1_general_ci NOT NULL default '',
  `Comment` text collate latin1_general_ci NOT NULL,
  PRIMARY KEY  (id),
  KEY fid_MSExperiment (fid_MSExperiment),
  KEY fid_MetaInfo (fid_MetaInfo)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- --------------------------------------------------------

--
-- Table structure for table `META_MSInstrument`
--

DROP TABLE IF EXISTS META_MSInstrument;

CREATE TABLE META_MSInstrument (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_MSExperiment bigint(20) unsigned NOT NULL default '0',
  fid_MetaInfo bigint(20) unsigned default NULL,
  Model varchar(30) collate latin1_general_ci NOT NULL default '',
  Vendor varchar(30) collate latin1_general_ci NOT NULL default '',
  Description text collate latin1_general_ci NOT NULL,
  IonOpticsType enum('UNKNOWN','MAGNETIC_DEFLECTION','DELAYED_EXTRACTION','COLLISION_QUADRUPOLE','SELECTED_ION_FLOW_TUBE','TIME_LAG_FOCUSING','REFLECTRON','EINZEL_LENS','ELECTRIC_FIELD_STRENGTH','FIRST_STABILITY_REGION','FRINGING_FIELD','KINETIC_ENERGY_ANALYZER','SPACE_CHARGE_EFFECT','STATIC_FIELD') collate latin1_general_ci NOT NULL default 'UNKNOWN',
  PRIMARY KEY  (id),
  KEY fid_MetaInfo (fid_MetaInfo),
  KEY fid_MSExperiment (fid_MSExperiment)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- --------------------------------------------------------

--
-- Table structure for table `META_IonDetector`
--

DROP TABLE IF EXISTS META_IonDetector;

CREATE TABLE META_IonDetector (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_MSInstrument bigint(20) unsigned NOT NULL default '0',
  fid_MetaInfo bigint(20) unsigned default NULL,
  AcquisitionMode enum('UNKNOWN','PULSECOUNTING','ADC','TDC','TRANSIENTRECORDER') NOT NULL default 'UNKNOWN',
  ADCSamplingFrequency float NOT NULL default '0',
  Resolution float NOT NULL default '0',
  `Type` enum('UNKNOWN','ELECTRONMULTIPLIER','PHOTOMULTIPLIER','FOCALPLANEARRAY','FARADAYCUP','CONVERSIONDYNODEELECTRONMULTIPLIER','CONVERSIONDYNODEPHOTOMULTIPLIER','MULTICOLLECTOR','CHANNELELECTRONMULTIPLIER') NOT NULL default 'UNKNOWN',
  InstrumentOrder tinyint default NULL,
  PRIMARY KEY  (id),
  KEY fid_MSInstrument (fid_MSInstrument),
  KEY fid_MetaInfo (fid_MetaInfo)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `META_IonSource`
--

DROP TABLE IF EXISTS META_IonSource;

CREATE TABLE META_IonSource (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_MSInstrument bigint(20) unsigned NOT NULL default '0',
  fid_MetaInfo bigint(20) unsigned default NULL,
  InletType enum('UNKNOWN','DIRECT','BATCH','CHROMATOGRAPHY','PARTICLEBEAM','MEMBRANESEPARATOR','OPENSPLIT','JETSEPARATOR','SEPTUM','RESERVOIR','MOVINGBELT','MOVINGWIRE','FLOWINJECTIONANALYSIS','ELECTROSPRAYINLET','THERMOSPRAYINLET','INFUSION','CONTINUOUSFLOWFASTATOMBOMBARDMENT','INDUCTIVELYCOUPLEDPLASMA') collate latin1_general_ci NOT NULL default 'UNKNOWN',
  IonizationMethod enum('UNKNOWN','ESI','EI','CI','FAB','TSP','LD','FD','FI','PD','SI','TI','API','ISI','CID','CAD','HN','APCI','APPI','ICP') collate latin1_general_ci NOT NULL default 'UNKNOWN',
  IonizationMode enum('UNKNOWN','POSITIVEIONMODE','NEGATIVEIONMODE') collate latin1_general_ci NOT NULL default 'UNKNOWN',
  InstrumentOrder tinyint default NULL,
  PRIMARY KEY  (id),
  KEY fid_MetaInfo (fid_MetaInfo),
  KEY fid_MSInstrument (fid_MSInstrument)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- --------------------------------------------------------

--
-- Table structure for table `META_MassAnalyzer`
--

DROP TABLE IF EXISTS META_MassAnalyzer;

CREATE TABLE META_MassAnalyzer (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_MSInstrument bigint(20) unsigned NOT NULL default '0',
  fid_MetaInfo bigint(20) unsigned default NULL,
  Accuracy float default NULL,
  FinalMSExponent int(11) default NULL,
  IsolationWidth float default NULL,
  MagneticFieldStrength float default NULL,
  ReflectronState enum('UNKNOWN','ON','OFF','NONE') NOT NULL default 'UNKNOWN',
  Resolution float default NULL,
  ResolutionMethod enum('UNKNOWN','FWHM','TENPERCENTVALLEY','BASELINE') NOT NULL default 'UNKNOWN',
  ResolutionType enum('UNKNOWN','CONSTANT','PROPORTIONAL') NOT NULL default 'UNKNOWN',
  ScanDirection enum('UNKNOWN','UP','DOWN') NOT NULL default 'UNKNOWN',
  ScanLaw enum('UNKNOWN','EXPONENTIAL','LINEAR','QUADRATIC') NOT NULL default 'UNKNOWN',
  ScanRate float default NULL,
  ScanTime float default NULL,
  TOFPathLength float default NULL,
  `Type` enum('UNKNOWN','QUADRUPOLE','IONTRAP','TOF','SECTOR','FOURIERTRANSFORM','IONSTORAGE') NOT NULL default 'UNKNOWN',
  InstrumentOrder tinyint default NULL,
  PRIMARY KEY  (id),
  KEY fid_MSInstrument (fid_MSInstrument),
  KEY fid_MetaInfo (fid_MetaInfo)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `META_Sample`
--

DROP TABLE IF EXISTS META_Sample;

CREATE TABLE META_Sample (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_MSExperiment bigint(20) unsigned NOT NULL default '0',
  fid_Sample bigint(20) unsigned default NULL,
  fid_MetaInfo bigint(20) unsigned default NULL,
  Name varchar(30) default NULL,
  SampleID varchar(30) default NULL,
  Mass float default NULL,
  Volume float default NULL,
  Concentration float default NULL,
  State enum('UNKNOWN','SOLID','LIQUID','GAS','SOLUTION','EULSION','SUSPENSION') NOT NULL default 'UNKNOWN',
  Organism varchar(40) default NULL,
  Description text character set latin1 collate latin1_general_ci NOT NULL,
  PRIMARY KEY  (id),
  KEY fid_MSExperiment (fid_MSExperiment),
  KEY fid_MetaInfo (fid_MetaInfo),
  KEY fid_Sample (fid_Sample)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `META_SampleTreatment`
--

DROP TABLE IF EXISTS META_SampleTreatment;

CREATE TABLE META_SampleTreatment (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_Sample bigint(20) unsigned NOT NULL default '0',
  fid_Digestion bigint(20) unsigned default NULL,
  fid_Modification bigint(20) unsigned default NULL,
  fid_MetaInfo bigint(20) unsigned default NULL,
  Description text collate latin1_general_ci NOT NULL,
  PRIMARY KEY  (id),
  KEY fid_Sample (fid_Sample),
  KEY fid_MetaInfo (fid_MetaInfo),
  KEY fid_Digestion (fid_Digestion),
  KEY fid_Modification (fid_Modification)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- --------------------------------------------------------


--
-- Table structure for table `META_Digestion`
--

DROP TABLE IF EXISTS META_Digestion;

CREATE TABLE META_Digestion (
  id bigint(20) unsigned NOT NULL auto_increment,
  Enzyme varchar(50) NOT NULL default '',
  DigestionTime float NOT NULL default '0',
  Ph float NOT NULL default '0',
  Temperature float NOT NULL default '0',
  PRIMARY KEY  (id)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `META_Modification`
--

DROP TABLE IF EXISTS META_Modification;

CREATE TABLE META_Modification (
  id bigint(20) unsigned NOT NULL auto_increment,
  ReagentName varchar(50) NOT NULL default '',
  AffectedAminoAcids varchar(30) NOT NULL default '',
  SpecificityType enum('AA','AA_AT_CTERM','AA_AT_NTERM','CTERM','NTERM') NOT NULL default 'AA',
  Mass float NOT NULL default '0',
  MassShift float default NULL,
  Variant enum('LIGHT','HEAVY') default NULL,
  PRIMARY KEY  (id)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `META_HPLC`
--

DROP TABLE IF EXISTS META_HPLC;

CREATE TABLE META_HPLC (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_MSExperiment bigint(20) unsigned NOT NULL default '0',
  InstrumentName varchar(60) collate latin1_general_ci NOT NULL default '',
  ColumnName varchar(40) collate latin1_general_ci NOT NULL default '',
  Description text collate latin1_general_ci NOT NULL,
  Flux int(11) NOT NULL default '0',
  Pressure int(11) NOT NULL default '0',
  Temperature int(11) NOT NULL default '0',
  PRIMARY KEY (id),
  KEY fid_MSExperiment (fid_MSExperiment)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- --------------------------------------------------------


--
-- Table structure for table `META_GradientEluent`
--

DROP TABLE IF EXISTS META_GradientEluent;

CREATE TABLE META_GradientEluent (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_HPLC bigint(20) unsigned NOT NULL default '0',
  Name varchar(40) collate latin1_general_ci NOT NULL default '',
  PRIMARY KEY  (id),
  KEY fid_Gradient (fid_HPLC)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- --------------------------------------------------------

--
-- Table structure for table `META_GradientTime`
--

DROP TABLE IF EXISTS META_GradientTime;

CREATE TABLE META_GradientTime (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_HPLC bigint(20) unsigned NOT NULL default '0',
  `Time` int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (id),
  KEY fid_Gradient (fid_HPLC)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- --------------------------------------------------------

--
-- Table structure for table `META_GradientPercentage`
--

DROP TABLE IF EXISTS META_GradientPercentage;

CREATE TABLE META_GradientPercentage (
  fid_GradientTime bigint(20) unsigned NOT NULL default '0',
  fid_GradientEluent bigint(20) unsigned NOT NULL default '0',
  Percentage int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (fid_GradientTime,fid_GradientEluent),
  KEY fid_GradientEluent (fid_GradientEluent)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- --------------------------------------------------------

--
-- Table structure for table `META_MetaInfo`
--

DROP TABLE IF EXISTS META_MetaInfo;

CREATE TABLE META_MetaInfo (
  id bigint(20) unsigned NOT NULL auto_increment,
  PRIMARY KEY  (id)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- --------------------------------------------------------

--
-- Table structure for table `META_Software`
--

DROP TABLE IF EXISTS META_Software;

CREATE TABLE META_Software (
  id bigint(20) unsigned NOT NULL auto_increment,
  fid_SoftwareApplicator bigint(20) unsigned NOT NULL default '0',
  fid_MetaInfo bigint(20) unsigned default NULL,
  SoftwareApplicator varchar(50) NOT NULL default '0',
  Name varchar(50) NOT NULL default '0',
  Version varchar(50) NOT NULL default '0',
  PRIMARY KEY  (id),
  KEY fid_MetaInfo (fid_MetaInfo),
  KEY fid_SoftwareApplicator (fid_SoftwareApplicator)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
-- --------------------------------------------------------

--
-- Table structure for table `META_TypeNameValue`
--

DROP TABLE IF EXISTS META_TypeNameValue;

CREATE TABLE META_TypeNameValue (
  fid_MetaInfo bigint(20) unsigned NOT NULL default '0',
  `Type` enum('string','double','int') collate latin1_general_ci NOT NULL default 'string',
  Name varchar(30) collate latin1_general_ci NOT NULL default '',
  `Value` text collate latin1_general_ci NOT NULL,
  PRIMARY KEY  (fid_MetaInfo,Name)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

--
-- Table structure for table `META_File`
--

DROP TABLE IF EXISTS META_File;

CREATE TABLE META_File (
  id bigint(20) unsigned NOT NULL auto_increment,
  FileName varchar(50) collate latin1_general_ci NOT NULL default '',
  FilePath varchar(80) collate latin1_general_ci NOT NULL default '',
  sha1 varchar(40) collate latin1_general_ci NOT NULL default '',
  Size float unsigned NOT NULL default '0',
  `Type` text collate latin1_general_ci NOT NULL,
  ChecksumType enum ('UNKNOWN', 'SHA1', 'MD5') collate latin1_general_ci NOT NULL default 'UNKNOWN',
  NativeIDType varchar(80) collate latin1_general_ci NOT NULL default '',
  PRIMARY KEY  (id)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- --------------------------------------------------------


--
-- Constraints for dumped tables
--

--
-- Constraints for table `META_MSExperiment`
--
ALTER TABLE `META_MSExperiment`
  ADD CONSTRAINT MSExperiment_ibfk_1 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL;

--
-- Constraints for table `ID_ProteinIdentification`
--
ALTER TABLE `ID_ProteinIdentification`
  ADD CONSTRAINT Identification_ibfk_2 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL,
  ADD CONSTRAINT ID_ProteinIdentification_ibfk_1 FOREIGN KEY (fid_File) REFERENCES META_File (id) ON DELETE SET NULL ON UPDATE SET NULL,
  ADD CONSTRAINT ID_ProteinIdentification_ibfk_2 FOREIGN KEY (fid_MSExperiment) REFERENCES META_MSExperiment (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `ID_ProteinHit`
--
ALTER TABLE `ID_ProteinHit`
  ADD CONSTRAINT ID_ProteinHit_ibfk_1 FOREIGN KEY (fid_ProteinIdentification) REFERENCES ID_ProteinIdentification (id) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT ID_ProteinHit_ibfk_2 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL;

--
-- Constraints for table `ID_SearchParameters`
--
ALTER TABLE `ID_SearchParameters`
  ADD CONSTRAINT ID_SearchParameters_ibfk_1 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL,
  ADD CONSTRAINT ID_SearchParameters_ibfk_2 FOREIGN KEY (fid_ProteinIdentification) REFERENCES ID_ProteinIdentification (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `ID_VariableModifications`
--
ALTER TABLE `ID_VariableModifications`
  ADD CONSTRAINT ID_SearchParameters_ibfk_3 FOREIGN KEY (fid_SearchParameters) REFERENCES ID_SearchParameters (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `ID_FixedModifications`
--
ALTER TABLE `ID_FixedModifications`
  ADD CONSTRAINT ID_SearchParameters_ibfk_4 FOREIGN KEY (fid_SearchParameters) REFERENCES ID_SearchParameters (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `DATA_Spectrum`
--
ALTER TABLE `DATA_Spectrum`
  ADD CONSTRAINT DATA_Spectrum_ibfk_1 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL,
  ADD CONSTRAINT DATA_Spectrum_ibfk_2 FOREIGN KEY (fid_File) REFERENCES META_File (id) ON DELETE SET NULL ON UPDATE SET NULL,
  ADD CONSTRAINT Spectrum_ibfk_2 FOREIGN KEY (fid_MSExperiment) REFERENCES META_MSExperiment (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `DATA_Peak`
--
ALTER TABLE `DATA_Peak`
  ADD CONSTRAINT DATA_Peak_ibfk_1 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL,
  ADD CONSTRAINT Peak_ibfk_1 FOREIGN KEY (fid_Spectrum) REFERENCES DATA_Spectrum (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `DATA_PeakMetaData`
--
ALTER TABLE `DATA_PeakMetaData`
  ADD CONSTRAINT DATA_PeakMetaData_ibfk_1 FOREIGN KEY (fid_Peak) REFERENCES DATA_Peak (id) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT DATA_PeakMetaData_ibfk_2 FOREIGN KEY (fid_MetaInfoDescription) REFERENCES META_MetaInfoDescription (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `META_MetaInfoDescription`
--
ALTER TABLE `META_MetaInfoDescription`
  ADD CONSTRAINT META_MetaInfoDescription_ibfk_1 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL,
  ADD CONSTRAINT META_MetaInfoDescription_ibfk_2 FOREIGN KEY (fid_Spectrum) REFERENCES DATA_Spectrum (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `DATA_Precursor`
--
ALTER TABLE `DATA_Precursor`
  ADD CONSTRAINT DATA_Precursor_ibfk_1 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL,
  ADD CONSTRAINT PrecursorInfo_ibfk_1 FOREIGN KEY (fid_Spectrum) REFERENCES DATA_Spectrum (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `DATA_PrecursorPCS`
--
ALTER TABLE `DATA_PrecursorPCS`
  ADD CONSTRAINT DATA_PrecursorPCS_ibfk_1 FOREIGN KEY (fid_Precursor) REFERENCES DATA_Precursor (id) ON DELETE CASCADE ON UPDATE CASCADE;

-- Constraints for table `DATA_PrecursorAM`
--
ALTER TABLE `DATA_PrecursorAM`
  ADD CONSTRAINT DATA_PrecursorAM_ibfk_1 FOREIGN KEY (fid_Precursor) REFERENCES DATA_Precursor (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `DATA_Products`
--
ALTER TABLE `DATA_Products`
  ADD CONSTRAINT DATA_Products_ibfk_1 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL,
  ADD CONSTRAINT ProductsInfo_ibfk_1 FOREIGN KEY (fid_Spectrum) REFERENCES DATA_Spectrum (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `ID_PeptideIdentification`
--
ALTER TABLE `ID_PeptideIdentification`
  ADD CONSTRAINT ID_PeptideIdentification_ibfk_1 FOREIGN KEY (fid_Spectrum) REFERENCES DATA_Spectrum (id) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT ID_PeptideIdentification_ibfk_2 FOREIGN KEY (fid_File) REFERENCES META_File (id) ON DELETE SET NULL ON UPDATE SET NULL,
  ADD CONSTRAINT ID_PeptideIdentification_ibfk_3 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL;

--
-- Constraints for table `ID_PeptideHit`
--
ALTER TABLE `ID_PeptideHit`
  ADD CONSTRAINT ID_PeptideHit_ibfk_1 FOREIGN KEY (fid_Identification) REFERENCES ID_PeptideIdentification (id) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT ID_PeptideHit_ibfk_2 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL;

--
-- Constraints for table `META_InstrumentSettings`
--
ALTER TABLE `META_InstrumentSettings`
  ADD CONSTRAINT META_InstrumentSettings_ibfk_1 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL,
  ADD CONSTRAINT META_InstrumentSettings_ibfk_2 FOREIGN KEY (fid_Spectrum) REFERENCES DATA_Spectrum (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `META_ScanWindows`
--
ALTER TABLE `META_ScanWindows`
  ADD CONSTRAINT META_ScanWindows_ibfk_1 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL,
  ADD CONSTRAINT META_ScanWindows_ibfk_2 FOREIGN KEY (fid_Spectrum) REFERENCES DATA_Spectrum (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `META_SpectrumQuality`
--
ALTER TABLE `META_SpectrumQuality`
  ADD CONSTRAINT SpectrumQuality_ibfk_1 FOREIGN KEY (fid_Spectrum) REFERENCES DATA_Spectrum (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `META_AcquisitionInfo`
--
ALTER TABLE `META_AcquisitionInfo`
  ADD CONSTRAINT META_AcquisitionInfo_ibfk_1 FOREIGN KEY (fid_Spectrum) REFERENCES DATA_Spectrum (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `META_Acquisition`
--
ALTER TABLE `META_Acquisition`
  ADD CONSTRAINT Acquisition_ibfk_1 FOREIGN KEY (fid_AcquisitionInfo) REFERENCES META_AcquisitionInfo (id) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT META_Acquisition_ibfk_1 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL;

--
-- Constraints for table `META_ContactPerson`
--
ALTER TABLE `META_ContactPerson`
  ADD CONSTRAINT META_ContactPerson_ibfk_1 FOREIGN KEY (fid_MSExperiment) REFERENCES META_MSExperiment (id) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT META_ContactPerson_ibfk_2 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL;

--
-- Constraints for table `META_MSInstrument`
--
ALTER TABLE `META_MSInstrument`
  ADD CONSTRAINT META_MSInstrument_ibfk_1 FOREIGN KEY (fid_MSExperiment) REFERENCES META_MSExperiment (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `META_IonDetector`
--
ALTER TABLE `META_IonDetector`
  ADD CONSTRAINT META_IonDetector_ibfk_1 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL,
  ADD CONSTRAINT META_IonDetector_ibfk_2 FOREIGN KEY (fid_MSInstrument) REFERENCES META_MSInstrument (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `META_IonSource`
--
ALTER TABLE `META_IonSource`
  ADD CONSTRAINT IonSource_ibfk_1 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL,
  ADD CONSTRAINT META_IonSource_ibfk_1 FOREIGN KEY (fid_MSInstrument) REFERENCES META_MSInstrument (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `META_MassAnalyzer`
--
ALTER TABLE `META_MassAnalyzer`
  ADD CONSTRAINT META_MassAnalyzer_ibfk_1 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL,
  ADD CONSTRAINT META_MassAnalyzer_ibfk_2 FOREIGN KEY (fid_MSInstrument) REFERENCES META_MSInstrument (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `META_Sample`
--
ALTER TABLE `META_Sample`
  ADD CONSTRAINT META_Sample_ibfk_1 FOREIGN KEY (fid_MSExperiment) REFERENCES META_MSExperiment (id) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT META_Sample_ibfk_2 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL,
  ADD CONSTRAINT Sample_ibfk_3 FOREIGN KEY (fid_Sample) REFERENCES META_Sample (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `META_SampleTreatment`
--
ALTER TABLE `META_SampleTreatment`
  ADD CONSTRAINT META_SampleTreatment_ibfk_1 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL,
  ADD CONSTRAINT SampleTreatment_ibfk_2 FOREIGN KEY (fid_Sample) REFERENCES META_Sample (id) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT SampleTreatment_ibfk_3 FOREIGN KEY (fid_Digestion) REFERENCES META_Digestion (id) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT SampleTreatment_ibfk_4 FOREIGN KEY (fid_Modification) REFERENCES META_Modification (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `META_Digestion`
--

-- ???

--
-- Constraints for table `META_Modification`
--

-- ???

--
-- Constraints for table `META_HPLC`
--
ALTER TABLE `META_HPLC`
  ADD CONSTRAINT META_HPLC_ibfk_1 FOREIGN KEY (fid_MSExperiment) REFERENCES META_MSExperiment (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `META_GradientEluent`
--
ALTER TABLE `META_GradientEluent`
  ADD CONSTRAINT META_GradientEluent_ibfk_1 FOREIGN KEY (fid_HPLC) REFERENCES META_HPLC (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `META_GradientPercentage`
--
ALTER TABLE `META_GradientPercentage`
  ADD CONSTRAINT META_GradientPercentage_ibfk_2 FOREIGN KEY (fid_GradientTime) REFERENCES META_GradientTime (id) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT META_GradientPercentage_ibfk_3 FOREIGN KEY (fid_GradientEluent) REFERENCES META_GradientEluent (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `META_GradientTime`
--
ALTER TABLE `META_GradientTime`
  ADD CONSTRAINT META_GradientTime_ibfk_1 FOREIGN KEY (fid_HPLC) REFERENCES META_HPLC (id) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `META_DataProcessing`
--
ALTER TABLE `META_DataProcessing`
  ADD CONSTRAINT META_DataProcessing_ibfk_1 FOREIGN KEY (fid_Spectrum) REFERENCES DATA_Spectrum (id) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT META_DataProcessing_ibfk_2 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL;
  --ADD CONSTRAINT META_DataProcessing_ibfk_3 FOREIGN KEY (fid_SoftwareMetaInfo) REFERENCES META_SoftwareMetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL,

--
-- Constraints for table `META_ProcessingActions`
--
ALTER TABLE `META_ProcessingActions`
  ADD CONSTRAINT META_ProcessingActions_ibfk_1 FOREIGN KEY (fid_DataProcessing) REFERENCES META_DataProcessing (id) ON DELETE CASCADE ON UPDATE CASCADE;


--
-- Constraints for table `META_Software`
--
ALTER TABLE `META_Software`
  ADD CONSTRAINT META_Software_ibfk_2 FOREIGN KEY (fid_MetaInfo) REFERENCES META_MetaInfo (id) ON DELETE SET NULL ON UPDATE SET NULL;

--
-- Constraints for table `META_TypeNameValue`
--

-- ???

--
-- Constraints for table `META_Info`
--

-- ???


-------------------------------------------------- DATA ----------------------------------------------------------

--
-- Dumping data for table `ADMIN_Version`
--

INSERT INTO `ADMIN_Version` (`version`) VALUES ('$Revision$');

SET FOREIGN_KEY_CHECKS=1;
