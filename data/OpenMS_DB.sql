SET FOREIGN_KEY_CHECKS=0;

-- 
-- Table structure for table `Acquisition`
-- 

CREATE TABLE `Acquisition` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `MetaInfo_id` bigint(20) unsigned default NULL,
  `Number` int(11) NOT NULL default '0',
  `AcquisitionInfo_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `AcquisitionInfo_id` (`AcquisitionInfo_id`),
  KEY `MetaInfo_id` (`MetaInfo_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `Acquisition`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `AcquisitionInfo`
-- 

CREATE TABLE `AcquisitionInfo` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Method_of_combination` varchar(30) collate latin1_general_ci NOT NULL default '',
  `Spectrum_type` varchar(30) collate latin1_general_ci NOT NULL default '',
  `SpectrumSettings_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `SpectrumSettings_id` (`SpectrumSettings_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `AcquisitionInfo`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `Analysis`
-- 

CREATE TABLE `Analysis` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `ClusterRun_id` bigint(20) unsigned default '0',
  `DataSet_id` bigint(20) unsigned default '0',
  `Configurable_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `ClusterRun_id` (`ClusterRun_id`),
  KEY `DataSet_id` (`DataSet_id`),
  KEY `Configurable_id` (`Configurable_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `Analysis`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `Cluster`
-- 

CREATE TABLE `Cluster` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Centroid` bigint(20) NOT NULL default '0',
  `Median` bigint(20) NOT NULL default '0',
  `DataSet_id` bigint(20) unsigned NOT NULL default '0',
  `Minmass` float NOT NULL default '0',
  `Maxmass` float NOT NULL default '0',
  `Sequence` varchar(50) collate latin1_general_ci default NULL,
  `Size` int(11) default '0',
  `Sequencecount` int(11) default '0',
  `ClusterRun_id` bigint(20) unsigned default '0',
  PRIMARY KEY  (`id`),
  KEY `ClusterRun_id` (`ClusterRun_id`),
  KEY `DataSet_id` (`DataSet_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `Cluster`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `ClusterRun`
-- 

CREATE TABLE `ClusterRun` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Configurable_id` bigint(20) unsigned NOT NULL default '0',
  `Binsize` float default '0',
  `Binspread` int(11) default '0',
  `Norm` enum('arithmetic','geometric') collate latin1_general_ci NOT NULL default 'arithmetic',
  PRIMARY KEY  (`id`),
  KEY `Configurable_id` (`Configurable_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `ClusterRun`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `Configurable`
-- 

CREATE TABLE `Configurable` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Name` varchar(50) collate latin1_general_ci NOT NULL default '',
  PRIMARY KEY  (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `Configurable`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `ConfigurableParam`
-- 

CREATE TABLE `ConfigurableParam` (
  `Configurable_id` bigint(20) unsigned NOT NULL default '0',
  `Name` varchar(50) collate latin1_general_ci NOT NULL default '',
  `Value` float NOT NULL default '0',
  KEY `configurable_id` (`Configurable_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `ConfigurableParam`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `ContactPerson`
-- 

CREATE TABLE `ContactPerson` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Contact_info` varchar(20) collate latin1_general_ci NOT NULL default '',
  `Email` varchar(40) collate latin1_general_ci NOT NULL default '',
  `Institution` varchar(50) collate latin1_general_ci NOT NULL default '',
  `MetaInfo_id` bigint(20) unsigned default NULL,
  `Name` varchar(50) collate latin1_general_ci NOT NULL default '',
  `ExperimentalSettings_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `ExperimentalSettings_id` (`ExperimentalSettings_id`),
  KEY `MetaInfo_id` (`MetaInfo_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `ContactPerson`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `DBSearch`
-- 

CREATE TABLE `DBSearch` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Date` date NOT NULL default '0000-00-00',
  `Charge` smallint(6) NOT NULL default '0',
  `Peptide_Significance_Threshold` float NOT NULL default '0',
  `Protein_Significance_Threshold` float NOT NULL default '0',
  `Spectrum_id` bigint(20) unsigned default NULL,
  `MetaInfo_id` bigint(20) unsigned default NULL,
  PRIMARY KEY  (`id`),
  KEY `Spectrum_id` (`Spectrum_id`),
  KEY `MetaInfo_id` (`MetaInfo_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `DBSearch`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `DBSearchParameters`
-- 

CREATE TABLE `DBSearchParameters` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Program` varchar(50) collate latin1_general_ci NOT NULL default '',
  `Database` varchar(50) collate latin1_general_ci default NULL,
  `Database_Version` varchar(20) collate latin1_general_ci default NULL,
  `Parameters_File` text collate latin1_general_ci,
  `Taxonomy` varchar(30) collate latin1_general_ci default NULL,
  `Fixed_Modifications` varchar(100) collate latin1_general_ci default NULL,
  `Variable_Modifications` varchar(100) collate latin1_general_ci default NULL,
  `Missed_Cleavages` smallint(6) default NULL,
  `Mass_Type` enum('monoisotopic','average') collate latin1_general_ci NOT NULL default 'monoisotopic',
  `Ion_Tolerance` float default NULL,
  `Peptide_Tolerance` float default NULL,
  `DBSearch_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `DBSearch_id` (`DBSearch_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `DBSearchParameters`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `DataSet`
-- 

CREATE TABLE `DataSet` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Name` varchar(100) collate latin1_general_ci NOT NULL default '',
  `Info` text collate latin1_general_ci,
  `MetaInfo_id` bigint(20) unsigned default NULL,
  PRIMARY KEY  (`id`),
  KEY `MetaInfo_id` (`MetaInfo_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `DataSet`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `DataSet_Member_Ref`
-- 

CREATE TABLE `DataSet_Member_Ref` (
  `DataSet_id` bigint(20) unsigned NOT NULL default '0',
  `Member_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`DataSet_id`,`Member_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `DataSet_Member_Ref`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `Digestion`
-- 

CREATE TABLE `Digestion` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Digestion_time` float NOT NULL default '0',
  `Enzyme` varchar(50) collate latin1_general_ci NOT NULL default '',
  `MetaInfo_id` bigint(20) unsigned default NULL,
  `Ph` float NOT NULL default '0',
  `Temperature` float NOT NULL default '0',
  `Type` varchar(50) collate latin1_general_ci NOT NULL default '',
  PRIMARY KEY  (`id`),
  KEY `MetaInfo_id` (`MetaInfo_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `Digestion`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `ExperimentalSettings`
-- 

CREATE TABLE `ExperimentalSettings` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Date` date NOT NULL default '0000-00-00',
  `Type` enum('UNKNOWN','MS','MS_MS','HPLC_MS','HPLC_MS_MS') collate latin1_general_ci default NULL,
  `MSExperiment_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `MSExperiment_id` (`MSExperiment_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `ExperimentalSettings`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `Gradient`
-- 

CREATE TABLE `Gradient` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `HPLC_id` bigint(20) unsigned NOT NULL default '0',
  `ExperimentalSettings_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `ExperimentalSettings_id` (`ExperimentalSettings_id`),
  KEY `HPLC_id` (`HPLC_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `Gradient`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `GradientEluent`
-- 

CREATE TABLE `GradientEluent` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Name` varchar(40) collate latin1_general_ci NOT NULL default '',
  `Gradient_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `GradientEluent`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `GradientPercentage`
-- 

CREATE TABLE `GradientPercentage` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Percentage` int(10) unsigned NOT NULL default '0',
  `GradientTime_id` bigint(20) unsigned NOT NULL default '0',
  `GradientEluent_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `GradientPercentage`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `GradientTime`
-- 

CREATE TABLE `GradientTime` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Time` int(10) unsigned NOT NULL default '0',
  `Gradient_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `GradientTime`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `HPLC`
-- 

CREATE TABLE `HPLC` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Column` varchar(40) collate latin1_general_ci NOT NULL default '',
  `Comment` text collate latin1_general_ci NOT NULL,
  `Flux` int(11) NOT NULL default '0',
  `Instrument` varchar(60) collate latin1_general_ci NOT NULL default '',
  `Pressure` int(11) NOT NULL default '0',
  `Temperature` int(11) NOT NULL default '0',
  `ExperimentalSettings_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `ExperimentalSettings_id` (`ExperimentalSettings_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `HPLC`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `Instrument`
-- 

CREATE TABLE `Instrument` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Customizations` text collate latin1_general_ci,
  `MetaInfo_id` bigint(20) unsigned default NULL,
  `Model` varchar(30) collate latin1_general_ci NOT NULL default '',
  `Name` varchar(30) collate latin1_general_ci NOT NULL default '',
  `Vendor` varchar(30) collate latin1_general_ci NOT NULL default '',
  `ExperimentalSettings_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `MetaInfo_id` (`MetaInfo_id`),
  KEY `ExperimentalSettings_id` (`ExperimentalSettings_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `Instrument`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `InstrumentSettings`
-- 

CREATE TABLE `InstrumentSettings` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Mz_range_start` float NOT NULL default '0',
  `Mz_range_stop` float NOT NULL default '0',
  `Polarity` enum('POLNULL','POSITIVE','NEGATIVE') collate latin1_general_ci default 'POLNULL',
  `Scan_mode` enum('SCANMODENULL','SELECTEDIONDETECTION','MASSSCAN') collate latin1_general_ci default 'SCANMODENULL',
  `MetaInfo_id` bigint(20) unsigned default NULL,
  `SpectrumSettings_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `MetaInfo_id` (`MetaInfo_id`),
  KEY `SpectrumSettings_id` (`SpectrumSettings_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `InstrumentSettings`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `IonDetector`
-- 

CREATE TABLE `IonDetector` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Acquisition_mode` enum('ACQMODENULL','PULSECOUNTING','ADC','TDC','TRANSIENTRECORDER') collate latin1_general_ci default 'ACQMODENULL',
  `ADC_sampling_frequency` float NOT NULL default '0',
  `Resolution` float NOT NULL default '0',
  `Type` enum('TYPENULL','ELECTRONMULTIPLIER','PHOTOMULTIPLIER','FOCALPLANEARRAY','FARADAYCUP','CONVERSIONDYNODEELECTRONMULTIPLIER','CONVERSIONDYNODEPHOTOMULTIPLIER','MULTICOLLECTOR','CHANNELELECTRONMULTIPLIER') collate latin1_general_ci default 'TYPENULL',
  `Instrument_id` bigint(20) unsigned NOT NULL default '0',
  `MetaInfo_id` bigint(20) unsigned default NULL,
  PRIMARY KEY  (`id`),
  KEY `Instrument_id` (`Instrument_id`),
  KEY `MetaInfo_id` (`MetaInfo_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `IonDetector`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `IonSource`
-- 

CREATE TABLE `IonSource` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Inlet_type` enum('INLETNULL','DIRECT','BATCH','CHROMATOGRAPHY','PARTICLEBEAM','MEMBRANESEPARATOR','OPENSPLIT','JETSEPARATOR','SEPTUM','RESERVOIR','MOVINGBELT','MOVINGWIRE','FLOWINJECTIONANALYSIS','ELECTROSPRAYINLET','THERMOSPRAYINLET','INFUSION','CONTINUOUSFLOWFASTATOMBOMBARDMENT','INDUCTIVELYCOUPLEDPLASMA') collate latin1_general_ci default 'INLETNULL',
  `Ionization_method` enum('IONMETHODNULL','ESI','EI','CI','FAB','TSP','LD','FD','FI','PD','SI','TI','API','ISI','CID','CAD','HN','APCI','APPI','ICP') collate latin1_general_ci default 'IONMETHODNULL',
  `Ionization_mode` enum('IONMODENULL','POSITIVEIONMODE','NEGATIVEIONMODE') collate latin1_general_ci default 'IONMODENULL',
  `MetaInfo_id` bigint(20) unsigned default NULL,
  `Instrument_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `MetaInfo_id` (`MetaInfo_id`),
  KEY `Instrument_id` (`Instrument_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `IonSource`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `MSExperiment`
-- 

CREATE TABLE `MSExperiment` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Description` varchar(100) collate latin1_general_ci default NULL,
  `Parameters_File` text collate latin1_general_ci,
  `MetaInfo_id` bigint(20) unsigned default NULL,
  PRIMARY KEY  (`id`),
  KEY `MetaInfo_id` (`MetaInfo_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `MSExperiment`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `MassAnalyzer`
-- 

CREATE TABLE `MassAnalyzer` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Accuracy` float NOT NULL default '0',
  `Final_MS_exponent` int(11) NOT NULL default '0',
  `Isolation_width` float NOT NULL default '0',
  `Magnetic_field_strength` float NOT NULL default '0',
  `Reflectron_state` enum('REFLSTATENULL','ON','OFF','NONE') collate latin1_general_ci default 'REFLSTATENULL',
  `Resolution` float NOT NULL default '0',
  `Resolution_method` enum('RESMETHNULL','FWHM','TENPERCENTVALLEY','BASELINE') collate latin1_general_ci default 'RESMETHNULL',
  `Resolution_type` enum('RESTYPENULL','CONSTANT','PROPORTIONAL') collate latin1_general_ci default 'RESTYPENULL',
  `Scan_direction` enum('SCANDIRNULL','UP','DOWN') collate latin1_general_ci default 'SCANDIRNULL',
  `Scan_function` enum('SCANFCTNULL','SELECTEDIONDETECTION','MASSSCAN') collate latin1_general_ci default 'SCANFCTNULL',
  `Scan_law` enum('SCANLAWNULL','EXPONENTIAL','LINEAR','QUADRATIC') collate latin1_general_ci default 'SCANLAWNULL',
  `Scan_rate` float NOT NULL default '0',
  `Scan_time` float NOT NULL default '0',
  `Tandem_scanning_method` enum('TANDEMNULL','PRODUCTIONSCAN','PRECURSORIONSCAN','CONSTANTNEUTRALLOSS','SINGLEREACTIONMONITORING','MULTIPLEREACTIONMONITORING','SINGLEIONMONITORING','MULTIPLEIONMONITORING') collate latin1_general_ci default 'TANDEMNULL',
  `TOF_total_path_length` float NOT NULL default '0',
  `Type` enum('ANALYZERNULL','QUADRUPOLE','PAULIONTRAP','RADIALEJECTIONLINEARIONTRAP','AXIALEJECTIONLINEARIONTRAP','TOF','SECTOR','FOURIERTRANSFORM','IONSTORAGE') collate latin1_general_ci default 'ANALYZERNULL',
  `Instrument_id` bigint(20) unsigned NOT NULL default '0',
  `MetaInfo_id` bigint(20) unsigned default NULL,
  PRIMARY KEY  (`id`),
  KEY `Instrument_id` (`Instrument_id`),
  KEY `MetaInfo_id` (`MetaInfo_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `MassAnalyzer`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `MetaInfo`
-- 

CREATE TABLE `MetaInfo` (
  `id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `MetaInfo`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `MetaInfoDescription`
-- 

CREATE TABLE `MetaInfoDescription` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Comment` text collate latin1_general_ci,
  `Name` varchar(30) collate latin1_general_ci NOT NULL default '',
  `SourceFile_id` bigint(20) unsigned NOT NULL default '0',
  `MetaInfo_id` bigint(20) unsigned default NULL,
  `SpectrumSettings_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `SpectrumSettings_id` (`SpectrumSettings_id`),
  KEY `MetaInfo_id` (`MetaInfo_id`),
  KEY `SourceFile_id` (`SourceFile_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `MetaInfoDescription`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `Modification`
-- 

CREATE TABLE `Modification` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Affected_amino_acids` varchar(30) collate latin1_general_ci NOT NULL default '',
  `Mass` float NOT NULL default '0',
  `MetaInfo_id` bigint(20) unsigned default NULL,
  `Reagent_name` varchar(50) collate latin1_general_ci NOT NULL default '',
  `Specificity_type` enum('AA','AA_AT_CTERM','AA_AT_NTERM','CTERM','NTERM') collate latin1_general_ci NOT NULL default 'AA',
  `Type` varchar(50) collate latin1_general_ci NOT NULL default '',
  `Mass_shift` float NOT NULL default '0',
  `Variant` enum('LIGHT','HEAVY') collate latin1_general_ci NOT NULL default 'LIGHT',
  PRIMARY KEY  (`id`),
  KEY `MetaInfo_id` (`MetaInfo_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `Modification`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `Peak`
-- 

CREATE TABLE `Peak` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `mz` float unsigned NOT NULL default '0',
  `Intensity` float unsigned default NULL,
  `Spectrum_id` bigint(20) unsigned NOT NULL default '0',
  `MetaInfo_id` bigint(20) unsigned default NULL,
  PRIMARY KEY  (`id`),
  KEY `peak_peak_list` (`Spectrum_id`),
  KEY `MetaInfo_id` (`MetaInfo_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `Peak`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `PeptideHit`
-- 

CREATE TABLE `PeptideHit` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Score` float NOT NULL default '0',
  `Score_Type` enum('Mascot','Sequest','Probability','CrossCorrelation') collate latin1_general_ci NOT NULL default 'Mascot',
  `Rank` smallint(5) unsigned NOT NULL default '0',
  `Sequence` varchar(100) collate latin1_general_ci NOT NULL default '',
  `DBSearch_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`,`DBSearch_id`),
  KEY `db_search` (`DBSearch_id`),
  KEY `Score` (`Score`),
  KEY `Rank` (`Rank`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `PeptideHit`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `PeptideHit_ProteinHit_Ref`
-- 

CREATE TABLE `PeptideHit_ProteinHit_Ref` (
  `PeptideHit_id` bigint(20) unsigned NOT NULL default '0',
  `ProteinHit_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`PeptideHit_id`,`ProteinHit_id`),
  KEY `ProteinHit_id` (`ProteinHit_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `PeptideHit_ProteinHit_Ref`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `Precursor`
-- 

CREATE TABLE `Precursor` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Activation_energy` float NOT NULL default '0',
  `Activation_energy_unit` enum('UNITSNULL','EV','PERCENT') collate latin1_general_ci default 'UNITSNULL',
  `Activation_method` enum('ACTMETHNULL','CID','PSD','PD','SID') collate latin1_general_ci default 'ACTMETHNULL',
  `SpectrumSettings_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `SpectrumSettings_id` (`SpectrumSettings_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `Precursor`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `PrecursorInfo`
-- 

CREATE TABLE `PrecursorInfo` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `mz` float unsigned NOT NULL default '0',
  `Charge` smallint(6) default NULL,
  `Intensity` float unsigned default NULL,
  `Spectrum_id` bigint(20) unsigned NOT NULL default '0',
  `MetaInfo_id` bigint(20) unsigned default NULL,
  PRIMARY KEY  (`id`),
  KEY `peak_list` (`Spectrum_id`),
  KEY `MetaInfo_id` (`MetaInfo_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `PrecursorInfo`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `ProcessingMethod`
-- 

CREATE TABLE `ProcessingMethod` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Charge_deconvolution` tinyint(1) NOT NULL default '0',
  `Deisotoping` tinyint(1) NOT NULL default '0',
  `Method` enum('UNKNOWN','PEAKS','RAWDATA') collate latin1_general_ci NOT NULL default 'UNKNOWN',
  `MetaInfo_id` bigint(20) unsigned default NULL,
  `ExperimentalSettings_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `MetaInfo_id` (`MetaInfo_id`),
  KEY `ExperimentalSettings_id` (`ExperimentalSettings_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `ProcessingMethod`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `ProteinHit`
-- 

CREATE TABLE `ProteinHit` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Score` float NOT NULL default '0',
  `Score_Type` enum('Mascot','Sequest','Probability','CrossCorrelation') collate latin1_general_ci NOT NULL default 'Mascot',
  `Rank` smallint(5) unsigned NOT NULL default '0',
  `Accession` varchar(20) collate latin1_general_ci NOT NULL default '0',
  `Accession_Type` enum('SwissProt','RefSeq') collate latin1_general_ci NOT NULL default 'SwissProt',
  `Sequence` varchar(100) collate latin1_general_ci NOT NULL default '',
  `DBSearch_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`,`DBSearch_id`),
  KEY `db_search` (`DBSearch_id`),
  KEY `Score` (`Score`),
  KEY `Rank` (`Rank`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `ProteinHit`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `Sample`
-- 

CREATE TABLE `Sample` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Comment` text collate latin1_general_ci,
  `Concentration` float NOT NULL default '0',
  `Mass` float NOT NULL default '0',
  `Name` varchar(30) collate latin1_general_ci default NULL,
  `Number` varchar(30) collate latin1_general_ci default NULL,
  `Organsim` varchar(40) collate latin1_general_ci default NULL,
  `State` enum('SAMPLENULL','SOLID','LIQUID','GAS','SOLUTION','EULSION','SUSPENSION') collate latin1_general_ci default 'SAMPLENULL',
  `Parent_Sample_id` bigint(20) unsigned NOT NULL default '0',
  `Volume` float NOT NULL default '0',
  `Experimental_settings_id` bigint(20) unsigned NOT NULL default '0',
  `MetaInfo_id` bigint(20) unsigned default NULL,
  PRIMARY KEY  (`id`),
  KEY `Experimental_settings_id` (`Experimental_settings_id`),
  KEY `MetaInfo_id` (`MetaInfo_id`),
  KEY `Parent_Sample_id` (`Parent_Sample_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `Sample`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `SampleTreatment`
-- 

CREATE TABLE `SampleTreatment` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `MetaInfo_id` bigint(20) unsigned default NULL,
  `Type` varchar(30) collate latin1_general_ci NOT NULL default '',
  `Sample_id` bigint(20) unsigned NOT NULL default '0',
  `Digestion_id` bigint(20) unsigned NOT NULL default '0',
  `Modification_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `Sample_id` (`Sample_id`),
  KEY `MetaInfo_id` (`MetaInfo_id`),
  KEY `Digestion_id` (`Digestion_id`),
  KEY `Modification_id` (`Modification_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `SampleTreatment`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `Software`
-- 

CREATE TABLE `Software` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Comment` text collate latin1_general_ci NOT NULL,
  `Completion_time` float NOT NULL default '0',
  `Name` varchar(50) collate latin1_general_ci NOT NULL default '',
  `Version` varchar(50) collate latin1_general_ci NOT NULL default '',
  `ExperimentalSettings_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `SoftwareRefsExperimentalSettings` (`ExperimentalSettings_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `Software`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `SourceFile`
-- 

CREATE TABLE `SourceFile` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `File_type` varchar(20) collate latin1_general_ci NOT NULL default '',
  `Name_of_file` varchar(50) collate latin1_general_ci NOT NULL default '',
  `Path_to_file` varchar(80) collate latin1_general_ci NOT NULL default '',
  `ExperimentalSettings_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `ExperimentalSettings_id` (`ExperimentalSettings_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `SourceFile`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `Spectrum`
-- 

CREATE TABLE `Spectrum` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `MS_Level` smallint(5) unsigned NOT NULL default '0',
  `Description` varchar(100) collate latin1_general_ci NOT NULL default '',
  `Mass_Type` enum('monoisotopic','average') collate latin1_general_ci NOT NULL default 'monoisotopic',
  `Retention_Start` float unsigned default NULL,
  `Retention_Stop` float unsigned default NULL,
  `Retention_Time` double unsigned default NULL,
  `Total_Ion_Current` int(11) unsigned default NULL,
  `Parent_Spectrum_id` bigint(20) unsigned default NULL,
  `MSExperiment_id` bigint(20) unsigned default NULL,
  `MetaInfo_id` bigint(20) unsigned default NULL,
  PRIMARY KEY  (`id`),
  KEY `Parent_Spectrum_id` (`Parent_Spectrum_id`),
  KEY `MSExperiment_id` (`MSExperiment_id`),
  KEY `MetaInfo_id` (`MetaInfo_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `Spectrum`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `SpectrumQuality`
-- 

CREATE TABLE `SpectrumQuality` (
  `Spectrum_id` bigint(20) unsigned NOT NULL default '0',
  `Score` float NOT NULL default '0',
  `Type` varchar(20) collate latin1_general_ci NOT NULL default '',
  `Info` varchar(20) collate latin1_general_ci default NULL,
  PRIMARY KEY  (`Spectrum_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `SpectrumQuality`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `SpectrumSettings`
-- 

CREATE TABLE `SpectrumSettings` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `Comment` text collate latin1_general_ci,
  `Type` enum('UNKNOWN','PEAKS','RAWDATA') collate latin1_general_ci default 'UNKNOWN',
  `Spectrum_id` bigint(20) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `Spectrum_id` (`Spectrum_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `SpectrumSettings`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `TypeNameValue`
-- 

CREATE TABLE `TypeNameValue` (
  `MetaInfo_id` bigint(20) unsigned NOT NULL default '0',
  `Type` enum('string','double','int') collate latin1_general_ci NOT NULL default 'string',
  `Name` varchar(30) collate latin1_general_ci NOT NULL default '',
  `Value` text collate latin1_general_ci NOT NULL,
  PRIMARY KEY  (`MetaInfo_id`,`Name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `TypeNameValue`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `_HighID`
-- 

CREATE TABLE `_HighID` (
  `id` int(10) unsigned NOT NULL auto_increment,
  PRIMARY KEY  (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `_HighID`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `_ID_Table`
-- 

CREATE TABLE `_ID_Table` (
  `id` bigint(20) unsigned NOT NULL default '0',
  `_Tables_id` int(11) unsigned NOT NULL default '0',
  PRIMARY KEY  (`id`),
  KEY `_tables_id` (`_Tables_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `_ID_Table`
-- 


-- --------------------------------------------------------

-- 
-- Table structure for table `_Tables`
-- 

CREATE TABLE `_Tables` (
  `id` int(11) unsigned NOT NULL auto_increment,
  `table_name` varchar(100) collate latin1_general_ci NOT NULL default '',
  PRIMARY KEY  (`id`),
  UNIQUE KEY `table` (`table_name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci;

-- 
-- Dumping data for table `_Tables`
-- 

INSERT INTO `_Tables` VALUES (1, 'MSExperiment');
INSERT INTO `_Tables` VALUES (2, 'Spectrum');
INSERT INTO `_Tables` VALUES (3, 'Peak');
INSERT INTO `_Tables` VALUES (4, 'Presursor');

-- 
-- Constraints for dumped tables
-- 

-- 
-- Constraints for table `Acquisition`
-- 
ALTER TABLE `Acquisition`
  ADD CONSTRAINT `Acquisition_ibfk_2` FOREIGN KEY (`MetaInfo_id`) REFERENCES `MetaInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Acquisition_ibfk_1` FOREIGN KEY (`AcquisitionInfo_id`) REFERENCES `AcquisitionInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `AcquisitionInfo`
-- 
ALTER TABLE `AcquisitionInfo`
  ADD CONSTRAINT `AcquisitionInfo_ibfk_1` FOREIGN KEY (`SpectrumSettings_id`) REFERENCES `SpectrumSettings` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `Analysis`
-- 
ALTER TABLE `Analysis`
  ADD CONSTRAINT `Analysis_ibfk_3` FOREIGN KEY (`Configurable_id`) REFERENCES `Configurable` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Analysis_ibfk_1` FOREIGN KEY (`ClusterRun_id`) REFERENCES `ClusterRun` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Analysis_ibfk_2` FOREIGN KEY (`DataSet_id`) REFERENCES `DataSet` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `Cluster`
-- 
ALTER TABLE `Cluster`
  ADD CONSTRAINT `Cluster_ibfk_1` FOREIGN KEY (`ClusterRun_id`) REFERENCES `ClusterRun` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Cluster_ibfk_2` FOREIGN KEY (`DataSet_id`) REFERENCES `DataSet` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `ClusterRun`
-- 
ALTER TABLE `ClusterRun`
  ADD CONSTRAINT `ClusterRun_ibfk_1` FOREIGN KEY (`Configurable_id`) REFERENCES `Configurable` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `ConfigurableParam`
-- 
ALTER TABLE `ConfigurableParam`
  ADD CONSTRAINT `ConfigurableParam_ibfk_1` FOREIGN KEY (`Configurable_id`) REFERENCES `Configurable` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `ContactPerson`
-- 
ALTER TABLE `ContactPerson`
  ADD CONSTRAINT `ContactPerson_ibfk_2` FOREIGN KEY (`ExperimentalSettings_id`) REFERENCES `ExperimentalSettings` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `ContactPerson_ibfk_1` FOREIGN KEY (`MetaInfo_id`) REFERENCES `MetaInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `DBSearch`
-- 
ALTER TABLE `DBSearch`
  ADD CONSTRAINT `DBSearch_ibfk_2` FOREIGN KEY (`MetaInfo_id`) REFERENCES `MetaInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `DBSearch_ibfk_1` FOREIGN KEY (`Spectrum_id`) REFERENCES `Spectrum` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `DBSearchParameters`
-- 
ALTER TABLE `DBSearchParameters`
  ADD CONSTRAINT `DBSearchParameters_ibfk_1` FOREIGN KEY (`DBSearch_id`) REFERENCES `DBSearch` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `DataSet`
-- 
ALTER TABLE `DataSet`
  ADD CONSTRAINT `DataSet_ibfk_1` FOREIGN KEY (`MetaInfo_id`) REFERENCES `MetaInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `DataSet_Member_Ref`
-- 
ALTER TABLE `DataSet_Member_Ref`
  ADD CONSTRAINT `DataSet_Member_Ref_ibfk_1` FOREIGN KEY (`DataSet_id`) REFERENCES `DataSet` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `Digestion`
-- 
ALTER TABLE `Digestion`
  ADD CONSTRAINT `Digestion_ibfk_1` FOREIGN KEY (`MetaInfo_id`) REFERENCES `MetaInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `ExperimentalSettings`
-- 
ALTER TABLE `ExperimentalSettings`
  ADD CONSTRAINT `ExperimentalSettings_ibfk_1` FOREIGN KEY (`MSExperiment_id`) REFERENCES `MSExperiment` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `Gradient`
-- 
ALTER TABLE `Gradient`
  ADD CONSTRAINT `Gradient_ibfk_2` FOREIGN KEY (`ExperimentalSettings_id`) REFERENCES `ExperimentalSettings` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Gradient_ibfk_1` FOREIGN KEY (`HPLC_id`) REFERENCES `HPLC` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `HPLC`
-- 
ALTER TABLE `HPLC`
  ADD CONSTRAINT `HPLC_ibfk_1` FOREIGN KEY (`ExperimentalSettings_id`) REFERENCES `ExperimentalSettings` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `Instrument`
-- 
ALTER TABLE `Instrument`
  ADD CONSTRAINT `Instrument_ibfk_1` FOREIGN KEY (`ExperimentalSettings_id`) REFERENCES `ExperimentalSettings` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `InstrumentSettings`
-- 
ALTER TABLE `InstrumentSettings`
  ADD CONSTRAINT `InstrumentSettings_ibfk_2` FOREIGN KEY (`SpectrumSettings_id`) REFERENCES `SpectrumSettings` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `InstrumentSettings_ibfk_1` FOREIGN KEY (`MetaInfo_id`) REFERENCES `MetaInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `IonDetector`
-- 
ALTER TABLE `IonDetector`
  ADD CONSTRAINT `IonDetector_ibfk_2` FOREIGN KEY (`MetaInfo_id`) REFERENCES `MetaInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `IonDetector_ibfk_1` FOREIGN KEY (`Instrument_id`) REFERENCES `Instrument` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `IonSource`
-- 
ALTER TABLE `IonSource`
  ADD CONSTRAINT `IonSource_ibfk_2` FOREIGN KEY (`Instrument_id`) REFERENCES `Instrument` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `IonSource_ibfk_1` FOREIGN KEY (`MetaInfo_id`) REFERENCES `MetaInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `MSExperiment`
-- 
ALTER TABLE `MSExperiment`
  ADD CONSTRAINT `MSExperiment_ibfk_1` FOREIGN KEY (`MetaInfo_id`) REFERENCES `MetaInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `MassAnalyzer`
-- 
ALTER TABLE `MassAnalyzer`
  ADD CONSTRAINT `MassAnalyzer_ibfk_2` FOREIGN KEY (`MetaInfo_id`) REFERENCES `MetaInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `MassAnalyzer_ibfk_1` FOREIGN KEY (`Instrument_id`) REFERENCES `Instrument` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `MetaInfoDescription`
-- 
ALTER TABLE `MetaInfoDescription`
  ADD CONSTRAINT `MetaInfoDescription_ibfk_3` FOREIGN KEY (`SpectrumSettings_id`) REFERENCES `SpectrumSettings` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `MetaInfoDescription_ibfk_1` FOREIGN KEY (`SourceFile_id`) REFERENCES `SourceFile` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `MetaInfoDescription_ibfk_2` FOREIGN KEY (`MetaInfo_id`) REFERENCES `MetaInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `Modification`
-- 
ALTER TABLE `Modification`
  ADD CONSTRAINT `Modification_ibfk_1` FOREIGN KEY (`MetaInfo_id`) REFERENCES `MetaInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `Peak`
-- 
ALTER TABLE `Peak`
  ADD CONSTRAINT `Peak_ibfk_1` FOREIGN KEY (`Spectrum_id`) REFERENCES `Spectrum` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Peak_ibfk_2` FOREIGN KEY (`MetaInfo_id`) REFERENCES `MetaInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `PeptideHit`
-- 
ALTER TABLE `PeptideHit`
  ADD CONSTRAINT `PeptideHit_ibfk_1` FOREIGN KEY (`DBSearch_id`) REFERENCES `DBSearch` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `PeptideHit_ProteinHit_Ref`
-- 
ALTER TABLE `PeptideHit_ProteinHit_Ref`
  ADD CONSTRAINT `PeptideHit_ProteinHit_Ref_ibfk_2` FOREIGN KEY (`ProteinHit_id`) REFERENCES `ProteinHit` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `PeptideHit_ProteinHit_Ref_ibfk_1` FOREIGN KEY (`PeptideHit_id`) REFERENCES `PeptideHit` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `Precursor`
-- 
ALTER TABLE `Precursor`
  ADD CONSTRAINT `Precursor_ibfk_1` FOREIGN KEY (`SpectrumSettings_id`) REFERENCES `SpectrumSettings` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `PrecursorInfo`
-- 
ALTER TABLE `PrecursorInfo`
  ADD CONSTRAINT `PrecursorInfo_ibfk_1` FOREIGN KEY (`Spectrum_id`) REFERENCES `Spectrum` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `PrecursorInfo_ibfk_2` FOREIGN KEY (`MetaInfo_id`) REFERENCES `MetaInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `ProcessingMethod`
-- 
ALTER TABLE `ProcessingMethod`
  ADD CONSTRAINT `ProcessingMethod_ibfk_2` FOREIGN KEY (`ExperimentalSettings_id`) REFERENCES `ExperimentalSettings` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `ProcessingMethod_ibfk_1` FOREIGN KEY (`MetaInfo_id`) REFERENCES `MetaInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `ProteinHit`
-- 
ALTER TABLE `ProteinHit`
  ADD CONSTRAINT `ProteinHit_ibfk_1` FOREIGN KEY (`DBSearch_id`) REFERENCES `DBSearch` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `Sample`
-- 
ALTER TABLE `Sample`
  ADD CONSTRAINT `Sample_ibfk_3` FOREIGN KEY (`Parent_Sample_id`) REFERENCES `Sample` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Sample_ibfk_1` FOREIGN KEY (`Experimental_settings_id`) REFERENCES `ExperimentalSettings` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Sample_ibfk_2` FOREIGN KEY (`MetaInfo_id`) REFERENCES `MetaInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `SampleTreatment`
-- 
ALTER TABLE `SampleTreatment`
  ADD CONSTRAINT `SampleTreatment_ibfk_4` FOREIGN KEY (`Modification_id`) REFERENCES `Modification` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `SampleTreatment_ibfk_1` FOREIGN KEY (`MetaInfo_id`) REFERENCES `MetaInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `SampleTreatment_ibfk_2` FOREIGN KEY (`Sample_id`) REFERENCES `Sample` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `SampleTreatment_ibfk_3` FOREIGN KEY (`Digestion_id`) REFERENCES `Digestion` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `Software`
-- 
ALTER TABLE `Software`
  ADD CONSTRAINT `SoftwareRefsExperimentalSettings` FOREIGN KEY (`ExperimentalSettings_id`) REFERENCES `ExperimentalSettings` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `SourceFile`
-- 
ALTER TABLE `SourceFile`
  ADD CONSTRAINT `SourceFile_ibfk_1` FOREIGN KEY (`ExperimentalSettings_id`) REFERENCES `ExperimentalSettings` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `Spectrum`
-- 
ALTER TABLE `Spectrum`
  ADD CONSTRAINT `Spectrum_ibfk_2` FOREIGN KEY (`MSExperiment_id`) REFERENCES `MSExperiment` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Spectrum_ibfk_3` FOREIGN KEY (`Parent_Spectrum_id`) REFERENCES `Spectrum` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Spectrum_ibfk_4` FOREIGN KEY (`MetaInfo_id`) REFERENCES `MetaInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `SpectrumQuality`
-- 
ALTER TABLE `SpectrumQuality`
  ADD CONSTRAINT `SpectrumQuality_ibfk_1` FOREIGN KEY (`Spectrum_id`) REFERENCES `Spectrum` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `SpectrumSettings`
-- 
ALTER TABLE `SpectrumSettings`
  ADD CONSTRAINT `SpectrumSettings_ibfk_1` FOREIGN KEY (`Spectrum_id`) REFERENCES `Spectrum` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

-- 
-- Constraints for table `TypeNameValue`
-- 
ALTER TABLE `TypeNameValue`
  ADD CONSTRAINT `TypeNameValue_ibfk_1` FOREIGN KEY (`MetaInfo_id`) REFERENCES `MetaInfo` (`id`) ON DELETE CASCADE ON UPDATE CASCADE;

SET FOREIGN_KEY_CHECKS=1;