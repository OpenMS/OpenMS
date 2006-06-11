CREATE DATABASE `specannotate`;
USE specannotate;

CREATE TABLE `aminoacid` (
  `aminoacid_ID` int(10) unsigned NOT NULL auto_increment,
  `aminoacid_name` varchar(128) binary NOT NULL default '',
  `one_letter_code` char(1) NOT NULL default '',
  `three_letter_code` char(3) NOT NULL default '',
  `middle_formula` varchar(128) NOT NULL default '',
  `single_formula` varchar(128) NOT NULL default '',
  `c_term_formula` varchar(128) NOT NULL default '',
  `n_term_formula` varchar(128) NOT NULL default '',
  `single_mono_mass` double unsigned NOT NULL default '0',
  `middle_mono_mass` double unsigned NOT NULL default '0',
  `c_term_mono_mass` double unsigned NOT NULL default '0',
  `n_term_mono_mass` double unsigned NOT NULL default '0',
  `single_average_mass` double unsigned NOT NULL default '0',
  `middle_average_mass` double unsigned NOT NULL default '0',
  `c_term_average_mass` double unsigned NOT NULL default '0',
  `n_term_average_mass` double unsigned NOT NULL default '0',
  PRIMARY KEY  (`aminoacid_ID`),
  KEY `name` (`aminoacid_name`),
  KEY `one_letter_code` (`one_letter_code`),
  KEY `three_letter_code` (`three_letter_code`)
) TYPE=MyISAM AUTO_INCREMENT=21 ;


INSERT INTO `aminoacid` (`aminoacid_ID`, `aminoacid_name`, `one_letter_code`, `three_letter_code`, `middle_formula`, `single_formula`, `c_term_formula`, `n_term_formula`, `single_mono_mass`, `middle_mono_mass`, `c_term_mono_mass`, `n_term_mono_mass`, `single_average_mass`, `middle_average_mass`, `c_term_average_mass`, `n_term_average_mass`) VALUES (1, 0x416c616e696e65, 'A', 'ALA', 'C3H5N1O1', 'C3H7N1O2', 'C3H6N1O2', 'C3H6N1O1', 89.047668, 71.037102, 88.039841, 72.04493, 89.09433, 71.078987, 88.086357, 72.08696),
(2, 0x417267696e696e65, 'R', 'ARG', 'C6H12N4O1', 'C6H14N4O2', 'C6H13N4O2', 'C6H13N4O1', 174.111649, 156.101089, 173.103821, 157.108917, 174.203339, 156.188004, 173.195374, 157.195969),
(3, 0x41737061726167696e65, 'N', 'ASN', 'C4H6N2O2', 'C4H8N4O3', 'C4H7N2O3', 'C4H7N2O2', 160.059616, 114.042908, 131.045639, 115.050735, 160.132919, 114.104095, 131.111465, 115.112068),
(4, 0x417370617274696341636964, 'D', 'ASP', 'C4H5N1O3', 'C4H7N1O4', 'C4H6N1O4', 'C4H6N1O3', 133.037491, 115.026924, 132.029663, 116.034752, 133.104126, 115.088791, 132.096161, 116.096756),
(5, 0x4379737465696e65, 'C', 'CYS', 'C3H5N1O1S1', 'C3H7N1O2S1', 'C3H6N1O2S1', 'C3H6N1O1S1', 121.019737, 103.009178, 120.011909, 104.016998, 121.160332, 103.144989, 120.152359, 104.152962),
(6, 0x476c7574616d696341636964, 'E', 'GLU', 'C5H7N1O3', 'C5H9N1O4', 'C5H8N1O4', 'C5H8N1O3', 147.053131, 129.042572, 146.045319, 130.0504, 147.131073, 129.115723, 146.123093, 130.123703),
(7, 0x476c7574616d696e65, 'Q', 'GLN', 'C5H8N2O2', 'C5H10N2O3', 'C5H9N2O3', 'C5H9N2O2', 146.069122, 128.058563, 145.061295, 129.066391, 146.146378, 128.131042, 145.138412, 129.139008),
(8, 0x476c7963696e65, 'G', 'GLY', 'C2H3N1O1', 'C2H5N1O2', 'C2H4N1O2', 'C2H4N1O1', 75.032013, 57.021454, 74.024193, 58.029282, 75.06739, 57.052048, 74.059418, 58.06002),
(9, 0x486973746964696e65, 'H', 'HIS', 'C6H7N3O1', 'C6H9N3O2', 'C6H8N3O2', 'C6H8N3O1', 155.069458, 137.058899, 154.06163, 138.066727, 155.156754, 137.141403, 154.148773, 138.149384),
(10, 0x49736f6c657563696e65, 'I', 'ILE', 'C6H11N1O1', 'C6H13N1O2', 'C6H12N1O2', 'C6H12N1O1', 131.09462, 113.084053, 130.086792, 114.091881, 131.17514, 113.159805, 130.167175, 114.167778),
(11, 0x4c657563696e65, 'L', 'LEU', 'C6H11N1O1', 'C6H13N1O2', 'C6H12N1O2', 'C6H12N1O1', 131.09462, 113.084053, 130.086792, 114.091881, 131.17514, 113.159805, 130.167175, 114.167778),
(12, 0x4c7973696e65, 'K', 'LYS', 'C6H12N2O1', 'C6H14N2O2', 'C6H13N2O2', 'C6H13N2O1', 146.105515, 128.094955, 145.097687, 129.102768, 146.18985, 128.174515, 145.181885, 129.18248),
(13, 0x4d657468696f6e696e65, 'M', 'MET', 'C5H9N1O1S1', 'C5H11N1O2S1', 'C5H10N1O2S1', 'C5H10N1O1S1', 149.051041, 131.040482, 148.043213, 132.048294, 149.214203, 131.198868, 148.206238, 132.206833),
(14, 0x5068656e796c616c616e696e65, 'F', 'PHE', 'C9H9N1O1', 'C9H11N1O2', 'C9H10N1O2', 'C9H10N1O1', 165.078964, 147.068405, 164.071136, 148.076233, 165.1922, 147.176865, 164.184235, 148.18483),
(15, 0x50726f6c696e65, 'P', 'PRO', 'C5H7N1O1', 'C5H9N1O2', 'C5H8N1O2', 'C5H8N1O1', 115.063316, 97.052757, 114.055489, 98.060577, 115.132271, 97.116928, 114.124298, 98.124901),
(16, 0x536572696e65, 'S', 'SER', 'C3H5N1O2', 'C3H7N1O3', 'C3H6N1O3', 'C3H6N1O2', 105.042572, 87.032013, 104.034752, 88.039841, 105.093727, 87.078392, 104.085762, 88.086357),
(17, 0x546872656f6e696e65, 'T', 'THR', 'C4H7N1O2', 'C4H9N1O3', 'C4H8N1O3', 'C4H8N1O2', 119.058228, 101.047668, 118.0504, 102.055489, 119.120667, 101.105331, 118.112701, 102.113297),
(18, 0x54727970746f7068616e, 'W', 'TRP', 'C11H10N2O1', 'C11H12N2O2', 'C11H11N2O2', 'C11H11N2O1', 204.089859, 186.0793, 203.082031, 187.087128, 204.228912, 186.213577, 203.220947, 187.221542),
(19, 0x5479726f73696e65, 'Y', 'TYR', 'C9H9N1O2', 'C9H11N1O3', 'C9H10N1O3', 'C9H10N1O2', 181.073883, 163.063309, 180.066055, 164.071136, 181.191605, 163.17627, 180.18364, 164.184235),
(20, 0x56616c696e65, 'V', 'VAL', 'C5H9N1O1', 'C5H11N1O2', 'C5H10N1O2', 'C5H10N1O1', 117.078964, 99.068405, 116.071144, 100.076233, 117.148209, 99.132866, 116.140236, 100.140839);


CREATE TABLE `annotation` (
  `annotation_ID` int(10) unsigned NOT NULL auto_increment,
  `sample_ID` int(10) unsigned NOT NULL default '0',
  `mass` double unsigned NOT NULL default '0',
  `digest_fragment_ID` int(10) NOT NULL default '-1',
  `realized_modification_ID` int(10) NOT NULL default '-1',
  `realized_modification_positionless_ID` int(11) NOT NULL default '-1',
  PRIMARY KEY  (`annotation_ID`),
  KEY `mass` (`mass`),
  KEY `sample_ID` (`sample_ID`,`mass`)
) TYPE=MyISAM AUTO_INCREMENT=1 ;


CREATE TABLE `digest_fragment` (
  `digest_fragment_ID` int(10) unsigned NOT NULL auto_increment,
  `protein_ID` int(10) unsigned NOT NULL default '0',
  `enzyme_ID` int(10) NOT NULL default '-1',
  `d_start_pos` int(10) unsigned NOT NULL default '0',
  `d_end_pos` int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (`digest_fragment_ID`),
  KEY `protein_ID` (`protein_ID`,`enzyme_ID`)
) TYPE=MyISAM AUTO_INCREMENT=212 ;


CREATE TABLE `enzyme` (
  `enzyme_ID` int(10) unsigned NOT NULL auto_increment,
  `enzyme_name` varchar(128) binary NOT NULL default '',
  `cleavage_sites` varchar(128) NOT NULL default '',
  `terminality` char(1) NOT NULL default '',
  PRIMARY KEY  (`enzyme_ID`)
) TYPE=MyISAM AUTO_INCREMENT=10 ;

INSERT INTO `enzyme` (`enzyme_ID`, `enzyme_name`, `cleavage_sites`, `terminality`) VALUES (1, 0x5472797073696e, 'KR', 'C'),
(2, 0x4368796d6f7472797073696e, 'FWY', 'C'),
(3, 0x5375626d6178696c6c61727573, 'R', 'C'),
(4, 0x53746170685f6175725f7638, 'DE', 'C'),
(5, 0x4173705f6e, 'DE', 'N'),
(6, 0x50657073696e, 'FWY', 'N'),
(7, 0x42726f6d6379616e, 'M', 'C');

CREATE TABLE `modification` (
  `modification_ID` int(10) unsigned NOT NULL auto_increment,
  `modification_name` varchar(128) binary NOT NULL default '',
  `plus_formula` varchar(128) NOT NULL default '',
  `minus_formula` varchar(128) default 'H2O1',
  `plus_mono_mass` double unsigned NOT NULL default '0',
  `minus_mono_mass` double unsigned NOT NULL default '0',
  `plus_average_mass` double unsigned NOT NULL default '0',
  `minus_average_mass` double unsigned NOT NULL default '18.01534',
  `modification_sites` varchar(128) NOT NULL default '',
  `note` text NOT NULL,
  PRIMARY KEY  (`modification_ID`),
  KEY `modification_name` (`modification_name`)
) TYPE=MyISAM AUTO_INCREMENT=46 ;

INSERT INTO `modification` (`modification_ID`, `modification_name`, `plus_formula`, `minus_formula`, `plus_mono_mass`, `minus_mono_mass`, `plus_average_mass`, `minus_average_mass`, `modification_sites`, `note`) VALUES (1, 0x4f7869646174696f6e, 'O', '', 15.99491, 0, 15.994, 0, 'M', 'Methionine Oxidation caused by Oxygen from the air'),
(2, 0x416c6b796c6174696f6e, 'C2H3O1N1', '', 57.0215, 0, 57.052, 0, 'C', 'Alkylation after using Iodacetamid for masking the free -SH groups of Cysteine after reductive cleavage of S-S disulfide bonds. (during preparation for digests)'),
(3, 0x4e315f616c706861, '', 'H2O1', 0, 0, 1438.311, 18.01534, 'N', ''),
(4, 0x4e325f616c706861, '', 'H2O1', 0, 0, 1276.16, 18.01534, 'N', ''),
(5, 0x4e335f616c706861, '', 'H2O1', 0, 0, 1787.65, 18.01534, 'N', ''),
(6, 0x4e345f616c706861, '', 'H2O1', 0, 0, 2152.98, 18.01534, 'N', ''),
(7, 0x4e355f616c706861, '', 'H201', 0, 0, 2152.98, 18.01534, 'N', ''),
(8, 0x4e365f616c706861, '', 'H2O1', 0, 0, 2518.32, 18.01524, 'N', ''),
(9, 0x4e315f62657461, '', 'H2O1', 0, 0, 2370.17, 18.01534, 'N', ''),
(10, 0x4e325f62657461, '', 'H2O1', 0, 0, 2078.91, 18.01534, 'N', ''),
(11, 0x4e335f62657461, '', 'H2O1', 0, 0, 2224.02, 18.01534, 'N', ''),
(12, 0x4e345f62657461, '', 'H2O1', 0, 0, 1932.77, 18.01534, 'N', ''),
(13, 0x4e355f62657461, '', 'H2O1', 0, 0, 1567.42, 18.01534, 'N', ''),
(14, 0x4e365f62657461, '', 'H2O1', 0, 0, 2880.62, 18.01534, 'N', ''),
(15, 0x4e375f62657461, '', 'H2O1', 0, 0, 2589.36, 18.01534, 'N', ''),
(16, 0x4e385f62657461, '', 'H2O1', 0, 0, 1787.65, 18.01534, 'N', ''),
(17, 0x4e395f62657461, '', 'H2O1', 0, 0, 1641.5, 18.01534, 'N', ''),
(18, 0x4e31305f62657461, '', 'H2O1', 0, 0, 1713.57, 18.01534, 'N', ''),
(19, 0x4f315f62657461, '', 'H2O1', 0, 0, 1331.21, 18.01534, 'S', ''),
(20, 0x4f325f62657461, '', 'H2O1', 0, 0, 965.87, 18.01534, 'S', ''),
(21, 0x4f335f62657461, '', 'H2O1', 0, 0, 674.61, 18.01534, 'S', ''),
(22, 0x4f345f62657461, '', 'H2O1', 0, 0, 512.47, 18.01534, 'S', ''),
(23, 0x4f355f62657461, '', 'H2O1', 0, 0, 748.69, 18.01534, 'S', ''),
(24, 0x4f365f62657461, '', 'H2O1', 0, 0, 383.35, 18.01534, 'S', ''),
(25, 0x756e6d6f646966696564, '', '', 0, 0, 0, 0, 'ARNDCEQGHILKMFPSTWYV', ''),
(26, 0x56696e796c7079726964696e5f416c6b, '', '', 0, 0, 105.1, 0, 'C', ''),
(27, 0x4e315f616c7068615f61, '', 'H2O1', 0, 0, 1600.45, 18.01534, 'N', ''),
(28, 0x4e315f616c7068615f62, '', 'H2O1', 0, 0, 1762.6, 18.01534, 'N', ''),
(29, 0x4e335f616c7068615f61, '', 'H2O1', 0, 0, 1641.51, 18.01534, 'N', ''),
(30, 0x4e345f616c7068615f61, '', 'H2O1', 0, 0, 2006.84, 18.01534, 'N', ''),
(31, 0x4e365f616c7068615f61, '', 'H2O1', 0, 0, 2372.18, 18.01534, 'N', ''),
(32, 0x4e375f616c706861, '', 'H2O1', 0, 0, 1990.84, 18.01534, 'N', ''),
(33, 0x4e375f616c7068615f61, '', 'H2O1', 0, 0, 1844.7, 18.01534, 'N', ''),
(35, 0x566572737563683032, '', 'H2O1', 0, 0, 1316.48, 18.01534, 'N', ''),
(34, 0x566572737563683031, '', 'H2O1', 0, 0, 1462.54, 18.01534, 'N', ''),
(36, 0x566572737563683033, '', 'H2O1', 0, 0, 1259.46, 18.01534, 'N', ''),
(37, 0x566572737563683034, '', 'H2O1', 0, 0, 1113.41, 18.01534, 'N', ''),
(38, 0x566572737563683035, '', 'H2O1', 0, 0, 1056.38, 18.01534, 'N', ''),
(39, 0x566572737563683036, '', 'H2O1', 0, 0, 910.33, 18.01534, 'N', ''),
(40, 0x566572737563683037, '', 'H2O1', 0, 0, 894.33, 18.01534, 'N', ''),
(41, 0x566572737563683038, '', 'H2O1', 0, 0, 748.27, 18.01534, 'N', ''),
(42, 0x566572737563683039, '', 'H2O1', 0, 0, 732.28, 18.01534, 'N', ''),
(43, 0x566572737563683130, '', 'H2O1', 0, 0, 568.20463, 18.01534, 'N', '');

CREATE TABLE `modification_combination` (
  `modification_combination_ID` int(10) unsigned NOT NULL auto_increment,
  `first_realized_modification_ID` int(10) unsigned NOT NULL default '0',
  `next_modification_combination_ID` int(10) unsigned default NULL,
  PRIMARY KEY  (`modification_combination_ID`)
) TYPE=MyISAM AUTO_INCREMENT=1 ;

CREATE TABLE `modification_combination_positionless` (
  `modification_combination_positionless_ID` int(10) unsigned NOT NULL auto_increment,
  `first_realized_modification_positionless_ID` int(10) unsigned NOT NULL default '0',
  `next_modification_combination_positionless_ID` int(10) unsigned default NULL,
  PRIMARY KEY  (`modification_combination_positionless_ID`)
) TYPE=MyISAM AUTO_INCREMENT=1 ;

CREATE TABLE `protein` (
  `protein_ID` int(10) unsigned NOT NULL auto_increment,
  `identifier` varchar(128) NOT NULL default '',
  `pdb_filename` varchar(128) binary NOT NULL default 'void',
  `fasta_filename` varchar(128) binary NOT NULL default 'void',
  `sequence_oneletter` text NOT NULL,
  `mono_mass` double unsigned NOT NULL default '0',
  `average_mass` double unsigned NOT NULL default '0',
  `no_of_aminoacids` int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (`protein_ID`),
  KEY `pdb_filename` (`pdb_filename`),
  KEY `fasta_filename` (`fasta_filename`),
  KEY `identifier` (`identifier`)
) TYPE=MyISAM AUTO_INCREMENT=4 ;

INSERT INTO `protein` (`protein_ID`, `identifier`, `pdb_filename`, `fasta_filename`, `sequence_oneletter`, `mono_mass`, `average_mass`, `no_of_aminoacids`) VALUES (1, 'hCG_alpha(pdb1e9j)', 0x2f686f6d652f6f726c616e646f2f4449504c4f4d4152424549542f446174656e2f70726f7465696e732f6843475f616c7068615f7064623165396a2e656e74, 0x2f686f6d652f6f726c616e646f2f4449504c4f4d4152424549542f446174656e2f70726f7465696e732f6843475f616c7068615f637573746f6d2e6661737461, 'APDVQDCPECTLQENPFFSQPGAPILQCMGCCFSRAYPTPLRSKKTMLVQKNVTSESTCCVAKSYNRVTVMGGFKVENHTACHCSTCYYHKS', 10198.666099, 10205.842246, 92),
(2, 'hCG_beta(P01233)', 0x766f6964, 0x2f686f6d652f6f726c616e646f2f4449504c4f4d4152424549542f446174656e2f70726f7465696e732f6843475f626574615f637573746f6d2e6661737461, 'SKEPLRPRCRPINATLAVEKEGCPVCITVNTTICAGYCPTMTRVLQGVLPALPQVVCNYRDVRFESIRLPGCPRGVNPVVSYAVALSCQCALCRRSTTDCGGPKDHPLTCDDPRFQDSSSSKAPPPSLPSPSRLPGPSDTPILPQ', 15521.734617, 15532.09232, 145),
(3, 'hcg_beta_core', 0x766f6964, 0x2f686f6d652f6f726c616e646f2f4449504c4f4d4152424549542f446174656e2f70726f7465696e732f6843475f626574615f636f72652e6661737461, 'RPRCRPINATLAVEKEGCPVCITVNTTICAGYCPTVVCNYRDVRFESIRLPGCPRGVNPVVSYAVALSCQCAL', 7883.907494, 7889.369377, 73);

CREATE TABLE `protein_modification_scenario` (
  `protein_modification_scenario_ID` int(10) unsigned NOT NULL auto_increment,
  `protein_ID` int(10) unsigned NOT NULL default '0',
  `overall_modifications` text,
  `partial_modifications` text,
  `annotation_method` varchar(30) default NULL,
  `modification_combination_ID` int(10) unsigned NOT NULL default '0',
  `modification_combination_positionless_ID` int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (`protein_modification_scenario_ID`)
) TYPE=MyISAM AUTO_INCREMENT=1 ;

CREATE TABLE `realized_modification` (
  `realized_modification_ID` int(10) unsigned NOT NULL auto_increment,
  `m_position` int(10) unsigned NOT NULL default '0',
  `modification_ID` int(10) unsigned NOT NULL default '0',
  `next_realized_modification_ID` int(10) unsigned default NULL,
  PRIMARY KEY  (`realized_modification_ID`)
) TYPE=MyISAM AUTO_INCREMENT=1 ;

CREATE TABLE `realized_modification_positionless` (
  `realized_modification_positionless_ID` int(10) unsigned NOT NULL auto_increment,
  `modification_ID` int(10) unsigned NOT NULL default '0',
  `no_of_occurrences` int(10) unsigned NOT NULL default '0',
  `next_realized_modification_positionless_ID` int(10) unsigned default NULL,
  PRIMARY KEY  (`realized_modification_positionless_ID`)
) TYPE=MyISAM AUTO_INCREMENT=1 ;

CREATE TABLE `sample` (
  `sample_ID` int(10) unsigned NOT NULL auto_increment,
  `enzyme_ID` int(10) NOT NULL default '-1',
  `protein_modification_scenario_ID` int(10) unsigned NOT NULL default '0',
  `annotation_method` varchar(30) default NULL,
  PRIMARY KEY  (`sample_ID`)
) TYPE=MyISAM AUTO_INCREMENT=1 ;

CREATE TABLE `sequence` (
  `protein_ID` int(10) unsigned NOT NULL default '0',
  `s_position` int(10) unsigned NOT NULL default '0',
  `aminoacid_ID` int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (`protein_ID`,`s_position`)
) TYPE=MyISAM;

