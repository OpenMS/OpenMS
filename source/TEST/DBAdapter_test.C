// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/DB/DBAdapter.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/Tagging.h>
#include <OpenMS/METADATA/Modification.h>
#include <OpenMS/METADATA/Digestion.h>
#include <OpenMS/KERNEL/StandardTypes.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DBAdapter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

/* This check does the following:
   - store experiment with 2 spectra to DB
   - load spectrum 1
   - load full experiment
   - modify 1st experiment
   - store 1st experiment (overwrite)
   - load experiment
   - store empty experiment
   
   So at the end of the test there should be 2 experiments stored in the database:
   - 1 full one, 2 spectra
   - 1 empty one.
*/

// check for credentials
// if they are not present, abort the test (successfully)
bool do_tests=true;
TextFile credentials;
try
{
	credentials.load("DB_credentials.txt",true);
}
catch(...)
{
	do_tests=false;
}

String db,host,user,password,port;

//read out connection data
for (TextFile::iterator it = credentials.begin(); it!= credentials.end(); ++it)
{
	//comments and empty lines
	if (it->hasPrefix('#') || *it == "")
	{
		continue;
	}
	
	//extract connection info
	if (it->hasPrefix("Host:")) host = it->suffix(':').trim();
	if (it->hasPrefix("Port:")) port = it->suffix(':').trim();
	if (it->hasPrefix("User:")) user = it->suffix(':').trim();
	if (it->hasPrefix("Password:")) password = it->suffix(':').trim();
	if (it->hasPrefix("DB:")) db = it->suffix(':').trim();
}

if (do_tests)
{

	//DB connection for DBAdapter
	DBConnection con;
	con.connect(db, user, password, host, port.toInt());
	DBConnection con2;
	con2.connect(db, user, password, host, port.toInt(), DB_PLUGIN, "alternateConnection");
	
	DBAdapter* ptr = 0;

	START_SECTION((DBAdapter(DBConnection& db_con)))
		ptr = new DBAdapter(con);
		TEST_NOT_EQUAL(ptr, 0)
END_SECTION
		
	START_SECTION((~DBAdapter()))
		delete ptr;
END_SECTION

	START_SECTION((void createDB()))
  	DBAdapter a(con);
		a.createDB();
		
		QSqlQuery result = con.executeQuery("SELECT id FROM META_MSExperiment");
	  TEST_EQUAL(result.size(),0)
END_SECTION

	//check if the DB is up-to-date
	bool db_up_to_date;
	START_SECTION((bool checkDBVersion(bool warning)))
		DBAdapter a(con);
		db_up_to_date = a.checkDBVersion(true);
		TEST_EQUAL(db_up_to_date,true)
END_SECTION
	
	if (db_up_to_date)
	{
	
		// create test data - one experiment containing 2 spectra.
		RichPeakMap exp_original;
		exp_original.setComment("bla");
		
		exp_original.getSample().setName("fruity loops");
		exp_original.getSample().setNumber("007");
		exp_original.getSample().setMass(30.1);
		vector<Sample> subsamples;
		Sample subsample;
		subsample.setVolume(60.1);
		subsample.setConcentration(101.1);
		Digestion digestion;
		digestion.setEnzyme("dhdh");
		digestion.setDigestionTime(36.6);
		digestion.setPh(7.2);
		digestion.setTemperature(37.7);
		subsample.addTreatment(digestion);
		subsamples.push_back(subsample);
		subsample = Sample();
		subsample.setState(Sample::GAS);
		subsample.setOrganism("isistius brasiliensis (cookiecutter shark, see www.isistius.de)");
		Modification modification;
		modification.setReagentName("reagent");
		modification.setAffectedAminoAcids("123");
		modification.setSpecificityType(Modification::CTERM);
		modification.setMass(12.3);
		subsample.addTreatment(modification);
		Tagging tagging;
		tagging.setReagentName("tagging");
		tagging.setMassShift(0.123);
		tagging.setVariant(Tagging::HEAVY);
		subsample.addTreatment(tagging);
		subsamples.push_back(subsample);
		subsample = Sample();
		subsample.setComment("nice");
		subsample.setMetaValue("label", String("pink"));
		subsamples.push_back(subsample);
		exp_original.getSample().setSubsamples(subsamples);

		// setting experiment's first protein identification (+ 2 proteine hits)
		ProteinIdentification pi;
		ProteinHit ph;
		pi.setSearchEngine("google");
		pi.setSearchEngineVersion("beta");
		DateTime datetime;
		// does not save time yet, DB schema must be changed from Date to DateTime
		datetime.set("2006-12-12 00:00:00");
		pi.setDateTime(datetime);
		pi.setScoreType("Type");
		pi.setHigherScoreBetter(true);
		pi.setSignificanceThreshold(3.456);
		ph.setAccession("0110110");
		std::vector<ProteinHit> vector_ph;
		vector_ph.push_back(ph);
		ph = ProteinHit();
		ph.setScore(4.567);
		ph.setAccession("1001001");
		ph.setSequence("ZXY");
		vector_ph.push_back(ph);
		pi.setHits(vector_ph);	
		exp_original.getProteinIdentifications().push_back(pi);

		// setting experiment's second protein identification (+ no proteine hits)
		pi = ProteinIdentification();
		pi.setHigherScoreBetter(false);
		exp_original.getProteinIdentifications().push_back(pi);
		
		ContactPerson contact;
		contact.setFirstName("Ferdinand");
		contact.setLastName("Piech");
		contact.setInstitution("aff");
		exp_original.getContacts().push_back(contact);
		contact = ContactPerson();	
		contact.setEmail("ferdi@porsche.de");
		contact.setContactInfo("ttss");
		contact.setMetaValue("label", String("polka-dotted"));
		exp_original.getContacts().push_back(contact);
		
		exp_original.getHPLC().setInstrument("guitar");
		exp_original.getHPLC().setColumn("bigone");
		exp_original.getHPLC().setComment("fhth");
		exp_original.getHPLC().setFlux(1);
		exp_original.getHPLC().setPressure(2);
		exp_original.getHPLC().setTemperature(3);
		
		exp_original.getHPLC().getGradient().addEluent("C2H5OH");
		exp_original.getHPLC().getGradient().addEluent("H2O");
		exp_original.getHPLC().getGradient().addTimepoint(1);
		exp_original.getHPLC().getGradient().addTimepoint(5);
		exp_original.getHPLC().getGradient().addTimepoint(7);
		exp_original.getHPLC().getGradient().setPercentage("C2H5OH", 1, 20);
		exp_original.getHPLC().getGradient().setPercentage("C2H5OH", 5, 40);
		exp_original.getHPLC().getGradient().setPercentage("C2H5OH", 7, 60);
		exp_original.getHPLC().getGradient().setPercentage("H2O", 1, 80);
		exp_original.getHPLC().getGradient().setPercentage("H2O", 5, 60);
		exp_original.getHPLC().getGradient().setPercentage("H2O", 7, 40);
		
		exp_original.getInstrument().setModel("Porsche 911");
		exp_original.getInstrument().setVendor("Porsche K.G. Zuffenhausen");
		exp_original.getInstrument().setCustomizations("340 PS");
		exp_original.getInstrument().setMetaValue("label", String("red"));
		exp_original.getInstrument().getIonDetectors().resize(1);
		exp_original.getInstrument().getIonDetectors()[0].setAcquisitionMode(IonDetector::PULSECOUNTING);
		exp_original.getInstrument().getIonDetectors()[0].setType(IonDetector::PHOTOMULTIPLIER);
		exp_original.getInstrument().getIonDetectors()[0].setResolution(6.7677);
		exp_original.getInstrument().getIonDetectors()[0].setADCSamplingFrequency(7.6766);
		exp_original.getInstrument().getIonDetectors()[0].setMetaValue("label", String("black"));
		exp_original.getInstrument().getIonSources().resize(1);
		exp_original.getInstrument().getIonSources()[0].setInletType(IonSource::DIRECT);
		exp_original.getInstrument().getIonSources()[0].setIonizationMethod(IonSource::ESI);
		exp_original.getInstrument().getIonSources()[0].setPolarity(IonSource::POSITIVE);
		exp_original.getInstrument().getIonSources()[0].setMetaValue("label", String("blue"));
			
		MassAnalyzer analyzer;
		analyzer.setAccuracy(1.2687);
		analyzer.setFinalMSExponent(8);
		analyzer.setIsolationWidth(8.456);
		analyzer.setMagneticFieldStrength(9.999);
		analyzer.setReflectronState(MassAnalyzer::NONE);
		analyzer.setResolution(7.444);
		analyzer.setResolutionMethod(MassAnalyzer::FWHM);
		analyzer.setResolutionType(MassAnalyzer::CONSTANT);
		exp_original.getInstrument().getMassAnalyzers().push_back(analyzer);	
		analyzer = MassAnalyzer();
		analyzer.setScanDirection(MassAnalyzer::UP);
		analyzer.setScanLaw(MassAnalyzer::LINEAR);
		analyzer.setScanRate(5.555);
		analyzer.setScanTime(6.666);
		analyzer.setTOFTotalPathLength(7.777);
		analyzer.setType(MassAnalyzer::TOF);
		analyzer.setMetaValue("label", String("pink"));
		exp_original.getInstrument().getMassAnalyzers().push_back(analyzer);
		
		// MS spectrum
		RichPeakMap::SpectrumType spec;
		RichPeakMap::SpectrumType::PeakType p;
		p.setIntensity(565);
		p.getPosition()[0] = 600.1;
		p.setMetaValue("label", String("peaklabel"));
		spec.push_back(p);
		p.setIntensity(620);
		p.getPosition()[0] = 700.1;
		p.removeMetaValue("label");
		spec.push_back(p);
		p.setIntensity(701);
		p.getPosition()[0] = 800.1;
		spec.push_back(p);
		spec.setRT(1.98);
		spec.setMSLevel(1);	
		
		InstrumentSettings settings;
		settings.getScanWindows().resize(1);
		settings.getScanWindows()[0].begin = 3.456;
		settings.getScanWindows()[0].end = 7.89;
		settings.setPolarity(IonSource::NEGATIVE);
		settings.setScanMode(InstrumentSettings::ZOOM);
		spec.setInstrumentSettings (settings);
		
		// set a spectrum source file
		SourceFile source_file;
		source_file.setNameOfFile("westberlin");
		source_file.setPathToFile("/osten/");	
		spec.setSourceFile(source_file);

		RichPeakSpectrum::MetaDataArray meta_data_array;
		meta_data_array.setName ("label");
		meta_data_array.setComment ("This represents some artful kind of label.");
		meta_data_array.setName ("icon");
		meta_data_array.setComment ("little icon with colors and stuff");
		meta_data_array.setMetaValue ("icon", String("an icon is an icon is an icon"));
		meta_data_array.push_back(3.14);
		meta_data_array.push_back(3.1);
		meta_data_array.push_back(3);
		source_file.setNameOfFile("this is the filename");
		source_file.setPathToFile("/slashdot/");
		source_file.setFileSize(1.234);
		source_file.setFileType("RAWDATA");
		source_file.setChecksum("6132b58967cf1ebc05062492c17145e5ee9f82a8",SourceFile::SHA1);
		meta_data_array.setSourceFile(source_file);
		spec.getMetaDataArrays().push_back(meta_data_array);
		
		// set acquisition info with 1 acquisition
		AcquisitionInfo info;
		info.setMethodOfCombination("combo");
		Acquisition acquisition;
		acquisition.setNumber(1);
		acquisition.setMetaValue ("icon", String("yet another icon"));
		info.push_back(acquisition);
		
		spec.setAcquisitionInfo(info);

		PeptideIdentification pei;
		PeptideHit peh;
		// first PeptideIdentification (+ 2 PeptideHits) for 1st Spectrum
		std::vector<PeptideIdentification> vec_pei;
		pei.setSignificanceThreshold(1.235);
		pei.setScoreType("ScoreType");
		pei.setHigherScoreBetter(true);		
//needs getter and setter methods first
//		source_file.setNameOfFile("testberlin");
//		source_file.setPathToFile("/testen/");	
//		pei.setSourceFile(source_file);
		std::vector<PeptideHit> vec_peh;
		peh.setScore(2.345);
		peh.setSequence("AACD");
		peh.setCharge(7);
		peh.setAABefore('b');
		peh.setAAAfter('c');
		vec_peh.push_back(peh);
		peh = PeptideHit();
		peh.setAABefore('d');
		peh.setAAAfter('e');
		vec_peh.push_back(peh);
		pei.setHits(vec_peh);
		vec_pei.push_back(pei);		

		// second PeptideIdentification (+ no PeptideHits) for 1st Spectrum
		pei = PeptideIdentification();
		pei.setHigherScoreBetter(false);		
		vec_pei.push_back(pei);
		spec.setPeptideIdentifications(vec_pei);

		exp_original.push_back(spec);
			
		//MSMS spectrum
		spec.clear();
		p.setIntensity(210);
		p.getPosition()[0] = 100.155;
		spec.push_back(p);
		p.setIntensity(101);
		p.getPosition()[0] = 150.25;
		spec.push_back(p);
		p.setIntensity(90);
		p.getPosition()[0] = 300.5;
		spec.push_back(p);
		spec.setRT(3.96);
		spec.setMSLevel(2);
		spec.getPrecursorPeak().getPosition()[0] = 600.1;
		spec.getPrecursorPeak().setIntensity(4711);
		spec.getPrecursorPeak().setCharge(2);
		spec.getPrecursor().setMetaValue("icon",String("Precursor"));
		spec.getPrecursor().setWindowSize(0.1456);
		spec.setComment("bla");

		spec.getMetaDataArrays().clear();
		
		// set empty AcquisitionInfo for spectrum 2
		spec.setAcquisitionInfo(AcquisitionInfo());
		
		exp_original.push_back(spec);
		
		//meta info
		exp_original.setMetaValue("label",5.55);
		exp_original.setMetaValue("icon",String("MSExperiment"));
		exp_original.setMetaValue("color",5);
		exp_original[0].setMetaValue("icon",String("Spectrum1"));
		exp_original[1].setMetaValue("icon",String("Spectrum2"));
		
		// to store the id of reading and writing
		UID tmp_id,spec_tmp_id,tmp_id2,spec_tmp_id2;
		
		// create a Peak1D experiment (raw data)
		// Peak1Ds are no MetaInfoInterfaces --> peak meta data should not be
		// tried to be stored in DB
		MSExperiment<Peak1D> exp_peak1d;
		MSSpectrum<Peak1D> spec_peak1d;
		Peak1D peak1d;
		
		peak1d.setIntensity(565);
		peak1d.getPosition()[0] = 600.1;
		spec_peak1d.push_back(peak1d);
		peak1d.setIntensity(620);
		peak1d.getPosition()[0] = 700.1;
		spec_peak1d.push_back(peak1d);
		peak1d.setIntensity(701);
		peak1d.getPosition()[0] = 800.1;
		spec_peak1d.push_back(peak1d);
		spec_peak1d.setRT(1.98);
		spec_peak1d.setMSLevel(1);	
		
		exp_peak1d.push_back(spec_peak1d);
		
		// save newly created experiments - should be added to database.
		// success is implicitly checked later when loading from database.
	START_SECTION((template<class ExperimentType> void storeExperiment(ExperimentType& exp)))
	  DBAdapter a(con);
	  a.storeExperiment(exp_original);
	  a.storeExperiment(exp_peak1d);
		tmp_id = exp_original.getPersistenceId();
		tmp_id2 = exp_peak1d.getPersistenceId();
		spec_tmp_id = exp_original[0].getPersistenceId();
		spec_tmp_id2 = exp_peak1d[0].getPersistenceId();
		QSqlQuery result = con.executeQuery("SELECT id FROM META_MSExperiment");
	  TEST_EQUAL(result.size(),2)
END_SECTION
	
	// add another experiment to the database (for TOPPView tests etc.)
	DBAdapter a(con);
	RichPeakMap exp_2;
	FileHandler fh;
	fh.loadExperiment("data/SimpleExtender_test.mzData", exp_2);
	a.storeExperiment(exp_2);
	
	
	// check if first spectrum of the first saved experiment can be loaded correctly
	START_SECTION((template <class SpectrumType> void loadSpectrum(UID id, SpectrumType &spec)))
	  	DBAdapter a(con);
	  	DBAdapter a2(con2);
		  
			RichPeakSpectrum spec;
			a.loadSpectrum(spec_tmp_id, spec);
						
		  TEST_EQUAL( spec.getRT() , exp_original.begin()->getRT() )
			TEST_EQUAL( spec.getMSLevel() , exp_original.begin()->getMSLevel() )
			TEST_EQUAL( spec.size() , exp_original.begin()->size() )
			TEST_EQUAL( spec.getInstrumentSettings().getScanWindows().size(),1)
			TEST_REAL_SIMILAR( spec.getInstrumentSettings().getScanWindows()[0].begin , exp_original.begin()->getInstrumentSettings().getScanWindows()[0].begin )
			TEST_REAL_SIMILAR( spec.getInstrumentSettings().getScanWindows()[0].end , exp_original.begin()->getInstrumentSettings().getScanWindows()[0].end )
			TEST_EQUAL( spec.getInstrumentSettings().getPolarity() , exp_original.begin()->getInstrumentSettings().getPolarity() )
			TEST_EQUAL( spec.getInstrumentSettings().getScanMode() , exp_original.begin()->getInstrumentSettings().getScanMode() )
			TEST_EQUAL( spec.getAcquisitionInfo().getMethodOfCombination(), "combo");
			// and how do we check	info.setSpectrumType("type"); ?
			TEST_EQUAL( spec.getAcquisitionInfo()[0].getNumber(), 1);
			TEST_EQUAL( spec.getAcquisitionInfo()[0].getMetaValue("icon"), "yet another icon");
	
			TEST_EQUAL( spec.getSourceFile().getNameOfFile() , exp_original.begin()->getSourceFile().getNameOfFile() )
			TEST_EQUAL( spec.getSourceFile().getPathToFile() , exp_original.begin()->getSourceFile().getPathToFile() )
			TEST_EQUAL( spec.getSourceFile().getChecksum() , exp_original.begin()->getSourceFile().getChecksum() )
			
			// make sure storing/loading of meta data works for RichPeaks
			TEST_EQUAL( spec[0].getMetaValue("label"), "peaklabel");
			
			RichPeakSpectrum::MetaDataArrays& meta_data_arrays = spec.getMetaDataArrays();
			TEST_EQUAL( meta_data_arrays[0].getComment(), "little icon with colors and stuff" )
			TEST_EQUAL( meta_data_arrays[0].getSourceFile().getNameOfFile(), "this is the filename" )
			TEST_EQUAL( meta_data_arrays[0].getSourceFile().getPathToFile(), "/slashdot/" )
			TEST_REAL_SIMILAR( meta_data_arrays[0].getSourceFile().getFileSize(), 1.234 )
			TEST_EQUAL( meta_data_arrays[0].getSourceFile().getFileType(), "RAWDATA" )
			TEST_EQUAL( meta_data_arrays[0].getMetaValue("icon"), "an icon is an icon is an icon" )
			TEST_REAL_SIMILAR( meta_data_arrays[0][0], 3.14 )
			TEST_REAL_SIMILAR( meta_data_arrays[0][1], 3.1 )
			TEST_REAL_SIMILAR( meta_data_arrays[0][2], 3 )

			
			TEST_EQUAL( spec.getSourceFile().getNameOfFile(), "westberlin" )
			TEST_EQUAL( spec.getSourceFile().getPathToFile(), "/osten/" )
			
			for (UInt i=0; i<3; ++i)
			{
				TEST_REAL_SIMILAR( spec[i].getIntensity() , exp_original.begin()->operator[](i).getIntensity() )
				TEST_REAL_SIMILAR( spec[i].getPosition()[0] , exp_original.begin()->operator[](i).getPosition()[0] )
			}

			PeakFileOptions options;
			options.setIntensityRange(DRange<1> (600, 1000));
			a.getOptions() = options;
			a.loadSpectrum(spec_tmp_id, spec);
			
			// check if the Intensity restriction worked - first peak (565) should have been skipped
			TEST_REAL_SIMILAR( spec[0].getIntensity() , 620 )
			TEST_REAL_SIMILAR( spec[1].getIntensity() , 701 )

			options = PeakFileOptions();
			options.setMZRange(DRange<1> (650, 1000));
			a.getOptions() = options;
			a.loadSpectrum(spec_tmp_id, spec);

			// check if the MZ restriction worked - first peak (600.1) should have been skipped
			TEST_REAL_SIMILAR( spec[0].getPosition()[0] , 700.1 )
			TEST_REAL_SIMILAR( spec[1].getPosition()[0] , 800.1 )

			// testing concurrent DB connections
			a2.loadSpectrum(spec_tmp_id, spec);
			TEST_REAL_SIMILAR( spec[0].getIntensity() , 565 )
END_SECTION
		
	// load first two experiments from database
	// (this implicitly checks if the new experiments were stored correctly)
	START_SECTION((template <class ExperimentType> void loadExperiment(UID id, ExperimentType &exp)))
		  DBAdapter a(con);
		  RichPeakMap exp_new;
		  std::map<String, MetaInfoDescription> descriptions;
			
			a.loadExperiment(tmp_id, exp_new);
			TEST_EQUAL(exp_new.getPersistenceId(), tmp_id)
			TEST_EQUAL(exp_new.getComment() , "bla" )
			
			TEST_EQUAL(exp_new.getSample().getName(), "fruity loops" )
			TEST_EQUAL(exp_new.getSample().getNumber(), "007" )
			TEST_REAL_SIMILAR(exp_new.getSample().getMass(), 30.1 )
			TEST_REAL_SIMILAR(exp_new.getSample().getSubsamples()[0].getVolume(), 60.1 )
			TEST_REAL_SIMILAR(exp_new.getSample().getSubsamples()[0].getConcentration(), 101.1 )
			const Digestion* digestion = dynamic_cast<const Digestion*>(&exp_new.getSample().getSubsamples()[0].getTreatment(0));
			TEST_EQUAL(digestion->getEnzyme(), "dhdh" )
			TEST_REAL_SIMILAR(digestion->getDigestionTime(), 36.6 )
			TEST_REAL_SIMILAR(digestion->getPh(), 7.2 )
			TEST_REAL_SIMILAR(digestion->getTemperature(), 37.7 )

			TEST_EQUAL(exp_new.getProteinIdentifications()[0].getSearchEngine(), "google" )
			TEST_EQUAL(exp_new.getProteinIdentifications()[0].getSearchEngineVersion(), "beta" )
			TEST_EQUAL(exp_new.getProteinIdentifications()[0].getDateTime().get(), "2006-12-12 00:00:00" )
			TEST_EQUAL(exp_new.getProteinIdentifications()[0].getScoreType(), "Type" )
			TEST_EQUAL(exp_new.getProteinIdentifications()[0].isHigherScoreBetter(), true )
			TEST_REAL_SIMILAR(exp_new.getProteinIdentifications()[0].getSignificanceThreshold(), 3.456 )
			TEST_EQUAL(exp_new.getProteinIdentifications()[0].getHits()[0].getAccession(), "0110110" )
			TEST_REAL_SIMILAR(exp_new.getProteinIdentifications()[0].getHits()[1].getScore(), 4.567 )
			TEST_EQUAL(exp_new.getProteinIdentifications()[0].getHits()[1].getAccession(), "1001001" )
			TEST_EQUAL(exp_new.getProteinIdentifications()[0].getHits()[1].getSequence(), "ZXY" )

			TEST_REAL_SIMILAR(exp_new[0].getPeptideIdentifications()[0].getSignificanceThreshold(), 1.235 )	
			TEST_EQUAL(exp_new[0].getPeptideIdentifications()[0].getScoreType(), "ScoreType" )
			TEST_EQUAL(exp_new[0].getPeptideIdentifications()[0].isHigherScoreBetter(), true )
// needs getter and setter methods first
//			TEST_EQUAL(exp_new[0].getPeptideIdentifications()[0].getSourceFile().getNameOfFile(), "testberlin" )
//			TEST_EQUAL(exp_new[0].getPeptideIdentifications()[0].getSourceFile().getPathToFile(), "/testen/" )
			TEST_EQUAL(exp_new[0].getPeptideIdentifications()[1].isHigherScoreBetter(), false )

			TEST_REAL_SIMILAR(exp_new[0].getPeptideIdentifications()[0].getHits()[0].getScore(), 2.345 )	
			TEST_EQUAL(exp_new[0].getPeptideIdentifications()[0].getHits()[0].getSequence(), "AACD" )	
			TEST_EQUAL(exp_new[0].getPeptideIdentifications()[0].getHits()[0].getCharge(), 7 )	
			TEST_EQUAL(exp_new[0].getPeptideIdentifications()[0].getHits()[0].getAABefore(), 'b' )	
			TEST_EQUAL(exp_new[0].getPeptideIdentifications()[0].getHits()[0].getAAAfter(), 'c' )	
			TEST_EQUAL(exp_new[0].getPeptideIdentifications()[0].getHits()[1].getAABefore(), 'd' )	
			TEST_EQUAL(exp_new[0].getPeptideIdentifications()[0].getHits()[1].getAAAfter(), 'e' )	
			
			TEST_EQUAL(exp_new.getSample().getSubsamples()[1].getState(), Sample::GAS )
			TEST_EQUAL(exp_new.getSample().getSubsamples()[1].getOrganism(), "isistius brasiliensis (cookiecutter shar" )
			const Modification* modification = dynamic_cast<const Modification*>(&exp_new.getSample().getSubsamples()[1].getTreatment(0));
			TEST_EQUAL(modification->getReagentName(), "reagent" )
			TEST_EQUAL(modification->getAffectedAminoAcids(), "123" )
			TEST_EQUAL(modification->getSpecificityType(), Modification::CTERM )
			TEST_REAL_SIMILAR(modification->getMass(), 12.3 )
			const Tagging* tagging = dynamic_cast<const Tagging*>(&exp_new.getSample().getSubsamples()[1].getTreatment(1));
			TEST_EQUAL(tagging->getReagentName(), "tagging" )
			TEST_REAL_SIMILAR(tagging->getMassShift(), 0.123 )
			TEST_EQUAL(tagging->getVariant(), Tagging::HEAVY )

			TEST_EQUAL(exp_new.getSample().getSubsamples()[2].getComment(), "nice" )
			TEST_EQUAL(exp_new.getSample().getSubsamples()[2].getMetaValue("label"), "pink" )
			
			TEST_EQUAL(exp_new.getContacts()[0].getFirstName() , "Ferdinand" )
			TEST_EQUAL(exp_new.getContacts()[0].getLastName() , "Piech" )
			TEST_EQUAL(exp_new.getContacts()[0].getInstitution() , "aff" )
			TEST_EQUAL(exp_new.getContacts()[1].getEmail() , "ferdi@porsche.de" )
			TEST_EQUAL(exp_new.getContacts()[1].getContactInfo() , "ttss" )
			TEST_EQUAL(exp_new.getContacts()[1].getMetaValue("label") , "polka-dotted" )
			
			TEST_EQUAL(exp_new.getHPLC().getInstrument() , "guitar" )
			TEST_EQUAL(exp_new.getHPLC().getColumn() , "bigone" )
			TEST_EQUAL(exp_new.getHPLC().getComment() , "fhth" )
			TEST_EQUAL(exp_new.getHPLC().getFlux() , 1 )
			TEST_EQUAL(exp_new.getHPLC().getPressure() , 2 )
			TEST_EQUAL(exp_new.getHPLC().getTemperature() , 3 )
			
			TEST_EQUAL(exp_new.getHPLC().getGradient().getPercentages()[0][0] , 20 )
			TEST_EQUAL(exp_new.getHPLC().getGradient().getPercentages()[0][1] , 40 )
			TEST_EQUAL(exp_new.getHPLC().getGradient().getPercentages()[0][2] , 60 )
			TEST_EQUAL(exp_new.getHPLC().getGradient().getPercentages()[1][0] , 80 )
			TEST_EQUAL(exp_new.getHPLC().getGradient().getPercentages()[1][1] , 60 )
			TEST_EQUAL(exp_new.getHPLC().getGradient().getPercentages()[1][2] , 40 )
			TEST_EQUAL(exp_new.getHPLC().getGradient().getEluents()[0] , "C2H5OH" )
			TEST_EQUAL(exp_new.getHPLC().getGradient().getEluents()[1] , "H2O" )
			TEST_EQUAL(exp_new.getHPLC().getGradient().getTimepoints()[0] , 1 )
			TEST_EQUAL(exp_new.getHPLC().getGradient().getTimepoints()[1] , 5 )
			TEST_EQUAL(exp_new.getHPLC().getGradient().getTimepoints()[2] , 7 )
			
			TEST_EQUAL(exp_new.getInstrument().getModel() , "Porsche 911" )
			TEST_EQUAL(exp_new.getInstrument().getVendor() , "Porsche K.G. Zuffenhausen" )
			TEST_EQUAL(exp_new.getInstrument().getCustomizations() , "340 PS" )
			TEST_EQUAL(exp_new.getInstrument().getMetaValue("label") , "red" )
			TEST_EQUAL(exp_new.getInstrument().getIonDetectors().size(),1)
			TEST_EQUAL(exp_new.getInstrument().getIonDetectors()[0].getType() , IonDetector::PHOTOMULTIPLIER )
			TEST_EQUAL(exp_new.getInstrument().getIonDetectors()[0].getAcquisitionMode() , IonDetector::PULSECOUNTING )
			TEST_REAL_SIMILAR(exp_new.getInstrument().getIonDetectors()[0].getResolution() , 6.7677 )
			TEST_REAL_SIMILAR(exp_new.getInstrument().getIonDetectors()[0].getADCSamplingFrequency() , 7.6766 )
			TEST_EQUAL(exp_new.getInstrument().getIonDetectors()[0].getMetaValue("label") , "black" )
			TEST_EQUAL(exp_new.getInstrument().getIonSources().size(),1)
			TEST_EQUAL(exp_new.getInstrument().getIonSources()[0].getInletType() , IonSource::DIRECT )
			TEST_EQUAL(exp_new.getInstrument().getIonSources()[0].getIonizationMethod() , IonSource::ESI )
			TEST_EQUAL(exp_new.getInstrument().getIonSources()[0].getPolarity() , IonSource::POSITIVE )
			TEST_EQUAL(exp_new.getInstrument().getIonSources()[0].getMetaValue("label") , "blue" )
			
			TEST_REAL_SIMILAR(exp_new.getInstrument().getMassAnalyzers()[0].getAccuracy() , 1.2687 )
			TEST_EQUAL(exp_new.getInstrument().getMassAnalyzers()[0].getFinalMSExponent() , 8 )
			TEST_REAL_SIMILAR(exp_new.getInstrument().getMassAnalyzers()[0].getIsolationWidth() , 8.456 )
			TEST_REAL_SIMILAR(exp_new.getInstrument().getMassAnalyzers()[0].getMagneticFieldStrength() , 9.999 )
			TEST_EQUAL(exp_new.getInstrument().getMassAnalyzers()[0].getReflectronState() , MassAnalyzer::NONE )
			TEST_REAL_SIMILAR(exp_new.getInstrument().getMassAnalyzers()[0].getResolution() , 7.444 )
			TEST_EQUAL(exp_new.getInstrument().getMassAnalyzers()[0].getResolutionMethod() , MassAnalyzer::FWHM )
			TEST_EQUAL(exp_new.getInstrument().getMassAnalyzers()[0].getResolutionType() , MassAnalyzer::CONSTANT )
			TEST_EQUAL(exp_new.getInstrument().getMassAnalyzers()[1].getScanDirection() , MassAnalyzer::UP )
			TEST_EQUAL(exp_new.getInstrument().getMassAnalyzers()[1].getScanLaw() , MassAnalyzer::LINEAR )
			TEST_REAL_SIMILAR(exp_new.getInstrument().getMassAnalyzers()[1].getScanRate() , 5.555 )
			TEST_REAL_SIMILAR(exp_new.getInstrument().getMassAnalyzers()[1].getScanTime() , 6.666 )
			TEST_REAL_SIMILAR(exp_new.getInstrument().getMassAnalyzers()[1].getTOFTotalPathLength() , 7.777 )
			TEST_EQUAL(exp_new.getInstrument().getMassAnalyzers()[1].getType() , MassAnalyzer::TOF )
			TEST_EQUAL(exp_new.getInstrument().getMassAnalyzers()[1].getMetaValue("label") , "pink" )
			
			//------ test if values are correct ------
			
			//SPECTRUM 1
			RichPeakMap::const_iterator itn(exp_new.begin());
			RichPeakMap::const_iterator ito(exp_original.begin());
				
		  TEST_EQUAL( itn->getRT() , ito->getRT() )
			TEST_EQUAL( itn->getMSLevel() , ito->getMSLevel() )
			TEST_EQUAL( itn->size() , ito->size() )
			for (UInt i=0; i<3; ++i)
			{
				TEST_REAL_SIMILAR( itn->operator[](i).getIntensity() , ito->operator[](i).getIntensity() )
				TEST_REAL_SIMILAR( itn->operator[](i).getPosition()[0] , ito->operator[](i).getPosition()[0] )
			}
		
			//SPECTRUM 2
			++itn;
			++ito;
				
		  TEST_EQUAL( itn->getRT() , ito->getRT() )
			TEST_EQUAL( itn->getMSLevel() , ito->getMSLevel() )
			TEST_EQUAL( itn->getPrecursorPeak().getPosition()[0] , ito->getPrecursorPeak().getPosition()[0] )
			TEST_EQUAL( itn->getPrecursorPeak().getIntensity() , ito->getPrecursorPeak().getIntensity() )
			TEST_EQUAL( itn->getPrecursorPeak().getCharge() , ito->getPrecursorPeak().getCharge() )
			TEST_EQUAL( itn->getPrecursor().getMetaValue("icon") , "Precursor" )
			TEST_REAL_SIMILAR( itn->getPrecursor().getWindowSize() , 0.1456)
	
			TEST_EQUAL( itn->getComment() , "bla" )
			TEST_EQUAL( itn->size() , ito->size() )
			for (UInt i=0; i<3; ++i)
			{
				TEST_REAL_SIMILAR( itn->operator[](i).getIntensity() , ito->operator[](i).getIntensity() )
				TEST_REAL_SIMILAR( itn->operator[](i).getPosition()[0] , ito->operator[](i).getPosition()[0] )
			}
			
			//META INFO
			TEST_REAL_SIMILAR((double)exp_new.getMetaValue("label"),5.55)
			TEST_EQUAL((string)exp_new.getMetaValue("icon"),"MSExperiment")
			TEST_EQUAL((int)exp_new.getMetaValue("color"),5)
			TEST_EQUAL((string)exp_new[0].getMetaValue("icon"),"Spectrum1")
			TEST_EQUAL((string)exp_new[1].getMetaValue("icon"),"Spectrum2")

			exp_new = RichPeakMap();
			PeakFileOptions options;
			options = PeakFileOptions();
			options.setRTRange(DRange<1> (2.5, 4.5));
			a.getOptions() = options;
			a.loadExperiment(tmp_id, exp_new);

			// check if the RT restriction worked - first spectrum should have been skipped
			TEST_REAL_SIMILAR( exp_new[0][0].getPosition()[0] , 100.155 )

			exp_new = RichPeakMap();
			options = PeakFileOptions();
			std::vector<int> levels;
			levels.push_back(2);
			options.setMSLevels(levels);
			a.getOptions() = options;
			a.loadExperiment(tmp_id, exp_new);

			// check if the MSLevel restriction worked - first spectrum should have been skipped
			TEST_REAL_SIMILAR( exp_new[0][0].getPosition()[0] , 100.155 )
END_SECTION
	
		// save modified version of already existing experiment - old records should be updated.
		// no checks are run, results are implicitly checked later when loading
		START_SECTION([EXTRA] updating of an existing dataset)
			exp_original.setComment("blubb");
	
			// modify first spectrum
			RichPeakMap::SpectrumType & modified_spec = exp_original[0];
			modified_spec[0].setIntensity(566);
			modified_spec[0].getPosition()[0] = 612.1;
			modified_spec[1].setIntensity(620);
			modified_spec[1].getPosition()[0] = 712.1;
			modified_spec[2].setIntensity(701);
			modified_spec[2].getPosition()[0] = 812.1;
			modified_spec.setRT(1.88);
			modified_spec.setMSLevel(1);
			modified_spec.getInstrumentSettings().getScanWindows()[0].begin = 3.567;
			modified_spec.getInstrumentSettings().getScanWindows()[0].end = 7.91;
			modified_spec.getInstrumentSettings().setPolarity(IonSource::POSITIVE);
			modified_spec.getInstrumentSettings().setScanMode(InstrumentSettings::ZOOM);
			modified_spec.getInstrumentSettings().setMetaValue("label", String("please bite here"));
			
			info.clear();
			acquisition.setNumber(1);
			acquisition.setMetaValue ("icon", String("one more icon"));
			info.push_back(acquisition);
			acquisition.setNumber(2);
			acquisition.setMetaValue ("label", String("yet another label"));
			info.push_back(acquisition);
			
			modified_spec.setAcquisitionInfo(info);
			// adding a meta data array
			modified_spec.getMetaDataArrays().clear();
			RichPeakSpectrum::MetaDataArray meta_data_array;
			meta_data_array.setName ("label");
			meta_data_array.setComment ("This represents some artful kind of label.");
			meta_data_array.setName ("icon");
			meta_data_array.setComment ("little icon with colors and stuff");
			meta_data_array.push_back(23);
			meta_data_array.push_back(42);
			meta_data_array.push_back(100.001);
			// setting a source file
			SourceFile source_file;
			source_file.setNameOfFile("this is the filename");
			source_file.setPathToFile("/slashdot/");
			source_file.setFileSize(1.234);
			source_file.setFileType("RAWDATA");
			meta_data_array.setSourceFile(source_file);
			
			modified_spec.getMetaDataArrays().push_back(meta_data_array);
			
			// modify 2nd spectrum
			exp_original[1].getPrecursor().setMetaValue("icon", String("NewPrecursor"));
	
		  DBAdapter a(con);
		  a.storeExperiment(exp_original);
			
			////////////PART 2 => LOADING
			
		  RichPeakMap exp_new;
			
			a.loadExperiment(tmp_id, exp_new);
			TEST_EQUAL(exp_new.getPersistenceId(), tmp_id)
			TEST_EQUAL(exp_new.getComment() , "blubb" )
			
			//------ test if values are correct ------
			
			//SPECTRUM 1
			RichPeakMap::const_iterator itn(exp_new.begin());
			RichPeakMap::const_iterator ito(exp_original.begin());
				
		  TEST_EQUAL( itn->getRT() , ito->getRT() )
			TEST_EQUAL( itn->getMSLevel() , ito->getMSLevel() )
			TEST_EQUAL( itn->size() , ito->size() )
			TEST_EQUAL( itn->getInstrumentSettings().getMetaValue("label") , "please bite here" )
			TEST_EQUAL( itn->getAcquisitionInfo()[0].getNumber(), 1);
			TEST_EQUAL( itn->getAcquisitionInfo()[0].getMetaValue("icon"), "one more icon");
			TEST_EQUAL( itn->getAcquisitionInfo()[1].getNumber(), 2);
			TEST_EQUAL( itn->getAcquisitionInfo()[1].getMetaValue("label"), "yet another label");
			for (UInt i=0; i<3; ++i)
			{
				TEST_REAL_SIMILAR( itn->operator[](i).getIntensity() , ito->operator[](i).getIntensity() )
				TEST_REAL_SIMILAR( itn->operator[](i).getPosition()[0] , ito->operator[](i).getPosition()[0] )
			}
		
			//SPECTRUM 2
			++itn;
			++ito;
				
		  TEST_EQUAL( itn->getRT() , ito->getRT() )
			TEST_EQUAL( itn->getMSLevel() , ito->getMSLevel() )
			TEST_EQUAL( itn->getPrecursorPeak().getPosition()[0] , ito->getPrecursorPeak().getPosition()[0] )
			TEST_EQUAL( itn->getPrecursorPeak().getIntensity() , ito->getPrecursorPeak().getIntensity() )
			TEST_EQUAL( itn->getPrecursorPeak().getCharge() , ito->getPrecursorPeak().getCharge() )
			TEST_EQUAL( itn->getPrecursor().getMetaValue("icon") , "NewPrecursor" )
			TEST_EQUAL( itn->getComment() , "bla" )
			TEST_EQUAL( itn->size() , ito->size() )
			for (UInt i=0; i<3; ++i)
			{
				TEST_REAL_SIMILAR( itn->operator[](i).getIntensity() , ito->operator[](i).getIntensity() )
				TEST_REAL_SIMILAR( itn->operator[](i).getPosition()[0] , ito->operator[](i).getPosition()[0] )
			}
			
			//META INFO
			TEST_REAL_SIMILAR((double)exp_new.getMetaValue("label"),5.55)
			TEST_EQUAL((string)exp_new.getMetaValue("icon"),"MSExperiment")
			TEST_EQUAL((int)exp_new.getMetaValue("color"),5)
			TEST_EQUAL((string)exp_new[0].getMetaValue("icon"),"Spectrum1")
			TEST_EQUAL((string)exp_new[1].getMetaValue("icon"),"Spectrum2")
			
			//load the Peak1D experiment
			//(peak meta data should not be tried to be loaded, because
			//Peak1D is no MetaInfoInterface)
			MSExperiment<Peak1D> exp2;
			a.loadExperiment(tmp_id2, exp2);
			TEST_EQUAL( exp2.size(), 1 );
			MSSpectrum<Peak1D>& spec2 = *exp2.begin();
			MSSpectrum<Peak1D>& spec2_original = *exp_peak1d.begin();
			TEST_EQUAL ( spec2.size(), 3 )
			//test if values are correct
			for(int i = 0; i < 3; i++)
			{
				TEST_REAL_SIMILAR( spec2[i].getIntensity(), spec2_original[i].getIntensity() )
				TEST_REAL_SIMILAR( spec2[i].getPosition()[0], spec2_original[i].getPosition()[0] )
			}
						
END_SECTION
	
		START_SECTION(([EXTRA] load and store of empty map))
			DBAdapter a(con);
		  RichPeakMap in, out;
		  a.storeExperiment(in);
			a.loadExperiment(in.getPersistenceId(),out);
			TEST_EQUAL(in==out, true)
END_SECTION

	START_SECTION((const PeakFileOptions& getOptions() const))
			DBAdapter a(con);
			TEST_EQUAL(a.getOptions().hasMSLevels(),false)
END_SECTION
		
	START_SECTION((PeakFileOptions& getOptions()))
			DBAdapter a(con);
			a.getOptions().addMSLevel(1);
			TEST_EQUAL(a.getOptions().hasMSLevels(),true);
END_SECTION

	//extra test with an empty spectrum
	START_SECTION(([EXTRA] template<class ExperimentType> void storeExperiment(ExperimentType& exp)))
		  RichPeakMap exp_tmp;
		  exp_tmp.resize(1);
		  DBAdapter a(con);
		  a.storeExperiment(exp_tmp);
		  TEST_NOT_EQUAL(exp_tmp[0].getPersistenceId(),0);
END_SECTION		

	} // DB up-to-date

}
else
{
	ADD_MESSAGE("skipped")
}
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



