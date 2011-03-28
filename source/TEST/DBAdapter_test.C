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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/DB/DBConnection.h>
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
	credentials.load(String(OPENMS_BINARY_PATH) + "/source/TEST/DB_credentials.txt",true);
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
	DBAdapter* nullPointer = 0;

	START_SECTION((DBAdapter(DBConnection& db_con)))
		ptr = new DBAdapter(con);
	TEST_NOT_EQUAL(ptr, nullPointer)
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
		ph.setRank(4u);
		vector_ph.push_back(ph);
		pi.setHits(vector_ph);
		pi.setMetaValue("label", String("proteinidentificationmetainfo"));

		ProteinIdentification::SearchParameters sp;
		sp.db = "register";
		sp.db_version = "0.0";
		sp.taxonomy = "bundesadler";
		sp.charges = "high";
		sp.mass_type = ProteinIdentification::AVERAGE;
		sp.enzyme = ProteinIdentification::TRYPSIN;
		sp.missed_cleavages = 6;
		sp.peak_mass_tolerance = 0.44;
		sp.precursor_tolerance = 0.55;
		sp.setMetaValue("label", String("searchparametersmetainfo"));

		std::vector<String> fm (3,"a");
		for(Size i = 0; i < fm.size(); ++i)
		{
			fm[i]+=i;
		}
		std::vector<String> vm (fm);
		for(Size i = 0; i < vm.size(); ++i)
		{
			vm[i]+=i;
		}

		sp.fixed_modifications = fm;
		sp.variable_modifications = vm;
		pi.setSearchParameters(sp);

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

		Software sw;
		sw.setName("tolle instrument-software");
		sw.setVersion("alpha");
		sw.setMetaValue("label", String("neu fuer msinstrument"));
		exp_original.getInstrument().setSoftware(sw);
		exp_original.getInstrument().setModel("Porsche 911");
		exp_original.getInstrument().setVendor("Porsche K.G. Zuffenhausen");
		exp_original.getInstrument().setCustomizations("340 PS");
		exp_original.getInstrument().setMetaValue("label", String("red"));
		exp_original.getInstrument().getIonDetectors().resize(1);
		exp_original.getInstrument().getIonDetectors()[0].setAcquisitionMode(IonDetector::PULSECOUNTING);
		exp_original.getInstrument().getIonDetectors()[0].setType(IonDetector::PHOTOMULTIPLIER);
		exp_original.getInstrument().getIonDetectors()[0].setResolution(6.7677);
		exp_original.getInstrument().getIonDetectors()[0].setADCSamplingFrequency(7.6766);
		exp_original.getInstrument().getIonDetectors()[0].setOrder(3);
		exp_original.getInstrument().getIonDetectors()[0].setMetaValue("label", String("black"));
		exp_original.getInstrument().getIonSources().resize(1);
		exp_original.getInstrument().getIonSources()[0].setInletType(IonSource::DIRECT);
		exp_original.getInstrument().getIonSources()[0].setIonizationMethod(IonSource::ESI);
		exp_original.getInstrument().getIonSources()[0].setPolarity(IonSource::POSITIVE);
		exp_original.getInstrument().getIonSources()[0].setMetaValue("label", String("blue"));
		exp_original.getInstrument().getIonSources()[0].setOrder(0);
		exp_original.getInstrument().setIonOptics(Instrument::FRINGING_FIELD);

		MassAnalyzer analyzer;
		analyzer.setAccuracy(1.2687);
		analyzer.setFinalMSExponent(8);
		analyzer.setIsolationWidth(8.456);
		analyzer.setMagneticFieldStrength(9.999);
		analyzer.setReflectronState(MassAnalyzer::NONE);
		analyzer.setResolution(7.444);
		analyzer.setResolutionMethod(MassAnalyzer::FWHM);
		analyzer.setResolutionType(MassAnalyzer::CONSTANT);
		analyzer.setOrder(1);
		exp_original.getInstrument().getMassAnalyzers().push_back(analyzer);
		analyzer = MassAnalyzer();
		analyzer.setScanDirection(MassAnalyzer::UP);
		analyzer.setScanLaw(MassAnalyzer::LINEAR);
		analyzer.setScanRate(5.555);
		analyzer.setScanTime(6.666);
		analyzer.setTOFTotalPathLength(7.777);
		analyzer.setType(MassAnalyzer::TOF);
		analyzer.setMetaValue("label", String("pink"));
		analyzer.setOrder(2);
		exp_original.getInstrument().getMassAnalyzers().push_back(analyzer);

		// MS spectrum
		RichPeakMap::SpectrumType spec;
		RichPeakMap::SpectrumType::PeakType p;
		p.setIntensity(565.0f);
		p.getPosition()[0] = 600.1;
		p.setMetaValue("label", String("peaklabel"));
		spec.push_back(p);
		p.setIntensity(620.0f);
		p.getPosition()[0] = 700.1;
		p.removeMetaValue("label");
		spec.push_back(p);
		p.setIntensity(701.0f);
		p.getPosition()[0] = 800.1;
		spec.push_back(p);
		spec.setRT(1.98);
		spec.setMSLevel(1);

		std::vector<Product> eier;
		Product ei;
		ei.setMZ(1.0);
		ei.setIsolationWindowLowerOffset(2.0);
		ei.setIsolationWindowUpperOffset(3.0);
		ei.setMetaValue("farbe", String("lilablassblau"));
		eier.push_back(ei);
		spec.setProducts(eier);
		ei.setMZ(4.0);
		ei.setIsolationWindowLowerOffset(5.0);
		ei.setIsolationWindowUpperOffset(6.0);
		ei.removeMetaValue("farbe");
		eier.push_back(ei);
		spec.setProducts(eier);

		InstrumentSettings settings;
		settings.getScanWindows().resize(1);
		settings.getScanWindows()[0].begin = 3.456;
		settings.getScanWindows()[0].end = 7.89;
		settings.getScanWindows()[0].setMetaValue("metavalue", String("info"));
		settings.setPolarity(IonSource::NEGATIVE);
		settings.setScanMode(InstrumentSettings::SIM);
		settings.setZoomScan(true);
		spec.setInstrumentSettings (settings);

		// set a spectrum source file
		SourceFile source_file;
		source_file.setNameOfFile("westberlin");
		source_file.setPathToFile("/osten/");
		source_file.setNativeIDType("Waters nativeID format");
		spec.setSourceFile(source_file);

		RichPeakSpectrum::FloatDataArray meta_data_array;
		meta_data_array.setName("icon");
		meta_data_array.setMetaValue ("icon", String("an icon is an icon is an icon"));
		meta_data_array.push_back(3.14f);
		meta_data_array.push_back(3.1f);
		meta_data_array.push_back(3.0f);
		spec.getFloatDataArrays().push_back(meta_data_array);

		// set acquisition info with 1 acquisition
		AcquisitionInfo info;
		info.setMethodOfCombination("combo");
		Acquisition acquisition;
		acquisition.setIdentifier("1");
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
		spec.clear(false);
		p.setIntensity(210.0f);
		p.getPosition()[0] = 100.155;
		spec.push_back(p);
		p.setIntensity(101.0f);
		p.getPosition()[0] = 150.25;
		spec.push_back(p);
		p.setIntensity(90.0f);
		p.getPosition()[0] = 300.5;
		spec.push_back(p);
		spec.setRT(3.96);
		spec.setMSLevel(2);

		DataProcessing dp;
		DateTime d;
		d.set("2006-12-12 00:00:00");
		dp.setCompletionTime(d);
		sw.setName("tolle software");
		sw.setVersion("alpha");
		sw.setMetaValue("label", String("echt"));
		dp.setSoftware(sw);
		spec.getDataProcessing().push_back(dp);
		dp.setCompletionTime(d);
		dp.getSoftware().setName("nicht so tolle software");
		dp.getSoftware().setVersion("alpha");
		dp.setMetaValue("label", String("prozessiert"));
		dp.getProcessingActions().insert(DataProcessing::ALIGNMENT);
		dp.getProcessingActions().insert(DataProcessing::SMOOTHING);
		spec.getDataProcessing().push_back(dp);

		//spectrum 2 gets 2 precursors
		spec.getPrecursors().resize(2);
		//1st precursor for spectrum 2
		spec.getPrecursors()[0].setMZ(600.1);
		spec.getPrecursors()[0].setIntensity(4711.0f);
		spec.getPrecursors()[0].setCharge(2);
		spec.getPrecursors()[0].setActivationEnergy(9.99);
		spec.getPrecursors()[0].setMetaValue("icon",String("Precursor1"));
		std::vector<Int> pcs;
		pcs.push_back(1);
		pcs.push_back(2);
		pcs.push_back(3);
		spec.getPrecursors()[0].setPossibleChargeStates(pcs);
		std::set<Precursor::ActivationMethod> am;
		am.insert(Precursor::LCID);
		am.insert(Precursor::CID);
		am.insert(Precursor::HCID);
		spec.getPrecursors()[0].setActivationMethods(am);
		spec.getFloatDataArrays().clear();
		//2nd precursor for spectrum 2
		spec.getPrecursors()[1].setMZ(600.1);
		spec.getPrecursors()[1].setIntensity(4711.0f);
		spec.getPrecursors()[1].setCharge(2);
		spec.getPrecursors()[1].setActivationEnergy(9.99);
		spec.getPrecursors()[1].setMetaValue("icon",String("Precursor2"));
		pcs[0]=(4);
		pcs[1]=(5);
		pcs[2]=(6);
		spec.getPrecursors()[1].setPossibleChargeStates(pcs);
		am.erase(Precursor::CID);
		spec.getPrecursors()[1].setActivationMethods(am);
		spec.setComment("bla");
		spec.getFloatDataArrays().clear();

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

		peak1d.setIntensity(565.0f);
		peak1d.getPosition()[0] = 600.1;
		spec_peak1d.push_back(peak1d);
		peak1d.setIntensity(620.0f);
		peak1d.getPosition()[0] = 700.1;
		spec_peak1d.push_back(peak1d);
		peak1d.setIntensity(701.0f);
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
	fh.loadExperiment(OPENMS_GET_TEST_DATA_PATH("SimpleExtender_test.mzData"), exp_2);
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
			TEST_EQUAL( spec.getInstrumentSettings().getZoomScan() , exp_original.begin()->getInstrumentSettings().getZoomScan() )
			TEST_EQUAL( spec.getInstrumentSettings().getScanWindows()[0].getMetaValue("metavalue") , exp_original.begin()->getInstrumentSettings().getScanWindows()[0].getMetaValue("metavalue"))
			for(Size ps = 0; ps < spec.getProducts().size(); ++ps)
			{
				TEST_EQUAL( spec.getProducts()[ps].getMZ() , exp_original.begin()->getProducts()[ps].getMZ() )
				TEST_EQUAL( spec.getProducts()[ps].getIsolationWindowLowerOffset() , exp_original.begin()->getProducts()[ps].getIsolationWindowLowerOffset() )
				TEST_EQUAL( spec.getProducts()[ps].getIsolationWindowUpperOffset() , exp_original.begin()->getProducts()[ps].getIsolationWindowUpperOffset() )
				TEST_EQUAL( spec.getProducts()[ps].getMetaValue("farbe") , exp_original.begin()->getProducts()[ps].getMetaValue("farbe"))
			}
			TEST_EQUAL( spec.getAcquisitionInfo().getMethodOfCombination(), "combo");
			// and how do we check	info.setSpectrumType("type"); ?
			TEST_EQUAL( spec.getAcquisitionInfo()[0].getIdentifier(), "1");
			TEST_EQUAL( spec.getAcquisitionInfo()[0].getMetaValue("icon"), "yet another icon");

			TEST_EQUAL( spec.getSourceFile().getNameOfFile() , exp_original.begin()->getSourceFile().getNameOfFile() )
			TEST_EQUAL( spec.getSourceFile().getPathToFile() , exp_original.begin()->getSourceFile().getPathToFile() )
			TEST_EQUAL( spec.getSourceFile().getNativeIDType() , exp_original.begin()->getSourceFile().getNativeIDType())
			TEST_EQUAL( spec.getSourceFile().getChecksum() , exp_original.begin()->getSourceFile().getChecksum() )

			// make sure storing/loading of meta data works for RichPeaks
			TEST_EQUAL( spec[0].getMetaValue("label"), "peaklabel");

			RichPeakSpectrum::FloatDataArrays& meta_data_arrays = spec.getFloatDataArrays();
			TEST_STRING_EQUAL( meta_data_arrays[0].getName(),"icon")
			TEST_EQUAL( meta_data_arrays[0].getMetaValue("icon"), "an icon is an icon is an icon" )
			TEST_REAL_SIMILAR( meta_data_arrays[0][0], 3.14 )
			TEST_REAL_SIMILAR( meta_data_arrays[0][1], 3.1 )
			TEST_REAL_SIMILAR( meta_data_arrays[0][2], 3 )


			TEST_EQUAL( spec.getSourceFile().getNameOfFile(), "westberlin" )
			TEST_EQUAL( spec.getSourceFile().getPathToFile(), "/osten/" )

			for (Size i=0; i<3; ++i)
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

			//test for products

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


			TEST_EQUAL(exp_new.getProteinIdentifications()[0].getMetaValue("label"), exp_original.getProteinIdentifications()[0].getMetaValue("label"))

		TEST_EQUAL(exp_new.getProteinIdentifications()[0].getSearchParameters().db, exp_original.getProteinIdentifications()[0].getSearchParameters().db)
		TEST_EQUAL(exp_new.getProteinIdentifications()[0].getSearchParameters().db_version, exp_original.getProteinIdentifications()[0].getSearchParameters().db_version)
		TEST_EQUAL(exp_new.getProteinIdentifications()[0].getSearchParameters().taxonomy, exp_original.getProteinIdentifications()[0].getSearchParameters().taxonomy)
		TEST_EQUAL(exp_new.getProteinIdentifications()[0].getSearchParameters().charges, exp_original.getProteinIdentifications()[0].getSearchParameters().charges)
		TEST_EQUAL(exp_new.getProteinIdentifications()[0].getSearchParameters().mass_type, exp_original.getProteinIdentifications()[0].getSearchParameters().mass_type)
		TEST_EQUAL(exp_new.getProteinIdentifications()[0].getSearchParameters().enzyme, exp_original.getProteinIdentifications()[0].getSearchParameters().enzyme)
		TEST_EQUAL(exp_new.getProteinIdentifications()[0].getSearchParameters().missed_cleavages, exp_original.getProteinIdentifications()[0].getSearchParameters().missed_cleavages)
		TEST_EQUAL(exp_new.getProteinIdentifications()[0].getSearchParameters().peak_mass_tolerance, exp_original.getProteinIdentifications()[0].getSearchParameters().peak_mass_tolerance)
		TEST_EQUAL(exp_new.getProteinIdentifications()[0].getSearchParameters().precursor_tolerance, exp_original.getProteinIdentifications()[0].getSearchParameters().precursor_tolerance)
		TEST_EQUAL(exp_new.getProteinIdentifications()[0].getSearchParameters().getMetaValue("label"), exp_original.getProteinIdentifications()[0].getSearchParameters().getMetaValue("label"))

		TEST_EQUAL(exp_new.getProteinIdentifications()[0].getSearchParameters().fixed_modifications.size(), exp_original.getProteinIdentifications()[0].getSearchParameters().fixed_modifications.size())
		for(Size i = 0; i < exp_original.getProteinIdentifications()[0].getSearchParameters().fixed_modifications.size(); ++i)
		{
			TEST_EQUAL(exp_new.getProteinIdentifications()[0].getSearchParameters().fixed_modifications[i], exp_original.getProteinIdentifications()[0].getSearchParameters().fixed_modifications[i])
		}
		TEST_EQUAL(exp_new.getProteinIdentifications()[0].getSearchParameters().variable_modifications.size(), exp_original.getProteinIdentifications()[0].getSearchParameters().variable_modifications.size())
		for(Size i = 0; i < exp_original.getProteinIdentifications()[0].getSearchParameters().variable_modifications.size(); ++i)
		{
			TEST_EQUAL(exp_new.getProteinIdentifications()[0].getSearchParameters().variable_modifications[i], exp_original.getProteinIdentifications()[0].getSearchParameters().variable_modifications[i])
		}


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

			TEST_EQUAL(exp_new.getInstrument().getSoftware().getName() , exp_original.getInstrument().getSoftware().getName() )
			TEST_EQUAL(exp_new.getInstrument().getSoftware().getVersion() , exp_original.getInstrument().getSoftware().getVersion() )
			TEST_EQUAL(exp_new.getInstrument().getSoftware().getMetaValue("label") , exp_original.getInstrument().getSoftware().getMetaValue("label") )

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
			TEST_EQUAL(exp_new.getInstrument().getIonDetectors()[0].getOrder(), 3)
			TEST_EQUAL(exp_new.getInstrument().getIonSources()[0].getOrder(), 0)
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
			TEST_EQUAL(exp_new.getInstrument().getMassAnalyzers()[0].getOrder() , 1 )
			TEST_EQUAL(exp_new.getInstrument().getMassAnalyzers()[1].getScanDirection() , MassAnalyzer::UP )
			TEST_EQUAL(exp_new.getInstrument().getMassAnalyzers()[1].getScanLaw() , MassAnalyzer::LINEAR )
			TEST_REAL_SIMILAR(exp_new.getInstrument().getMassAnalyzers()[1].getScanRate() , 5.555 )
			TEST_REAL_SIMILAR(exp_new.getInstrument().getMassAnalyzers()[1].getScanTime() , 6.666 )
			TEST_REAL_SIMILAR(exp_new.getInstrument().getMassAnalyzers()[1].getTOFTotalPathLength() , 7.777 )
			TEST_EQUAL(exp_new.getInstrument().getMassAnalyzers()[1].getType() , MassAnalyzer::TOF )
			TEST_EQUAL(exp_new.getInstrument().getMassAnalyzers()[1].getMetaValue("label") , "pink" )
			TEST_EQUAL(exp_new.getInstrument().getMassAnalyzers()[1].getOrder() , 2 )

			TEST_EQUAL(exp_new.getInstrument().getIonOptics() , Instrument::FRINGING_FIELD )

			//------ test if values are correct ------

			//SPECTRUM 1
			RichPeakMap::const_iterator itn(exp_new.begin());
			RichPeakMap::const_iterator ito(exp_original.begin());

		  TEST_EQUAL( itn->getRT() , ito->getRT() )
			TEST_EQUAL( itn->getMSLevel() , ito->getMSLevel() )
			TEST_EQUAL( itn->size() , ito->size() )
			for (Size i=0; i<3; ++i)
			{
				TEST_REAL_SIMILAR( itn->operator[](i).getIntensity() , ito->operator[](i).getIntensity() )
				TEST_REAL_SIMILAR( itn->operator[](i).getPosition()[0] , ito->operator[](i).getPosition()[0] )
			}
			//~ TEST_EQUAL( itn->getPrecursors()[0].getMetaValue("icon") , "Precursor" )



			TEST_EQUAL(itn->getDataProcessing().size(), ito->getDataProcessing().size())
			for (Size i=0; i<itn->getDataProcessing().size(); ++i)
			{
				TEST_EQUAL(itn->getDataProcessing()[i].getSoftware().getName(), ito->getDataProcessing()[i].getSoftware().getName())
				TEST_EQUAL(itn->getDataProcessing()[i].getSoftware().getVersion(), ito->getDataProcessing()[i].getSoftware().getVersion())
				TEST_EQUAL(itn->getDataProcessing()[i].getCompletionTime().get(), ito->getDataProcessing()[i].getCompletionTime().get())
				TEST_EQUAL(itn->getDataProcessing()[i].getMetaValue("label"), ito->getDataProcessing()[i].getMetaValue("label"))
				TEST_EQUAL(itn->getDataProcessing()[i].getProcessingActions().size(), ito->getDataProcessing()[i].getProcessingActions().size())
				std::set< DataProcessing::ProcessingAction >::const_iterator set_it_original, set_it_new;
				set_it_new = (itn->getDataProcessing()[i].getProcessingActions()).begin();
				set_it_original = ito->getDataProcessing()[i].getProcessingActions().begin();
				for (set_it_new=set_it_new; set_it_new!=itn->getDataProcessing()[i].getProcessingActions().end() && set_it_original!=ito->getDataProcessing()[i].getProcessingActions().end() ; (++set_it_new), (++set_it_original))
				{
					TEST_EQUAL(*set_it_new, *set_it_original)
				}
			}

			//SPECTRUM 2
			++itn;
			++ito;

		  TEST_EQUAL( itn->getRT() , ito->getRT() )
			TEST_EQUAL( itn->getMSLevel() , ito->getMSLevel() )

			TEST_EQUAL( itn->getPrecursors().size() , ito->getPrecursors().size() )
			for(Size i = 0; i < itn->getPrecursors().size(); ++i)
			{
				TEST_EQUAL( itn->getPrecursors()[i].getMZ() , ito->getPrecursors()[i].getMZ() )
				TEST_EQUAL( itn->getPrecursors()[i].getIntensity() , ito->getPrecursors()[i].getIntensity() )
				TEST_EQUAL( itn->getPrecursors()[i].getCharge() , ito->getPrecursors()[i].getCharge() )
				TEST_EQUAL( itn->getPrecursors()[i].getActivationEnergy() , ito->getPrecursors()[i].getActivationEnergy() )
				TEST_EQUAL( itn->getPrecursors()[i].getMetaValue("icon") , ito->getPrecursors()[i].getMetaValue("icon") )
				TEST_EQUAL( itn->getPrecursors()[i].getPossibleChargeStates().size() , ito->getPrecursors()[i].getPossibleChargeStates().size() )
				for (Size j=0; j<itn->getPrecursors()[i].getPossibleChargeStates().size(); ++j)
				{
					TEST_EQUAL( itn->getPrecursors()[i].getPossibleChargeStates()[j] , ito->getPrecursors()[i].getPossibleChargeStates()[j] )
				}
				TEST_EQUAL( itn->getPrecursors()[i].getActivationMethods().size() , ito->getPrecursors()[i].getActivationMethods().size() )
				for ( std::set<Precursor::ActivationMethod>::const_iterator amn_it(itn->getPrecursors()[i].getActivationMethods().begin()), amo_it(ito->getPrecursors()[i].getActivationMethods().begin()); amn_it != itn->getPrecursors()[i].getActivationMethods().end(); ++amn_it,++amo_it)
				{
					TEST_EQUAL( *amn_it , *amo_it )
				}
			}

			TEST_EQUAL( itn->getComment() , "bla" )
			TEST_EQUAL( itn->size() , ito->size() )
			for (Size i=0; i<3; ++i)
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
			modified_spec[0].setIntensity(566.0f);
			modified_spec[0].getPosition()[0] = 612.1;
			modified_spec[1].setIntensity(620.0f);
			modified_spec[1].getPosition()[0] = 712.1;
			modified_spec[2].setIntensity(701.0f);
			modified_spec[2].getPosition()[0] = 812.1;
			modified_spec.setRT(1.88);
			modified_spec.setMSLevel(1);
			modified_spec.getInstrumentSettings().getScanWindows()[0].begin = 3.567;
			modified_spec.getInstrumentSettings().getScanWindows()[0].end = 7.91;
			modified_spec.getInstrumentSettings().setPolarity(IonSource::POSITIVE);
			modified_spec.getInstrumentSettings().setScanMode(InstrumentSettings::SIM);
			modified_spec.getInstrumentSettings().setZoomScan(false);
			modified_spec.getInstrumentSettings().setMetaValue("label", String("please bite here"));

			modified_spec.getProducts()[1].setMZ(5);
			modified_spec.getProducts()[1].setIsolationWindowLowerOffset(6);
			modified_spec.getProducts()[1].setIsolationWindowUpperOffset(7);
			modified_spec.getProducts()[1].setMetaValue("farbe", String("erbrochengruengelb"));

			info.clear();
			acquisition.setIdentifier("1");
			acquisition.setMetaValue ("icon", String("one more icon"));
			info.push_back(acquisition);
			acquisition.setIdentifier("2");
			acquisition.setMetaValue ("label", String("yet another label"));
			info.push_back(acquisition);

			modified_spec.setAcquisitionInfo(info);
			// adding a meta data array
			modified_spec.getFloatDataArrays().clear();
			RichPeakSpectrum::FloatDataArray meta_data_array;
			meta_data_array.setName("icon");
			meta_data_array.push_back(23.0f);
			meta_data_array.push_back(42.0f);
			meta_data_array.push_back(100.001f);

			modified_spec.getFloatDataArrays().push_back(meta_data_array);


			// modify 2nd spectrum
			exp_original[1].getPrecursors()[0].setMetaValue("icon", String("NewPrecursor"));

			// update others
			exp_original.getProteinIdentifications()[0].getHits()[1].setRank( 5u );
			exp_original.getInstrument().getMassAnalyzers()[0].setOrder(2);
			exp_original.getInstrument().getMassAnalyzers()[1].setOrder(3);
			exp_original.getInstrument().getIonDetectors()[0].setOrder(4);
			exp_original.getInstrument().getIonSources()[0].setOrder(1);
			exp_original.getInstrument().setIonOptics(Instrument::EINZEL_LENS);
			ProteinIdentification::SearchParameters s = exp_original.getProteinIdentifications()[0].getSearchParameters();
			s.missed_cleavages = 66;
			exp_original.getProteinIdentifications()[0].setSearchParameters(s);

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
			TEST_EQUAL( itn->getAcquisitionInfo()[0].getIdentifier(), "1");
			TEST_EQUAL( itn->getAcquisitionInfo()[0].getMetaValue("icon"), "one more icon");
			TEST_EQUAL( itn->getAcquisitionInfo()[1].getIdentifier(), "2");
			TEST_EQUAL( itn->getAcquisitionInfo()[1].getMetaValue("label"), "yet another label");
			for (Size i=0; i<3; ++i)
			{
				TEST_REAL_SIMILAR( itn->operator[](i).getIntensity() , ito->operator[](i).getIntensity() )
				TEST_REAL_SIMILAR( itn->operator[](i).getPosition()[0] , ito->operator[](i).getPosition()[0] )
			}
			TEST_EQUAL( itn->getInstrumentSettings().getZoomScan() , ito->getInstrumentSettings().getZoomScan() )
			for(Size ps = 0; ps < itn->getProducts().size(); ++ps)
			{
				TEST_EQUAL( itn->getProducts()[ps].getMZ() , ito->getProducts()[ps].getMZ() )
				TEST_EQUAL( itn->getProducts()[ps].getIsolationWindowLowerOffset() , ito->getProducts()[ps].getIsolationWindowLowerOffset() )
				TEST_EQUAL( itn->getProducts()[ps].getIsolationWindowUpperOffset() , ito->getProducts()[ps].getIsolationWindowUpperOffset() )
				TEST_EQUAL( itn->getProducts()[ps].getMetaValue("farbe") , ito->getProducts()[ps].getMetaValue("farbe"))
			}


			//SPECTRUM 2
			++itn;
			++ito;

		  TEST_EQUAL( itn->getRT() , ito->getRT() )
			TEST_EQUAL( itn->getMSLevel() , ito->getMSLevel() )
			TEST_EQUAL( itn->getPrecursors().size() , ito->getPrecursors().size() )
			TEST_EQUAL( itn->getPrecursors()[0].getMZ() , ito->getPrecursors()[0].getMZ() )
			TEST_EQUAL( itn->getPrecursors()[0].getIntensity() , ito->getPrecursors()[0].getIntensity() )
			TEST_EQUAL( itn->getPrecursors()[0].getCharge() , ito->getPrecursors()[0].getCharge() )
			TEST_EQUAL( itn->getPrecursors()[0].getMetaValue("icon") , "NewPrecursor" )
			TEST_EQUAL( itn->getComment() , "bla" )
			TEST_EQUAL( itn->size() , ito->size() )
			for (Size i=0; i<3; ++i)
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

			//test update of others
			TEST_EQUAL(exp_new.getProteinIdentifications()[0].getHits()[1].getRank(), 5u )
			TEST_EQUAL(exp_new[0].getSourceFile().getChecksumType(), SourceFile::UNKNOWN_CHECKSUM )
			TEST_STRING_EQUAL(exp_new[0].getSourceFile().getNativeIDType(), "Waters nativeID format" )
			TEST_EQUAL(exp_new.getInstrument().getMassAnalyzers()[0].getOrder() , 2 )
			TEST_EQUAL(exp_new.getInstrument().getMassAnalyzers()[1].getOrder() , 3 )
			TEST_EQUAL(exp_new.getInstrument().getIonDetectors()[0].getOrder(), 4)
			TEST_EQUAL(exp_new.getInstrument().getIonSources()[0].getOrder(), 1)
			TEST_EQUAL(exp_new.getInstrument().getIonOptics() , Instrument::EINZEL_LENS )

			TEST_EQUAL(exp_new.getProteinIdentifications()[0].getSearchParameters().missed_cleavages, 66)
			TEST_EQUAL(exp_new.getProteinIdentifications()[0].getSearchParameters().peak_mass_tolerance, exp_original.getProteinIdentifications()[0].getSearchParameters().peak_mass_tolerance)


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



