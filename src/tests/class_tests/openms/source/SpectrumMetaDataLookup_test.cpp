// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/SpectrumMetaDataLookup.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SpectrumMetaDataLookup, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SpectrumMetaDataLookup* ptr = nullptr;
SpectrumMetaDataLookup* null_ptr = nullptr;

START_SECTION((SpectrumMetaDataLookup()))
{
	ptr = new SpectrumMetaDataLookup();
	TEST_NOT_EQUAL(ptr, null_ptr);
  TEST_EQUAL(ptr->empty(), true);
}
END_SECTION

START_SECTION((~SpectrumMetaDataLookup()))
{
	delete ptr;
}
END_SECTION

vector<MSSpectrum> spectra;
MSSpectrum spectrum;
spectrum.setNativeID("spectrum=0");
spectrum.setRT(1.0);
spectrum.setMSLevel(1);
spectra.push_back(spectrum);
spectrum.setNativeID("spectrum=1");
spectrum.setRT(2.0);
spectrum.setMSLevel(2);
Precursor prec;
prec.setMZ(1000.0);
prec.setCharge(2);
spectrum.getPrecursors().push_back(prec);
spectra.push_back(spectrum);
spectrum.setNativeID("spectrum=2");
spectrum.setRT(3.0);
spectrum.setMSLevel(2);
prec.setMZ(500.0);
prec.setCharge(3);
spectrum.getPrecursors()[0] = prec;
spectra.push_back(spectrum);

SpectrumMetaDataLookup lookup;

START_SECTION((template <typename SpectrumContainer> void readSpectra(const SpectrumContainer&, const String&, bool)))
{
  lookup.readSpectra(spectra, SpectrumLookup::default_scan_regexp, true);
  TEST_EQUAL(lookup.empty(), false);
}
END_SECTION

START_SECTION((void getSpectrumMetaData(Size, SpectrumMetaData&) const))
{
  SpectrumMetaDataLookup::SpectrumMetaData meta;
  lookup.getSpectrumMetaData(0, meta);
  TEST_EQUAL(meta.rt, 1.0);
  TEST_EQUAL(meta.ms_level, 1);
  TEST_EQUAL(meta.native_id, "spectrum=0");
  TEST_EQUAL(meta.scan_number, 0);

  lookup.getSpectrumMetaData(1, meta);
  TEST_EQUAL(meta.rt, 2.0);
  TEST_EQUAL(meta.precursor_rt, 1.0);
  TEST_EQUAL(meta.precursor_mz, 1000.0);
  TEST_EQUAL(meta.precursor_charge, 2);
  TEST_EQUAL(meta.ms_level, 2);
  TEST_EQUAL(meta.native_id, "spectrum=1");
  TEST_EQUAL(meta.scan_number, 1);
}
END_SECTION

START_SECTION((static void getSpectrumMetaData(const MSSpectrum&, SpectrumMetaData&, const boost::regex&, const map<Size, double>&)))
{
  SpectrumMetaDataLookup::SpectrumMetaData meta;
  SpectrumMetaDataLookup::getSpectrumMetaData(spectrum, meta);
  TEST_EQUAL(meta.rt, 3.0);
  TEST_EQUAL(meta.precursor_mz, 500.0);
  TEST_EQUAL(meta.precursor_charge, 3);
  TEST_EQUAL(meta.ms_level, 2);
  TEST_EQUAL(meta.native_id, "spectrum=2");
  TEST_EQUAL(meta.scan_number, -1); // not extracted

  map<Size, double> precursor_rts;
  precursor_rts[1] = 1.0;
  boost::regex scan_regexp("=(?<SCAN>\\d+)$");
  SpectrumMetaDataLookup::getSpectrumMetaData(spectrum, meta, scan_regexp,
                                              precursor_rts);
  TEST_EQUAL(meta.precursor_rt, 1.0);
  TEST_EQUAL(meta.scan_number, 2);
}
END_SECTION

START_SECTION((void getSpectrumMetaData(const String&, SpectrumMetaData&, MetaDataFlags) const))
{
  SpectrumMetaDataLookup::SpectrumMetaData meta;
  lookup.addReferenceFormat(SpectrumLookup::default_scan_regexp);
  lookup.getSpectrumMetaData("scan_number=1", meta);
  TEST_EQUAL(meta.rt, 2.0);
  TEST_EQUAL(meta.native_id, "spectrum=1");

  lookup.addReferenceFormat(R"(rt=(?<RT>\d+(\.\d+)?),mz=(?<MZ>\d+(\.\d+)?))");
  SpectrumMetaDataLookup::SpectrumMetaData meta2;
  SpectrumMetaDataLookup::MetaDataFlags flags =
    (SpectrumMetaDataLookup::MDF_RT | SpectrumMetaDataLookup::MDF_PRECURSORMZ);
  // no actual look-up of the spectrum necessary:
  lookup.getSpectrumMetaData("rt=5.0,mz=1000.0", meta2, flags);
  TEST_EQUAL(meta2.rt, 5.0);
  TEST_EQUAL(meta2.precursor_mz, 1000.0);
  TEST_EQUAL(meta2.precursor_charge, 0);
  TEST_EQUAL(meta2.native_id, "");

  // look-up of the spectrum necessary:
  SpectrumMetaDataLookup::SpectrumMetaData meta3;
  lookup.getSpectrumMetaData("rt=2.0,mz=1000.0", meta3);
  TEST_EQUAL(meta3.rt, 2.0);
  TEST_EQUAL(meta3.precursor_mz, 1000.0);
  TEST_EQUAL(meta3.precursor_charge, 2);
  TEST_EQUAL(meta3.native_id, "spectrum=1");

  TEST_EXCEPTION(Exception::ElementNotFound, lookup.getSpectrumMetaData("rt=5.0,mz=1000.0", meta3));
}
END_SECTION


START_SECTION((bool addMissingRTsToPeptideIDs(vector<PeptideIdentification>& peptides, const String& filename, bool stop_on_error)))
{
  vector<PeptideIdentification> peptides(1);
  peptides[0].setRT(1.0);
  String filename = "this_file_does_not_exist.mzML";
  // no missing RTs -> no attempt to load mzML file:
  SpectrumMetaDataLookup::addMissingRTsToPeptideIDs(peptides, filename);
  TEST_EQUAL(peptides[0].getRT(), 1.0);

  peptides.resize(2);
  peptides[0].setSpectrumReference( "index=0");
  peptides[1].setSpectrumReference( "index=2");
  filename = OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML");
  SpectrumMetaDataLookup::addMissingRTsToPeptideIDs(peptides, filename);
  TEST_EQUAL(peptides[0].getRT(), 1.0); // this doesn't get overwritten
  TEST_REAL_SIMILAR(peptides[1].getRT(), 5.3);
}
END_SECTION


START_SECTION((bool addMissingSpectrumReferences(vector<PeptideIdentification>& peptides, 
  const String& filename, 
  bool stop_on_error, 
  bool override_spectra_data, 
  bool override_spectra_references, 
  vector<ProteinIdentification> proteins)))
{
  vector<PeptideIdentification> peptides(1);
  peptides[0].setRT(5.1);
  peptides[0].setSpectrumReference( "index=666");
  String filename = "this_file_does_not_exist.mzML";
  SpectrumMetaDataLookup lookup;
  // missing file -> exception, no non-effective executions
  TEST_EXCEPTION(Exception::FileNotFound, SpectrumMetaDataLookup::addMissingSpectrumReferences(
    peptides, filename, false, false));
  // no lookup, no spectrum_references
  TEST_EQUAL(peptides[0].getSpectrumReference(), "index=666");

  peptides.resize(2);
  peptides[1].setRT(5.3);
  filename = OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML");

  SpectrumMetaDataLookup::addMissingSpectrumReferences(peptides, filename, false, false, false);

  TEST_EQUAL(peptides[0].getSpectrumReference(), "index=666"); // no overwrite
  TEST_EQUAL(peptides[1].getSpectrumReference(), "index=2");

  SpectrumMetaDataLookup::addMissingSpectrumReferences(peptides, filename, false, true, true);

  TEST_EQUAL(peptides[0].getSpectrumReference(), "index=0"); // gets updated
  TEST_EQUAL(peptides[1].getSpectrumReference(), "index=2");
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
