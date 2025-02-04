// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SpectrumSettings, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SpectrumSettings* ptr = nullptr;
SpectrumSettings* nullPointer = nullptr;

START_SECTION((SpectrumSettings()))
	ptr = new SpectrumSettings();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~SpectrumSettings()))
	delete ptr;
END_SECTION

START_SECTION((const String& getNativeID() const))
	SpectrumSettings tmp;
	TEST_STRING_EQUAL(tmp.getNativeID(),"");
END_SECTION

START_SECTION((void setNativeID(const String& native_id)))
	SpectrumSettings tmp;
	tmp.setNativeID("nid");
	TEST_STRING_EQUAL(tmp.getNativeID(),"nid");
END_SECTION

START_SECTION((const std::vector<DataProcessing>& getDataProcessing() const))
  SpectrumSettings tmp;
  TEST_EQUAL(tmp.getDataProcessing().size(),0);
END_SECTION

START_SECTION((void setDataProcessing(const std::vector< DataProcessing > &data_processing)))
  SpectrumSettings tmp;
  std::vector<DataProcessingPtr > dummy;
  dummy.resize(1);
  tmp.setDataProcessing(dummy);
  TEST_EQUAL(tmp.getDataProcessing().size(),1);
END_SECTION

START_SECTION((std::vector<DataProcessing>& getDataProcessing()))
  SpectrumSettings tmp;
  tmp.getDataProcessing().resize(1);
  TEST_EQUAL(tmp.getDataProcessing().size(),1);
END_SECTION

START_SECTION((AcquisitionInfo& getAcquisitionInfo()))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getAcquisitionInfo().empty(), true);
END_SECTION

START_SECTION((void setAcquisitionInfo(const AcquisitionInfo& acquisition_info)))
	SpectrumSettings tmp;
	AcquisitionInfo ai;
	ai.setMethodOfCombination("test");
	tmp.setAcquisitionInfo(ai);
	TEST_EQUAL(tmp.getAcquisitionInfo() == AcquisitionInfo(), false);
END_SECTION

START_SECTION((const AcquisitionInfo& getAcquisitionInfo() const))
	SpectrumSettings tmp;
	tmp.getAcquisitionInfo().setMethodOfCombination("test");
	TEST_EQUAL(tmp.getAcquisitionInfo() == AcquisitionInfo(), false);
END_SECTION

START_SECTION((SourceFile& getSourceFile()))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getSourceFile()==SourceFile(), true);
END_SECTION

START_SECTION((void setSourceFile(const SourceFile& source_file)))
	SpectrumSettings tmp;
	SourceFile sf;
	sf.setNameOfFile("test");
	tmp.setSourceFile(sf);
	TEST_EQUAL(tmp.getSourceFile()==SourceFile(), false);
END_SECTION

START_SECTION((const SourceFile& getSourceFile() const))
	SpectrumSettings tmp;
	tmp.getSourceFile().setNameOfFile("test");
	TEST_EQUAL(tmp.getSourceFile()==SourceFile(), false);
END_SECTION

START_SECTION((const InstrumentSettings& getInstrumentSettings() const))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getInstrumentSettings()==InstrumentSettings(), true);	  
END_SECTION

START_SECTION((void setInstrumentSettings(const InstrumentSettings& instrument_settings)))
	SpectrumSettings tmp;
	InstrumentSettings is;
	is.getScanWindows().resize(1);
	tmp.setInstrumentSettings(is);
	TEST_EQUAL(tmp.getInstrumentSettings()==InstrumentSettings(), false);
END_SECTION

START_SECTION((InstrumentSettings& getInstrumentSettings()))
	SpectrumSettings tmp;
	tmp.getInstrumentSettings().getScanWindows().resize(1);
	TEST_EQUAL(tmp.getInstrumentSettings()==InstrumentSettings(), false);	
END_SECTION

START_SECTION((const std::vector<Precursor>& getPrecursors() const))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getPrecursors().size(),0);	  
END_SECTION

START_SECTION((void setPrecursors(const std::vector<Precursor>& precursors)))
	SpectrumSettings tmp;
	tmp.setPrecursors(vector<Precursor>(2));
	TEST_EQUAL(tmp.getPrecursors().size(), 2);
END_SECTION

START_SECTION((std::vector<Precursor>& getPrecursors()))
	SpectrumSettings tmp;
	tmp.getPrecursors().resize(4);
	TEST_EQUAL(tmp.getPrecursors().size(), 4);	
END_SECTION

START_SECTION((const std::vector<Product>& getProducts() const))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getProducts().size(),0);	  
END_SECTION

START_SECTION((void setProducts(const std::vector<Product>& products)))
	SpectrumSettings tmp;
	tmp.setProducts(vector<Product>(2));
	TEST_EQUAL(tmp.getProducts().size(), 2);
END_SECTION

START_SECTION((std::vector<Product>& getProducts()))
	SpectrumSettings tmp;
	tmp.getProducts().resize(4);
	TEST_EQUAL(tmp.getProducts().size(), 4);	
END_SECTION

START_SECTION((SpectrumType getType() const))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getType(), SpectrumSettings::UNKNOWN);	  
END_SECTION

START_SECTION((void setType(SpectrumType type)))
	SpectrumSettings tmp;
	tmp.setType(SpectrumSettings::CENTROID);
	TEST_EQUAL(tmp.getType(), SpectrumSettings::CENTROID);
END_SECTION

START_SECTION((const String& getComment() const))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getComment(), "");
END_SECTION

START_SECTION((void setComment(const String& comment)))
	SpectrumSettings tmp;
	tmp.setComment("bla");
	TEST_EQUAL(tmp.getComment(), "bla");
END_SECTION

START_SECTION((const std::vector<PeptideIdentification>& getPeptideIdentifications() const))
	SpectrumSettings tmp;
	vector<PeptideIdentification> vec(tmp.getPeptideIdentifications());
	TEST_EQUAL(vec.size(),0);
END_SECTION

START_SECTION((void setPeptideIdentifications(const std::vector<PeptideIdentification>& identifications)))
	SpectrumSettings tmp;
	vector<PeptideIdentification> vec;
	
	tmp.setPeptideIdentifications(vec);
	TEST_EQUAL(tmp.getPeptideIdentifications().size(),0);
	
	PeptideIdentification dbs;
	vec.push_back(dbs);
	tmp.setPeptideIdentifications(vec);
	TEST_EQUAL(tmp.getPeptideIdentifications().size(),1);
END_SECTION

START_SECTION((std::vector<PeptideIdentification>& getPeptideIdentifications()))
	SpectrumSettings tmp;
	vector<PeptideIdentification> vec;
	
	tmp.getPeptideIdentifications().resize(1);
	TEST_EQUAL(tmp.getPeptideIdentifications().size(),1);
END_SECTION

START_SECTION((SpectrumSettings& operator= (const SpectrumSettings& source)))
  SpectrumSettings tmp;
  tmp.setMetaValue("bla","bluff");
	tmp.getAcquisitionInfo().setMethodOfCombination("test");
	tmp.getInstrumentSettings().getScanWindows().resize(1);
	tmp.getPrecursors().resize(1);
	tmp.getProducts().resize(1);
	tmp.getPeptideIdentifications().resize(1);
	tmp.setType(SpectrumSettings::CENTROID);
	tmp.setComment("bla");
	tmp.setNativeID("nid");
	tmp.getDataProcessing().resize(1);
	
	SpectrumSettings tmp2(tmp);
	TEST_EQUAL(tmp2.getComment(), "bla");
	TEST_EQUAL(tmp2.getType(), SpectrumSettings::CENTROID);
	TEST_EQUAL(tmp2.getPeptideIdentifications().size(), 1);	
	TEST_EQUAL(tmp2.getPrecursors().size(),1);	
	TEST_EQUAL(tmp2.getProducts().size(),1);	
	TEST_EQUAL(tmp2.getInstrumentSettings()==InstrumentSettings(), false);
	TEST_EQUAL(tmp2.getAcquisitionInfo()==AcquisitionInfo(), false);
	TEST_EQUAL(tmp2.getAcquisitionInfo().empty(), true);
	TEST_STRING_EQUAL(tmp2.getNativeID(),"nid");
	TEST_EQUAL(tmp2.getDataProcessing().size(),1);
	TEST_EQUAL(tmp2.getMetaValue("bla")=="bluff",true);
END_SECTION

START_SECTION((SpectrumSettings(const SpectrumSettings& source)))
  SpectrumSettings tmp;
	tmp.getAcquisitionInfo().setMethodOfCombination("test");
	tmp.getInstrumentSettings().getScanWindows().resize(1);
	tmp.getPrecursors().resize(1);
	tmp.getProducts().resize(1);
	tmp.setType(SpectrumSettings::CENTROID);
	tmp.setComment("bla");
	tmp.getPeptideIdentifications().resize(1);
	tmp.setNativeID("nid");
	tmp.getDataProcessing().resize(1);
	tmp.setMetaValue("bla","bluff");
	
	SpectrumSettings tmp2;
	tmp2 = tmp;
	TEST_EQUAL(tmp2.getComment(), "bla");
	TEST_EQUAL(tmp2.getType(), SpectrumSettings::CENTROID);
	TEST_EQUAL(tmp2.getPrecursors().size(), 1);
	TEST_EQUAL(tmp2.getProducts().size(), 1)
	TEST_EQUAL(tmp2.getInstrumentSettings()==InstrumentSettings(), false);	
	TEST_EQUAL(tmp2.getAcquisitionInfo().empty(), true);
	TEST_EQUAL(tmp2.getAcquisitionInfo()==AcquisitionInfo(), false);
	TEST_EQUAL(tmp2.getPeptideIdentifications().size(), 1);	
	TEST_STRING_EQUAL(tmp2.getNativeID(),"nid");
	TEST_EQUAL(tmp2.getDataProcessing().size(),1);
	TEST_STRING_EQUAL(tmp2.getMetaValue("bla"),"bluff");


	tmp2 = SpectrumSettings();
	TEST_EQUAL(tmp2.getComment(), "");
	TEST_EQUAL(tmp2.getType(), SpectrumSettings::UNKNOWN);
	TEST_EQUAL(tmp2.getPrecursors().size(),0);	
	TEST_EQUAL(tmp2.getProducts().size(),0);	
	TEST_EQUAL(tmp2.getInstrumentSettings()==InstrumentSettings(), true);	
	TEST_EQUAL(tmp2.getAcquisitionInfo().empty(), true);
	TEST_EQUAL(tmp2.getPeptideIdentifications().size(), 0);	
	TEST_STRING_EQUAL(tmp2.getNativeID(),"");
	TEST_EQUAL(tmp2.getDataProcessing().size(),0);
	TEST_EQUAL(tmp2.metaValueExists("bla"),false);

END_SECTION

START_SECTION((bool operator== (const SpectrumSettings& rhs) const))
  SpectrumSettings edit, empty;
  
  TEST_TRUE(edit == empty);
  
	edit.getAcquisitionInfo().setMethodOfCombination("test");
	TEST_EQUAL(edit==empty, false);
	
	edit = empty;
	edit.setNativeID("nid");
	TEST_EQUAL(edit==empty, false);
	
	edit = empty;
	edit.getInstrumentSettings().getScanWindows().resize(1);
	TEST_EQUAL(edit==empty, false);
	
	edit = empty;
	edit.getPrecursors().resize(1);
	TEST_EQUAL(edit==empty, false);
	
	edit = empty;
	edit.setType(SpectrumSettings::CENTROID);
	TEST_EQUAL(edit==empty, false);
	
	edit = empty;
	edit.setComment("bla");
	TEST_EQUAL(edit==empty, false);

	edit = empty;
	edit.getPrecursors().resize(1);
	TEST_EQUAL(edit==empty, false);

	edit = empty;
	edit.getProducts().resize(1);
	TEST_EQUAL(edit==empty, false);
	
	edit = empty;
	edit.getPeptideIdentifications().resize(1);
	TEST_EQUAL(edit==empty, false);

	edit = empty;
    DataProcessingPtr dp = boost::shared_ptr<DataProcessing>(new DataProcessing); 
	edit.getDataProcessing().push_back(dp);
	TEST_EQUAL(edit==empty, false);

	edit = empty;
	edit.setMetaValue("bla","bluff");
	TEST_EQUAL(edit==empty, false);

END_SECTION

START_SECTION((bool operator!= (const SpectrumSettings& rhs) const))
  SpectrumSettings edit, empty;
  
  TEST_EQUAL(edit!=empty, false);
  
	edit.getAcquisitionInfo().setMethodOfCombination("test");
	TEST_FALSE(edit == empty);
	
	edit = empty;
	edit.setNativeID("nid");
	TEST_FALSE(edit == empty);
	
	edit = empty;
	edit.getInstrumentSettings().getScanWindows().resize(1);
	TEST_FALSE(edit == empty);
	
	edit = empty;
	edit.getPrecursors().resize(1);
	TEST_FALSE(edit == empty);
	
	edit = empty;
	edit.setType(SpectrumSettings::CENTROID);
	TEST_FALSE(edit == empty);
	
	edit = empty;
	edit.setComment("bla");
	TEST_FALSE(edit == empty);

	edit = empty;
	edit.getPrecursors().resize(1);
	TEST_FALSE(edit == empty);

	edit = empty;
	edit.getProducts().resize(1);
	TEST_FALSE(edit == empty);

	edit = empty;
	edit.getPeptideIdentifications().resize(1);
	TEST_FALSE(edit == empty);

	edit = empty;
    DataProcessingPtr dp = boost::shared_ptr<DataProcessing>(new DataProcessing); 
	edit.getDataProcessing().push_back(dp);
	TEST_FALSE(edit == empty);

	edit = empty;
	edit.setMetaValue("bla","bluff");
	TEST_FALSE(edit == empty);


END_SECTION

START_SECTION((void unify(const SpectrumSettings &rhs)))
{
  SpectrumSettings org, appended;

  // MetaValues
  org.setMetaValue(1, "will be gone");
  org.setMetaValue(2, "will be still present");
  appended.setMetaValue(1, "will overwrite org comment");

  // Comments
  org.setComment("Original Comment");
  appended.setComment("Appended to org Commment");

  // Precursors
  Precursor org_precursor;
  org_precursor.setMZ(1.0);
  org.getPrecursors().push_back(org_precursor);

  Precursor appended_precursor;
  appended_precursor.setMZ(2.0);
  appended.getPrecursors().push_back(appended_precursor);

  // type
  org.setType(SpectrumSettings::PROFILE);
  appended.setType(SpectrumSettings::PROFILE);

  // Products
  Product org_product;
  org_product.setMZ(1.0);
  org.getProducts().push_back(org_product);

  Product appended_product;
  appended_product.setMZ(2.0);
  appended.getProducts().push_back(appended_product);

  // Identifications
  PeptideIdentification org_ident;
  org_ident.setIdentifier("org_ident");
  org.getPeptideIdentifications().push_back(org_ident);

  PeptideIdentification appended_ident;
  appended_ident.setIdentifier("appended_ident");
  appended.getPeptideIdentifications().push_back(appended_ident);

  // DataProcessings
  DataProcessingPtr org_processing = boost::shared_ptr<DataProcessing>(new DataProcessing);
  Software org_software;
  org_software.setName("org_software");
  org_processing->setSoftware(org_software);
  org.getDataProcessing().push_back(org_processing);

  DataProcessingPtr appended_processing = boost::shared_ptr<DataProcessing>(new DataProcessing);
  Software appended_software;
  appended_software.setName("appended_software");
  appended_processing->setSoftware(appended_software);
  appended.getDataProcessing().push_back(appended_processing);

  org.unify(appended);

  // MetaValues
  TEST_EQUAL(org.getMetaValue(1), "will overwrite org comment")
  TEST_EQUAL(org.getMetaValue(2), "will be still present")

  // Comments
  TEST_EQUAL(org.getComment(), "Original CommentAppended to org Commment")

  // Precursors
  TEST_EQUAL(org.getPrecursors().size(), 2)
  ABORT_IF(org.getPrecursors().size()!=2)

  TEST_EQUAL(org.getPrecursors()[0].getMZ(), 1.0)
  TEST_EQUAL(org.getPrecursors()[1].getMZ(), 2.0)

  // type
  TEST_EQUAL(org.getType(), SpectrumSettings::PROFILE)

  // Products
  TEST_EQUAL(org.getProducts().size(), 2)
  ABORT_IF(org.getProducts().size()!=2)

  TEST_EQUAL(org.getProducts()[0].getMZ(), 1.0)
  TEST_EQUAL(org.getProducts()[1].getMZ(), 2.0)

  // Identifications
  TEST_EQUAL(org.getPeptideIdentifications().size(), 2)
  ABORT_IF(org.getPeptideIdentifications().size()!=2)

  TEST_EQUAL(org.getPeptideIdentifications()[0].getIdentifier(), "org_ident")
  TEST_EQUAL(org.getPeptideIdentifications()[1].getIdentifier(), "appended_ident")

  // Identifications
  TEST_EQUAL(org.getDataProcessing().size(), 2)
  ABORT_IF(org.getDataProcessing().size()!=2)

  TEST_EQUAL(org.getDataProcessing()[0]->getSoftware().getName(), "org_software")
  TEST_EQUAL(org.getDataProcessing()[1]->getSoftware().getName(), "appended_software")

  // unify should set Type to unknown in case of type mismatch
  SpectrumSettings empty;
  empty.setType(SpectrumSettings::CENTROID);
  org.unify(empty);

  TEST_EQUAL(org.getType(), SpectrumSettings::UNKNOWN)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



