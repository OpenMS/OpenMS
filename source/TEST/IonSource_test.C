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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/IonSource.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(IonSource, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IonSource* ptr = 0;
CHECK((IonSource()))
	ptr = new IonSource();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~IonSource()))
	delete ptr;
RESULT

    enum InletType{INLETNULL,DIRECT, BATCH, CHROMATOGRAPHY, PARTICLEBEAM, MEMBRANESEPARATOR, OPENSPLIT, JETSEPARATOR, 
										 SEPTUM, RESERVOIR, MOVINGBELT, MOVINGWIRE, FLOWINJECTIONANALYSIS, ELECTROSPRAYINLET, 
										 THERMOSPRAYINLET, INFUSION, CONTINUOUSFLOWFASTATOMBOMBARDMENT, INDUCTIVELYCOUPLEDPLASMA};
	    /// ionization method
	    enum IonizationMethod{IONMETHODNULL,ESI, EI, CI, FAB, TSP, LD, FD, FI, PD, SI, TI, API,
													ISI, CID, CAD, HN, APCI, APPI, ICP};
      /// polarity of the ion source
      enum Polarity{POLNULL,POSITIVE, NEGATIVE};

CHECK((InletType getInletType() const))
  IonSource tmp;
  TEST_EQUAL(tmp.getInletType(),IonSource::INLETNULL);
RESULT

CHECK((void setInletType(InletType inlet_type)))
  IonSource tmp;
  tmp.setInletType(IonSource::DIRECT);
  TEST_EQUAL(tmp.getInletType(),IonSource::DIRECT);
RESULT

CHECK((IonizationMethod getIonizationMethod() const))
  IonSource tmp;
  TEST_EQUAL(tmp.getIonizationMethod(),IonSource::IONMETHODNULL);
RESULT

CHECK((void setIonizationMethod(IonizationMethod ionization_type)))
  IonSource tmp;
  tmp.setIonizationMethod(IonSource::ESI);
  TEST_EQUAL(tmp.getIonizationMethod(),IonSource::ESI);
RESULT

CHECK((Polarity getPolarity() const))
  IonSource tmp;
  TEST_EQUAL(tmp.getPolarity(),IonSource::POLNULL);
RESULT

CHECK((void setPolarity(Polarity polarity)))
	IonSource tmp;
  tmp.setPolarity(IonSource::POSITIVE);
  TEST_EQUAL(tmp.getPolarity(),IonSource::POSITIVE);
RESULT

CHECK((IonSource(const IonSource& source)))
  IonSource tmp;
  tmp.setInletType(IonSource::DIRECT);
  tmp.setIonizationMethod(IonSource::ESI);
  tmp.setPolarity(IonSource::POSITIVE);
  tmp.setMetaValue("label",String("label"));
  	
  IonSource tmp2(tmp);
  TEST_EQUAL(tmp2.getPolarity(),IonSource::POSITIVE);
  TEST_EQUAL(tmp2.getInletType(),IonSource::DIRECT);
  TEST_EQUAL(tmp2.getIonizationMethod(),IonSource::ESI);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
RESULT

CHECK((IonSource& operator= (const IonSource& source)))
  IonSource tmp;
  tmp.setInletType(IonSource::DIRECT);
  tmp.setIonizationMethod(IonSource::ESI);
  tmp.setPolarity(IonSource::POSITIVE);
  tmp.setMetaValue("label",String("label"));
  
  IonSource tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getPolarity(),IonSource::POSITIVE);
  TEST_EQUAL(tmp2.getInletType(),IonSource::DIRECT);
  TEST_EQUAL(tmp2.getIonizationMethod(),IonSource::ESI);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
  
  tmp2 = IonSource();
  TEST_EQUAL(tmp2.getPolarity(),IonSource::POLNULL);
  TEST_EQUAL(tmp2.getInletType(),IonSource::INLETNULL);
  TEST_EQUAL(tmp2.getIonizationMethod(),IonSource::IONMETHODNULL);
  TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
RESULT

CHECK((bool operator== (const IonSource& rhs) const))
  IonSource edit,empty;
  
  TEST_EQUAL(edit==empty,true);
  
  edit = empty;
  edit.setInletType(IonSource::DIRECT);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setIonizationMethod(IonSource::ESI);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setPolarity(IonSource::POSITIVE);
	TEST_EQUAL(edit==empty,false);
	
	edit = empty;
	edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit==empty,false);
RESULT

CHECK((bool operator!= (const IonSource& rhs) const))
  IonSource edit,empty;
  
  TEST_EQUAL(edit!=empty,false);
  
  edit = empty;
  edit.setInletType(IonSource::DIRECT);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setIonizationMethod(IonSource::ESI);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setPolarity(IonSource::POSITIVE);
	TEST_EQUAL(edit!=empty,true);

	edit = empty;
	edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit!=empty,true);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



