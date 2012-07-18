// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>

using namespace OpenMS;
using namespace std;

START_TEST(TraMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


TraMLFile * ptr = 0;
TraMLFile * nullPointer = 0;

START_SECTION((TraMLFile()))
{
  ptr = new TraMLFile;
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~TraMLFile()))
{
  delete ptr;
}
END_SECTION

START_SECTION((void load(const String &filename, TargetedExperiment & id)))
{
  TraMLFile file;
  TargetedExperiment exp;
  file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.TraML"), exp);
}
END_SECTION

START_SECTION((void store(const String &filename, const TargetedExperiment &id) const))
{
  TraMLFile file;

  //load map
  TargetedExperiment exp_original;
  file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.TraML"), exp_original);

  //store map
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  file.store(tmp_filename, exp_original);
  //file.store("include.traML", TargetedExperiment());
  //load written map
  TargetedExperiment exp;
  file.load(tmp_filename, exp);

  //test if everything worked
  TEST_EQUAL(exp == exp_original, true)
}
END_SECTION

START_SECTION((void equal()))
{
  TraMLFile file;

  TargetedExperiment exp_original;
  TargetedExperiment exp_second;
  file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.TraML"), exp_original);
  file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.TraML"), exp_second);

  TEST_EQUAL(exp_second == exp_original, true)
}
END_SECTION

START_SECTION((void assign()))
{
  TraMLFile file;

  //load map
  TargetedExperiment exp_original;
  TargetedExperiment exp_added;
  file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.TraML"), exp_original);

  //store map
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  exp_added = exp_original;
  TEST_EQUAL(exp_original.getTargetCVTerms().getCVTerms().size(), 1)
  TEST_EQUAL(exp_added.getTargetCVTerms().getCVTerms().size(), 1)

  TEST_EQUAL(exp_added == exp_original, true)
}
END_SECTION

START_SECTION((void add()))
{
  TraMLFile file;

  //load map
  TargetedExperiment exp_original;
  TargetedExperiment exp_added;
  file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.TraML"), exp_original);

  //store map
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  exp_added += exp_original;

  TEST_EQUAL(exp_added == exp_original, true)
}
END_SECTION

START_SECTION([EXTRA] bool isValid(const String & filename, std::ostream & os = std::cerr))
{
  std::string tmp_filename;
  TraMLFile file;
  TargetedExperiment e;

//written empty file
  NEW_TMP_FILE(tmp_filename);
  file.store(tmp_filename, e);
  TEST_EQUAL(file.isValid(tmp_filename), true);

//written filled file
  NEW_TMP_FILE(tmp_filename);
  file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.TraML"), e);
  file.store(tmp_filename, e);
  TEST_EQUAL(file.isValid(tmp_filename), true);
}
END_SECTION

START_SECTION(bool isSemanticallyValid(const String & filename, StringList & errors, StringList & warnings))
{
  std::string tmp_filename;
  TraMLFile file;
  StringList errors, warnings;
  TargetedExperiment e;

//written empty file
  NEW_TMP_FILE(tmp_filename);
  file.store(tmp_filename, e);
  TEST_EQUAL(file.isSemanticallyValid(tmp_filename, errors, warnings), true);
  TEST_EQUAL(errors.size(), 0)
  TEST_EQUAL(warnings.size(), 0)

//written filled file
  NEW_TMP_FILE(tmp_filename);
  file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.TraML"), e);
  file.store(tmp_filename, e);
//TEST_EQUAL(file.isSemanticallyValid(tmp_filename, errors, warnings),true);
  TEST_EQUAL(errors.size(), 0)
  TEST_EQUAL(warnings.size(), 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
