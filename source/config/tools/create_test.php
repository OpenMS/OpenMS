<?

	########################usage###############################
	if ($argc<2 OR $argc>3)
	{
?>
Usage: create_test.php <header_file> [class]
  <header_file> -- the header of the class to create the test from
	[class]       -- if class is given, only methods of that class are tested
  
<?
	exit;	
	}
	
	// if class set -> use it
  if ($argc==3) $class = $argv[2];
	// else parse it from filename
	else $class = substr($argv[1],strrpos($argv[1],'/')+1,-2);

?>
// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-<? print date("Y"); ?> -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/TODO/<? print $class; ?>.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(<? print $class; ?>, "$<? print "Id"; ?>$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

<? print $class; ?>* ptr = 0;
CHECK(<? print $class; ?>())
	ptr = new <? print $class; ?>();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~<? print $class; ?>())
	delete ptr;
RESULT

<?
passthru('./check_test '.$argv[1].' '.$argv[2].' -c');
?>

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



