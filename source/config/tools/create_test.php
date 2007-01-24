<?
# -*- Mode: C++; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# --------------------------------------------------------------------------
# $Maintainer: Marc Sturm $
# --------------------------------------------------------------------------

	error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);
	
	include "common_functions.php";

	########################usage###############################
	if ($argc!=3)
	{
		print "Usage: create_test.php <Absolut path to OpenMS> <Absolut path to header>\n";
  				"  <header_file> -- the header of the class to create the test from.\n";
 		exit;
 	}
	
	$path = $argv[1];
	$file = $argv[2];
	$class = substr(basename($file),0,-2);

	#load file info
	$class_info = getClassInfo($path,$file,0);	
	
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

START_TEST(<? print $class; ?>, "$Id$")

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
foreach ($class_info["public-long"] as $c)
{
	if (trim($c) != $class_info["classname"]."()")
	{
		print "CHECK($c)\n";
		print "  // TODO\n";
		print "RESULT\n";
		print "\n";
	}
}
?>

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



