<?
# -*- Mode: C++; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
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

if ($argc<2 OR $argc>2)
{
	print "\n\nUsage: tests_missing.php <Path to OpenMS>\n\n";
	exit;	
}

$ignore = array(
	"/FORMAT/HANDLERS/",
	"/VISUAL/",
	"/CONCEPT/Types.h",
	"/CONCEPT/Macros.h",
	"/CONCEPT/Exception.h",
	"/CONCEPT/Benchmark.h",
	"/CONCEPT/Constants.h",
	"/config.h"
	);

$path = $argv[1];
exec("find $path/include/ -name \"*.h\"", $headers);
exec("find $path/source/TEST/ -name \"*_test.C\"", $tests);

foreach ($tests as $test)
{
	$test_names[] = trim(substr($test, strrpos($test,"/")+1,-7));
}
sort($test_names);

foreach ($headers as $header)
{
	//ignore handlersstuff that needs no tests
	foreach ($ignore as $i)
	{
		if (strpos($header,$i)!==FALSE)
		{
			continue(2);
		}		
	}

	$maintainer_line = "";
	$handle = fopen($header, "r");
	while (!feof($handle)) 
	{
	  $line = fgets($handle);
		if (strpos($line,"\$Maintainer")!==FALSE)
		{
			$maintainer_line = $line;
			break;
		}
	}
	fclose($handle);
	$begin = strpos($maintainer_line,":");
	$end = strrpos($maintainer_line,"\$");
	$people[trim(substr($maintainer_line, $begin+1,$end-$begin-1))][] = $header;
}

foreach ($people as $name => $headers)
{
	$missing_tests = false;
	foreach ($headers as $header)
	{
		$class = substr($header, strrpos($header,"/")+1,-2);
		
		if (!in_array($class,$test_names))
		{
			if ($missing_tests == false)
			{
				if ($name=="")
				{
					$name = "NO MAINTAINER";
				}
				print "\n\n-- MISSING TESTS: $name --\n\n";
				$missing_tests = true;
			}
			print "$class\n";
		}
	}
}

//----------------------------------------Which tests are called -------------------------------------------

print "\n\n-- TESTS NOT IN MAKEFILE --\n\n";

$makefile = file("$path/source/TEST/Makefile");
foreach($makefile as $line)
{
	$line = trim($line);
	if (strpos($line,"_test")!== FALSE && strpos($line,":")=== FALSE && $line[0]!="#" && $line[0]!="@")
	{
		$called_tests[] = strtr($line,array("\\"=>""," "=>"", "	"=>""));
	}
}

foreach($test_names as $name)
{
	if (!in_array($name."_test",$called_tests))
	{
		print $name."\n";
	}
}

//----------------------------------------Which tests do nothing -------------------------------------------

print "\n\n-- TESTS THAT ARE PROBABLY NOT COMPLETE --\n\n";
print "Name CHECKs TESTS\n";

foreach($test_names as $name)
{
	$file = file("$path/source/TEST/$name"."_test.C");
	$check = 0;
	$test = 0;
	foreach($file as $line)
	{
		$line = trim($line);
		
		if ($line[0].$line[1]!="//")
		{
			if (substr($line,0,5)=="CHECK")
			{
				$check++;
			}
			if (substr($line,0,15)=="TEST_REAL_EQUAL" || substr($line,0,10)=="TEST_EQUAL" || substr($line,0,14)=="TEST_NOT_EQUAL" || substr($line,0,14)=="TEST_EXCEPTION" || substr($line,0,9)=="TEST_FILE" )
			{
				$test++;
			}
		}
	}
	if ($check<3 || $test < $check)
	print "$name $check $test\n";
}

?>
