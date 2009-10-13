<?php
# -*- mode: C++; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
# $Maintainer: $
# $Authors: Marc Sturm $
# --------------------------------------------------------------------------

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

include "common_functions.php";

if ( $argc<3 || $argc>4 )
{
	print "\n\nUsage: correct_test.php <Absolute path to OpenMS> <Absolute path to header> [-v]\n\n";
	exit;
}

########################## auxilary functions ##################################

/// penalizes differing method names with 100
function penalize_name($f1, $f2)
{
	# extract name (between whitepace and bracket
	$tmp = trim(substr($f1,0,strpos($f1,"(")));
	$tmp = strtr($tmp,array("operator "=>"operator","operator\t"=>"operator"));
	$n1 = trim(substr($tmp,max(strrpos($tmp," "), strrpos($tmp,"\t"))));

	# extract name (between whitepace and bracket
	$tmp = trim(substr($f2,0,strpos($f2,"(")));
	$tmp = strtr($tmp,array("operator "=>"operator","operator\t"=>"operator"));
	$n2 = trim(substr($tmp,max(strrpos($tmp," "), strrpos($tmp,"\t"))));

	if ($n1==$n2)
	{
		return 0;
	}
	return 100;
}

######################## parameter handling ####################################
$verbose = false;
if (in_array("-v",$argv))
{
	$verbose = true;
}

$path = $argv[1];
$header = $argv[2];
$basename = basename($header);
$test_name = "$path/source/TEST/".substr($basename,0,-2)."_test.C";

######################## determine tested methods ##############################
$tmp = parseTestFile($test_name);
$tests = $tmp["tests"];

if ($verbose)
{
	print "\n\nTests:\n";
	foreach ($tests as $t)
	{
		print "  '$t'\n";
	}
}

######################## determine declared methods ##########################
$class_info = getClassInfo($path,$header,0);
$methods = $class_info["public-long"];

// print methods in verbose mode
if ($verbose)
{
	print "\n\Methods to correct:\n";
	foreach ($methods as $m)
	{
		print "  '$m'\n";
	}
}


######################deter#################

if (count($tests)==0 || count($methods)==0)
{
	print "Nothing to do (no tests or methods)\n";
	exit;
}

$replace_strings = array(
	"\t" => "",
	" " => "",
	"std::" => "",
	"OpenMS::" => "",
	);

$dists = array();
for($i=0; $i<count($tests); ++$i)
{
	for($j=0; $j<count($methods); ++$j)
	{
	  //print "\nTest: ".strtr($tests[$i],$replace_strings)."\nMethod: ".strtr($methods[$j],$replace_strings)."\n";
		$penalty = levenshtein ( strtr($tests[$i],$replace_strings), strtr($methods[$j],$replace_strings), 1, 10, 10 );
//		if ($penalty!=0)
//		{
			$penalty += penalize_name($tests[$i],$methods[$j]);
//		}
		$dists[$i][$j] = $penalty;
		//print "-- '$tests[$i]'\n";
		//print "-- '$methods[$j]'\n";
		//print "-------> '".$dists[$i][$j]."'\n";
		//print "done\n";
	}
}

$fp=fopen("php://stdin","r");

for($i=0; $i<count($tests); ++$i)
{
	$array = $dists[$i];
	asort($array);
	// abort if exact match
	if (current($array)==0)
	{
		$replace[] = $tests[$i];
		continue;
	}

	print "\n\nTest:     ".$tests[$i]."\n\n";
	$j=0;
	foreach ($array as $index => $score)
	{
		print "$j) ".str_pad($score, 4, " ", STR_PAD_LEFT)." - ".$methods[$index]."\n";

		//abort after 10
		++$j;
		if ($j==10)
		{
			break;
		}
	}
	print "\n[enter]  => 0\n";
	print   "[i]      => ignore this test\n";
	print   "[x]      => make [EXTRA] test (is ignored by checker.php)\n";
	print   "[CTRL+C] => abort\n";
	@ob_flush();
	flush();

	//read in choise
	do
	{
		$line = trim(fgets($fp));
	}
	while($line!="" AND !ereg("^[0-9]$",$line) AND $line!="i" AND $line!="x");

	if ($line == "i")
	{
		$replace[] = $tests[$i];
	}
	else if ($line == "x")
	{
		$replace[] = "[EXTRA]".$tests[$i];
	}
	else
	{
		if ($line == "")
		{
			$line = 0;
		}
		$tmp = array_keys($array);
		$replace[] = $methods[$tmp[$line]];
	}
}

fclose($fp);

//backup original test
exec("mv $test_name $test_name.bak");

//write test
$fp=fopen($test_name,"w");
$i=0;
foreach(file($test_name.".bak") as $line)
{
	if (substr(trim($line),0,13)=="START_SECTION" && strpos($line,"[EXTRA]")===FALSE)
	{
		fwrite($fp,"START_SECTION((".$replace[$i]."))\n");
		++$i;
	}
	else
	{
		fwrite($fp,"$line");
	}
}

fclose($fp);

?>
