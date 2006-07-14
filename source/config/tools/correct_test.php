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

if ( $argc<3 || $argc>4 || ( $argc==4 &&  $argv[$argc-1]!="-verbose") )
{
	print "\n\nUsage: correct_test.php <Absolut path to OpenMS> <Absolut path to header> [-verbose]\n\n";
	exit;	
}

$verbose = false;
if ($argc==4 && $argv[$argc-1]=="-verbose")
{
	$verbose = true;
}

function penalize_name($f1, $f2)
{
	$p = 100;
	eregi("[ 	]+([^ 	]+)[ 	]?\(", " ".$f1, $r1);
	eregi("[ 	]+([^ 	]+)[ 	]?\(", " ".$f2, $r2);
	
	if ($r1[1]==$r2[1])
	{
		$p = 0;
	}
	//print_r($r1);
	//print_r($r2);
	//print "\nF1: $f1\nR1: $r1[1]\nF2: $f2\nR2: $r2[1]\nScore: $p\n";
	return $p;
}

//determine methods
$check_test = $argv[1]."/source/config/tools/check_test";
exec("$check_test $argv[2] -a",$out, $return);
if (!file_exists($check_test))
{
	print "Tool $check_test not present => Aborting!";
	exit(1);
}

for ($i=1; $i< count($out); ++$i)
{
	$methods[] = substr(trim($out[$i]),1,-1);
}

// print methods in verbose mode
if ($verbose)
{
	print "\n\nDecalared methods:\n";
	foreach ($methods as $m) print "  $m\n";
}

//parse test
$test_name = $argv[1]."/source/TEST/".substr($argv[2], strrpos($argv[2],"/")+1,-2)."_test.C";

if (!file_exists($test_name))
{
	print "Test $test_name not present => Aborting!";
	exit(1);
}

$test = file($test_name);
foreach($test as $line)
{
	$line = trim($line);
	if (substr($line,0,5)=="CHECK")
	{
		$tests[] = trim(substr($line,strpos($line,"(")+1,-1));
	}
}

//calculate diff
$replace_whitespaces = array("\t"=>""," "=>"");

for($i=0; $i<count($tests); ++$i)
{
	for($j=0; $j<count($methods); ++$j)
	{
		$dists[$i][$j] = penalize_name($tests[$i],$methods[$j]) + levenshtein ( strtr($tests[$i],$replace_whitespaces), strtr($methods[$j],$replace_whitespaces), 1, 10, 10 );
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
		$replace[] = $methods[key($array)];
		continue;
	}
	
	if (strtoupper(substr($tests[$i],0,7))=="[EXTRA]")
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
	print "\n[enter] = 0\n";
	print   "[x]     = do not change\n";
	
	//read in choise
	do
	{
		$line = trim(fgets($fp));
	}
	while($line!="" AND !ereg("^[0-9]$",$line) AND $line!="x");
	
	if ($line == "x")
	{
		$replace[] = $tests[$i];
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
foreach($test as $line)
{
	if (substr(trim($line),0,5)=="CHECK")
	{
		fwrite($fp,"CHECK(".$replace[$i].")\n");
		++$i;
	}
	else
	{
		fwrite($fp,"$line");
	}
}

fclose($fp);

?>
