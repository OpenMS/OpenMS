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
#
# --------------------------------------------------------------------------
# $Id: checker.php,v 1.4 2006/03/28 12:53:15 marc_sturm Exp $
# $Author: marc_sturm $
# $Maintainer: Marc Sturm $
# --------------------------------------------------------------------------

include "common_functions.php";

	########################usage###############################
	if ($argc<2 OR $argc>4)
	{
?>
Usage: checker.php <Path to OpenMS> [<check option>] [-v]

If no option is given, all tests are executed. Otherwise only the given specfic test is executed.

check options:
  -g  check if header guards present and correct
  -i  check if include guards correct (if present)
  -h  check headers (tab settings for editors, cvs headers, maintainer)
  -m  check headers for empty maintainer field
  -u  check for unneeded includes (does not check references or pointers yet)
other options:
  -v  verbose mode
<?
	exit;	
	}
	
	######################command line parsing###########################
	if (endsWith($argv[1],"/"))
	{
		$path = substr($argv[1],0,-1);
	}
	else
	{
		$path = $argv[1];
	}
	
	$do_all = false;
	if (!isset($argv[2]) OR ($argv[2]!="-g" AND $argv[2]!="-i" AND $argv[2]!="-h" AND $argv[2]!="-u" AND $argv[2]!="-m"))
	{
		$do_all = true;
	}
	
	$verbose = false;
	if (in_array("-v",$argv))
	{
		$verbose = true;
	}
	
	print "\n\n";
	#################header guards###############################
	if ($do_all OR $argv[2]=="-g")
	{
		
		$dont_report = array("TypeNameIdStringMiscellanyDefs.h");
		
		$files=array();
		exec("find $path/include/ -name \"*.h\" ! -name \"*Template.h\"", $files);
		foreach ($files as $header)
		{
			if ($verbose) print "##file: $header\n";
			$file = file($header);
			for ($i=0;$i<count($file);$i++)
			{
				$line = trim($file[$i]);
				if (beginsWith($line,"#ifndef"))
				{
					$guard = trim(substr($line,8));
					$nextline = trim($file[$i+1]);		
					//header guards
					if (beginsWith($nextline,"#define") AND trim(substr($nextline,8))==$guard)
					{
						$right_guard = includeToGuard(suffix($header,strlen($guard)));
						if ($right_guard!=$guard OR !beginsWith($guard,"OPENMS_"))
						{
							print ">>WRONG HEADER GUARD\n";
							print "  guard   : $guard\n";
							print "  filename: $header\n";
						}					
						break;
					}
				}
				
				$class = trim(substr($header,strrpos($header,"/")+1));
				if ($i==count($file)-1 AND !in_array($class,$dont_report))
				{
							print ">>MISSING HEADER GUARD\n";
							print "  filename: $header\n";
				}
			}	
		}
	}

	#################include guards###############################
	if ($do_all OR $argv[2]=="-i")
	{
		$files=array();
		exec("find $path/include/ -name \"*.h\"", $files);
		foreach ($files as $header)
		{
			$file = file($header);
			$out = ">>WRONG INCLUDE GUARD in $header\n";
			$ok = true;
			for ($i=0;$i<count($file);$i++)
			{
				$line = trim($file[$i]);
				if (beginsWith($line,"#ifndef"))
				{
					$guard = trim(substr($line,8));
					$nextline = trim($file[$i+1]);
					$nextnextline = trim($file[$i+2]);
					//header guards
					if (isIncludeLine($nextline,$include) AND beginsWith($nextnextline,"#endif"))
					{
						$right_guard = includeToGuard($include);
						if ($guard!=$right_guard)
						{
							$out .= "  line    : $i\n";
							$out .= "  guard   : $guard\n";
							$out .= "  inlcude : $include\n";
							$ok = false;
						}
					}
				}
			}
			if (!$ok)
			{
				print $out;
			}	
		}
	}

	#################headers###############################
	if ($do_all OR $argv[2]=="-h")
	{
		$files=array();
		exec("find $path/include/ -name \"*.h\" ! -name \"*Template.h\"", $files);
		exec("find $path/source/ -name \"*.C\" ! -name \"*_moc.C\" ! -name \"moc_*.C\" ! -name \"*Template.C\"", $files);
		exec("find $path/source/ -name \"Makefile\"", $files);
		foreach ($files as $infile)
		{
			if ($verbose) print "##file: $infile\n";
			$file = file($infile);
			$out = ">>HEADER ERRORS in $infile\n";
			$ok = true;
			$maintainer_count = 0;
			$tab_count = 0;
			$cvs_id = 0;
			$cvs_author = 0;
			for ($i=0;$i<min(30,count($file));$i++)
			{
				$line = trim($file[$i]);
				if (strpos($line,"\$Maintainer")!==false )
				{
					$maintainer_count++;
				}
				if (strpos($line,"vi: set ts=2:")!==false )
				{
					$tab_count++;
				}
				if (strpos($line,"\$Id")!==false )
				{
					$cvs_id++;
				}
				if (strpos($line,"\$Author")!==false )
				{
					$cvs_author++;
				}
			}
			if ($maintainer_count!=1)
			{
				$out .="  Maintainer lines: $maintainer_count\n";
				$ok=false;
			}
			if ($tab_count!=1)
			{
				$out .="  Editor tab settings lines: $tab_count\n";
				$ok=false;
			}
			if ($cvs_id!=1)
			{
				$out .="  CVS id lines: $cvs_id\n";
				$ok=false;
			}
			if ($cvs_author!=1)
			{
				$out .="  CVS author lines: $cvs_author\n";
				$ok=false;
			}
			if ($verbose)
			{ 
				print "  Maintainer lines: $maintainer_count\n";	
				print "  Editor tab settings lines: $tab_count\n";	
				print "  CVS id: $cvs_id\n";	
				print "  CVS author: $cvs_author\n";
			}		
			if (!$ok)
			{
				print $out;
			}	
		}
	}

	#################empty maintainer lines###############################
	if ($do_all OR $argv[2]=="-m")
	{
		$files=array();
		exec("find $path/include/ -name \"*.h\"", $files);
		exec("find $path/source/ -name \"*.C\" ! -name \"*_moc.C\" ! -name \"moc_*.C\"", $files);
		foreach ($files as $infile)
		{
			if ($verbose) print "##file: $infile\n";
			$file = file($infile);
			
			for ($i=0;$i<min(30,count($file));$i++)
			{
				if (strpos($file[$i],"\$Maintainer")!==false )
				{
					if (ereg("Maintainer[:]?[ 	\$]*$",trim($file[$i])))
					{
						print "Empty maintainer in $infile\n";
					}
				}
			}
		}
	}

	#################unneeded includes###############################
	if ($do_all OR $argv[2]=="-u")
	{
		$dont_report = array("config","Types","helper","TypeNameIdStringMiscellanyDefs","Macros","Exception","Base64","ComparatorUtils", "GaussFitUtils","PreprocessorFlags","MathFunctions","TimeStamp");
		//todo multi-class headers: TimeStamp
		
		$files=array();
		exec("find $path/include/ -name \"*.h\"", $files);
		foreach ($files as $header)
		{
			$count = array();
			$file = file($header);
			for ($i=0;$i<count($file);$i++)
			{
				$line = trim($file[$i]);
				if (isIncludeLine($line,$include))
				{
					if (beginsWith($include,"OpenMS/"))
					{
						$class = substr($include,strrpos($include,"/")+1,-2);
						if (!in_array($class,$dont_report))
						{
							$count["$class"] = 0;
						}
					}
				}
				else
				{
					foreach ($count as $class => $number)
					{
						if (strpos($line,$class)!== FALSE)
						{
							$count[$class]++;
						}
					}
				}
			}
			
			$out = ">>UNNEEDED INCLUDE in $header\n";
			foreach ($count as $class => $number)
			{
				if ($number == 0)
				{
					$out .= "  $class\n";
				}
			}
			if ($out != ">>UNNEEDED INCLUDE in $header\n")
			{
				print $out;
			}
		}		
	}

?>
