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

include "common_functions.php";

	########################usage###############################
	if ($argc<2 OR $argc>4)
	{
?>
Usage: checker.php <Path to OpenMS> [<check option>] [-v]

If no option is given, all tests are executed. Otherwise only the given specfic test is executed.

check options:
  -g  check if header guards present and correct
  -h  check headers (tab settings for editors, cvs headers, maintainer)
  -m  check headers for empty maintainer field
  -u  check for unneeded includes
  -c  check for unneeded includes in C-files (results are not always reliable because of forward declarations)
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
	if (!isset($argv[2]) OR ($argv[2]!="-g" AND $argv[2]!="-h" AND $argv[2]!="-u" AND $argv[2]!="-m" AND $argv[2]!="-c"))
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
			if ($verbose)
			{ 
				print "  Maintainer lines: $maintainer_count\n";	
				print "  Editor tab settings lines: $tab_count\n";	
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
	//hide
	$dont_report = array("config","TypeNameIdStringMiscellanyDefs","Constants","Benchmark","helper");
	//multiple classes per file (filename => classes)
	$multiclass = array(
									"Macros" => array("OPENMS_PRECONDITION","OPENMS_POSTCONDITION"),
									"ComparatorUtils" => array("PointerComparator","ReverseComparator","LexicographicComparator","pointerComparator","reverseComparator","lexicographicComparator"),
	 								"MathFunctions" => array("ceil_decimal","round_decimal","intervalTransformation","linear2log","log2linear","isOdd"),
	 								"HashFunction" => array("hashPointer","hashPJWString","hashElfString","getNextPrime"),
									"RangeUtils" => array("RTRange","MSLevelRange","ScanModePredicate","SpectrumEmptyPredicate","MzRange","IntensityRange"),
									"StandardTypes" => array("RawDataPoint","RawDataPoint2D","RawSpectrum","RawMap","Peak","Peak2D","PeakSpectrum","PeakMap","Feature","FeatureMap"),
									"Types" => array("Distance","Handle","Index","SignedInt","UnsignedInt","Size","Time","HashIndex","Position","Real","DoubleReal","Property","ErrorCode","Byte","PointerSizeUInt","PointerSizeInt","UID"),
									"Exception" => array("throw","Base"),
									"DPickedPeak" => array("PeakShapeType","DPickedPeak"),
									"TimeStamp" => array("PreciseTime","TimeStamp"),
									"Identification" => array("Identification", "IdentificationData")
									
	 							);

	###unneeded includes for header files###
	if ($do_all OR $argv[2]=="-u")
	{
		####unneeded declarations or forward declartion possible###
		$files=array();
		exec("find $path/include/ -name \"*.h\" ! -name \"*Template.h\"", $files);

		foreach ($files as $filename)
		{
			$count = array();
			$forward = array();
			$file = file($filename);
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
							$forward["$class"] = 0;
						}
					}
				}
				else
				{
					$tokens = tokenize($line."\\");
					foreach ($count as $class => $number)
					{
						$classes = array();
						//multiple classes per file
						if (isset($multiclass[$class]))
						{
							$classes = $multiclass[$class];
						}
						// class name == file name
						else
						{
							$classes[0] = $class;
						}
						
						// scan tokens for class name
						for ($j = 0; $j != sizeof($tokens); $j++)
						{
							if (in_array($tokens[$j], $classes))
							{
								$count[$class]++;
								// check next token: if it is '*' or '&', the include might be unneeded
								if (isset($tokens[$j + 1]) && ( $tokens[$j + 1] == "*" || $tokens[$j + 1] == "&"))
								{
									$forward[$class]++;
								}
							}
						}
					}
				}
			}
			
			$out = ">>UNNEEDED INCLUDE in $filename\n";
			foreach ($count as $class => $number)
			{
				if ($number == 0)
				{
					$out .= "  $class\n";
				}
				else if (isset($forward[$class]) && $forward[$class] == $number)
				{
					$out .= "  $class (forward)\n";
				}
			}
			if ($out != ">>UNNEEDED INCLUDE in $filename\n")
			{
				print $out;
			}
		}
		
		###check for duplicate includes###
		$files=array();
		exec("find $path/source/ -name \"*.C\" ! -name \"*_test.C\" ! -wholename \"*/EXAMPLES/*\" ! -wholename \"*/TEST/*\"", $files);

		foreach ($files as $filename)
		{
			$includes = array();
			$duplicates = array();
			$file = file($filename);
			for ($i=0;$i<count($file);$i++)
			{
				$line = trim($file[$i]);
				if (isIncludeLine($line,$include))
				{
					if (beginsWith($include,"OpenMS/"))
					{
						$includes[] = substr($include,strrpos($include,"/")+1,-2);
					}
				}
			}
			//header for source file
			$header = substr($filename,0,strpos($filename,"/source/"))."/include/OpenMS/".substr($filename,strpos($filename,"/source/")+8,-2).".h";
			if (file_exists($header))
			{
				$headerfile = file($header);
			}
			else
			{
				$headerfile = array();
			}
			for ($i=0;$i<count($headerfile);$i++)
			{
				$line = trim($headerfile[$i]);
				if (isIncludeLine($line,$include))
				{
					if (beginsWith($include,"OpenMS/"))
					{
						$class = substr($include,strrpos($include,"/")+1,-2);
						if (in_array($class,$includes))
						{
							$duplicates[] = $class;
						}
					}
				}
			}
			if (count($duplicates)!=0)
			{
				print ">>UNNEEDED INCLUDE in $filename\n";
				foreach ($duplicates as $d)
				{
					print "  $d\n";
				}
			}
		}
	}

	###unneeded includes for source files###
	if ($argv[2]=="-c")
	{
		$files=array();
		exec("find $path/source/ -name \"*.C\" ! -name \"*_test.C\" ! -wholename \"*/EXAMPLES/*\" ! -wholename \"*/TEST/*\"", $files);

		foreach ($files as $filename)
		{
			$count = array();
			$file = file($filename);
			for ($i=0;$i<count($file);$i++)
			{
				$line = trim($file[$i]);
				if (isIncludeLine($line,$include))
				{
					if (beginsWith($include,"OpenMS/"))
					{
						$class = substr($include,strrpos($include,"/")+1,-2);
						//Do not report header of source file
						if (!in_array($class,$dont_report) && $class != substr($filename,strrpos($filename,"/")+1,-2))
						{
							$count["$class"] = 0;
						}
					}
				}
				else
				{
					foreach ($count as $class => $number)
					{
						//multiple classes per file
						if (isset($multiclass[$class]))
						{
							foreach ($multiclass[$class] as $subclass)
							{
								if (strpos($line,$subclass)!== FALSE)
								{
									$count[$class]++;
								}
							}
						}
						//class name == file
						else if (strpos($line,$class)!== FALSE)
						{
							$count[$class]++;
						}
					}
				}
			}

			// find header for source files
		  $header = substr($filename,0,strpos($filename,"/source/"))."/include/OpenMS/".substr($filename,strpos($filename,"/source/")+8,-2).".h";
			if (file_exists($header))
			{
				$headerfile = file($header);
			}
			else
			{
				$headerfile = array();
			}
			
			$out = ">>UNNEEDED INCLUDE in $filename\n";
			foreach ($count as $class => $number)
			{
				if ($number == 0)
				{
					//check for forward declaration
					$forward_declared = false;
					foreach ($headerfile as $line)
					{
						if (trim($line)=="class $class;")
						{
							$forward_declared = true;
							break;
						}
					}
					if (!$forward_declared) $out .= "  $class ?\n";
				}
			}
			if ($out != ">>UNNEEDED INCLUDE in $filename\n")
			{
				print $out;
			}
		}		
	}

?>
