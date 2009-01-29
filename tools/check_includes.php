<?
# -*- mode: C++; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
Usage: checker.php <Path to OpenMS> <check option>

check options:
  -u  check for unneeded includes
  -c  check for unneeded includes in C-files (results are not always reliable because of forward declarations)
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
	
	if ($argv[2]!="-u" AND $argv[2]!="-c")
	{
		print "No check option given -> ABORTING!\n\n";
		exit;
	}
	
	print "\n\n";
	
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
									"Identification" => array("Identification", "IdentificationData"),
									
	 							);

	###unneeded includes for header files###
	if ($argv[2]=="-u")
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
		exec("find $path/source/ -name \"*.C\" ! -name \"*_test.C\" ! -regex \"*/EXAMPLES/*\" ! -regex \"*/TEST/*\"", $files);

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
		exec("find $path/source/ -name \"*.C\" ! -name \"*_test.C\" ! -regex \"*/EXAMPLES/*\" ! -regex \"*/TEST/*\"", $files);

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
