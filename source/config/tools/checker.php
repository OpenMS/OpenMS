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
	
	#TODO: add check for self assignment
	
	######################## helper functions ###############################
	function printUsage()
	{
		print "Usage: checker.php <Absolut path to OpenMS> [-u \"user name\"] [-t test] [options]\n";
		print "\n";
		print "If no user name is given, the tests are performed for all users.\n";
		print "\n";
		print "tests:\n";
		foreach ($GLOBALS["tests"] as $name => $desc)
		{
				print "  $name -- $desc\n";
		}
		print "\n";
		print "options:\n";
		print "  -v         verbose mode\n";
		print "  -d <level> debug mode\n";
		print "  --help     shows this help\n";
	}

	function realOutput($text,$user,$verbose,$filename)
	{
		if ($verbose)
		{
			print "------> ";
		}
		print $text;
		if ($user=="all" && $filename!="")
		{
			print " (".$GLOBALS["file_maintainers"][$filename].")";
		}
		print "\n";
		#error count
		if ($filename!="")
		{
			foreach (explode(", ",$GLOBALS["file_maintainers"][$filename]) as $u)
			{
				$GLOBALS["maintainer_info"][trim($u)]["errors"]++;
			}
		}
	}
	
	######################## declarations ###############################
	$GLOBALS["tests"] = array(
															"all"				     => "performs all tests (default)",
															"guards"         => "check if header guards present and correct",
															"tab"            => "check tab settings for editors",
															"maintainers"    => "check if maintainers are consistent in header, source and test file",
															"missing_tests"  => "check for missing tests",
															"brief"          => "check for doxygen tag @brief",
															"doxygen_errors" => "check for errors in dogygen-error.log",
															"check_test"     => "check the class test for completeness",
															"test_output"    => "check test output for warnings and errors",
															"topp_output"    => "check TOPP test output for warnings and errors",
															"svn_keywords"   => "check if the SVN keywords are set for tests",
															"coding"				 => "check if coding convention is followed"
														);
	
	$options = array("-u","-t","-d");
	
	$flags = array("-v");
	
	######################## parameter handling ###############################
	
	#no parameters
	if ($argc==1 || in_array("--help",$argv))
	{
		printUsage();
		exit;
	}
	
	#wrong parameter count
	for ($i=1; $i<count($argv); ++$i)
	{
		#option
		if (beginsWith($argv[$i],"-"))
		{
			#no path given
			if ($i==1)
			{
				print "\nError: No path given!\n\n";
				printUsage();
				exit;
			}
			#not registered option or flag
			if (!in_array($argv[$i],array_merge($options,$flags)))
			{
				print "\nError: Unregistered option '".$argv[$i]."'!\n\n";
				printUsage();
				exit;
			}
			# no argument to an option
			if (in_array($argv[$i],$options) && ( !isset($argv[$i+1]) || beginsWith($argv[$i+1],"-")))
			{
				print "\nError: No argument to option '".$argv[$i]."'!\n\n";
				printUsage();
				exit;
			}
		}
	}
	
	#remove slash from path
	if (endsWith($argv[1],"/"))
	{
		$path = substr($argv[1],0,-1);
	}
	else
	{
		$path = $argv[1];
	}
	
	# verbose
	$verbose = false;
	if (in_array("-v",$argv))
	{
		$verbose = true;
	}

	#debug
	$debug = 0; 
	if (in_array("-d",$argv))
	{
		$debug = $argv[array_search("-d",$argv)+1];
	}
	
	#user
	$user = "all"; 
	if (in_array("-u",$argv))
	{
		$user = $argv[array_search("-u",$argv)+1];
	}
	
	#test
	$test = "all";
	if (in_array("-t",$argv))
	{
		$test = $argv[array_search("-t",$argv)+1];
		if (!in_array($test, array_keys($GLOBALS["tests"])))
		{
			print "\nError: Unknown test '$test'!\n\n";
			printUsage();
			exit;
		}
	}
	
	if ($verbose)
	{
		print "Path : '$path'\n";
		print "User : '$user'\n";
		print "Test : '$test'\n";
		print "Debug: '$debug'\n";
	}

	########################### NEEDED FILES ###############################
	
	$abort = false;
	$out = array();
	exec("cd $path/doc/xml/ && ls -a *.xml | wc -l",$out);
	if (trim($out[0]) < 100)
	{
		print "Error: For this script, doxygen XML output is needed!\n";
		print "       Please execute 'make idoc' in '$path/doc/'.\n";
		$abort = true;
	}
	if ($test == "all" || $test == "doxygen_errors")
	{
		if (!file_exists("$path/doc/doxygen-error.log"))
		{
			print "Error: For the 'doxygen_errors' test, the file '$path/doc/doxygen-error.log' is needed!\n";
			print "       Please execute 'make idoc' in '$path/doc/'.\n";
			$abort = true;
		}
	}
	if ($test == "all" || $test == "test_output")
	{
		$out = array();
		exec("cd $path/source/TEST/ && ls -a *.output 2> /dev/null | wc -l",$out);
		if (trim($out[0]) < 100)
		{
			print "Error: For the 'test_output' test, test output files are needed!\n";
			print "       Please execute 'make test' in '$path/source/'.\n";
			$abort = true;
		}
	}
	if ($test == "all" || $test == "topp_output")
	{
		$out = array();
		exec("cd $path/source/TEST/TOPP/ && ls -a *.output 2> /dev/null | wc -l",$out);
		if (trim($out[0]) < 10)
		{
			print "Error: For the 'topp_output' test, TOPP test output files are needed!\n";
			print "       Please execute 'make TOPPtest' in '$path/source/'.\n";
			$abort = true;
		}
	}
	if ($abort)
	{
		exit;
	}
	
	########################### MAINTAINERS ################################
	
	$files=array();
	exec("cd $path && find include/OpenMS/ -name \"*.h\" ! -name \"*Template.h\"", $files);
	exec("cd $path && find source/ -name \"*.C\" ! -regex \".*/EXAMPLES/.*\" ! -regex \".*/tools/.*\" ! -name \"*_moc.C\" ! -name \"moc_*.C\" ! -name \"*Template.C\"", $files);
	
	//look up Maintainer in first 40 lines of files
	
	$GLOBALS["maintainer_info"] = array();
	$files_todo = array();
	$GLOBALS["file_maintainers"] = array();
	
	foreach ($files as $f)
	{
		$maintainerline = array();
		exec("cd $path && head -40 $f | grep Maintainer ", $maintainerline);
		if (count($maintainerline) == 0)
		{
			if ($user=="all")
			{
				print "No Maintainer lines in '$f'\n";
				$GLOBALS["file_maintainers"][$f] = "";
			}
		}
		else if (count($maintainerline) > 1)
		{
			if ($user=="all")
			{
				print "Several maintainer lines in '$f'\n";
				$GLOBALS["file_maintainers"][$f] = "";
			}
		}
		else
		{
			$maintainers = parseMaintainerLine($maintainerline[0]);
			if (count($maintainers) == 0)
			{
				if ($user=="all")
				{
					print "No Maintainer for '$f'\n";
					$GLOBALS["file_maintainers"][$f] = "";
				}
			}
			else
			{
				$GLOBALS["file_maintainers"][$f] = implode(", ",$maintainers);
				
				#count files per maintainer
				foreach ($maintainers as $m)
				{
					if (!isset($GLOBALS["maintainer_info"][$m]))
					{
						$GLOBALS["maintainer_info"][$m]["files"] = 0;
						$GLOBALS["maintainer_info"][$m]["errors"] = 0;
					}
					$GLOBALS["maintainer_info"][$m]["files"]++;
				}
				#check for misspelled maintainers
				if ($user!="all")
				{
					foreach ($maintainers as $m)
					{
						$dist = levenshtein($m,$user);
						if ($dist==0)
						{
							$files_todo[] = $f;
						}
						else if($dist<=4)
						{
							print "Possibly misspelled maintainer '$user'<-$dist->'$m' in '$f'\n";
						}
					}
				}
			}
		}
	}
	//files to parse
	if ($user == "all")
	{
		$files_todo = $files;
	}
	
	########################### auxilary files #############################
	$called_tests = array();
	$makefile = file("$path/source/TEST/Makefile");
	foreach($makefile as $line)
	{
		$line = trim($line);
		if (strpos($line,"_test")!== FALSE && strpos($line,":")===FALSE && $line[0]!="#" && $line[0]!="@")
		{
			$called_tests[] = "source/TEST/".strtr($line,array("\\"=>""," "=>"", "	"=>"")).".C";
		}
	}

	if ($test == "all" || $test == "doxygen_errors")
	{
		$doxygen_errors = array();
		$errorfile = file("$path/doc/doxygen-error.log");
		foreach($errorfile as $line)
		{
			if (ereg("(.*/[a-zA-Z0-9_]+\.[hC]):[0-9]+:",$line,$parts))
			{
				//skip warning where doxygen cannot resolve members
				if (strpos($line,"no uniquely matching class member")===FALSE && strpos($line,"no matching class member")===FALSE)
				{
					$pos = strpos($parts[1],"source/");
					if($pos!==FALSE)
					{
						$doxygen_errors[] = substr($parts[1],$pos);
					}
					$pos = strpos($parts[1],"include/OpenMS/");
					if($pos!==FALSE)
					{
						$doxygen_errors[] = substr($parts[1],$pos);
					}
				}
			}
		}
		$doxygen_errors = array_unique($doxygen_errors);
	}
	########################################################################
	########################### TESTS ######################################
	########################################################################
	foreach ($files_todo as $f)
	{
		if ($debug>0)
		{
			print "File name: '$f'\n";
		}
		
		//file name (without path)
		$basename = basename($f);
		//class name (for source and header files)
		$classname = substr($basename,0,-2);
		//test name (for source and header files) 
		$testname = "source/TEST/".$classname."_test.C";
		
		// file content
		$file = file($path."/".$f);
		
		######################### load class info ################################
		$dont_load = array(
			"RangeUtils.h",
			"ComparatorUtils.h",
			"KernelTraits.h",
			"StandardTypes.h",
			"DimensionDescription.h",
			"MathFunctions.h",
			"ClassTest.h",
			"LayerData.h",
			"config.h",
			"XMLSchemes.h",
			"Serialization.h",
			"Exception.h",
			"Types.h",
			"Macros.h",
			"Benchmark.h",
			"Constants.h",
			);

		if (!endsWith($f,"_registerChildren.h") && endsWith($f,".h") && !in_array($basename,$dont_load))
		{
			$class_info = getClassInfo($path,$f,$debug);
		}
		else
		{
			unset($class_info);
		}

		########################### guards ######################################
		if ($test == "all" || $test == "guards")
		{
			$dont_report = array("TypeNameIdStringMiscellanyDefs.h");
			
			if (endsWith($f,".h"))
			{
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
							$right_guard = includeToGuard(suffix($f,strlen($guard)));
							if ($right_guard!=$guard OR !beginsWith($guard,"OPENMS_"))
							{
								realOutput("Wrong header guard '$guard' in '$f'",$user,$verbose,$f);
							}					
							break;
						}
					}
					
					$class = trim(substr($f,strrpos($f,"/")+1));
					if ($i==count($file)-1 AND !in_array($class,$dont_report))
					{
						realOutput("Missing header guard in '$f' ",$user,$verbose,$f);
					}
				}	
			}
		}
	
		########################### tab settings #####################################
		if ($test == "all" || $test == "tab")
		{
			$tab_count = 0;
			for ($i=0;$i<min(30,count($file));$i++)
			{
				$line = trim($file[$i]);
				if (strpos($line,"vi: set ts=2:")!==false )
				{
					$tab_count++;
				}
			}
			if ($tab_count!=1)
			{
				realOutput("Missing tab settings in '$f'",$user,$verbose,$f);
			}
		}

		########################### maintainers  #####################################
		if ($test == "all" || $test == "maintainers")
		{
			if (endsWith($f,".h") )
			{
				# maintainer of test file
				if (in_array($testname,$files))
				{
					if ($file_maintainers[$testname] != $file_maintainers[$f])
					{
						realOutput("Inconsistent maintainers in '$f' and '$testname'",$user,$verbose,$f);
						if ($verbose)
						{
							print "  '$file_maintainers[$testname]'<->'$file_maintainers[$f]'\n";
						}
					}
				}
				# maintainer of source file
				$source_name = "source/".substr($f,15,-2).".C";
				if (in_array($source_name,$files))
				{
					if ($file_maintainers[$source_name] != $file_maintainers[$f])
					{
						realOutput("Inconsistent maintainers in '$f' and '$source_name'",$user,$verbose,$f);
						if ($verbose)
						{
							print "  '$file_maintainers[$source_name]'<->'$file_maintainers[$f]'\n";
						}
					}
				}
			}
		}

		########################### missing tests  #####################################
		if ($test == "all" || $test == "missing_tests")
		{
			$dont_report = array(
				"/FORMAT/HANDLERS/",
				"/VISUAL/",
				"/CONCEPT/Types.h",
				"/CONCEPT/Macros.h",
				"/CONCEPT/Exception.h",
				"/CONCEPT/Benchmark.h",
				"/CONCEPT/Constants.h",
				"/config.h",
				"include/OpenMS/APPLICATIONS/TOPPViewBase.h",
				"_registerChildren.h",
				"DataReducer.h",
				"SchemaFile.h",
				"Serialization.h",
				);

			if (endsWith($f,".h") )
			{
				$ignore = false;
				foreach ($dont_report as $i)
				{
					if (strpos($f,$i)!==FALSE)
					{
						$ignore = true;
					}		
				}
				
				if (!$ignore)
				{
					if (!in_array($testname,$files))
					{
						realOutput("Missing test for '$f'",$user,$verbose,$f);
					}
					else if (!in_array($testname,$called_tests))
					{
						realOutput("Test not in Makefile for '$f'",$user,$verbose,$f);
					}
				}
			}
		}
		
		########################### @brief  #####################################
		if ($test == "all" || $test == "brief")
		{
			if (endsWith($f,".h") )
			{
				for ($i=0; $i<count($file); ++$i)
				{
					#match a class declaration
					if (preg_match("/(class|struct)[\s]+[\w]+[\s]*({|:[^:]|$)/i",$file[$i]))
					{
						#take the second line too in case : or { come in the next line
						if (isset($file[$i+1]) && preg_match("/(class|struct)[\s]+([\w]+)[\s]*({|:[^:])/i",$file[$i].$file[$i+1],$parts))
						{
							$brief = false;
							#search for /// comment
							if (strpos($file[$i-1],"///")!==FALSE)
							{
								$brief = true;
							}
							# backward search for @brief until comment closes 
							for ($j=$i; $j!=0; --$j)
							{
								if (strpos($file[$j],"@brief")!==FALSE)
								{
									$brief = true;
									break;
								}
								if (strpos($file[$j],"/**")!==FALSE)
								{
									break;
								}
							}
							if (!$brief)
							{
								realOutput("No @brief description for '$parts[2]' in '$f'",$user,$verbose,$f);
							}
						}
					}
				}
			}
		}
		
		########################### doxygen errors  ####################################
		if ($test == "all" || $test == "doxygen_errors")
		{
			if (in_array($f,$doxygen_errors))
			{
				realOutput("Doxygen errors in '$f'",$user,$verbose,$f);
			}
		}
		
		########################### test errors  #####################################
		
		if (isset($class_info) && ($test == "all" || $test == "check_test"))
		{
			if (in_array($testname,$files))
			{
				
 				#parse test
 				parseTestFile("$path/$testname",$tests);
				
				#compare declarations and tests
				$out = compareDeclarationsAndTests($class_info["public-long"],$tests);
				
				#output
				if (count($out["missing"])!=0 || count($out["unknown"])!=0)
				{
					realOutput("Test errors in '$f'",$user,$verbose,$testname);
					if ($verbose)
					{
						if (count($out["unknown"])!=0)
						{
							print "  Tests of unknown methods:\n";
							foreach ($out["unknown"] as $u)
							{
								print "    - '$u'\n";
							}
						}
						if (count($out["missing"])!=0)
						{
							print "  Missing tests:\n";
							foreach ($out["missing"] as $m)
							{
								print "    - '$m'\n";
							}
						}
					}
				}
			}
		}
		
		
		
		############################## coding ##########################################
		if (isset($class_info) && ($test == "all" || $test == "coding"))
		{
			$out = array();
			#variables
			foreach($class_info["variables"] as $tmp)
			{
				if (!endswith($tmp,'_'))
				{
					$out[] = "  - invalid non-public variable name '$tmp'\n";
				}
			}
			# non-public member
			foreach($class_info["non-public"] as $tmp)
			{
				# constructor, destructor, serialize methods and QT events are allowed
				if ( endswith($tmp,'Event') || endsWith($tmp,'load')  || endsWith($tmp,'save') || endsWith($tmp,'serialize') || $tmp==$class_info["classname"] || $tmp=='~'.$class_info["classname"] )
				{
					continue;
				}
				if (!endswith($tmp,'_'))
				{
					$out[] = "  - invalid non-public method name '$tmp'\n";
				}
				else if (strpos(substr($tmp,0,-1),'_')!==FALSE)
				{
					$out[] = "  - invalid non-public method name '$tmp'\n";
				}
			}
			#output
			if (count($out)!=0)
			{
				realOutput("Coding convention violation in '$f'",$user,$verbose,$f);
				if ($verbose)
				{
					foreach ($out as $o)
					{
						print $o;
					}
				}
			}
		}
		
		
		########################### warnings test  #####################################
		if ($test == "all" || $test == "test_output")
		{
			$outputfile = substr("$path/$testname",0,-2).".output";
			if (endsWith($f,".h") && in_array($testname,$files) && file_exists($outputfile))
			{
				$testfile = file($outputfile);
				$errors = array();
				foreach ($testfile as $line)
				{
					if (stripos($line,"warning")!==FALSE || stripos($line,"error")!==FALSE)
					{
						$errors[] = trim($line);
					}
				}
				if (count($errors)!=0)
				{
					realOutput("Error/warnings in test output of '$testname'",$user,$verbose,$testname);
					if ($verbose)
					{
						foreach ($errors as $e)
						{
							print "  '$e'\n";
						}
					}
				}
			}
		}

		######################### 'Id' keyword in tests  ###############################
		if ($test == "all" || $test == "svn_keywords")
		{
			if (endsWith($f,".h") )
			{
				if (in_array($testname,$files))
				{
					$out = array();
					exec("svn proplist -v $path/$testname",$out);
					$kw = false;
					foreach ($out as $line)
					{
						if (strpos($line,"svn:keywords")!==FALSE)
						{
							$kw = true;
							if (strpos($line,"Id")===FALSE)
							{
								realOutput("svn:keyword 'Id' not set for '$testname'",$user,$verbose,$testname);
							}
						}
					}
					if (!$kw)
					{
						realOutput("svn:keyword 'Id' not set for '$testname'",$user,$verbose,$testname);
					}
				}	
			}
		}
	}

	################### doxygen errors in .doxygen-files  ##########################
	if ($user == "all" && ($test == "all" || $test == "doxygen_errors"))
	{
		$file = file("$path/doc/doxygen-error.log");
		foreach ($file as $line)
		{
			$line = trim($line);
			if (ereg("(.*/[a-zA-Z0-9_]+\.doxygen):[0-9]+:",$line,$parts))
			{
				realOutput("Doxygen errors in '".$parts[1]."'",$user,$verbose,"");
			}
		}
	}

	
	########################### warnings TOPPtest  #################################
	if ($test == "all" || $test == "topp_output")
	{
		$out = array();
		exec("cd $path/source/TEST/TOPP/ && ls -a *.output",$out);
		foreach ($out as $file)
		{
			$file = trim($file);
			$topp_file = "source/APPLICATIONS/TOPP/".substr($file,0,-6)."C";
			if (in_array($topp_file,$files_todo))
			{
				$testfile = file("$path/source/TEST/TOPP/$file");
				$errors = array();
				foreach ($testfile as $line)
				{
					if (stripos($line,"warning")!==FALSE || stripos($line,"error")!==FALSE)
					{
						$errors[] = trim($line);
					}
				}
				if (count($errors)!=0)
				{
					realOutput("Error/warnings in TOPP test output 'source/TEST/TOPP/$file'",$user,$verbose,$topp_file);
					if ($verbose)
					{
						foreach ($errors as $e)
						{
							print "  '$e'\n";
						}
					}
				}
			}
		}
	}
	
	########################### maintainer summary #################################
	if ($user == "all")
	{
		print "\nMaintainers:\n";
		foreach ($GLOBALS["maintainer_info"] as $m => $info)
		{
			print "  $m (Files: ".$info["files"]."  Errors: ".$info["errors"].")\n";
		}
	}

?>
