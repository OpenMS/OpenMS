<?php
# -*- mode: php; tab-width: 2; -*-
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
# $Maintainer:$
# $Authors: Marc Sturm $
# --------------------------------------------------------------------------
	
	error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);
	ini_set('memory_limit', '1000M');
	
	include "common_functions.php";
	
	######################## helper functions ###############################
	function printUsage()
	{
		print "Usage: checker.php <OpenMS src path> <OpenMS build path> [-u \"user name\"] [-t test] [options]\n";
		print "\n";
		print "This script works only if an OpenMS copy is used, where\n";
		print "- the internal documentation was built (doc_internal),\n";
		print "- the XML documentation was built (doc_xml),\n";
		print "- all tests were executed.\n";
		print "\n";
		print "The user name given with the '-u' option is matched against the maintainers in the\n";
		print "file headers. If no user name is given, the tests are performed for all users.\n";
		print "\n";
		print "tests:\n";
		foreach ($GLOBALS["all_tests"] as $name => $desc)
		{
				print "  $name -- $desc\n";
		}
		print "\n";
		print "options:\n";
		print "  -x         do not rebuild doxygen xml output\n";
		print "  -s         comma-sparated list of tests to skip\n";
		print "  -d <level> debug mode\n";
		print "  --help     shows this help\n";
	}

	function realOutput($text,$user,$filename)
	{
		print "------> ";
		print $text;
		if ($user=="all" && $filename!="")
		{
			print " (".$GLOBALS["file_maintainers"][$filename].")";
		}
		print "\n";
		#error count
		if ($filename!="")
		{
			if (isset($GLOBALS["file_maintainers"][$filename]))
			{
				foreach (explode(", ",$GLOBALS["file_maintainers"][$filename]) as $u)
				{
					$GLOBALS["maintainer_info"][trim($u)]["errors"]++;
				}
			}
		}
	}
	
	######################## declarations ###############################
	$GLOBALS["all_tests"] = array(
															"(none)"		     => "performs all tests",
															"guards"         => "check if header guards present and correct",
															"tab"            => "check tab settings for editors",
															"maintainers"    => "check if maintainers are consistent in header, source and test file",
															"missing_tests"  => "check for missing tests",
															"old_files"			 => "check for unneeded .C files",
															"brief"          => "check for doxygen tag @brief",
															"doxygen_errors" => "check for errors in dogygen-error.log",
															"check_test"     => "check the class test for completeness",
															"test_output"    => "check test output for warnings and errors",
															"topp_output"    => "check TOPP test output for warnings and errors",
															"svn_keywords"   => "check if the SVN keywords are set for tests",
															"coding"				 => "check if coding convention is followed",
															"defaults"		   => "check if DefautParamHandler classes all have a linked parameters section"
														);
	
	$options = array("-u","-t","-d","-s");
	
	$flags = array("-x");
	
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
			if ($i<3)
			{
				print "\nError: No paths given!\n\n";
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
		$src_path = substr($argv[1],0,-1);
	}
	else
	{
		$src_path = $argv[1];
	}

	#remove slash from path
	if (endsWith($argv[2],"/"))
	{
		$bin_path = substr($argv[2],0,-1);
	}
	else
	{
		$bin_path = $argv[2];
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
	$tests = array_keys($GLOBALS["all_tests"]);
	unset($tests[0]);
	if (in_array("-t",$argv))
	{
		$test = $argv[array_search("-t",$argv)+1];
		if (!in_array($test, $tests))
		{
			print "\nError: Unknown test '$test'!\n\n";
			printUsage();
			exit;
		}
		$tests = array($test);
	}
	
	#skip
	if (in_array("-s",$argv))
	{
		$skip = explode(',',$argv[array_search("-s",$argv)+1]);
		//trim all
		array_map("trim", $skip);
		foreach ($skip as $test)
		{
			if (!in_array($test, $tests))
			{
				print "\nError: Unknown test '$test'!\n\n";
				printUsage();
				exit;
			}
		}
		$tests = array_diff($tests,$skip);
	}
	$tests = array_values($tests);

	#doxygen XML output
	$rebuilt_xml = true;
	if (in_array("-x",$argv))
	{
		$rebuilt_xml = false;
	}
	
	######################## doxygen XML output ############################
	if ($rebuilt_xml)
	{
		if ($debug>0)
		{
			print "Rebuilding doxygen XML output\n";
		}
		exec("cd $bin_path && make doc_xml");
		if ($debug>0)
		{
			print "Done\n";
		}
	}

	########################### NEEDED FILES ###############################
	
	$abort = false;
	$out = array();
	exec("cd $bin_path/doc/xml_output/ && ls -a *.xml | wc -l",$out);
	if (trim($out[0]) < 100)
	{
		print "Error: For this script, doxygen XML output is needed!\n";
		print "       Please execute 'make doc_xml' first!.\n";
		$abort = true;
	}
	if (in_array("doxygen_errors",$tests))
	{
		if (!file_exists("$bin_path/doc/doxygen/doxygen-error.log"))
		{
			print "Error: For the 'doxygen_errors' test, the file '$bin_path/doc/doxygen/doxygen-error.log' is needed!\n";
			print "       Please execute 'make doc_internal' first!'.\n";
			$abort = true;
		}
	}
	$test_log = array(); //Array from test name to warnings/errors
	if (in_array("test_output",$tests) || in_array("topp_output",$tests))
	{
		if (!file_exists("$bin_path/Testing/Temporary/LastTest.log"))
		{
			print "Error: For the tests 'test_output' and 'topp_output', the test log is needed!\n";
			print "       Please execute 'make test'.\n";
			$abort = true;
		}
		else
		{
			$current_test_name = "";
			$log_file = file("$bin_path/Testing/Temporary/LastTest.log");
			foreach ($log_file as $line)
			{
				if (ereg("[0-9]+/[0-9]+ Testing: (.*)",$line,$parts))
				{
					$current_test_name = trim($parts[1]);
				}
				if (beginsWith($line,"warning") || beginsWith($line,"Warning") || beginsWith($line,"error") || beginsWith($line,"Error"))
				{
					$test_log[$current_test_name][] = trim($line);
				}
			}
		}
	}
	if ($abort)
	{
		exit;
	}
	
	########################### MAINTAINERS ################################
	
	$files=array();
	exec("cd $src_path && find include/OpenMS -name \"*.h\" ! -name \"ui_*.h\" ! -name \"nnls.h\"", $files);
	exec("cd $src_path && find source -name \"*.C\" ! -regex \".*/EXAMPLES/.*\" ! -regex \".*/tools/.*\" ! -name \"*_moc.C\" ! -name \"moc_*.C\" ! -name \"*Template.C\"", $files);
	
	//look up Maintainer in first 40 lines of files
	
	$GLOBALS["maintainer_info"] = array();
	$files_todo = array();
	$GLOBALS["file_maintainers"] = array();
	
	foreach ($files as $f)
	{
		$maintainerline = array();
		exec("cd $src_path && head -40 $f | grep Maintainer ", $maintainerline);
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
	$makefile = file("$src_path/source/TEST/executables.cmake");
	foreach($makefile as $line)
	{
		$line = trim($line);
		if (strpos($line,"_test")!== FALSE)
		{
			//Appended conditional test
			if (strpos($line,"APPEND")!== FALSE)
			{
				$line = substr($line,strrpos($line,' ')+1,-1);
				$called_tests[] = "source/TEST/".$line.".C";
			}
			else
			{
				$called_tests[] = "source/TEST/".strtr($line,array("_1"=>"","_2"=>"")).".C";
			}
		}
	}
	if (in_array("doxygen_errors",$tests))
	{
		$doxygen_errors = array();
		$errorfile = file("$bin_path/doc/doxygen/doxygen-error.log");
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
		$file = file($src_path."/".$f);
		
		######################### load class info ################################
		$dont_load = array(
			"IsotopeCluster.h",
			"RangeUtils.h",
			"ComparatorUtils.h",
			"StatisticFunctions.h",
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
			"IsotopeWaveletConstants.h",
			"IsotopeWaveletCudaKernel.h",
			"IsotopeWaveletParallelFor.h",
			"openms_svn_revision.h",
			"openms_package_version.h"
			);

		if (!endsWith($f,"_impl.h") && endsWith($f,".h") && !in_array($basename,$dont_load))
		{
			$class_info = getClassInfo($bin_path,$f,$debug);
		}
		else
		{
			unset($class_info);
		}

		########################### guards ######################################
		if (in_array("guards",$tests))
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
								realOutput("Wrong header guard '$guard' in '$f' should be '$right_guard'",$user,$f);
							}					
							break;
						}
					}
					
					$class = trim(substr($f,strrpos($f,"/")+1));
					if ($i==count($file)-1 AND !in_array($class,$dont_report))
					{
						realOutput("Missing header guard in '$f' ",$user,$f);
					}
				}	
			}
		}
	
		########################### tab settings #####################################
		if (in_array("tab",$tests))
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
				realOutput("Missing tab settings in '$f'",$user,$f);
			}
		}

		########################### maintainers  #####################################
		if (in_array("maintainers",$tests))
		{
			if (endsWith($f,".h") )
			{
				# maintainer of test file
				if (in_array($testname,$files))
				{
					if (isset($file_maintainers[$testname]) && $file_maintainers[$testname] != $file_maintainers[$f])
					{
						realOutput("Inconsistent maintainers in '$f' and '$testname'",$user,$f);
						print "  '$file_maintainers[$testname]'<->'$file_maintainers[$f]'\n";
					}
				}
				# maintainer of source file
				$source_name = "source/".substr($f,15,-2).".C";
				if (in_array($source_name,$files))
				{
					if ($file_maintainers[$source_name] != $file_maintainers[$f])
					{
						realOutput("Inconsistent maintainers in '$f' and '$source_name'",$user,$f);
						print "  '$file_maintainers[$source_name]'<->'$file_maintainers[$f]'\n";
					}
				}
			}
		}

		########################### missing tests  #####################################
		if (in_array("missing_tests",$tests))
		{
			$dont_report = array(
				"/FORMAT/HANDLERS/",
				"/VISUAL/",
				"/CONCEPT/Types.h",
				"/FORMAT/FileTypes.h",
				"/CONCEPT/Macros.h",
				"/CONCEPT/Exception.h",
				"/CONCEPT/Benchmark.h",
				"/CONCEPT/Constants.h",
				"CONCEPT/ProgressLogger.h",
				"/config.h",
				"include/OpenMS/APPLICATIONS/TOPPViewBase.h",
				"include/OpenMS/APPLICATIONS/TOPPASBase.h",
				"include/OpenMS/APPLICATIONS/INIFileEditorWindow.h",
				"include/OpenMS/DATASTRUCTURES/SeqanIncludeWrapper.h",
				"include/OpenMS/ANALYSIS/DENOVO/CompNovoIdentificationBase.h",
				"include/OpenMS/ANALYSIS/DENOVO/CompNovoIonScoringBase.h",
				"_registerChildren.h",
				"DataReducer.h",
				"SchemaFile.h",
				"Serialization.h",
				"IsotopeCluster.h",
				"Param.h",
				"IsotopeWaveletCudaKernel.h",
				"IsotopeWaveletConstants.h",
				"IsotopeWaveletParallelFor.h",
				"include/OpenMS/openms_svn_revision.h",
				"include/OpenMS/openms_package_version.h"
				);

			if (endsWith($f,".h") && !endsWith($f,"_impl.h"))
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
						realOutput("Missing test for '$f'",$user,$f);
					}
					else if (!in_array($testname,$called_tests))
					{
						realOutput("Test not in executables.cmake for '$f'",$user,$f);
					}
				}
			}
		}

		########################### guards ######################################
		if (in_array("old_files",$tests))
		{
			$ignore = array(
				"SampleTreatment_test.C",
				"Exception_Base_test.C",
				"NumericDiff.C",
			);
			
			if (!in_array($basename,$ignore)  && !beginsWith($f,"source/APPLICATIONS/TOPP/")  && !beginsWith($f,"source/APPLICATIONS/UTILS/"))
			{
				if (endsWith($f,"_test.C"))
				{
					$hits = array();
					foreach($file as $line)
					{
						if (isIncludeLine($line,$include) && strpos($line,substr($basename,0,-7))!==FALSE )
						{
							$hits[] = "include/".$include;
						}
					}
					if (count($hits)==1)
					{
						//print "$f -> $hits[0]\n";
						if (!file_exists($src_path."/".$hits[0]))
						{
							realOutput("Outdated test file '$f'",$user,$f);
						}
					}
				}
				//source file -> look for header
				else if (endsWith($f,".C"))
				{
					if (!file_exists($src_path."/include/OpenMS/".substr($f,7,-2).".h"))
					{
						realOutput("Outdated source file '$f'",$user,$f);
					}
				}
			}
		}

		########################### @brief  #####################################
		if (in_array("brief",$tests))
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
								realOutput("No @brief description for '$parts[2]' in '$f'",$user,$f);
							}
						}
					}
				}
			}
		}
		
		########################### doxygen errors  ####################################
		if (in_array("doxygen_errors",$tests))
		{
			if (in_array($f,$doxygen_errors))
			{
				realOutput("Doxygen errors in '$f'",$user,$f);
				print "  See 'OpenMS/doc/doxygen/doxygen-error.log'\n";
			}
		}
		
		########################### test errors  #####################################
		
		if (isset($class_info) && in_array("check_test",$tests))
		{
			if (in_array($testname,$files))
			{
 				#parse test
 				$tmp = parseTestFile("$src_path/$testname");
 				$todo_tests = $tmp["todo"];
				$tests2 = $tmp["tests"];
								
				#compare declarations and tests
				$out = compareDeclarationsAndTests($class_info["public-long"],$tests2);
				
				# remove methods that can be tested although they are not defined
				$new_unknown = array();
				foreach ($out["unknown"] as $m)
				{
					//print ">>> '".$m."'  '$classname()'\n";
					if ($m != "$classname()" && $m != "~$classname()" && !beginsWith($m,"const $classname& operator=(") && !beginsWith($m, "$classname(const $classname&") )
					{
						$new_unknown[] = $m;
					}
				}
				$out["unknown"] = $new_unknown;
				
				#output
				if (count($out["missing"])!=0 || count($out["unknown"])!=0 || count($out["double"])!=0 || count($todo_tests)!=0)
				{
					realOutput("Test errors in '$f'",$user,$testname);
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
					if (count($out["double"])!=0)
					{
						$out["double"] = array_unique($out["double"]);
						print "  Methods tested several times:\n";
						foreach ($out["double"] as $m)
						{
							print "    - '$m'\n";
						}
					}
					if (count($todo_tests)!=0)
					{
						print "  Tests that contain 'TODO' or '????':\n";
						foreach ($todo_tests as $m)
						{
							print "    - '$m'\n";
						}
					}
				}
			}
		}
		
		
		
		############################## coding ##########################################
		if (isset($class_info) && in_array("coding",$tests) )
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
				//print "NP: '".$tmp."'\n";
				//Allow special methods to conflict with the coding convention
				if ( endswith($tmp,'Event') || //Qt events
				     endsWith($tmp,'load')  || endsWith($tmp,'save') || endsWith($tmp,'serialize') || //serialize methods
				     $tmp==$class_info["classname"] || $tmp=='~'.$class_info["classname"] || //constructor or destructor
				     $tmp=="operator=" || //assignment
				     $tmp=="startElement" || $tmp=="endElement" || $tmp=="characters" || //Xerces data handler
				     $tmp=="warning" || $tmp=="error" || $tmp=="fatalError" || $tmp=="resetErrors" //Xerces error handler
				   )
				{
					continue;
				}
				
				if (!endswith($tmp,'_'))
				{
					$out[] = "  - invalid non-public method name '$tmp'\n";
				}
				//check if there is a underscore in the middle, that must not be there
				else if (strpos(substr($tmp,0,-1),'_')!==FALSE)
				{
					$out[] = "  - invalid non-public method name '$tmp'\n";
				}
			}
			#output
			if (count($out)!=0)
			{
				realOutput("Coding convention violation in '$f'",$user,$f);
				foreach ($out as $o)
				{
					print $o;
				}
			}
		}
		
		
		########################### warnings test  #####################################
		if (in_array("test_output",$tests) )
		{
			if (endsWith($f,".h") && in_array($testname,$files))
			{
				$errors = array();
				if (isset($test_log[$classname."_test"])) $errors = array_unique($test_log[$classname."_test"]);
				if (count($errors)!=0)
				{
					realOutput("Error/warnings in test output of '$testname'",$user,$testname);
					foreach ($errors as $e)
					{
						print "  '$e'\n";
					}
				}
			}
		}

		######################### 'Id' keyword in tests  ###############################
		if (in_array("svn_keywords",$tests))
		{
			if (endsWith($f,".h") )
			{
				if (in_array($testname,$files))
				{
					$out = array();
					exec("svn proplist -v $src_path/$testname",$out);
					$kw = false;
					foreach ($out as $line)
					{
						if (strpos($line,"svn:keywords")!==FALSE)
						{
							$kw = true;
							if (strpos($line,"Id")===FALSE)
							{
								realOutput("svn:keyword 'Id' not set for '$testname'",$user,$testname);
							}
						}
					}
					if (!$kw)
					{
						realOutput("svn:keyword 'Id' not set for '$testname'",$user,$testname);
					}
				}	
			}
		}

		########################### DefaultParamHandler  #################################
		if (in_array("defaults",$tests))
		{
			if (endsWith($f,".h") && !endsWith($f,"_impl.h"))
			{
				//check if defaults are set in .h file
				$is_dph = false;
				$output = array();
				exec("grep -l defaults_.setValue $src_path/$f", $output);
				if ($output!=array())
				{
					$is_dph = true;
				}
				//check if defaults are set in .C file
				else
				{
					$c_file = "$src_path/source/".substr($f,15,-2).".C";
					if (file_exists($c_file))
					{
						exec("grep -l defaults_.setValue $c_file", $output);
						if ($output!=array())
						{
							$is_dph = true;
						}
					}
				}
				//check if reference to paramters docu page is present
				if ($is_dph)
				{
					$output = array();
					exec("grep -l OpenMS_".$classname.".parameters $src_path/$f", $output);
					if ($output==array())
					{
						realOutput("Missing reference to parameters page in '$f' unless abstract base class",$user,$f);
					}
				}
			}
		}

	}//End of files loop

	################### doxygen errors in .doxygen-files  ##########################
	if ($user == "all" && in_array("doxygen_errors",$tests) )
	{
		$file = file("$bin_path/doc/doxygen/doxygen-error.log");
		foreach ($file as $line)
		{
			$line = trim($line);
			if (ereg("(.*/[a-zA-Z0-9_]+\.doxygen):[0-9]+:",$line,$parts))
			{
				realOutput("Doxygen errors in '".$parts[1]."'",$user,"");
				print "  See 'OpenMS/doc/doxygen/doxygen-error.log'\n";
			}
		}
	}

	
	########################### warnings TOPP test  #################################
	if (in_array("topp_output",$tests))
	{
		$file_warnings = array();
		foreach ($test_log as $name => $warnings)
		{
			if (beginsWith($name,"TOPP_"))
			{
				$name = substr($name,5);
				$name = substr($name,0,strpos($name,'_'));
				$topp_file = "source/APPLICATIONS/TOPP/".$name.".C";
				if (in_array($topp_file,$files_todo))
				{
					if (!isset($file_warnings[$topp_file])) $file_warnings[$topp_file] = array();
					$file_warnings[$topp_file] = array_merge($file_warnings[$topp_file],$warnings);
				}
			}
		}
		//print errors/warnings bundled for each TOPP tool
		foreach($file_warnings as $file => $warnings)
		{
			realOutput("Error/warnings in TOPP tool test of '$file'",$user,$file);
			$warnings = array_unique($warnings);
			foreach ($warnings as $e)
			{
				print "  '$e'\n";
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
