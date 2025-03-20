<?php
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer:$
# $Authors: Marc Sturm $
# --------------------------------------------------------------------------
error_reporting(E_ERROR| E_WARNING| E_PARSE| E_NOTICE);
ini_set('memory_limit', '1000M');
date_default_timezone_set('UTC');

include "common_functions.php";

######################## helper functions ###############################
function printUsage() {
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
  print "  -r         report as Testing.xml for ctest submission\n";
  print "  --help     shows this help\n";
}

function realOutput($text, $user, $filename) {
  print "------> ";
  print $text;
  if ($user == "all" && $filename != "")
  {
    print " (".$GLOBALS["file_maintainers"][$filename].")";
  }
  print "\n";
  #error count
  if ($filename != "")
  {
    if (isset($GLOBALS["file_maintainers"][$filename]))
    {
      foreach (explode(", ", $GLOBALS["file_maintainers"][$filename]) as $u)
      {
        $GLOBALS["maintainer_info"][trim($u)]["errors"]++;
      }
    }
  }
}

/**
  Store results of tests for Ctest-Reporting
*/
function reportTestResult($text, $user, $testname, $file, $result) {
  # check if we have the correct user
  if ($user == "all" && $file != "")
  {
    $user = $GLOBALS["file_maintainers"][$file];
  }
  $reportAs                                  = $file."_".$testname;
  $GLOBALS["TestList"][$reportAs]            = array();
  $GLOBALS["TestList"][$reportAs]["message"] = $text;
  $GLOBALS["TestList"][$reportAs]["result"]  = $result;
  $GLOBALS["TestList"][$reportAs]["user"]    = $user;
}

function xmlentities($string) {
  return str_replace(array("&", "<", ">", "\"", "'"), array("&amp;", "&lt;", "&gt;", "&quot;", "&apos;"), $string);
}



######################## declarations ###############################
$GLOBALS["all_tests"] = array(
  "(none)"               => "performs all tests",
  "maintainers"          => "check if maintainers are consistent in header, source and test file",
  "missing_tests"        => "check for missing tests",
  "old_files"            => "check for unneeded .cpp files",
  "brief"                => "check for doxygen tag @brief",
  "doxygen_errors"       => "check for errors in dogygen-error.log",
  "check_test"           => "check the class test for completeness",
  "test_output"          => "check test output for warnings and errors",
  "topp_output"          => "check TOPP test output for warnings and errors",
  "coding"               => "check if coding convention is followed",
  "defaults"             => "check if DefautParamHandler classes all have a linked parameters section",
  "license"              => "check if the license header is correctly included in the file",
  "inactive-maintainers" => "checks if the found maintainers are listed in the ACTIVE_MAINTAINERS list",
);

$options = array(
  "-u",
  "-t",
  "-d",
  "-s",
);

$flags = array(
  "-x",
  "-r",
);

######################## parameter handling ###############################
#no parameters
if ($argc == 1 || in_array("--help", $argv))
{
  printUsage();
  exit;
}

#wrong parameter count
for ($i = 1;$i < count($argv);++$i)
{
  #option
  if (beginsWith($argv[$i], "-"))
  {
    #no path given
    if ($i < 3)
    {
      print "\nError: No paths given!\n\n";
      printUsage();
      exit;
    }
    #not registered option or flag
    if (!in_array($argv[$i], array_merge($options, $flags)))
    {
      print "\nError: Unregistered option '".$argv[$i]."'!\n\n";
      printUsage();
      exit;
    }
    # no argument to an option
    if (in_array($argv[$i], $options) && (!isset($argv[$i+1]) || beginsWith($argv[$i+1], "-")))
    {
      print "\nError: No argument to option '".$argv[$i]."'!\n\n";
      printUsage();
      exit;
    }
  }
}

#remove slash from path
if (endsWith($argv[1], "/"))
{
  $src_path = substr($argv[1], 0,-1);
}
else
{
  $src_path = $argv[1];
}

#remove slash from path
if (endsWith($argv[2], "/"))
{
  $bin_path = substr($argv[2], 0,-1);
}
else
{
  $bin_path = $argv[2];
}

#debug
$debug = 0;
if (in_array("-d", $argv))
{
  $debug = $argv[array_search("-d", $argv)+1];
}

#user
$user = "all";
if (in_array("-u", $argv))
{
  $user = $argv[array_search("-u", $argv)+1];
}

#test
$tests = array_keys($GLOBALS["all_tests"]);
unset($tests[0]);
if (in_array("-t", $argv))
{
  $test = $argv[array_search("-t", $argv)+1];
  if (!in_array($test, $tests))
  {
    print "\nError: Unknown test '$test'!\n\n";
    printUsage();
    exit;
  }
  $tests = array(
    $test,
  );
}

#skip
if (in_array("-s", $argv))
{
  $skip = explode(',', $argv[array_search("-s", $argv)+1]);
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
  $tests = array_diff($tests, $skip);
}
$tests = array_values($tests);

#ctestReporting
$ctestReporting = false;
$ctestReportingPath = "";
if (in_array("-r", $argv))
{
  $ctestReporting = true;
}

#doxygen XML output
$rebuilt_xml = true;
if (in_array("-x", $argv))
{
  $rebuilt_xml = false;
}

######################## doxygen XML output ############################
if ($rebuilt_xml)
{
  if ($debug > 0)
  {
    print "Rebuilding doxygen XML output\n";
  }
  exec("cd $bin_path && cmake --build $bin_path --target doc_xml");
  if ($debug > 0)
  {
    print "Done\n";
  }
}

########################### NEEDED FILES ###############################
$abort = false;
$out = array();
exec("cd $bin_path/doc/xml_output/ && ls -a *.xml | wc -l", $out);
if (trim($out[0]) < 100)
{
  print "Error: For this script, doxygen XML output is needed!\n";
  print "       Please execute 'make doc_xml' first!.\n";
  $abort = true;
}
if (in_array("doxygen_errors", $tests))
{
  if (!file_exists("$bin_path/doc/doxygen/doxygen-error.log"))
  {
    print "Error: For the 'doxygen_errors' test, the file '$bin_path/doc/doxygen/doxygen-error.log' is needed!\n";
    print "       Please execute 'make doc_internal' first!'.\n";
    $abort = true;
  }
}
$test_log = array();
//Array from test name to warnings/errors
if (in_array("test_output", $tests) || in_array("topp_output", $tests))
{
  # try to find a possible test log from ctest
  #  * "$bin_path/Testing/Temporary/LastTest.log"
  #  * "$bin_path/Testing/TAG ->
  $lastTestFile = "$bin_path/Testing/Temporary/LastTest.log";
  if (!file_exists($lastTestFile))
  {
    if (file_exists("$bin_path/Testing/TAG"))
    {
      $taginfo = array();
      exec("cat $bin_path/Testing/TAG", $taginfo);
      if (file_exists("$bin_path/Testing/Temporary/LastTest_".$taginfo[0].".log"))
      {
        $lastTestFile = "$bin_path/Testing/Temporary/LastTest_".$taginfo[0].".log";
      }
      else
      {
        print "Error: For the tests 'test_output' and 'topp_output', the test log is needed!\n";
        print "       Please execute 'make test'.\n";
        $abort = true;
      }
    }
  }

  if (!$abort)
  {
    $current_test_name = "";
    $log_file = file($lastTestFile);
    foreach ($log_file as $line)
    {
      if ( preg_match('/[0-9]+\/[0-9]+ Testing: (.*)/', $line, $parts))
      {
        $current_test_name = trim($parts[1]);
      }
      if (beginsWith($line, "warning") || beginsWith($line, "Warning") || beginsWith($line, "error") || beginsWith($line, "Error"))
      {
        $test_log[$current_test_name][] = trim($line);
      }
    }
  }
}
## check if all files required to report in CDash are present
if ($ctestReporting)
{
  # check if we find the last run log
  if (!file_exists("$bin_path/Testing/TAG"))
  {
    print "Error: Missing nightly test information $bin_path/Testing/TAG!\n";
    $abort = true;
  }
  else
  {
    # find the time information
    $taginfo = array();
    exec("cat $bin_path/Testing/TAG", $taginfo);
    $ctestReportingPath = "$bin_path/Testing/$taginfo[0]";
    if (!file_exists("$bin_path/Testing/$taginfo[0]"))
    {
      $abort = true;
    }
  }
}

# read active maintainers file
$activemaintainers = array();
if (!file_exists("$src_path/tools/ACTIVE_MAINTAINERS"))
{
  print "Error: Missing $src_path/tools/ACTIVE_MAINTAINERS file!\n";
  $abort = true;
}
else
{
  exec("cat $src_path/tools/ACTIVE_MAINTAINERS", $activemaintainers);
}

if ($abort)
{
  exit;
}

########################### MAINTAINERS ################################
$files = array();

$includePaths = array("src/openms/include/OpenMS",
                      "src/openms_gui/include/OpenMS",
                      "src/openswathalgo/include/OpenMS",
                      "src/tests/class_tests/openms/include/OpenMS");

$sourcePaths = array("src/openms/source",
                     "src/openms_gui/source",
                     "src/openswathalgo/source",
                     "src/tests/class_tests/openms/source",
                     "src/tests/class_tests/openms_gui/source",
                     "src/tests/class_tests/openswathalgo/",
                     "src/topp"
                    );

exec("cd $src_path && find ".implode(" ", $includePaths)." -name \"*.h\" ! -name \"ui_*.h\" ! -name \"nnls.h\" ! -name \"MSNumpress*.h\"", $files);
exec("cd $src_path && find ".implode(" ", $sourcePaths)." -name \"*.cpp\" ! -regex \".*/EXAMPLES/.*\" ! -regex \".*/tools/.*\" ! -name \"*_moc.cpp\" ! -name \"moc_*.cpp\" ! -name \"*Template.cpp\" ! -name \"MSNumpress*.cpp\"", $files);

//look up Maintainer in first 40 lines of files
$GLOBALS["maintainer_info"]  = array();
$files_todo                  = array();
$GLOBALS["file_maintainers"] = array();
$GLOBALS["TestList"]         = array();

foreach ($files as $f)
{
  $maintainerline = array();
  exec("cd $src_path && head -40 $f | grep Maintainer ", $maintainerline);
  if (count($maintainerline) == 0)
  {
    if ($user == "all")
    {
      print "No Maintainer lines in '$f'\n";
      $GLOBALS["file_maintainers"][$f] = "";
    }
  }
  elseif (count($maintainerline) > 1)
  {
    if ($user == "all")
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
      if ($user == "all")
      {
        print "No Maintainer for '$f'\n";
        $GLOBALS["file_maintainers"][$f] = "";
      }
    }
    else
    {
      $GLOBALS["file_maintainers"][$f] = implode(", ", $maintainers);

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
      if ($user != "all")
      {
        foreach ($maintainers as $m)
        {
          $dist = levenshtein($m, $user);
          if ($dist == 0)
          {
            $files_todo[] = $f;
          }
          elseif ($dist <= 4)
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
$makefile = file("$src_path/src/tests/class_tests/openms/executables.cmake");
foreach ($makefile as $line)
{
  $line = trim($line);
  if (strpos($line, "_test") !== FALSE)
  {
    //Appended conditional test
    if (strpos($line, "APPEND") !== FALSE)
    {
      $line = substr($line, strrpos($line, ' ')+1,-1);
      $called_tests[] = "src/tests/class_tests/openms/source/".$line.".cpp";
    }
    else
    {
      $called_tests[] = "src/tests/class_tests/openms/source/".strtr($line, array("_1" => "", "_2" => "")).".cpp";
    }
  }
}
if (in_array("doxygen_errors", $tests))
{
  $doxygen_errors = array();
  $errorfile = file("$bin_path/doc/doxygen/doxygen-error.log");
  foreach ($errorfile as $line)
  {
    if ( preg_match('/(.*\/[a-zA-Z0-9_]+\.[hC]):[0-9]+:/', $line, $parts))
    {
      //skip warning where doxygen cannot resolve members
      if (strpos($line, "no uniquely matching class member") === FALSE && strpos($line, "no matching class member") === FALSE)
      {
        $pos = strpos($parts[1], "source/");
        if ($pos !== FALSE)
        {
          $doxygen_errors[] = substr($parts[1], $pos);
        }
        $pos = strpos($parts[1], "include/OpenMS/");
        if ($pos !== FALSE)
        {
          $doxygen_errors[] = substr($parts[1], $pos);
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
  if ($debug > 0)
  {
    print "File name: '$f'\n";
  }

  //file name (without path)
  $basename = basename($f);
  //class name (for source and header files)
  $classname = substr($basename, 0,-2);
  //test name (for source and header files)
  $testnames = array(
    "src/tests/class_tests/openms/source/".$classname."_test.cpp",
    "src/tests/class_tests/openms_gui/source/".$classname."_test.cpp"
                    );

  // file content
  $file = file($src_path."/".$f);

  ######################### load class info ################################
  $dont_load = array(
    "IsotopeCluster.h",
    "RangeUtils.h",
    "StatisticFunctions.h",
    "KernelTraits.h",
    "StandardTypes.h",
    "DimensionDescription.h",
    "MathFunctions.h",
    "ClassTest.h",
    "LayerData.h",
    "config.h",
    "OpenSwathAlgoConfig.h",
    "XMLSchemes.h",
    "Serialization.h",
    "Exception.h",
    "GlobalExceptionHandler.h",
    "Types.h",
    "TypeAsName.h",
    "PrecisionWrapper.h",
    "Macros.h",
    "Benchmark.h",
    "Constants.h",
    "openms_package_version.h",
  );

  if (!endsWith($f, "_impl.h") && endsWith($f, ".h") && !in_array($basename, $dont_load))
  {
    $class_info = getClassInfo($bin_path, $f, $debug);
  }
  else
  {
    unset($class_info);
  }

  ########################### maintainers #####################################
  if (in_array("maintainers", $tests))
  {
    if (endsWith($f, ".h"))
    {
      # maintainer of test file
      $testname = current(array_intersect($testnames, $files));
      if ($testname == FALSE)
      {
        if (isset($file_maintainers[$testname]) && $file_maintainers[$testname] != $file_maintainers[$f])
        {
          $message = "Inconsistent maintainers in '$f' and '$testname'";
          realOutput($message, $user, $f);
          print "  '$file_maintainers[$testname]'<->'$file_maintainers[$f]'\n";
          reportTestResult($message, $user, "maintainers_test", $f, false);
        }
        else
        {
          #report test result to ctest
          reportTestResult("", $user, "maintainers_test", $f, true);
        }
      }

      # maintainer of source file
      $source_name = "source/".substr($f, 15,-2).".cpp";
      if (in_array($source_name, $files))
      {
        if ($file_maintainers[$source_name] != $file_maintainers[$f])
        {
          $message = "Inconsistent maintainers in '$f' and '$source_name'";
          realOutput("Inconsistent maintainers in '$f' and '$source_name'", $user, $f);
          print "  '$file_maintainers[$source_name]'<->'$file_maintainers[$f]'\n";
          reportTestResult($message, $user, "maintainers_source", $f, false);
        }
        else
        {
          reportTestResult("", $user, "maintainers_source", $f, true);
        }
      }
    }
  }
  ########################### (in)active maintainer #####################################
  if (in_array("inactive-maintainers", $tests))
  {
    if (endsWith($f, ".h"))
    {
      # check if the maintainer is set
      if (isset($file_maintainers[$f]) && $file_maintainers[$f] != '')
      {
        $all_active = true;
        $maintainers = explode(',', $file_maintainers[$f]);
        foreach ($maintainers as $m)
        {
          $realMaintainer = trim($m);
          if (!in_array($realMaintainer, $activemaintainers))
          {
            reportTestResult("Inactive maintainer ($realMaintainer) given in '$f'.", $user, "active_maintainer", $f, false);
            realOutput("Inactive maintainer ($realMaintainer) given in '$f'.", $user, $f);
            $all_active = false;
          }
        }
        if ($all_active)
        {
          reportTestResult("All maintainers given in '$f' are active maintainers.", $user, "active_maintainer", $f, true);
        }
      }
      else
      {
        reportTestResult("No maintainer given in '$f'.", $user, "active_maintainer", $f, false);
      }
    }
  }


  ########################### missing tests  #####################################
  if (in_array("missing_tests", $tests))
  {
    $dont_report = array(
      "/FORMAT/HANDLERS/",
      "/VISUAL/",
      "/CONCEPT/Types.h",
      "/CONCEPT/PrecisionWrapper.h",
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
      "_registerChildren.h",
      "DataReducer.h",
      "SchemaFile.h",
      "Serialization.h",
      "IsotopeCluster.h",
      "Param.h",
      "include/OpenMS/openms_package_version.h",
      "include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h",
      "include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h",
      "include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h",
      "include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSCached.h",
    );

    if (endsWith($f, ".h") && !endsWith($f, "_impl.h"))
    {
      $ignore = false;
      foreach ($dont_report as $i)
      {
        if (strpos($f, $i) !== FALSE)
        {
          $ignore = true;
        }
      }

      // Exclude all OpenSwathAlgo from the tests since they are not properly recognized
      if (preg_match("/OPENSWATHALGO/i", $f))
      {
        $ignore = true;
      }

      if (!$ignore)
      {
        # check if test exists
        $testname = current(array_intersect($testnames, $files));
        if (!$testname)
        {
          realOutput("Missing test for '$f'", $user, $f);
          reportTestResult("Missing test for '$f'", $user, "missing_tests", $f, false);
        }
        elseif (!in_array($testname, $called_tests))
        # check if test was executed
        {
          realOutput("Test not in executables.cmake for '$f'", $user, $f);
          reportTestResult("Test not in executables.cmake for '$f'", $user, "missing_tests", $f, false);
        }
        else
        {
          reportTestResult("", $user, "missing_tests", $f, true);
        }
      }
    }
  }

  ########################### guards ######################################
  if (in_array("old_files", $tests))
  {
    $ignore = array(
      "SampleTreatment_test.cpp",
      "Exception_Base_test.cpp",
      "NumericDiff.cpp",
      "ExampleLibraryFile.cpp",
      "TestExternalCode.cpp",
    );

    if (!in_array($basename, $ignore) && !beginsWith($f, "src/topp/") && !beginsWith($f, "src/openms_gui/source/VISUAL/APPLICATIONS/GUITOOLS/"))
    {
      $message = "";
      $result = true;
      if (endsWith($f, "_test.cpp"))
      {
        $hits = array();
        foreach ($file as $line)
        {
          if (isIncludeLine($line, $include) && strpos($line, substr("/".$basename.".h", 0,-7)) !== FALSE)
          {
            # check all potential include locations
            $incl_path = substr($include, 7);
            foreach ($includePaths as $incl)
            {
              array_push($hits, $incl."/".$incl_path);
            }
          }
        }
        if (count($hits) > 1)
        {
          $file_exists = false;
          foreach($hits as $h)
          {
            if(endsWith($f, "String_test.cpp"))
            {
              echo "String_test -> ".$src_path."/".$h."\n";
            }
            $file_exists |= file_exists($src_path."/".$h);
          }

          if (!$file_exists)
          {
            $message = "Outdated test file '$f'";
            $result = false;
            realOutput($message, $user, $f);
          }

        }
      }
      //source file -> look for header
      elseif (endsWith($f, ".cpp"))
      {
        # $h_file = substr($f, 7,-2).".h";
        $h_file = $src_path."/".str_replace("/source", "/include/OpenMS", $f);
        $h_file = preg_replace("/cpp$/", "h", $h_file);
        if (!file_exists($h_file))
        {
          $message = "Outdated source file '$f' (no header -> $h_file)";
          $result = false;
          realOutput($message, $user, $f);
        }
      }
      reportTestResult($message, $user, "old_files", $f, $result);
    }
  }

  ########################### @brief  #####################################
  if (in_array("brief", $tests))
  {
    if (endsWith($f, ".h"))
    {
      for ($i = 0;$i < count($file);++$i)
      {
        #match a class declaration
        if (preg_match("/(class|struct)[\s]+[\w]+[\s]*({|:[^:]|$)/i", $file[$i]))
        {
          #take the second line too in case : or { come in the next line
          if (isset($file[$i+1]) && preg_match("/(class|struct)[\s]+([\w]+)[\s]*({|:[^:])/i", $file[$i].$file[$i+1], $parts))
          {
            $brief = false;
            #search for /// comment
            if (strpos($file[$i-1], "///") !== FALSE)
            {
              $brief = true;
            }
            # backward search for @brief until comment closes
            for ($j = $i;$j != 0;--$j)
            {
              if (strpos($file[$j], "@brief") !== FALSE)
              {
                $brief = true;
                break;
              }
              if (strpos($file[$j], "/**") !== FALSE)
              {
                break;
              }
            }
            if (!$brief)
            {
              realOutput("No @brief description for '$parts[2]' in '$f'", $user, $f);
              reportTestResult("No @brief description for '$parts[2]' in '$f'", $user, "brief", $f, false);
            }
            else
            {
              reportTestResult("", $user, "brief", $f, true);
            }
          }
        }
      }
    }
  }

  ########################### doxygen errors  ####################################
  if (in_array("doxygen_errors", $tests))
  {
    if (in_array($f, $doxygen_errors))
    {
      realOutput("Doxygen errors in '$f'", $user, $f);
      print "  See 'OpenMS/doc/doxygen/doxygen-error.log'\n";
    }
  }

  ########################### test errors  #####################################
  if (isset($class_info) && in_array("check_test", $tests))
  {
    # get the actual test name -> if it exists
    $testname = current(array_intersect($testnames, $files));
    if (!$testname)
    {
      #parse test
      $tmp        = parseTestFile("$src_path/$testname");
      $todo_tests = $tmp["todo"];
      $tests2     = $tmp["tests"];

      #compare declarations and tests
      $out = compareDeclarationsAndTests($class_info["test-name"], $tests2);

      # remove methods that can be tested although they are not defined
      $new_unknown = array();
      foreach ($out["unknown"] as $m)
      {
        //print ">>> '".$m."'  '$classname()'\n";
        if ($m != "$classname()" && $m != "~$classname()" && !beginsWith($m, "const $classname& operator=(") && !beginsWith($m, "$classname(const $classname&"))
        {
          $new_unknown[] = $m;
        }
      }
      $out["unknown"] = $new_unknown;

      #output
      if (count($out["missing"]) != 0 || count($out["unknown"]) != 0 || count($out["double"]) != 0 || count($todo_tests) != 0)
      {
        realOutput("Test errors in '$f'", $user, $testname);
        $message = "Test errors in '$f'\n";
        if (count($out["unknown"]) != 0)
        {
          print "  Tests of unknown methods:\n";
          $message .= "  Tests of unknown methods:\n";
          foreach ($out["unknown"] as $u)
          {
            print "    - '$u'\n";
            $message .= "    - '$u'\n";
          }
        }
        if (count($out["missing"]) != 0)
        {
          print "  Missing tests:\n";
          $message .= "  Missing tests:\n";
          foreach ($out["missing"] as $m)
          {
            print "    - '$m'\n";
            $message .= "    - '$m'\n";
          }
        }
        if (count($out["double"]) != 0)
        {
          $out["double"] = array_unique($out["double"]);
          print "  Methods tested several times:\n";
          $message .= "  Methods tested several times:\n";
          foreach ($out["double"] as $m)
          {
            print "    - '$m'\n";
            $message .= "    - '$m'\n";
          }
        }
        if (count($todo_tests) != 0)
        {
          print "  Tests that contain 'TODO' or '????':\n";
          $message .= "  Tests that contain 'TODO' or '????':\n";
          foreach ($todo_tests as $m)
          {
            print "    - '$m'\n";
            $message .= "    - '$m'\n";
          }
        }
        reportTestResult($message, $user, "check_test", $f, false);
      }
      else
        # Test is OK
        {
          reportTestResult("", $user, "check_test", $f, true);
      }
    }
  }

  ############################## coding ##########################################
  if (isset($class_info) && in_array("coding", $tests))
  {
    $out = array();
    #variables
    foreach ($class_info["variables"] as $tmp)
    {
      if (!endswith($tmp, '_'))
      {
        $out[] = "  - invalid non-public variable name '$tmp'\n";
      }
    }
    # non-public member
    foreach ($class_info["non-public"] as $tmp)
    {
      //print "NP: '".$tmp."'\n";
      //Allow special methods to conflict with the coding convention
      if (endswith($tmp, 'Event') || //Qt events
      endsWith($tmp, 'load') || endsWith($tmp, 'save') || endsWith($tmp, 'serialize') || //serialize methods
      $tmp == $class_info["classname"] || $tmp == '~'.$class_info["classname"] || //constructor or destructor
      $tmp == "operator=" || //assignment
      $tmp == "startElement" || $tmp == "endElement" || $tmp == "characters" || //Xerces data handler
      $tmp == "warning" || $tmp == "error" || $tmp == "fatalError" || $tmp == "resetErrors")
        //Xerces error handler
        {
          continue;
      }

      if (!endswith($tmp, '_'))
      {
        $out[] = "  - invalid non-public method name '$tmp'\n";
      }
      //check if there is a underscore in the middle, that must not be there
      elseif (strpos(substr($tmp, 0,-1), '_') !== FALSE)
      {
        $out[] = "  - invalid non-public method name '$tmp'\n";
      }
    }
    #output
    if (count($out) != 0)
    {
      $message = "Coding convention violation in '$f'\n";
      realOutput("Coding convention violation in '$f'", $user, $f);
      foreach ($out as $o)
      {
        print $o;
        $message .= $o."\n";
      }
      reportTestResult($message, $user, "coding", $f, false);
    }
    else
    {
      reportTestResult("", $user, "coding", $f, true);
    }
  }

  ########################### warnings test  #####################################
  if (in_array("test_output", $tests))
  {
    # get the actual test name
    $testname = current(array_intersect($testnames, $files));

    if (endsWith($f, ".h") && !$testname)
    {
      $errors = array();
      if (isset($test_log[$classname."_test"]))
        $errors = array_unique($test_log[$classname."_test"]);
      if (count($errors) != 0)
      {
        realOutput("Error/warnings in test output of '$testname'", $user, $testname);
        $message = "Error/warnings in test output of '$testname'\n";
        foreach ($errors as $e)
        {
          print "  '$e'\n";
          $message .= "  '$e'\n";
        }
        reportTestResult($message, $user, "test_output", $f, false);
      }
      else
      {
        reportTestResult("", $user, "test_output", $f, true);
      }
    }
  }

  ########################### DefaultParamHandler  #################################
  if (in_array("defaults", $tests))
  {
    // don't report e.g. abstract base classes
    $dont_report = array(
      "src/openms/include/OpenMS/VISUAL/PlotCanvas.h",
      "src/openms/include/OpenMS/APPLICATIONS/TOPPViewBase.h",
      "src/openms/include/OpenMS/APPLICATIONS/TOPPASBase.h",
      "src/openms/include/OpenMS/FEATUREFINDER/BaseModel.h",
      "src/openms/include/OpenMS/FEATUREFINDER/LevMarqFitter1D.h",
    );

    $ignore = false;
    foreach ($dont_report as $i)
    {
      if (strpos($f, $i) !== FALSE)
      {
        $ignore = true;
      }
    }

    if (endsWith($f, ".h") && !endsWith($f, "_impl.h") && !$ignore)
    {
      //check if defaults are set in .h file
      $is_dph = false;
      $output = array();
      exec("grep -l defaults_.setValue $src_path/$f", $output);
      if ($output != array())
      {
        $is_dph = true;
      }
      //check if defaults are set in .cpp file
      else
      {
        $c_file = "$src_path/source/".substr($f, 15,-2).".cpp";
        if (file_exists($c_file))
        {
          exec("grep -l defaults_.setValue $c_file", $output);
          if ($output != array())
          {
            $is_dph = true;
          }
        }
      }
      //check if reference to parameters docu page is present
      if ($is_dph)
      {
        $output = array();
        exec("grep -l OpenMS_".$classname.".parameters $src_path/$f", $output);
        if ($output == array())
        {
          realOutput("Missing reference to parameters page in '$f' unless abstract base class", $user, $f);
          reportTestResult("Missing reference to parameters page in '$f' unless abstract base class", $user, "defaults", $f, false);
        }
        else
        {
          reportTestResult("", $user, "defaults", $f, true);
        }
      }
    }
  }
  // if (in_array("defaults",$tests))
  ########################### Check license header  ################################
  if (in_array("license", $tests))
  {
    $LICENSE_FILE = file($src_path."/LICENSE");
    $LICENSE = array();
    foreach ($LICENSE_FILE as $line)
    {
      array_push($LICENSE, rtrim($line));
    }
    $LICENSE_LEN = count($LICENSE);

    $isEqual       = True;
    $offendingLine = "";
    $message       = "";
    # line breaks to remove
    $breaks = array(
      "\r\n",
      "\n",
      "\r",
    );

    # every file should contain at least as much lines as the LICENSE header
    if (count($file) >= $LICENSE_LEN)
    {
      foreach ($LICENSE as $index => $line)
      {
        $licenseLine = rtrim("// ".str_replace($breaks, "", $line));
        $fileLine = str_replace($breaks, "", rtrim($file[$index]));
        if (strcmp($licenseLine, $fileLine) != 0)
        {
          $isEqual       = False;
          $offendingLine = "\tshould be: ".$licenseLine."\n\tbut is:    ".$fileLine;
          $message       = "Invalid LICENSE header in '$f'";
          break;
        }
      }
    }
    else
    {
      $isEqual = false;
      $message = "Missing LICENSE header in '$f'";
    }

		if ($debug > 0 && $offendingLine != "")
		{
			$message .= "\nOffending line: \n".$offendingLine;
		}

    # report result
    if (!$isEqual)
    {
      realOutput($message, $user, $f);
      reportTestResult($message, $user, "license", $f, false);
    }
    else
    {
      reportTestResult($message, $user, "license", $f, true);
    }
  }
}
//End of files loop
################### doxygen errors in .doxygen-files  ##########################
if ($user == "all" && in_array("doxygen_errors", $tests))
{
  $file = file("$bin_path/doc/doxygen/doxygen-error.log");
  foreach ($file as $line)
  {
    $line = trim($line);
    if ( preg_match('/(.*\/[a-zA-Z0-9_]+\.doxygen):[0-9]+:/', $line, $parts))
    {
      realOutput("Doxygen errors in '".$parts[1]."'", $user, "");
      print "  See 'OpenMS/doc/doxygen/doxygen-error.log'\n";
    }
  }
}

########################### warnings TOPP test  #################################
if (in_array("topp_output", $tests))
{
  $file_warnings = array();
  foreach ($test_log as $name => $warnings)
  {
    if (beginsWith($name, "TOPP_"))
    {
      $name      = substr($name, 5);
      $name      = substr($name, 0, strpos($name, '_'));
      $topp_file = "src/topp/".$name.".cpp";
      if (in_array($topp_file, $files_todo))
      {
        if (!isset($file_warnings[$topp_file]))
          $file_warnings[$topp_file] = array();
        $file_warnings[$topp_file] = array_merge($file_warnings[$topp_file], $warnings);
      }
    }
  }
  //print errors/warnings bundled for each TOPP tool
  foreach ($file_warnings as $file => $warnings)
  {
    realOutput("Error/warnings in TOPP tool test of '$file'", $user, $file);
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

if ($ctestReporting)
{
  # we assume that ctest was already executed so we can use the original file as
  # template
  print "Report on ".count($GLOBALS["TestList"])." tests.";

  # load template head
  $template = file($ctestReportingPath."/Test.xml");
  $newTestFile = array();
  foreach ($template as $line)
  {
    array_push($newTestFile, $line);
    if (trim($line) == "<Testing>")
    {
      break;
    }
  }

  array_push($newTestFile, "<StartDateTime>".date("M j G:i:s T")."</StartDateTime>\n");
  array_push($newTestFile, "<StartTestTime>".time()."</StartTestTime>\n");

  # define all executed tests
  array_push($newTestFile, "<TestList>\n");
  foreach ($GLOBALS["TestList"] as $test => $testdata)
  {
    array_push($newTestFile, "    <Test>".$test."</Test>\n");
  }
  array_push($newTestFile, "</TestList>\n");

  # and now the actual results
  foreach ($GLOBALS["TestList"] as $test => $testdata)
  {
    /*
    $GLOBALS["TestList"][$reportAs]["message"] = $text;
    $GLOBALS["TestList"][$reportAs]["result"]  = $result;
    $GLOBALS["TestList"][$reportAs]["user"]    = $user;
    */
    $testName = $testdata["user"]."-".$test;

    array_push($newTestFile, "    <Test Status=\"".($testdata["result"]?"passed":"failed")."\">\n");
    # failed, passed
    array_push($newTestFile, "      <Name>".$testName."</Name>\n");
    array_push($newTestFile, "      <Path>./tools/</Path>\n");
    array_push($newTestFile, "      <FullName>".$testName."</FullName>\n");
    array_push($newTestFile, "      <FullCommandLine>checker ".$test."</FullCommandLine>\n");
    array_push($newTestFile, "      <Results>\n");
    array_push($newTestFile, "              <NamedMeasurement type=\"numeric/double\" name=\"Execution Time\"><Value>0.001</Value></NamedMeasurement>\n");
    array_push($newTestFile, "              <NamedMeasurement type=\"text/string\" name=\"Completion Status\"><Value>Completed</Value></NamedMeasurement>\n");
    array_push($newTestFile, "              <NamedMeasurement type=\"text/string\" name=\"Command Line\">\n");
    array_push($newTestFile, "                <Value>checker.php</Value>\n");
    array_push($newTestFile, "              </NamedMeasurement>\n");
    array_push($newTestFile, "              <Measurement>\n");
    array_push($newTestFile, "                <Value>\n");
    array_push($newTestFile, xmlentities($testdata["message"]));
    array_push($newTestFile, "                </Value>\n");
    array_push($newTestFile, "              </Measurement>\n");
    array_push($newTestFile, "      </Results>\n");
    array_push($newTestFile, "    </Test>\n");
  }

  array_push($newTestFile, "<EndDateTime>".date("M j G:i:s T")."</EndDateTime>\n");
  array_push($newTestFile, "<EndTestTime>".time()."</EndTestTime>\n");
  array_push($newTestFile, "<ElapsedMinutes></ElapsedMinutes>\n");
  array_push($newTestFile, "</Testing>\n");
  array_push($newTestFile, "</Site>\n");

  rename($ctestReportingPath."/Test.xml", $ctestReportingPath."/RegularTest.xml");
  file_put_contents($ctestReportingPath."/Test.xml", $newTestFile);

  /*
  <?xml version="1.0" encoding="UTF-8"?>
  <Site BuildName="Darwin-clang++"
          BuildStamp="20121021-2300-Nightly"
          Name="laphroaig.imp.fu-berlin.de"
          Generator="ctest-2.8.9"
          CompilerName="/usr/bin/clang++"
          OSName="Mac OS X"
          Hostname="laphroaig.imp.fu-berlin.de"
          OSRelease="10.7.5"
          OSVersion="11G63"
          OSPlatform="x86_64"
          Is64Bits="1"
          VendorString="GenuineIntel"
          VendorID="Intel Corporation"
          FamilyID="6"
          ModelID="37"
          ProcessorCacheSize="32768"
          NumberOfLogicalCPU="4"
          NumberOfPhysicalCPU="2"
          TotalVirtualMemory="512"
          TotalPhysicalMemory="8192"
          LogicalProcessorsPerPhysical="8"
          ProcessorClockFrequency="2660"
  >
  <Testing>
  <StartDateTime>Oct 22 18:36 CEST</StartDateTime>
  <StartTestTime>1350923805</StartTestTime>
  <TestList>
    <Test></Test>
  </TestList>
  <Test Status="passed">
    <Name>SomeTool_test</Name>
    <Path>./source/TEST</Path>
    <FullName>./source/TEST/SomeTool_test</FullName>
    <FullCommandLine>/Users/aiche/dev/openms/openms-src/build/ninja/source/TEST/bin/SomeTool_test</FullCommandLine>
    <Results>
            <NamedMeasurement type="numeric/double" name="Execution Time"><Value>0.469694</Value></NamedMeasurement>
            <NamedMeasurement type="text/string" name="Completion Status"><Value>Completed</Value></NamedMeasurement>
            <NamedMeasurement type="text/string" name="Command Line"><Value>/Users/aiche/dev/openms/openms-src/build/ninja/source/TEST/bin/SomeTool_test</Value></NamedMeasurement>
            <Measurement>
              <Value>
              </Value>
            </Measurement>
    </Results>
  </Test>
  	<EndDateTime>Oct 22 18:43 CEST</EndDateTime>
  	<EndTestTime>1350924239</EndTestTime>
  <ElapsedMinutes>7.2</ElapsedMinutes></Testing>
</Site>
  */
}

print "\nchecker.php finished: ".date("c")."\n";

?>
