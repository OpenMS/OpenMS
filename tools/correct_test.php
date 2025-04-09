<?php
# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright OpenMS Inc. -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-present.
#
# This software is released under a three-clause BSD license:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of any author or any participating institution
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
# For a full list of authors, refer to the file AUTHORS.
# --------------------------------------------------------------------------
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
# INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# --------------------------------------------------------------------------
# $Maintainer: $
# $Authors: Marc Sturm $
# --------------------------------------------------------------------------
error_reporting(E_ERROR| E_WARNING| E_PARSE| E_NOTICE);

include "common_functions.php";

if($argc < 3 || $argc > 4)
{
  print "\n\nUsage: correct_test.php <Absolute path to OpenMS> <Absolute path to header> [-v]\n\n";
  exit;
}

########################## auxilary functions ##################################
/// penalizes differing method names with 100
function penalize_name($f1, $f2) {
  # extract name (between whitepace and bracket
  $tmp = trim(substr($f1, 0, strpos($f1, "(")));
  $tmp = strtr($tmp, array("operator " => "operator", "operator\t" => "operator"));
  $n1  = trim(substr($tmp, max(strrpos($tmp, " "), strrpos($tmp, "\t"))));

  # extract name (between whitepace and bracket
  $tmp = trim(substr($f2, 0, strpos($f2, "(")));
  $tmp = strtr($tmp, array("operator " => "operator", "operator\t" => "operator"));
  $n2  = trim(substr($tmp, max(strrpos($tmp, " "), strrpos($tmp, "\t"))));

  if($n1 == $n2)
  {
    return 0;
  }
  return 100;
}

######################## parameter handling ####################################
$verbose = false;
if(in_array("-v", $argv))
{
  $verbose = true;
}

$path      = $argv[1];
$header    = $argv[2];
$basename  = basename($header);
$test_name = "$path/source/TEST/".substr($basename, 0,-2)."_test.cpp";

######################## determine tested methods ##############################
$tmp = parseTestFile($test_name);
$tests = $tmp["tests"];

if($verbose)
{
  print "\n\nTests:\n";
  foreach($tests as $t)
  {
    print "  '$t'\n";
  }
}

######################## determine declared methods ##########################
$class_info = getClassInfo($path, $header, 0);
$methods = $class_info["public-long"];

// print methods in verbose mode
if($verbose)
{
  print "\n\Methods to correct:\n";
  foreach($methods as $m)
  {
    print "  '$m'\n";
  }
}


######################deter#################
if(count($tests) == 0 || count($methods) == 0)
{
  print "Nothing to do (no tests or methods)\n";
  exit;
}

$replace_strings = array(
  "\t"       => "",
  " "        => "",
  "std::"    => "",
  "OpenMS::" => "",
);

$dists = array();
for($i = 0;$i < count($tests);++$i)
{
  for($j = 0;
  $j < count($methods);
  ++$j)
  {
    //print "\nTest: ".strtr($tests[$i],$replace_strings)."\nMethod: ".strtr($methods[$j],$replace_strings)."\n";
    $penalty = levenshtein(strtr($tests[$i], $replace_strings), strtr($methods[$j], $replace_strings), 1, 10, 10);
    //		if ($penalty!=0)
    //		{
    $penalty += penalize_name($tests[$i], $methods[$j]);
    //		}
    $dists[$i][$j] = $penalty;
    //print "-- '$tests[$i]'\n";
    //print "-- '$methods[$j]'\n";
    //print "-------> '".$dists[$i][$j]."'\n";
    //print "done\n";
  }
}

$fp = fopen("php://stdin", "r");

for($i = 0;$i < count($tests);++$i)
{
  $array = $dists[$i];
  asort($array);
  // abort if exact match
  if(current($array) == 0)
  {
    $replace[] = $tests[$i];
    continue;
  }

  print "\n\nTest:     ".$tests[$i]."\n\n";
  $j = 0;
  foreach($array as $index => $score)
  {
    print "$j) ".str_pad($score, 4, " ", STR_PAD_LEFT)." - ".$methods[$index]."\n";

    //abort after 10
    ++$j;
    if($j == 10)
    {
      break;
    }
  }
  print "\n[enter]  => 0\n";
  print "[i]      => ignore this test\n";
  print "[x]      => make [EXTRA] test (is ignored by checker.php)\n";
  print "[CTRL+C] => abort\n";
  @ob_flush();
  flush();

  //read in choise
  do
  {
    $line = trim(fgets($fp));
  }
  while($line != "" AND !ereg("^[0-9]$", $line) AND $line != "i" AND $line != "x");

  if($line == "i")
  {
    $replace[] = $tests[$i];
  }
  elseif($line == "x")
  {
    $replace[] = "[EXTRA]".$tests[$i];
  }
  else
  {
    if($line == "")
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
$fp = fopen($test_name, "w");
$i = 0;
foreach(file($test_name.".bak") as $line)
{
  if(substr(trim($line), 0, 13) == "START_SECTION" && strpos($line, "[EXTRA]") === FALSE)
  {
    fwrite($fp, "START_SECTION((".$replace[$i]."))\n");
    ++$i;
  }
  else
  {
    fwrite($fp, "$line");
  }
}

fclose($fp);

?>
