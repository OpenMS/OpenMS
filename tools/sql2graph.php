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
# $Maintainer:$
# $Authors: Marc Sturm $
# --------------------------------------------------------------------------
////////////////////////////////METHODS/////////////////////////////////////
include("common_functions.php");

function createUMLNode($name, $attributes, $methods) {
  global $out;
  
  //make sure arrays and attributes are arrays
  if(!is_array($attributes)) 
    $attributes = array(
      $attributes,
    );
  if(!is_array($methods)) 
    $methods = array(
      $methods,
    );
  
  //calculate width and hight
  $height  = 35+14*(count($attributes)+count($methods));
  $max_att = array_reduce(array_map("strlen", $attributes), "max", 0);
  $width   = 10+max(8*strlen($name), 7*$max_att);
  
  //begin
  $out[] = "    <node id=\"$name\">\n";
  $out[] = "      <data key=\"d0\" >\n";
  $out[] = "        <y:UMLClassNode >\n";
  $out[] = "          <y:Geometry width=\"$width\" height=\"$height\"/>\n";
  $out[] = "          <y:Fill hasColor=\"false\" />\n";
  $out[] = "          <y:NodeLabel fontStyle=\"bold\">$name</y:NodeLabel>\n";
  $out[] = "          <y:UML use3DEffect=\"false\">\n";
  //attributes
  $out[] = "            <y:AttributeLabel>";
  foreach($attributes as $a)
  {
    $out[] = "$a\n";
  }
  $out[] = "            </y:AttributeLabel>\n";
  //methods
  $out[] = "            <y:MethodLabel>";
  foreach($methods as $m)
  {
    $out[] = "$m\n";
  }
  $out[] = "            </y:MethodLabel>\n";
  //end
  $out[] = "          </y:UML>\n";
  $out[] = "        </y:UMLClassNode>\n";
  $out[] = "      </data>\n";
  $out[] = "    </node>\n";
}

function createUMLEdge($from, $to) {
  global $out;
  
  static $count = 0;
  $out[] = "    <edge id=\"e$count\" source=\"$from\" target=\"$to\" />\n";
  ++$count;
}

function rp($string) {
  return strtr($string, array(")" => "", "(" => "", " " => "", "	" => "", "," => ""));
}

function parseSQLFile($filename, &$tables, &$edges) {
  //init
  $tables = array();
  $edges = array();
  
  //parse file
  $file = file($filename);
  for($i = 0;$i < count($file);++$i)
  {
    $line = trim($file[$i]);
    //node
    if(beginsWith($line, "CREATE TABLE"))
    {
      $parts                  = explode(" ", $line);
      $current_table          = $parts[2];
      $tables[$current_table] = array();
      $j                      = $i+1;
      $line                   = trim($file[$j]);
      //table
      while(!beginsWith($line, "CREATE TABLE") && !beginsWith($line, "ALTER TABLE") && $line != "")
      {
        //primary key
        if(beginsWith($line, "PRIMARY KEY"))
        {
          $parts = explode(" ", $line);
          $tables[$current_table][rp($parts[3])][] = "P";
        }
        // index
        elseif(beginsWith($line, "KEY"))
        {
          $parts = explode(" ", $line);
          $tables[$current_table][rp($parts[2])][] = "I";
        }
        // normal field
        else
        {
          $parts = explode(" ", $line);
          if($parts[0] != ")") 
            $tables[$current_table][$parts[0]] = array();
        }
        ++$j;
        $line = trim($file[$j]);
      }
    }
    //edges
    elseif(beginsWith($line, "ALTER TABLE"))
    {
      $parts = explode("`", $line);
      $from  = $parts[1];
      $j     = $i+1;
      while(beginsWith(trim($file[$j]), "ADD CONSTRAINT"))
      {
        $parts                          = explode(" ", $file[$j]);
        $to                             = $parts[9];
        $tables[$from][rp($parts[7])][] = "F";
        $edges[] = array(
          $from,
          $to,
        );
        ++$j;
      }
    }
  }
}

////////////////////////////////MAIN/////////////////////////////////////
if($argc != 3)
{
  print "Usage: sql2graph.php <input> <output>\n\n";
  exit;
}

//begin
$out[] = "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns/graphml\" xmlns:y=\"http://yworks.com/xml/graphml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns/graphml http://www.yworks.com/xml/schema/graphml/1.0/ygraphml.xsd\">\n";
$out[] = "<key id=\"d0\" for=\"node\" yfiles.type=\"nodegraphics\"/>\n";
$out[] = "<graph id=\"G\" edgedefault=\"undirected\">\n";

//parsing and graph creation
parseSQLFile($argv[1], $tables, $edges);

foreach($tables as $t => $fields)
{
  $text = array();
  foreach($fields as $field => $prop)
  {
    if(in_array("P", $prop))
    {
      $field .= "**";
    }
    elseif(in_array("I", $prop))
    {
      $field .= "*";
    }
    if(in_array("F", $prop))
    {
      $field .= " (FK)";
    }
    $text[] = $field;
  }
  createUMLNode($t, $text, "");
}
foreach($edges as $e)
{
  if($e[1] != "META_MetaInfo" && $e[1] != "META_File") 
    createUMLEdge($e[0], $e[1]);
}

//end
$out[] = "</graph>\n";
$out[] = "</graphml>\n";

//write file
$fp = fopen($argv[2], "w");
foreach($out as $line)
{
  fwrite($fp, $line);
}
fclose($fp);

?>
