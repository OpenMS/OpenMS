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

function addBody(&$array)
{
	$array[] = "	{";
	$array[] = "	  //????";
	$array[] = "	}";
	$array[] = "";
}

	########################usage###############################
	if ($argc<2 OR $argc>3)
	{
?>
Usage: create_methods.php <header> [-v]

Create source file and adds member access method declarations to the header file.

options:
  -v verbose output

<?
	exit;	
	}
	
	######################command line parsing###########################
	$header = $argv[1];

	if ($argv[2] == "-v")
	{
		$verbose = true;
	}
	else
	{
		$verbose = false;
	}

	##############################variables###############################
	$vars = array();
	
	$h_out = array();
	$s_out = array();
		
	//types that are passed by value (enums are added after parsing)
	$base_simple_types = array("bool","float","double","SignedInt","UnsignedInt","Size","Index");
	//types that do not need a non-const get-method
	$normal_types = array("String");
	
	$namespaces = array();
	$enums = array();
	
	##############################parse input#############################
	if ($verbose) print ">> $header\n";

	$file = file($header);
	
	### parse inlude file path determine source file path
	$pos = strpos($header,"include/OpenMS");
	if ($pos!==FALSE)
	{
		$class_path = substr($header,$pos+15);
		$source_file_path = substr($header,0,$pos)."source/".substr($class_path,0,-2).".C";
		if ($verbose)
		{
			print "  class path: $class_path\n";
			print "  source file: $source_file_path\n";
		}
	}
	else
	{
		print "Error: start this script from outside the include folder of OpenMS!\n";
		print "       The header file path must contain 'include/OpenMS'!\n";
		exit();
	}
	
	### header lines
	$i=0;
	while(beginsWith(trim($file[$i]),"//"))
	{
		$s_out[] = trim($file[$i]);
		++$i;
	}
	if ($verbose) print "  header lines: ".count($s_out)."\n";
	
	## inlcude in header
	$s_out[] = "";
	$s_out[] = "#include <OpenMS/$class_path>";
	$s_out[] = "";
  $s_out[] = "using namespace std;";
  $s_out[] = "";
  $s_out[] = "namespace OpenMS";
  $s_out[] = "{";
  
	for (;$i<count($file);$i++)
	{
		$line = trim($file[$i]);
		### namespaces
		if(beginsWith($line,"namespace") AND ereg("^[namespac]{9}[ 	]+([A-Za-z0-9]+)",$line,$parts) AND $parts[1]!="OpenMS")
		{
			$namespaces[] = $parts[1];
		}		
		
		
		### class name
		if(beginsWith($line,"class") AND ereg("^[clas]{5}[ 	]+([A-Za-z0-9_]+)",$line,$parts))
		{
			$class = $parts[1];
			if ($verbose) print "  class: $class\n";
		}
		
		### public
		if (beginsWith($line,"public:"))
		{
			$public_line = $i;
			if ($verbose) print "  public line: $public_line\n";
		}
		### public
		if (beginsWith($line,"enum") AND ereg("^[enum]{4}[ 	]+([A-Za-z0-9_]+)",$line,$parts))
		{
			$enums[] = $parts[1];
			$public_line = $i;
			if ($verbose) print "  enum line: $public_line ($parts[1])\n";
		}
				
		###member definitions
		if(ereg("^(.*)[ 	]+([a-zA-Z0-9_]+);$",$line,$parts))
		{
			
			//extract type and name
			$type = $parts[1];
			$name = $parts[2];
			
			//derive method name from name
			$method="";
			$parts = explode("_",$name);
			foreach ($parts as $part)
			{
				if (in_array($part,array("db")))
				{
					$method .= strtoupper($part);
				}
				else
				{
					$method .= ucfirst($part);
				}
			}
			if ($verbose) print "  definition: $name ($type) ($method)\n"; 
			$vars[]=array("n"=>$name,"t"=>$type,"m"=>$method);
		}
		
	}
	
	###namespaces
	$namespaces = implode($namespaces,"::");
	if ($verbose) print "  namespaces: $namespaces";
	
	### add constructors, destroctor and assignment operator to header
  $i_out[] = "      $class();";
  $i_out[] = "      $class(const $class& source);";
  $i_out[] = "      ~$class();";
  $i_out[] = "      $class& operator = (const $class& source);";
  $i_out[] = "      bool operator == (const $class& source) const;";
  $i_out[] = "      bool operator != (const $class& source) const;";
  $i_out[] = "";	
  
  ### add constructors, destroctor and assignment operator to source
  $s_out[] = "	$class::$class() //????:";
  addBody($s_out);
  $s_out[] = "	$class::$class(const $class& source):";
  foreach ($vars as $v)
	{
		$s_out[] = "	  ".$v["n"]."(source.".$v["n"]."),";  
	}
	$s_out[count($s_out)-1] = substr($s_out[count($s_out)-1],0,-1);
  addBody($s_out);
  $s_out[] = "	$class::~$class()";
  addBody($s_out);
  $s_out[] = "	$class& $class::operator = (const $class& source)";  
	$s_out[] = "	{";
	$s_out[] = "	  if (&source != this)";
	$s_out[] = " 	 {";
	foreach ($vars as $v)
	{
		$s_out[] = " 	   ".$v["n"]." = source.".$v["n"].";";  
	}
	$s_out[] = "	  }";
	$s_out[] = "	  return *this;";
	$s_out[] = "	}";
	$s_out[] = "	";

  $s_out[] = "	bool $class::operator == (const $class& source) const";  
	$s_out[] = "	{";
	$s_out[] = "	  return";
	for ($i=0; $i< count($vars)-1; ++$i)
	{
		$s_out[] = " 	   ".$vars[$i]["n"]." == source.".$vars[$i]["n"]." &&";  
	}
	$s_out[] = " 	   ".$vars[count($vars)-1]["n"]." == source.".$vars[count($vars)-1]["n"]."";  
	$s_out[] = "	  ;";
	$s_out[] = "	}";
	$s_out[] = "	";

  $s_out[] = "	bool $class::operator != (const $class& source) const";  
	$s_out[] = "	{";
	$s_out[] = "	  return ! operator==(source);";
	$s_out[] = "	}";
	$s_out[] = "	";
 
 	$simple_types = array_merge($base_simple_types, $enums);
 	
	### add member accessors to headers
	// header output
	foreach ($vars as $v)
	{
		$name = $v["n"];
		$type = $v["t"];
		$method = $v["m"];
		
		//get method
		if (in_array($type,$simple_types))
		{
			$i_out[] = "      $type get$method() const;";
		}
		else
		{
			$i_out[] = "      const $type& get$method() const;";
			// add non-const get method of complex type
			if (! in_array($type,$normal_types))
			{
				$i_out[] = "      $type& get$method();";
			}
		}
		
		//set method
		if (in_array($type,$simple_types))
		{
			$i_out[] = "      void set$method($type ".substr($name,0,-1).");";
		}
		else
		{
			$i_out[] = "      void set$method(const $type& ".substr($name,0,-1).");";
		}
		
		$i_out[] = "";		
	}	
	
	### add member accessors to source file
	foreach ($vars as $v)
	{
		$name = $v["n"];
		$type = $v["t"];
		$method = $v["m"];

		
		if (in_array($type,$enums))
		{
			$real_type = "	$class::$type";
		}
		else
		{
			$real_type = $type;	
		}
				
		//get method
		if (in_array($type,$simple_types))
		{
			$s_out[] = "	$real_type $class::get$method() const ";
		}
		else
		{
			$s_out[] = "	const $real_type& $class::get$method() const ";
			// add non-const get method of complex type
			if (! in_array($type,$normal_types))
			{
				$s_out[] = "	{";
				$s_out[] = "	  return $name; ";
				$s_out[] = "	}";
				$s_out[] = "	";
				$s_out[] = "	$real_type&  $class::get$method()";
			}
		}
		
		$s_out[] = "	{";
		$s_out[] = "	  return $name; ";
		$s_out[] = "	}";
		
		$s_out[] = "	";
		
		//set method
		if (in_array($type,$simple_types))
		{
			$s_out[] = "	void $class::set$method($real_type ".substr($name,0,-1).")";
		}
		else
		{
			$s_out[] = "	void $class::set$method(const $real_type& ".substr($name,0,-1).")";
		}
	
		$s_out[] = "	{";
		$s_out[] = "	  $name = ".substr($name,0,-1)."; ";
		$s_out[] = "	}";
		
		$s_out[] = "	";		
	}
	
	$s_out[] = "} //namespace";
	$s_out[] = "";
	
	### write files
	array_splice($file,$public_line+1,0,$i_out);
	writeFile($header,$file);
	writeFile($source_file_path,$s_out);
?>
