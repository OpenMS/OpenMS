<?
# -*- Mode: C++; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

# This script calls svn status and rearranges the output to make it human readable
#
# Call the following command in the root directory of an OpenMS SVN working copy:
# > php -qC source/config/tools/svn_status.php

	error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);


	//different types of files
	$patched = array();
	$modified = array();
	$conflict = array();
	$unknown = array();
	$added = array();
	$removed = array();
	

	// array with all variables that store different types of files (for output)
	$vars = array("unknown","patched","modified","added","removed","conflict");

	$templates_to_ignore=array(
	"^*.log$",
	"^*_test$",
	"^*_moc.C$",
	"^*.Dependencies$",
	"^include/OpenMS/VISUAL/UIC*.h$",
	"^include/OpenMS/VISUAL/DIALOGS/UIC*.h$",
	"^source/VISUAL/DIALOGS/*Template.C$",
	"^source/VISUAL/*Template.C$",
	"^source/APPLICATIONS/TOPP/!$",
	"^source/TEST/feature*$",
	"^source/TEST/*.vgr$",
	"^source/TEST/*.chk$",
	"^source/TEST/*.output$",
	"^source/TEST/model*$",
	"^source/TEST/TOPP/feature*$",
	"^source/TEST/TOPP/model*$",
	"^source/TEST/TOPP/region*$",
	"^source/TEST/TOPP/*.tmp$",
	);

	$templates_to_ignore = str_replace("*","[/.//A-Z0-9a-z_-]*",$templates_to_ignore);
	$templates_to_ignore = str_replace("!","[A-Z0-9a-z_-]*",$templates_to_ignore);
	
	$files_to_ignore=array(
		"bin",
		"lib",
		"doc/html",
		"doc/Doxyfile",
		"doc/misc/coding_style.pdf",
		"include/OpenMS/CONFIG",
		"include/OpenMS/config.h",
		"source/Makefile",
		"source/common.mak",
		"source/config.lic",
		"source/config.mak",
		"source/config.status",
		"source/config_defs.mak",
		"source/libOpenMS.objects",
		"source/libOpenMS_GUI.objects",
		"source/autom4te.cache",
		"source/conf.diag.tar",
		"source/config.h",
		"source/TEST/TOPP/NumericDiff",
		"source/TEST/DB_credentials.txt",
		);
	
	exec("svn status",$lines);
	print "\n";
	foreach ($lines as $line)
	{
		$line = trim($line);
		// P
		if (substr($line,0,2)=="P ") $patched[]=substr($line,2);
		// U
		else if (substr($line,0,2)=="U ") $patched[]=substr($line,2);
		// M
		else if (substr($line,0,2)=="M " OR substr($line,0,2)=="MM")  $modified[]=substr($line,2);
		// C
		else if (substr($line,0,2)=="C ") $conflict[]=substr($line,2);
		// A
		else if (substr($line,0,2)=="A ") $added[]=substr($line,2);
		// D
		else if (substr($line,0,2)=="D ") $removed[]=substr($line,2);
		// ?
		else if (substr($line,0,2)=="? ")
		{
			//files to ignore
			if (!in_array(trim(substr($line,2)),$files_to_ignore))
			{
				foreach ($templates_to_ignore as $template)
				{
					if (ereg($template,trim(substr($line,2))))
					{
						//print substr($line,2)." --> ".$template."\n";
						continue 2;
					}
				}
				$unknown[]=substr($line,2);
			}
		}
		else
		{
			if ($line!="") print $line;
		}
	}

	foreach ($vars as $var)
	{
		if (count($$var)!=0)
		{
			print "\n>>> ".strtoupper($var)." <<<\n";
			foreach ($$var as $u)
			{
				print trim($u)."\n";
			}
		}
	}

	print "\n";

?>
