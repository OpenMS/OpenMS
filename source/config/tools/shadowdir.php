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
# $Id: shadowdir.php,v 1.1 2006/04/05 07:26:45 marc_sturm Exp $
# $Author: marc_sturm $
# $Maintainer: Marc Sturm $
# --------------------------------------------------------------------------

include "common_functions.php";

	########################usage###############################
	if ($argc!=3)
	{
?>

Usage: shadowdir.php <source directory> <shadow directory> 

Shadows the source OpenMS directory in the shadow directory recursively.
I.e. the directroy structure is copied and a symbolic link from the 
source directory to the shadow directroy is created for each file. 

<?
	exit;	
	}
	
	//rename arguments
	$from = $argv[1];
	$to = $argv[2];
	
	//copy directory structure
	$files=array();
	exec("find $from -type d", $dirs);
	foreach ($dirs as $dir)
	{
		passthru("mkdir -p $to/".substr($dir,strlen($from))."\n",$return);
		if ($return!=0)
		{
			print "Error: could not create directory!\n";
			exit;
		}
	}

	//copy files
	$files=array();
	exec("find $from -not -type d", $files);
	foreach ($files as $file)
	{
		if (!(
			endsWith($file,".so") OR 
			endsWith($file,".o") OR
			endsWith($file,".objects") OR
			endsWith($file,"common.mak") OR 
			endsWith($file,"config.mak") OR 
			endsWith($file,"config_defs.mak") OR 
			endsWith($file,"configure") OR 
			endsWith($file,"source/Makefile") OR 
			endsWith($file,"autom4te.cache") 
			))
		{
			passthru("ln -f $file $to/".substr($file,strlen($from))."\n",$return);
			if ($return!=0)
			{
				print "Error: could not create link!\n";
				exit;
			}
		}
	}

?>
