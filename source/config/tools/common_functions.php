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

//write array to file line by line
function writeFile($name,$array)
{
	$fp = fopen($name,"w");
	foreach ($array as $line){
		fwrite($fp,rtrim($line,"\r\n")."\n");
		}
	fclose($fp);
}

// returns true if $string begins with $begin, return false otherwise
function beginsWith($string,$begin)
{
	if (substr($string,0,strlen($begin)) == $begin )
	{
		return true;
	}
	return false;
}

// returns true if $string ends with $end, return false otherwise
function endsWith($string,$end)
{
	if (substr($string,-1*strlen($end)) == $end )
	{
		return true;
	}
	return false;
}

// returns the prefix of length $length of string $string
function prefix($string,$length)
{
	return substr($string,$length);
}

// returns the suffix of length $length of string $string
function suffix($string,$length)
{
	return substr($string,-1*$length);
}

// returns the include guard for the included file $include
// give only the path and filename i.e. the part between "#include <" and ">"
function includeToGuard($include)
{
	return strtoupper(str_replace(".","_",str_replace("/","_",$include)));
}

//return true if the given line is a include line, returns false otherwise
//the fielname and path of the include are written into $include if it is given 
function isIncludeLine($line,&$include)
{
	if (ereg ("^#[ \t]*include[ \t]*<(.*)>", ltrim($line),$parts))
	{
		$include = $parts[1];
		return true;
	}
	return false;
}

function isAlpha($char)
{
	$c = ord($char);
	return ((($c >= 64) && ($c <= 90)) || (($c >= 97) && ($c <= 122)));
}

function isSpace($char)
{
	return ($char == " ") || ($char == "\n") || ($char == "\t");
}

function isDigit($char)
{
	return ('0' <= $char) && ($char <= '9');
}

function isAlnum($char)
{
	return isAlpha($char) || isDigit($char);
}

function isIdentifierChar($char)
{
	return isAlnum($char) || $char == '_';
}

// tokenize a line into something like C++-Tokens (not exactly, but good enough)
function tokenize($line)
{
	$result = array();
	$start = 0;
	// states: start => whitespace
	//         ident => identifier / number etc.
	$state = 0;
	for ($i = 0; $i != strlen($line); $i++)
	{
		$char = $line[$i];
		if ($state == 0)
		{
			// start state
			if (isIdentifierChar($char))
			{
				// identifier / number etc.
				$start = $i;
				$state = 1;
			}
			else if (!isSpace($char))
			{
				// operator
				$result[] = $char;
			}
		}
		else
		{
			// identifier state
			if (!isIdentifierChar($char))
			{
				// end of identifier
				$state = 0;
				$result[] = substr($line, $start, $i - $start);
				
				if (!isSpace($char))
				{
					// additionally, the current char belongs to an operator
					$result[] = $char;
				}
			}
		}
	}
	
/*	print("tokenize: ");
	
	for ($i = 0; $i != count($result); $i++) {
		print("\"" . $result[$i] . "\" ");
	}
	print("\n");*/
	return $result;
}

?>