// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------


/// Small tool to update the svn revision number in a tiny header file.
/// The header file is only modified if the revision number within the header
/// and the one from a call to svnversion(.exe) actually differs


#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

/// print usage of this tool
string usage(const string toolname)
{
	return string("Usage: \n  ") + toolname + string(" <project source dir> <header file>\n");
}

/// grab output of given command
/// returns 1 on error, 0 on success
int getStdoutFromCommand(string cmd, string& data)
{
	data = "";
	FILE *stream;
	const int MAX_BUFFER = 256;
	char buffer[MAX_BUFFER];
	cmd.append(" 2>&1");
#if WIN32
	stream = _popen(cmd.c_str(), "r");
#else
	stream =  popen(cmd.c_str(), "r");
#endif	
	if (!stream) 
	{
		data="svnversion call failed!";
	}
	else
	{
		while (!feof(stream))
		{
			if (fgets(buffer, MAX_BUFFER, stream) != NULL) data.append(buffer);
		}
	}
#if WIN32
	_pclose(stream);
#else
	pclose(stream);
#endif		
	return 0;
}

/// splits a string
void splitString(const string& s, const char splitter, vector<string>& substrings)
{
	string::const_iterator begin = s.begin();		
	string::const_iterator end = s.begin();
	
	for (; end != s.end() && *end != '/'; ++end)
	{
		if (*end == splitter)
		{
			substrings.push_back(string(begin, end));
			begin = end + 1;
		}
	}
	substrings.push_back(string(begin, s.end()));
	
	return;
}

int main( int argc, const char* argv[] )
{
	if (argc!=3) 
	{
		std::cerr << "Error. Invalid number of arguments!\n" << usage(string(argv[0])) << "\n";
		return 1;
	}
	
	string svn_dir         = string(argv[1]);
	string svn_header_file = string(argv[2]);
	
	// use svnversion command to get the current svn revision
	string svn_revision;
	getStdoutFromCommand(string("svnversion \"") + svn_dir + string("\" -n"), svn_revision);
	svn_revision = string("\"") + svn_revision + string("\"");
	//std::cout << "DEBUG - got svn revision: " << svn_revision << endl;
	
	// extract svn revision from header and see if it needs updating
	ifstream hfile (svn_header_file.c_str());
  string line;
  if (hfile.is_open())
  {
    if (! hfile.eof() )
    {
      getline (hfile,line);
      //std::cout << "DEBUG - got header line: " << line << endl;
    }
    hfile.close();
  }
  
  // dissect line
  vector<string> substrings;
  splitString(line, ' ', substrings);
	
	if (substrings.size() != 4)
	{ 
		std::cerr << "Input file " << svn_header_file << " not formatted as expected: got " << substrings.size() << " substrings, expected 4\n";
		for (unsigned int i=0;i<substrings.size();++i)	
		{
			std::cerr << " " << substrings[i] << "\n";
		}		
		return 1;
	}
	
	//std::cout << "DEBUG - comparing: '" << substrings[2] << "' vs '" << svn_revision << "'" << endl;
	// compare the two revisions:
	if (substrings[2] != svn_revision)
	{
		substrings[2] = svn_revision; // replace with new revision
		ofstream hfile;
		hfile.open (svn_header_file.c_str());
		for (unsigned int i=0;i<substrings.size();++i)	
		{
			if (i!=0) hfile << " ";
			hfile << substrings[i];
		}
		hfile << "\n";
		hfile.close();	
	}
	else
	{
		// nothing changed.
		//std::cout << "DEBUG - nothing to be done. Header file is up-to-date." << endl;
	}
	
}
