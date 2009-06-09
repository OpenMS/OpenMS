// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                    AutoExecute file support
// --------------------------------------------------------------------------
//  Copyright (C) 2009 -- Guillaume Belz (guillaume.belz@chu-lyon.fr)
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
// $Maintainer: Guillaume Belz
// $Authors: Guillaume Belz
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/AutoExecuteFile.h>

using namespace std;

namespace OpenMS
{

	AutoExecuteFile::AutoExecuteFile()
	{
    PosOnScout_ = 0;
    SpectrumDirectory_ = 0;
    SpectrumFilename_ = 0;
    ChipOnScout_ = 0;
	}
	
	
	StringList AutoExecuteFile::getFileList(
	      const String & filename,  
	      const bool isAutoExecute,
	      const unsigned int begin,
	      const unsigned int end,
	      const String src_dir)
	{
    std::ifstream is(filename.c_str());
    if (!is)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
    
    //temporary variables
    StringList list;
    String line;

    if(isAutoExecute)
    {
cout << "isAutoExecute" << endl;    
      // read headers
      readAutoExecuteHeader(is);
    }
cout << "readAutoExecuteHeader ok" << endl;
    // skip first lines
    unsigned int index = 0;
    for(; index<begin; index++)
    {
cout << "1 line skip" << endl;
      getline(is, line, '\n');
    }
    
    //read lines
    while(getline(is, line, '\n') && index<end)
    {
cout << "line " << index << " : " << line << endl;    
      if(isAutoExecute)
      {
        String file = autoExecuteToFilename(line, src_dir);
        list.push_back(file);
      }
      else
        list.push_back(line);
      
      index++; 
    }
    
    return list;
	}
	
	
	String AutoExecuteFile::autoExecuteToFilename(const String & line, const String & src_dir)
	{
	  vector<String> cells;
	  String filename;
	  
	  line.split('\t', cells);
	  
	  if(src_dir == String(""))
	    filename = cells[SpectrumDirectory_];
	  else
	    filename = src_dir;
	    
	  filename += cells[SpectrumFilename_] + String("/");
	  
	  // '1' read only first spectrum
	  // '1SLin' for TOF linear mode
	  filename += cells[ChipOnScout_] + String("_") + cells[PosOnScout_].remove(':') + String("/1/1SLin/fid");
	  
	  cells.clear();
cout << "autoExecuteToFilename : " << filename << endl; 	  
	  return filename;
	}
	
	
	void AutoExecuteFile::readAutoExecuteHeader(std::ifstream & is)
	{
	  String line;
	  vector<String> headers;
	  
	  // read headers
	  getline(is, line, '\n');
	  line.split('\t', headers);
    for(unsigned int index=0; index<headers.size(); index++)
    {
      if(headers[index] == String("Pos_on_Scout"))
        PosOnScout_ = index;
      else if(headers[index] == String("Spectrum_Directory"))
        SpectrumDirectory_ = index;
      else if(headers[index] == String("Spectrum_Filename"))
        SpectrumFilename_ = index;
      else if(headers[index] == String("Chip_on_Scout"))
        ChipOnScout_ = index;
    }
    headers.clear();
cout << "headers : " << line << endl;
cout << "PosOnScout : " << PosOnScout_ << " ; SpectrumDirectory : " << SpectrumDirectory_ << " ; SpectrumFilename : " << SpectrumFilename_ << " ; ChipOnScout : " << ChipOnScout_ << endl;
	  
	  // read '************' line
	  getline(is, line, '\n');
	}
} // namespace OpenMS

