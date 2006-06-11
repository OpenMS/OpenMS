// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: SequestOutfile.C,v 1.3 2006/06/10 06:40:18 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/SequestOutfile.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

#include <cassert>
#include <string>
using namespace std;

namespace OpenMS {

  SequestOutfile::SequestOutfile(const std::string& filename)
    : filename_(filename) ,phfinished_(false), phstartingpos_(0)
  {
  }

  SequestOutfile::SequestOutfile(const SequestOutfile& source)
    : filename_(source.filename_) , phfinished_(source.phfinished_), phstartingpos_(source.phstartingpos_)
  {
  }

  SequestOutfile::~SequestOutfile()
  {
  }
  
  SequestOutfile& SequestOutfile::operator>>(Identification& identification)
  {
  	DateTime date;
    std::ifstream is(filename_.c_str());
    	
    if (!is)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename_);
    }
    identification.clear();
    
    std::string str;
    for (int i = 0; i < 6; i++ ) 
    {
      std::getline(is,str,'\n');
    }
    // corrupt outfile
    if (is.eof() || str.size() < 11 )
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename_);
    }
//    identification.id_date().month = atoi(&str[1]);
//    identification.id_date().day = atoi(&str[4]);
//    identification.id_date().year = atoi(&str[7]);
    
    char datetemp[11];
    datetemp[0] = str[7];
    datetemp[1] = str[8];
    datetemp[2] = str[9];
    datetemp[3] = str[10];
    datetemp[4] = 0;
    datetemp[5] = str[1];
    datetemp[6] = str[2];
    datetemp[7] = 0;
    datetemp[8] = str[4];
    datetemp[9] = str[5];
    datetemp[10] = 0;
    
    // if date is not readable, use default
    if (!(atoi(&datetemp[0]) && atoi(&datetemp[5]) && atoi(&datetemp[8]) ) )
    {
      for(uint i = 0; i < 10; ++i)
      {
        datetemp[i] = '0';
      }
    }
    datetemp[4] = datetemp[7] = 0;
		
		date.set(datetemp);
		
    identification.setDateTime(date);
    
    PeptideHit temphit= PeptideHit();
    *this >> temphit;
    while (temphit.getScore()) 
    {
     if ( temphit.getSequence() != "" ) 
     {
       identification.insertPeptideHit(temphit);
     }
      *this >> temphit;
    }

    is.close();
    return *this;
  }
  
  SequestOutfile& SequestOutfile::operator>>(PeptideHit& peptidehit)
  {
    if (phfinished_) {
      peptidehit.setScore(0);
      return *this;
    }
    peptidehit.clear();
    std::ifstream is(filename_.c_str());
    if (!is)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename_);
    }
    
    std::string str;
    if (!phstartingpos_) {
      for (int i = 0; i < 15; i++ ) 
      {
       std::getline(is,str,'\n');
      }
    }
    else is.seekg(phstartingpos_);
    char* sep = " ";
    char* token ;
    //string token = "";
    int ignorelines = 0;
    uint rank = 0;
    uint sprank = 0;
    double XCorr, deltaCn, mh, Sp;
    std::string ions, reference, peptide;
    std::getline(is,str,'\n');
    while (!containlength_(str.c_str())) 
    {
      if ( is.eof() ) {
        peptidehit.setScore(0);
        phfinished_ = true;
        is.close();
        return *this;
      }
      std::getline(is,str,'\n');
    }
    char line[400];
    strcpy(line,str.c_str());
    //token = stringtoken(str,sep);
    token = strtok(line,sep);
    int i = 0;
    while ((token = strtok(0,sep))) 
    //while ((token = stringtoken(str,sep)).size()) 
    {
      switch (i) 
      {
        case 0:
        if (!atoi(token) && length_(str.c_str()) ) {
          peptidehit.setScore(0);
          phfinished_ = true;
          is.close();
          return *this;
        }
        rank = atoi(token);
        i++;
        break;
        case 1:
        if (*token == '/')
        {
          if ( length_(token) > 1 )
          {
            sprank = atoi(token+1);
            ++i;
          }
        }
        else 
        {
          sprank = atoi(token); 
          ++i;
        }
        break;
        case 2:
        mh = atof(token);
        i++;
        break;
        case 3:
        deltaCn = atof(token);
        i++;
        break;
        case 4:
        XCorr = atof(token);
        i++;
        break;
        case 5:
        Sp = atof(token);
        i++;
        break;
        case 6:
        ions = token;
        i++;
        break;
        case 7:
        reference = token;
        replaceapostrophe_(reference);
        i++;
        break;
        case 8:
        if (!ignorelines ) ignorelines = atoi(token); 
        peptide = token;
        break;
        default:
        cerr << "input too long " << std::endl;
        break;
      } 
    }
    peptidehit.setScore(XCorr);
    peptidehit.setSequence(peptide);
    peptidehit.setRank(rank);
    for (int i = 0; i < ignorelines; i++) {
     std::getline(is,str,'\n');
     if ( is.eof() ) {
        phfinished_ = true;
      }
    }
    phstartingpos_ = is.tellg();
    is.close();
    return *this;
  }
  
  int SequestOutfile::length_(const char* charp)
  {
    int i = 0;
    while (charp[i]) i++;
    return i;
  }
  
  int SequestOutfile::containlength_(const char* charp) 
  {
    int i = 0;
    for (int j = 0; j < length_(charp); j++ )
    {
      if ( charp[j]!=' ' ) i++;
    }
    return i;
  }

  void SequestOutfile::replaceapostrophe_(string& str)
  {
    for(uint i = 0; i < str.size(); ++i)
    {
      if (str[i] == '\'')
      {
        str.insert(i,"\\");
        ++i;
      }
    }
  }
  
} //namespace OpenMS
