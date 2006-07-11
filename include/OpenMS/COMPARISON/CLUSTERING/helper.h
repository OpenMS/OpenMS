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
// $Maintainer:  $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_CLUSTERING_HELPER_H
#define OPENMS_COMPARISON_CLUSTERING_HELPER_H

#include <streambuf>
#include <algorithm>
#include <iterator>
#include <ostream>
#include <cassert>



#define DEBUG_ cout << __LINE__ << endl;

//some helpful templates like maximum, minimum and absolute
//and something to automatically indent output (used for writing the xml 
//in ClusterExperiment)

namespace OpenMS
{

/*  template<typename T>
  int numdigits(T a,int base = 10)
  {
    int res = 0;
    while  ( a >= 1) {
      res++;
      a = a/base;
    }
    return res;
  }
*/

  //to make nice xml files without DOM Level 3
  struct indentbuf: public std::streambuf { 
    indentbuf(std::streambuf* sbuf):
      m_sbuf(sbuf), m_indent(0), m_need(true) {}
    
    int indent() const { return m_indent; }
    void indent(int i) { m_indent = 2*i; } //indentation in steps of 2, while using ++/--
    std::streambuf* sbuf() const { return m_sbuf; }

  private:
    int overflow(int c) ;

    std::streambuf* m_sbuf;
    int             m_indent;
    bool            m_need;
  };

  struct indent {
    indent(int i): m_indent(i) {}
    int m_indent;
  };
    
  std::ostream& operator<< (std::ostream& out, indent const& ind) ;

}
#endif //OPENMS_COMPARISON_CLUSTERING_HELPER_H
