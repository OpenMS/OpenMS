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

#include <OpenMS/COMPARISON/CLUSTERING/helper.h>

namespace OpenMS
{

  int indentbuf::overflow(int c)
  {
    if (traits_type::eq_int_type(c, traits_type::eof()))
      return m_sbuf->sputc(c);

    if (m_need)
    {
      std::fill_n(std::ostreambuf_iterator<char>(m_sbuf), m_indent, ' ');
      m_need = false;
    }
    if (traits_type::eq_int_type(m_sbuf->sputc(c),
          traits_type::eof()))
      return traits_type::eof();
    if (traits_type::eq_int_type(c, traits_type::to_char_type('\n')))
      m_need = true;
    return traits_type::not_eof(c);
  }

  static int const index = std::ios_base::xalloc();
  
  void callback(std::ios_base::event ev, std::ios_base& ios, int) {
    if (ev == std::ios_base::erase_event) {
      std::ostream* out = static_cast<std::ostream*>(ios.pword(index));
      assert(out != 0);
      indentbuf* sbuf = dynamic_cast<indentbuf*>(out->rdbuf());
      assert(sbuf != 0);
      out->rdbuf(sbuf->sbuf());

      delete sbuf;
    }
  }
  
  std::ostream& operator<< (std::ostream& out, indent const& ind) {
    indentbuf* sbuf = dynamic_cast<indentbuf*>(out.rdbuf());
    if (!sbuf)
    {
      sbuf = new indentbuf(out.rdbuf());
      out.pword(index) = &out;
      out.register_callback(callback, 0);
      out.rdbuf(sbuf);
    }

   sbuf->indent(ind.m_indent);
   return out;
  }
}
