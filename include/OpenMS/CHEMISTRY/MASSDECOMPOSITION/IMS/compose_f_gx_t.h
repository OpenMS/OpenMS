// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_COMPOSE_F_GX_T_H
#define OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_COMPOSE_F_GX_T_H

namespace OpenMS {

namespace ims {

#include <functional>
/** 
 * Represents the function object unary adapter.
 * This simplest and most fundamental compose function adapter uses the result of 
 * a unary operation as input to another unary operation. For more details, please
 * refer to the book "The C++ Standard Library" by Nicolay Josuttis.
 * This class is implemented in order to reduce (substitute) the usage of user defined functors.
 * You can use it for example to combine several operation in one function call. 
 * Here is an example of combining "minus 5 and multiply by 10": 
 * @c compose_f_gx(bind2nd(multiplies<int>(),10),bind2nd(minus<int>(), 5))
 * 
 * Another example of @class compose_f_gx_t usage is along with @c mem_ref_fun function 
 * to search for elements in container based on the certain element's property.
 * Here is an example:
 * @code
include <string>;
include <vector>;
using namespace std;
class Element {
 public:
  Element(const string& a) : _a(a) {}
  const string& getA() const { return _a; }
 private:
  string _a;
};
int main(int argc, char** argv) {

 Element a("4"), b("3"), c("2"), d("1");

 vector<Element> elements;
 elements.push_back(a);
 elements.push_back(b);
 elements.push_back(c);
 elements.push_back(d);

 vector<Element>::const_iterator it = find_if(elements.begin(), elements.end(),
  compose_f_gx(bind2nd(equal_to<string>(), "3"), mem_fun_ref(&Element::getA)));
	
 if (it != elements.end()) {
  cout << "element 3 has index " << it-elements.begin() << endl;
 } else {
  cout << "element 3 is not in the container." << endl;
 }
 return 0;
}
@endcode
 * The Example has the following output:
@code
element 3 has index 1
@endcode
 *
 *
 * Copyright 1999 by Addison Wesley Longman, Inc. and Nicolai M. Josuttis.
 * All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software for personal and
 * educational use is hereby granted without fee, provided that the above copyright
 * notice appears in all copies and that both that copyright notice and this
 * permission notice appear in supporting documentation, and that the names of Addison
 * Wesley Longman or the author are not used in advertising or publicity pertaining
 * to distribution of the software without specific, written prior permission. Addison
 * Wesley Longman and the author make no representations about the suitability of this
 * software for any purpose. It is provided "as is" without express or implied warranty.
 *
 * ADDISON WESLEY LONGMAN AND THE AUTHOR DISCLAIM ALL WARRANTIES WITH REGARD TO THIS
 * SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT
 * SHALL ADDISON WESLEY LONGMAN OR THE AUTHOR BE LIABLE FOR ANY SPECIAL, INDIRECT OR
 * CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
 * PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING
 * OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 */
template <class OP1, class OP2>
class compose_f_gx_t
    : public std::unary_function<typename OP2::argument_type, typename OP1::result_type>
{
private:
  OP1 op1;    // process: op1(op2(x))
  OP2 op2;
public:
  /// Constructor
  compose_f_gx_t(const OP1& o1, const OP2& o2)
    : op1(o1), op2(o2)
  {
  }

  /// function call
  typename OP1::result_type
  operator()(const typename OP2::argument_type& x) const
  {
    return op1(op2(x));
  }
};

/**
 * @brief convenience function for the compose_f_gx adapter
 */
template <class OP1, class OP2>
inline compose_f_gx_t<OP1,OP2>
compose_f_gx (const OP1& o1, const OP2& o2)
{
  return compose_f_gx_t<OP1,OP2>(o1,o2);
}

} // namespace ims
} // namespace OpenMS

#endif // IMS_COMPOSE_F_GX_T_H
