// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CONCEPT_UNARYCOMPOSEFUNCTIONADAPTER_H
#define OPENMS_CONCEPT_UNARYCOMPOSEFUNCTIONADAPTER_H

#include <functional>

namespace OpenMS
{

  /**
    @brief Represents the function object unary adapter.

     This simplest and most fundamental compose function adapter uses the result of
     a unary operation as input to another unary operation. For more details, please
     refer to the book "The C++ Standard Library" by Nicolay Josuttis.
     This class is implemented in order to reduce (substitute) the usage of user defined functors.
     You can use it for example to combine several operation in one function call.
     Here is an example of combining "minus 5 and multiply by 10":
     @c compose_f_gx(bind2nd(multiplies<int>(),10),bind2nd(minus<int>(), 5))

     Another example of @c compose_f_gx_t usage is along with @c mem_ref_fun function
     to search for elements in container based on the certain element's property.
     Here is an example:
     @code
#include <string>
#include <vector>
#include <utility>
#include <algorithm>

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
  unaryCompose(bind2nd(equal_to<string>(), "3"), mem_fun_ref(&Element::getA)));

 if (it != elements.end())
 {
  cout << "element 3 has index " << it-elements.begin() << endl;
 }
 else
 {
  cout << "element 3 is not in the container." << endl;
 }
 return 0;
}
    @endcode
     The Example has the following output:
    @code
    element 3 has index 1
    @endcode

     Copyright 1999 by Addison Wesley Longman, Inc. and Nicolai M. Josuttis.
     All rights resuntitlederved.

     Permission to use, copy, modify and distribute this software for personal and
     educational use is hereby granted without fee, provided that the above copyright
     notice appears in all copies and that both that copyright notice and this
     permission notice appear in supporting documentation, and that the names of Addison
     Wesley Longman or the author are not used in advertising or publicity pertaining
     to distribution of the software without specific, written prior permission. Addison
     Wesley Longman and the author make no representations about the suitability of this
     software for any purpose. It is provided "as is" without express or implied warranty.

     ADDISON WESLEY LONGMAN AND THE AUTHOR DISCLAIM ALL WARRANTIES WITH REGARD TO THIS
     SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT
     SHALL ADDISON WESLEY LONGMAN OR THE AUTHOR BE LIABLE FOR ANY SPECIAL, INDIRECT OR
     CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
     PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING
     OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
  */
  template <class OP1, class OP2>
  class UnaryComposeFunctionAdapter :
    public std::unary_function<typename OP2::argument_type, typename OP1::result_type>
  {
private:
    OP1 op1_;    // process: op1(op2(x))
    OP2 op2_;
public:
    /// Constructor
    UnaryComposeFunctionAdapter(const OP1 & o1, const OP2 & o2) :
      op1_(o1), op2_(o2)
    {
    }

    /// function call
    typename OP1::result_type
    operator()(const typename OP2::argument_type & x) const
    {
      return op1_(op2_(x));
    }

  };

  /**
    @brief convenience function for the @c UnaryComposeFunctionAdapter adapter
  */
  template <class OP1, class OP2>
  inline UnaryComposeFunctionAdapter<OP1, OP2>
  unaryCompose(const OP1 & o1, const OP2 & o2)
  {
    return UnaryComposeFunctionAdapter<OP1, OP2>(o1, o2);
  }

} // namespace OpenMS

#endif // OPENMS_CONCEPT_UNARYCOMPOSEFUNCTIONADAPTER_H
