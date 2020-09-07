// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/TypeAsString.h>
#include <OpenMS/KERNEL/DPeak.h>
#include <iostream>

double someFunc(int, int *)
{
  return 9.9;
}

template <typename, int = 5>
struct WOW
{
  static float staticMemberFunc(float &) {return 0; }
  double memberFunc(const double &) {return 0; }
};

template <typename X, int I, unsigned J = 666, template <class, int> class WauWau = WOW>
struct Oink {};

int main()
{
  using namespace OpenMS;
  using namespace std;

  cout << typeAsString(1) << endl;
  cout << typeAsString(2u) << endl;
  cout << typeAsString(3.) << endl;
  cout << typeAsString(4.f) << endl;
  cout << endl;
  cout << typeAsString<Int>() << endl;
  cout << typeAsString<PointerSizeUInt>() << endl;
  cout << endl;
  cout << typeAsString(Peak1D()) << endl;
  cout << typeAsString(DPeak<1>::Type()) << endl;
  cout << typeAsString(DPeak<1>::Type::PositionType()) << endl;
  cout << typeAsString<DPeak<1>::Type::CoordinateType>() << endl;
  cout << typeAsString<DPeak<1>::Type::IntensityType>() << endl;
  cout << endl;
  cout << typeAsString(&someFunc) << endl;
  cout << typeAsString<WOW<char const * const *** const &> >() << endl;
  cout << typeAsString<Oink<double, 55> >() << endl;
  cout << typeAsString(&WOW<string, 8>::staticMemberFunc) << endl;
  cout << typeAsString(&WOW<char, 8>::memberFunc) << endl;

  return 0;
} // end of Tutorial_typeAsString.cpp
