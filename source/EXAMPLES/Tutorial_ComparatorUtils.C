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

#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/KERNEL/StandardTypes.h>
using OpenMS::Int;
using OpenMS::UInt;
using OpenMS::Size;
using OpenMS::Real;
using OpenMS::String;

// Simple class with three data members
class IntRealString
{
public:
  IntRealString(Int i, Real r, String s) :
    i_(i), r_(r), s_(s) {}
  IntRealString(const IntRealString & rhs) :
    i_(rhs.i_), r_(rhs.r_), s_(rhs.s_) {}

  void print() const
  {
    std::cout << "(" << i_ << ", " << r_ << ", " << s_ << ")" << std::endl;
  }

  Int i_;
  Real r_;
  String s_;
};

// A vector of IntRealString objects
class IntRealStringVector :
  public std::vector<IntRealString>
{
public:
  void print() const
  {
    for (Size i = 0; i < size(); ++i) (*this)[i].print();
    std::cout << std::endl;
  }

};

// Comparison function
bool lessByInt(IntRealString left, IntRealString right)
{
  return left.i_ < right.i_;
}

// Comparator class
struct LessByInt :
  std::binary_function<IntRealString, IntRealString, bool>
{
  bool operator()(IntRealString left, IntRealString right) const { return left.i_ < right.i_; }
};

// Comparator class
struct LessByReal :
  std::binary_function<IntRealString, IntRealString, bool>
{
  bool operator()(IntRealString left, IntRealString right) const { return left.r_ < right.r_; }
};

// Comparator class
struct LessByString :
  std::binary_function<IntRealString, IntRealString, bool>
{
  bool operator()(IntRealString left, IntRealString right) const { return left.s_ < right.s_; }
};

Int main()
{
  IntRealStringVector vec;
  vec.push_back(IntRealString(1, 4.5f, "paul"));
  vec.push_back(IntRealString(2, 4.5f, "josie"));
  vec.push_back(IntRealString(1, 4.5f, "john"));
  vec.push_back(IntRealString(2, 3.9f, "kim"));

  std::cout << "After initialization:" << std::endl;
  vec.print();

  std::cout << "Sorted using lessByInt function:" << std::endl;
  std::sort(vec.begin(), vec.end(), lessByInt);
  vec.print();

  std::cout << "Sorted using LessByInt comparator class:" << std::endl;
  std::sort(vec.begin(), vec.end(), LessByInt());
  vec.print();

  std::cout << "Sorted using reversed LessByInt comparator class:" << std::endl;
  std::sort(vec.begin(), vec.end(), OpenMS::reverseComparator(LessByInt()));
  vec.print();

  std::cout << "Sorted using lexicographic order: 1. LessByInt, 2. LessByReal" << std::endl;
  std::sort(vec.begin(), vec.end(), OpenMS::lexicographicComparator(LessByInt(), LessByReal()));
  vec.print();

  std::cout << "Sorted using lexicographic order: 1. reversed LessByInt, 2. LessByReal, 3. LessByString" << std::endl;
  std::sort(vec.begin(), vec.end(),
            OpenMS::lexicographicComparator
              (OpenMS::lexicographicComparator
              (
                OpenMS::reverseComparator(LessByInt()),
                LessByReal()
              ),
              LessByString()
              )
            );
  vec.print();

  // vector of pointers into vec
  std::vector<const IntRealString *> ptr_vec;
  for (Size i = 0; i < vec.size(); ++i)
  {
    ptr_vec.push_back(&vec[i]);
  }

  std::cout << "ptr_vec before sorting" << std::endl;
  for (Size i = 0; i < ptr_vec.size(); ++i)
    ptr_vec[i]->print();
  std::cout << std::endl;

  std::sort(ptr_vec.begin(), ptr_vec.end(), OpenMS::pointerComparator(LessByString()));

  std::cout << "ptr_vec after sorting with pointerComparator(LessByString())" << std::endl;
  for (Size i = 0; i < ptr_vec.size(); ++i)
    ptr_vec[i]->print();
  std::cout << std::endl;

  return 0;
} //end of main
