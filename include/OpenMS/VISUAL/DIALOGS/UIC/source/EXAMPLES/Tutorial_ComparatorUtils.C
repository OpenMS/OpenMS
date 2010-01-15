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
  IntRealString(Int i, Real r, String s) : i_(i), r_(r), s_(s) {}
  IntRealString(const IntRealString& rhs) : i_(rhs.i_), r_(rhs.r_), s_(rhs.s_) {}

  void print() const
  {
    std::cout << "(" << i_ << ", " << r_ << ", " << s_ << ")" << std::endl;
  }

  Int i_;
  Real r_;
  String s_;
};

// A vector of IntRealString objects
class IntRealStringVector : public std::vector<IntRealString>
{
 public:
  void print() const
  {
    for (Size i = 0; i < size(); ++i) (*this)[i].print();
    std::cout << std::endl;
  }
};

// Comparison function
bool lessByInt ( IntRealString left, IntRealString right ) { return left.i_ < right.i_; }

// Comparator class
struct LessByInt
  : std::binary_function<IntRealString,IntRealString,bool>
{
  bool operator() (IntRealString left, IntRealString right ) const { return left.i_ < right.i_; }
};

// Comparator class
struct LessByReal
  : std::binary_function<IntRealString,IntRealString,bool>
{
  bool operator() (IntRealString left, IntRealString right ) const { return left.r_ < right.r_; }
};

// Comparator class
struct LessByString
  : std::binary_function<IntRealString,IntRealString,bool>
{
  bool operator() (IntRealString left, IntRealString right ) const { return left.s_ < right.s_; }
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
  std::sort(vec.begin(),vec.end(),lessByInt);
  vec.print();

  std::cout << "Sorted using LessByInt comparator class:" << std::endl;
  std::sort(vec.begin(),vec.end(),LessByInt());
  vec.print();

  std::cout << "Sorted using reversed LessByInt comparator class:" << std::endl;
  std::sort(vec.begin(),vec.end(),OpenMS::reverseComparator(LessByInt()));
  vec.print();

  std::cout << "Sorted using lexicographic order: 1. LessByInt, 2. LessByReal" << std::endl;
  std::sort(vec.begin(),vec.end(),OpenMS::lexicographicComparator(LessByInt(),LessByReal()));
  vec.print();

  std::cout << "Sorted using lexicographic order: 1. reversed LessByInt, 2. LessByReal, 3. LessByString" << std::endl;
  std::sort(vec.begin(),vec.end(),
      OpenMS::lexicographicComparator
      ( OpenMS::lexicographicComparator
        (
         OpenMS::reverseComparator(LessByInt()),
         LessByReal()
        ),
        LessByString()
      )
     );
  vec.print();

  // vector of pointers into vec
  std::vector<const IntRealString*> ptr_vec;
  for (Size i = 0; i < vec.size(); ++i)
  {
    ptr_vec.push_back(&vec[i]);
  }

  std::cout << "ptr_vec before sorting" << std::endl;
  for (Size i = 0; i < ptr_vec.size(); ++i) ptr_vec[i]->print();
  std::cout << std::endl;
  
  std::sort(ptr_vec.begin(),ptr_vec.end(),OpenMS::pointerComparator(LessByString()));

  std::cout << "ptr_vec after sorting with pointerComparator(LessByString())" << std::endl;
  for (Size i = 0; i < ptr_vec.size(); ++i) ptr_vec[i]->print();
  std::cout << std::endl;
  
  return 0;
} //end of main
