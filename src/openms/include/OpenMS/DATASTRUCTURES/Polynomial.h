// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Nikos Patikas $
// $Authors: Nikos Patikas $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Types.h>
#include <vector>
#include <deque>

namespace OpenMS
{
  typedef std::pair<Size, Size> Range;

  // TODO(Nikos) add openMP support
  // TODO(Nikos) add documentation
  class OPENMS_DLLAPI CounterSet
  {
  public:
    class RangeCounter;
    typedef std::deque<Size> ContainerType;
    typedef std::vector<RangeCounter> Ranges;
    typedef ContainerType::reverse_iterator ReverseIterator;
    typedef ContainerType::iterator Iterator;
    typedef ContainerType::const_iterator ConstIterator;
    
    class RangeCounter
    {
   private:
      Size min_;
      Size max_;
      Size& value;

   public:
      RangeCounter(Size min, Size max, Size& value);
      RangeCounter(const RangeCounter& other);
      RangeCounter& operator++();
      RangeCounter& operator--();
      UInt operator+=(const UInt&);
      Size& getValue() const;
      bool wasReset() const;
      void reset();
      const Size& min() const;
      const Size& max() const;
    };

    struct initCounter
    {
      ContainerType& container_;
      initCounter(ContainerType&);
      RangeCounter operator()(Range&);
      void operator()(RangeCounter& c);
    };
    
    CounterSet(UInt, std::vector<Range>);
    CounterSet(UInt);
    const ContainerType& getCounters() const;
    Size& operator[](const Size& index);
    CounterSet& operator++();
    void addCounter(Size min, Size max);
    void reset();
    ConstIterator begin() const;
    ConstIterator end() const;
    const bool& notLast();
    UInt sum() const;
  private:
    UInt N;
    UInt min_sum;
    struct initCounter initializer;
    Ranges::reverse_iterator count_it;
    void finished();
    UInt maxAllowedValue(RangeCounter&) const;
    void print(char*);
    //struct counterIterator generator;
    std::vector<RangeCounter> range_counters;
    ContainerType counters;
    bool hasNext;
    
  };


}
