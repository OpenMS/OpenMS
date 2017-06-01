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
#include <numeric>
#include <vector>
#include <deque>

namespace OpenMS
{
  // TODO(Nikos) add openMP support
  // TODO(Nikos) add documentation
  class OPENMS_DLLAPI CounterSet
  {
  public:
    class RangeCounter;
    typedef std::deque<UInt> ContainerType;
    typedef std::vector<RangeCounter> Ranges;

    class RangeCounter
    {
   private:
      UInt min_;
      UInt max_;
      UInt max_allowed_;
      UInt& value;

   public:
      RangeCounter(UInt min, UInt max, UInt& value);
      RangeCounter& operator++();
      UInt operator+=(const UInt&);
      void setMaxAllowedValue(UInt);
      inline UInt& getValue() const { return value; }
      inline void reset() { value = min_; }
      inline const UInt& min() const { return min_; }
      inline const UInt& max() const { return max_; }
      inline const UInt& maxAllowed() const { return max_allowed_; }
    };
    
    CounterSet(UInt);
    const ContainerType& getCounters() const {return counters;}
    CounterSet& operator++();
    void addCounter(UInt, UInt);
    void reset();
    inline UInt sum() const { return accumulate(counters.begin(), counters.end(), 0); }
    const bool& hasNext() const { return has_next; }
    void print(char*);
  
  private:
    UInt N;
    UInt min_sum;
    std::vector<RangeCounter> range_counters;
    ContainerType counters;
    bool has_next;
    Ranges::reverse_iterator count_it;

    void prepare();


  };

}
