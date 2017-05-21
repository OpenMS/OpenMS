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

#include <OpenMS/DATASTRUCTURES/Polynomial.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <functional>

using namespace std;

namespace OpenMS
{ 
  
  CounterSet::RangeCounter::RangeCounter(Size min, Size max, Size& val): 
    min_(min), max_(max), value(val)
  {
    value = min_;
  }

  CounterSet::RangeCounter& CounterSet::RangeCounter::operator++()
  {
    value = min_ + (this->value + 1) % (max_ - min_);
    return *this;
  }
  
  Size& CounterSet::RangeCounter::getValue()
  {
    return value;
  }
  
  void CounterSet::RangeCounter::reset()
  {
    value = min_;
  }

  bool CounterSet::RangeCounter::wasReset() const
  {
    return value == min_;
  }

  const Size& CounterSet::RangeCounter::min() const
  {
    return min_;
  }
  
  const Size& CounterSet::RangeCounter::max() const
  {
    return max_;
  }

  CounterSet::initCounter::initCounter(CounterSet::ContainerType& container):container_(container)
  {
    
  }
  CounterSet::RangeCounter CounterSet::initCounter::operator()(Range& r)
  {
    container_.push_back(0);
    return CounterSet::RangeCounter(r.first, r.second, container_.back());
  }

  void CounterSet::initCounter::operator()(RangeCounter& c)
  {
    c.reset();
  }

  CounterSet::CounterSet(UInt n, vector<Range> ranges):N(n), initializer(this->counters)
  {
    transform(ranges.begin(), ranges.end(), back_inserter(range_counters), initializer);
  }

  CounterSet::CounterSet(UInt n): N(n), initializer(this->counters)
  {
  }

  void CounterSet::reset()
  {
    for_each(range_counters.begin(), range_counters.end(), initializer);
    hasNext = true;
  }
  
  Size& CounterSet::operator[](const Size& index)
  {
      return counters[index];
  }

  CounterSet& CounterSet::operator++()
  {
    if(!hasNext)
    {
      LOG_WARN << "Overflow in counter set" << endl;
      return *this;
    }

    for(Ranges::reverse_iterator it = range_counters.rbegin(); it != range_counters.rend(); ++it)
    {
      UInt sum = 0;
      accumulate(counters.begin(), counters.end(), sum, plus<UInt>());
      if(sum > this->N)
      {
        it->reset();
        continue;
      }
      else if(sum < this->N)
      {
        ++(*it);
        it = it == range_counters.rbegin() ? it : --it;
      }
      else
      {
        return *this;
      }
      //if(!it->wasReset())
      //{
      //  return *this;
      //}
    }
    //hasNext=false;
    return *this;
  }
  
  const CounterSet::ContainerType& CounterSet::getCounters()
  {
    return counters;
  }
        
  CounterSet::ConstIterator CounterSet::begin() const
  {
    return counters.begin();
  }
  
  CounterSet::ConstIterator CounterSet::end() const
  {
    return counters.end();
  }
  
  void CounterSet::addCounter(Size min, Size max)
  {
    counters.push_back(0);
    range_counters.push_back(RangeCounter(min, max, counters.back()));
  }
  

  const bool& CounterSet::notLast() const
  {
    return hasNext;
  }

}
