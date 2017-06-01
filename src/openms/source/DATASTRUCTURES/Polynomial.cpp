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
#include <boost/utility.hpp>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <functional>
#include <iterator>

using namespace std;

namespace OpenMS
{ 
  
  CounterSet::RangeCounter::RangeCounter(Size min, Size max, Size& val): 
    min_(min), max_(max), value(val)
  {
    value = min_;
    
  }

  CounterSet::RangeCounter::RangeCounter(const CounterSet::RangeCounter& obj):
    min_(obj.min()), max_(obj.max()), value(obj.getValue())
  {
    
  }

  UInt CounterSet::RangeCounter::operator+=(const UInt& num)
  {
    UInt diff = max_ - value;
    if(num > diff)
    {
      value = max_;
      return num - diff;
    }
    else
    {
      value = min_ + ((this->value - min_ + num) % (max_ + 1 - min_));
      return 0;
    }
  }

  CounterSet::RangeCounter& CounterSet::RangeCounter::operator++()
  {
    value = min_ + ((this->value - min_ + 1) % (max_ + 1 - min_));
    return *this;
  }

  CounterSet::RangeCounter& CounterSet::RangeCounter::operator--()
  {
    value = min_ + ((this->value - min_ - 1) % (max_ + 1 - min_));
    return *this;
  }
  
  Size& CounterSet::RangeCounter::getValue() const
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
    count_it = range_counters.rend();
    hasNext = true;
    reset();
  }

  CounterSet::CounterSet(UInt n): N(n), initializer(this->counters)
  {
    hasNext = true;
  }


  void CounterSet::reset()
  {
    for_each(range_counters.begin(), range_counters.end(), initializer);
    hasNext = true;
    min_sum = accumulate(counters.begin(), counters.end(), 0);

    for(Ranges::reverse_iterator it = range_counters.rbegin(); it != range_counters.rend(); ++it)
    {
      
      RangeCounter& range = *it;
      UInt remain = (range += (N - sum()));

      if(remain == 0)
      {
        count_it = it;
        break;
      }
    }


  }
  
  Size& CounterSet::operator[](const Size& index)
  {
      return counters[index];
  }

  UInt CounterSet::sum() const
  {
    return accumulate(counters.begin(), counters.end(), 0);
  }

  CounterSet& CounterSet::operator++()
  {
    if(!hasNext)
    {
      LOG_WARN << "Overflow in counter set" << endl;
      return *this;
    }

    
    for(Ranges::reverse_iterator it = range_counters.rbegin(); it != count_it; ++it)
    {
      
      RangeCounter& range = *it;
      UInt remain = (range += (N - sum()));
      

      if(remain == 0)
      {
        
        count_it = it;
        if(&range_counters.back() == &range || count_it->getValue() == maxAllowedValue(*it))
        {
          for(Ranges::reverse_iterator it2 = count_it; it2 != range_counters.rend(); ++it2)
          {
            if(count_it->getValue() < maxAllowedValue(*it2))
            {
              count_it = it2;
              break;
            }
          }
        }
        return *this;
      }
 
    }
    return *this;
  }
  
  const CounterSet::ContainerType& CounterSet::getCounters() const
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
    //Add a dummy value, this will be changed to min at the constructor
    counters.push_back(0);
    range_counters.push_back(RangeCounter(min, max, counters.back()));
    reset();
  }
  

  const bool& CounterSet::notLast()
  {
    finished();
    return hasNext;
  }

  void CounterSet::finished()
  {

    for(Ranges::reverse_iterator it = count_it; it != range_counters.rend(); ++it)
    {
      // TODO(Nikos) Optimize here maxAllowedValue as max
      if(it->getValue() < maxAllowedValue(*it) && &range_counters.back() != &(*it))
      {
        ++(*it);
        count_it = it;
        for(Ranges::reverse_iterator it2 = range_counters.rbegin(); it2 != count_it; ++it2)
        {
          it2->reset();
        }
        
        if(sum() > N)
        {
          if(boost::next(it) != range_counters.rend())
          {
            it->reset();
            ++it;
            ++(*it);
            count_it = it;
            if(sum() > N)
            {
              hasNext = false;
            }
          }
          else
          {
            hasNext = false;
          }
          return;
        }
        else
        {
          return;
        }
      }
    }
   
    hasNext = false;
   
  }

  UInt CounterSet::maxAllowedValue(RangeCounter& counter) const
  {
    if(N < counter.max() + min_sum - counter.min())
    {
      return N - min_sum + counter.min();
    }
    else
    {
      return counter.max();
    }
  }
  
  void CounterSet::print(char* fn)
  {
    cout<<"@@---"<<fn<<endl;
    copy(counters.begin(), counters.end(), ostream_iterator<UInt>(cout, " "));
    cout << endl;
    cout<<"@@--END"<<endl;
  }

}
