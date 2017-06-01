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
#include <algorithm>
#include <functional>

using namespace std;

namespace OpenMS
{ 
  
  CounterSet::RangeCounter::RangeCounter(UInt min, UInt max, UInt& val): 
    min_(min), max_(max), value(val)
  {
    reset();
  }

  UInt CounterSet::RangeCounter::operator+=(const UInt& num)
  {
    UInt diff = max_allowed_ - value;
    if(num > diff)
    {
      value = max_allowed_;
      return num - diff;
    }
    else
    {
      value = min_ + ((this->value - min_ + num) % (max_allowed_ + 1 - min_));
      return 0;
    }
  }

  CounterSet::RangeCounter& CounterSet::RangeCounter::operator++()
  {
    (*this)+=1;
    return *this;
  }

  void CounterSet::RangeCounter::setMaxAllowedValue(UInt counters_range)
  {
    max_ = counters_range < max() - min() ? counters_range + min() : max();
  }

  CounterSet::CounterSet(UInt n): N(n), has_next(true)
  {
  }

  void CounterSet::reset()
  {
    has_next = true;
    for(Ranges::iterator it = range_counters.begin(); it != range_counters.end(); ++it)
    {
      it->reset();
    }
    min_sum = accumulate(counters.begin(), counters.end(), 0);

    for(Ranges::reverse_iterator it = range_counters.rbegin(); it != range_counters.rend(); ++it)
    {
      UInt remain = ((*it) += (N - sum()));

      if(remain == 0)
      {
        count_it = it;
        break;
      }
    }
  }

  CounterSet& CounterSet::operator++()
  {

    if(!has_next)
    {
      LOG_WARN << "Overflow in counter set" << endl;
      return *this;
    }

    prepare();
    
    for(Ranges::reverse_iterator it = range_counters.rbegin(); it != count_it; ++it)
    {
      
      RangeCounter& range = *it;
      UInt remain = (range += (N - sum()));
      if(remain == 0)
      {
        count_it = it;
        if(&range_counters.back() == &range || count_it->getValue() == it->maxAllowed())
        {
          for(Ranges::reverse_iterator it2 = count_it; it2 != range_counters.rend(); ++it2)
          {
            if(count_it->getValue() < it->maxAllowed())
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
  
  void CounterSet::addCounter(UInt min, UInt max)
  {
    //Add a dummy value, this will be changed to min at the constructor
    counters.push_back(0);
    range_counters.push_back(RangeCounter(min, max, counters.back()));
    reset();
  }

  void CounterSet::prepare()
  {

    for(Ranges::reverse_iterator it = count_it; it != range_counters.rend(); ++it)
    {
      if(it->getValue() < it->maxAllowed() && &range_counters.back() != &(*it))
      {
        ++(*it);
        count_it = it;
        for(Ranges::reverse_iterator it2 = range_counters.rbegin(); it2 != count_it; ++it2)
        {
          it2->reset();
        }
        
        if(sum() > N)
        {
          while(&(*it) != &range_counters.front())
          {
            it->reset();
            ++it;
            UInt val = it->getValue();
            ++(*it);
            count_it = it; 
            if(it->getValue() != val && sum() <= N)
            {
              break;
            }
          }
          if(sum() > N)
          {
            has_next = false;
          }
          
        }
        return;
      }
    }
    has_next = false;
   
  }

}
