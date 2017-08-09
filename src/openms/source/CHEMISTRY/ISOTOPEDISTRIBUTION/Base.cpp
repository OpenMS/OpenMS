// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Chris Bielow $
// $Authors: Clemens Groepl, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------
//


#include <cmath>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <functional>
#include <numeric>
#include <fstream>


#include <boost/utility.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/range/adaptor/reversed.hpp>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/Base.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/DATASTRUCTURES/Polynomial.h>

using namespace std;

namespace OpenMS
{

  IsotopeDistribution::IsotopeDistribution() :
    max_isotope_(0)
  {
    distribution_.push_back(IsotopeDistribution::MassAbundance(0, 1));
  }

  IsotopeDistribution::IsotopeDistribution(Size max_isotope) :
    max_isotope_(max_isotope)
  {
    distribution_.push_back(IsotopeDistribution::MassAbundance(0, 1));
  }

  IsotopeDistribution::IsotopeDistribution(const IsotopeDistribution & isotope_distribution) :
    max_isotope_(isotope_distribution.max_isotope_),
    distribution_(isotope_distribution.distribution_)
  {
  }

  IsotopeDistribution::~IsotopeDistribution()
  {
  }

  void IsotopeDistribution::setMaxIsotope(Size max_isotope)
  {
    max_isotope_ = max_isotope;
  }

  Size IsotopeDistribution::getMaxIsotope() const
  {
    return max_isotope_;
  }

  IsotopeDistribution & IsotopeDistribution::operator=(const IsotopeDistribution & iso)
  {
    if (this != &iso)
    {
      distribution_ = iso.distribution_;
      max_isotope_ = iso.max_isotope_;
    }
    return *this;
  }


  void IsotopeDistribution::set(const ContainerType & distribution)
  {
    distribution_ = distribution;
  }

  const IsotopeDistribution::ContainerType & IsotopeDistribution::getContainer() const
  {
    return distribution_;
  }

  Size IsotopeDistribution::getMax() const
  {
    if (distribution_.empty())
    {
      return 0;
    }
    return distribution_[distribution_.size() - 1].first;
  }

  Size IsotopeDistribution::getMin() const
  {
    if (distribution_.empty())
    {
      return 0;
    }
    return distribution_[0].first;
  }


  Size IsotopeDistribution::size() const
  {
    return distribution_.size();
  }

  void IsotopeDistribution::clear()
  {
    distribution_.clear();
    max_isotope_ = 0;
  }

  void IsotopeDistribution::trimIntensities(double cutoff)
  {
    distribution_.erase(
      remove_if(distribution_.begin(),
                distribution_.end(),
                [&cutoff](MassAbundance& sample)
                {
                  return sample.second > cutoff;
                }));
  }

  void IsotopeDistribution::sort_(function<bool(pair<Mass, Abundance>&, pair<Mass, Abundance>&)> sorter)
  {
    sort(distribution_.begin(), distribution_.end(), sorter);
  }

  void IsotopeDistribution::sortByIntensity()
  {
    sort_([](MassAbundance& p1, MassAbundance& p2){
        return p1.second > p2.second;
      });
  }
  
  void IsotopeDistribution::transform_(function<void(MassAbundance&)> lambda)
  {
    for_each(distribution_.begin(), distribution_.end(), lambda);
  }


  IsotopeDistribution IsotopeDistribution::operator+(const IsotopeDistribution & iso) const
  {
    throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    return *this;
  }

  IsotopeDistribution & IsotopeDistribution::operator+=(const IsotopeDistribution & iso)
  {
    throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    return *this;
  }

  IsotopeDistribution & IsotopeDistribution::operator*=(Size factor)
  {
    throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    return *this;
  }

  IsotopeDistribution IsotopeDistribution::operator*(Size factor) const
  {
    throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    return *this;
  }

  bool IsotopeDistribution::operator==(const IsotopeDistribution & isotope_distribution) const
  {
    return max_isotope_ == isotope_distribution.max_isotope_ &&
           distribution_ == isotope_distribution.distribution_;
  }

  bool IsotopeDistribution::operator!=(const IsotopeDistribution & isotope_distribution) const
  {
    return !(isotope_distribution == *this);
  }

  void IsotopeDistribution::renormalize()
  {
    if (distribution_.size() != 0)
    {
      double sum(0);
      // loop backwards as most distributions contains a lot of small values at the end
      for (ContainerType::const_reverse_iterator it = distribution_.rbegin(); it != distribution_.rend(); ++it)
      {
        sum += it->second;
      }

      for (Iterator it = distribution_.begin(); it != distribution_.end(); ++it)
      {
        it->second /= sum;
      }
    }
    return;
  }

  void IsotopeDistribution::trimRight(double cutoff)
  {
    ContainerType::reverse_iterator riter = distribution_.rbegin();

    // loop from right to left until an entry is larger than the cutoff
    for (; riter != distribution_.rend(); ++riter)
    {
      if (riter->second >= cutoff)
        break;
    }
    // trim the container
    distribution_.resize(riter.base() - distribution_.begin());
  }

  void IsotopeDistribution::trimLeft(double cutoff)
  {
    for (ContainerType::iterator iter = distribution_.begin(); iter != distribution_.end(); ++iter)
    {
      if (iter->second >= cutoff)
      {
        distribution_.erase(distribution_.begin(), iter);
        break;
      }
    }
  }

  inline bool desc_prob(const struct MIDAsPolynomialID::PMember& p0, const struct MIDAsPolynomialID::PMember& p)
  {
    return p0.probability > p.probability;
  }

  inline bool by_power(const struct MIDAsPolynomialID::PMember& p0, const struct MIDAsPolynomialID::PMember& p)
  {
    return p0.power < p.power;
  }

  inline bool zero_prob(const MIDAsPolynomialID::PMember& m)
  {
    return m.probability == 0;
  }

  inline bool zero_power(const MIDAsPolynomialID::PMember& m)
  {
    return m.power == 0;
  }

  inline bool lightest(const IsotopeDistribution::MassAbundance& a, const IsotopeDistribution::MassAbundance& b)
  {
    return a.first < b.first;
  }

  double lightest_element(const Element& el)
  {
    return min_element(el.getIsotopeDistribution().begin(), el.getIsotopeDistribution().end(), lightest)->first;
  }


  /* Start of the midas interface */
  MIDAs::MIDAs(EmpiricalFormula& formula, double resolution, UInt N_):
    IsotopeDistribution(),
    min_prob(1e-16),
    formula_(formula),
    resolution_(resolution),
    N(N_)
  {
    
  }
  

  MIDAs::MIDAs() : IsotopeDistribution()
  {
  }

  MIDAs::MIDAs(const IsotopeDistribution& isotope_distribution) : 
    IsotopeDistribution(isotope_distribution)
  {
  }

  void MIDAs::merge(Polynomial& raw, double resolution)
  {
    //raw must be ordered to work correctly ascending order on power field

    UInt output_size = ceil((raw.back().power - raw.front().power)/resolution);
    LOG_INFO << "output size " << output_size << endl;
    LOG_INFO << "raw size " << raw.size() <<endl;

    distribution_.clear();
    distribution_.resize(output_size, make_pair<double, double>(0, 0));

    for(auto& p : raw)
    {
      // Is this the case?
      
      UInt index = round((p.power - raw.front().power)/resolution);
      if(index >= distribution_.size()){

        LOG_INFO << index <<endl;
        
      }
      distribution_[index].first = distribution_[index].first == 0 ?
                                    raw.front().power * index :
                                    distribution_[index].first;
      distribution_[index].second += p.probability;
    }

  }

  void MIDAs::dumpIDToFile(String file)
  {
    ofstream out(file.c_str());
    for(auto& sample : distribution_)
    {
      out << sample.first << sample.second << endl;
    }
    
    out.close();
  }


  /* Start of the midas polynomial method */


  MIDAsPolynomialID::MIDAsPolynomialID(EmpiricalFormula& formula, double resolution):
    MIDAs(formula, resolution, 10),
    lighter_isotope(0)
  {
    for(EmpiricalFormula::const_iterator el = formula_.begin(); el != formula_.end(); ++el)
    {
      lighter_isotope += lightest_element(*(el->first)) * (el->second);
    }

    LOG_INFO << "Fine resolution: " << resolution_ << endl;

    mw_resolution = 1e-12;
    


  }

  inline double MIDAsPolynomialID::fact_ln(UInt x)
  {
    return boost::math::lgamma(x+1);
  }

  void MIDAsPolynomialID::run()
  {
    vector<Polynomial> el_dist;
    for(EmpiricalFormula::ConstIterator element = formula_.begin(); element != formula_.end(); ++element)
    {
      el_dist.push_back(generatePolynomial(*(element->first), element->second));
      LOG_INFO << element->first->getName() <<" has " << el_dist.back().size() << " data points " << endl;
    }
    Polynomial& T = *(el_dist.begin());

    for(vector<Polynomial>::iterator pol = boost::next(el_dist.begin()); pol != el_dist.end(); ++pol)
    {
      multiplyPolynomials(T, *pol);
    }

    LOG_INFO << "T after multiplication has " << T.size() <<" elements" << endl;

    //BZ2_bzWrite ( NULL, NULL, (void*)NULL, NULL );

    LOG_INFO << "RESULTS---------------" << endl;
    double probability = 0;
    for(Polynomial::const_iterator it = T.begin(); it != T.end(); ++it)
    {
      //LOG_INFO << it->power*mw_resolution << " " << it->probability << endl;
      probability += it->probability;
    }
    LOG_INFO << "probability sum " << probability <<endl;

    sort(T.begin(), T.end(), by_power);
    for(auto& pmember : T)
    {
      pmember.power *= mw_resolution;
    }

    merge(T, 0.0001);

    LOG_INFO << "Lightest theoretical element " << lighter_isotope << endl;

    trimRight(0.0001);
    trimLeft(0.0001);
    LOG_INFO << "Final distribution has " << distribution_.size() <<endl;
    for(ContainerType::const_iterator it = distribution_.begin(); it != distribution_.end(); ++it)
    {

    }


    LOG_INFO << "Isotope Distribution of " << formula_.toString() << " successfully computed " << endl;
    LOG_INFO << "Isotope Distribution has " << T.size() << " data points " << endl;

  }

  

  void addCounter(CounterSet& c, const double& abundance, const UInt& size, const UInt& N)
  {
    double expectation = size * abundance;
    double var = size * abundance *(1 - abundance);
    UInt U = expectation + (N * sqrt(1 + var));
    UInt B = expectation > (N * sqrt(1 + var)) ? ceil(expectation - (N * sqrt(1 + var))) : 0;
    //LOG_INFO << "Added counter with values " << B << " " << U <<endl;
    c.addCounter(B, U);
  }


  MIDAsPolynomialID::Polynomial MIDAsPolynomialID::generatePolynomial(const Element& p, const SignedSize size)
  {
    std::vector<unsigned long> base_power;
    std::vector<double> log_prob;
    const IsotopeDistribution::ContainerType& isotope = p.getIsotopeDistribution().getContainer();
    CounterSet c(size);
    Polynomial pol;

    for(IsotopeDistribution::ConstIterator iso_it = isotope.begin(); iso_it != isotope.end(); ++iso_it)
    {
      if(iso_it->second == 0)
      {
        continue;
      }
      addCounter(c, iso_it->second, size, N);
      base_power.push_back(round(iso_it->first / mw_resolution));
      log_prob.push_back(log(iso_it->second));
    }

    for(const CounterSet::ContainerType& counters = c.getCounters(); c.hasNext(); ++c)
    {
      MIDAsPolynomialID::PMember member;
      // UInt s = 0;
      member.power = 0;
      member.probability = fact_ln(size);
      UInt index = 0;
      for(CounterSet::ContainerType::const_iterator iso_count = counters.begin(); iso_count != counters.end(); ++iso_count, ++index)
      {
        member.probability += ((*iso_count) * log_prob[index]) - fact_ln((*iso_count));
      }

      member.probability = exp(member.probability);

      if(member.probability < min_prob)
      {
        continue;
      }
      // check if it is faster having another iteration
      index = 0;
      for(CounterSet::ContainerType::const_iterator iso_count = counters.begin(); iso_count != counters.end(); ++iso_count, ++index)
      {
        //LOG_INFO << *iso_count <<" "<<base_power[index]<<endl;
        member.power += (*iso_count)*base_power[index];
      }
//      LOG_INFO << member.power <<" "<<member.probability<<endl;

      pol.push_back(member);

    }

    return pol;
  }


  void MIDAsPolynomialID::multiplyPolynomials(Polynomial& f, Polynomial& g)
  {

    LOG_INFO << "Sorting polynomial" << f.size() << " and " << g.size() <<endl;
    sort(f.begin(), f.end(), desc_prob);
    sort(g.begin(), g.end(), desc_prob);

    LOG_INFO << "Multiplying polynomial" << f.size() << " and " << g.size() <<endl;
    double min_mass = min_element(g.begin(), g.end(), by_power)->power
                      + min_element(f.begin(), f.end(), by_power)->power;
    double max_mass = max_element(g.begin(), g.end(), by_power)->power
                      + max_element(f.begin(), f.end(), by_power)->power;
    double delta_mass = resolution_/mw_resolution;
    UInt size = round((max_mass-min_mass)/delta_mass);
    Polynomial fgid(size, PMember());

    for(Polynomial::iterator g_it = g.begin(); g_it != g.end(); ++g_it)
    {
      for(Polynomial::iterator f_it = f.begin(); f_it != f.end(); ++f_it)
      {
        double prob = f_it->probability*g_it->probability;
        if(prob > min_prob)
        {
          double mass = f_it->power + g_it->power;
          UInt bin = round((mass - min_mass) / delta_mass);
          fgid[bin].probability += prob;
          fgid[bin].power += mass * prob;
        }
        else
        {
          // Polynomials are sorted based on probability so we can safely break
          break;
        }
      }
    }

    fgid.erase(remove_if(fgid.begin(), fgid.end(), zero_prob), fgid.end());
    for(Polynomial::iterator f_it = fgid.begin(); f_it != fgid.end(); ++f_it)
    {
      f_it->power /= f_it->probability;
    }
    f = fgid ;
  }

  MIDAsFFTID::MIDAsFFTID(EmpiricalFormula& formula, double resolution):
    MIDAs(formula, resolution, 15),
    cutoff_amplitude_factor_(2)
  {
    UInt sample_size,k=0;
    Stats coarse(formulaMeanAndVariance()), fine(formulaMeanAndVariance(resolution_));
    double sigma;
    
    sigma = coarse.variance;
    double range = N * sqrt(1 + sigma);
    mass_range_ = pow(2, ceil(log2(ceil(range))));
    LOG_INFO << "Resolution " << resolution_ << endl;
    LOG_INFO << "Coarse Average mass " << coarse.mean << " Variance: " << coarse.variance << endl;
    LOG_INFO << "Fine Average mass " << fine.mean << " Variance: " << fine.variance << endl;

    do
    {
      resolution_ = resolution/pow(2,k);  
      LOG_INFO << "Mass range " << mass_range_ << endl;
      sample_size = pow(2, ceil(log2(mass_range_ / resolution_)));
      delta_ = 1.0 / sample_size;
      resolution_ = mass_range_/sample_size;
      LOG_INFO << "-Resolution " << resolution << endl;
      k++;
    }while(resolution_ > resolution);
    
    average_mass_ = round(fine.mean) / resolution_;
    LOG_INFO <<"Mass Range " << mass_range_ <<endl;

    fft_complex s = {0,0};
    input_.resize(sample_size, s);
    output_.resize(sample_size, s);
    
    
    LOG_INFO << "Sample size: " << input_.size() << "   " << output_.size() << endl;
    init();

  }

  void MIDAsFFTID::init()
  {

    LOG_INFO <<"Average mass "<< average_mass_ <<endl;
    LOG_INFO << "Resolution " << resolution_ << endl;
    Int k = 0;//-input_.size()/2;
    for(auto& sample : input_)
    {
      Int j = k > input_.size() / 2?  k++ - input_.size() : k++;
  
      double phi = 0, angle = 0, radius = 1, freq = j * delta_;
      double phase = (2 * Constants::PI * average_mass_ * freq);

      //LOG_INFO << "Delta:" << freq << endl;
      for(const auto& element : formula_)
      {
        //Perform temporary calculations on sample data structure
        auto& atoms = element.second;
        
        sample.r = sample.i = 0;
        for(const auto& iso : element.first->getIsotopeDistribution())
        {
          auto mass = round(iso.first / resolution_);
          auto& prob = iso.second;
          if(!(prob > 0))
          {
            continue;
          }
          phi = 2 * Constants::PI * mass * freq;
          sample.r += prob * cos(phi);
          sample.i += prob * sin(phi);
          
        }
        //LOG_INFO<<"x,y " << sample.r << " " <<sample.i << endl;
        radius *= pow(hypot(sample.r, sample.i), atoms);
        angle += atoms * atan2(sample.i, sample.r);
        //LOG_INFO<<"radius,angle " << radius << " " << angle << endl;
        
      }
      
      //After looping assign the value
      
      //LOG_INFO << " radius " << radius << " angle "<< angle << endl;
      sample.r = radius * cos(angle - phase);
      sample.i = radius * sin(angle - phase);
      //LOG_INFO << sample.r << " " << sample.i << endl;
    }

    //input_.back().r = input_.back().i = 0;
    LOG_INFO << "End of initialization" << endl;
    
  }

  void four1(double *Data, int nn, int isign)
  {
    unsigned long i, j, m, n, mmax, istep;
    double wr, wpr, wpi, wi, theta;
    double wtemp, tempr, tempi;
    double one_pi = acos(-1);
    double two_pi= 2*one_pi;

    /* Perform bit reversal of Data[] */
    n = nn << 1;
    j=1;
    for (i=1; i<n; i+=2)
      {
        if (j > i)
          {
            wtemp = Data[i];
            Data[i] = Data[j];
            Data[j] = wtemp;
            wtemp = Data[i+1];
            Data[i+1] = Data[j+1];
            Data[j+1] = wtemp;
          }
        m = n >> 1;
        while (m >= 2 && j > m)
          {
            j -= m;
            m >>= 1;
          }
        j += m;
      }
 
    /* Perform Danielson-Lanczos section of FFT */
    n = nn << 1;
    mmax = 2;
    while (n > mmax)  /* Loop executed log(2)nn times */
      {
	istep = mmax << 1;
	theta = isign * (two_pi/mmax);  
	wtemp = sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m=1; m<mmax; m+=2)
	  {
	    for (i=m; i<=n; i+=istep)
	      {
		j = i+mmax;                      

		tempr = wr*Data[j]-wi*Data[j+1];
		tempi = wr*Data[j+1]+wi*Data[j];
		Data[j] = Data[i]-tempr;
		Data[j+1] = Data[i+1]-tempi;
		Data[i] += tempr;
		Data[i+1] += tempi;
	      }
	    wr = (wtemp=wr)*wpr-wi*wpi+wr;
	    wi = wi*wpr+wtemp*wpi+wi;
	  }
	mmax = istep;
      }
  }

  void MIDAsFFTID::run()
  {
    
    kiss_fft_cfg cfg = kiss_fft_alloc(input_.size(), INVERSE, NULL, NULL);
    kiss_fft(cfg, &(*input_.begin()), &(*output_.begin()));
    kiss_fft_cleanup();

    bool fft_midas_bypass = true;

    if(!fft_midas_bypass)
    {
      LOG_INFO << "Using midas native fft method " << endl;
      UInt size = output_.size()*2 + 1;
      double *input = new double[ size ];
      //double *output = new double[ size ];
      UInt k = 1;
      for(auto& sample : input_)
      {
        if(2*k >= size)
        {
          continue;
        }
        input[2*k-1] = sample.r;
        input[2*k] = sample.i;
        k++;
      }

      four1(input, size/2, -1);
    
      k = 1;
      for(auto& sample : output_)
      {
        if(2*k >= size)
        {
          continue;
        }
        sample.r = input[2*k-1];
        sample.i = input[2*k];
        k++;
      }
    }

    // Resume normal operation
    output_.resize(output_.size()/2);

    LOG_INFO << "IFFT done " <<" Sample size: "<< output_.size() <<endl;
    
    double min_prob = -cutoff_amplitude_factor_ * 
      min_element(output_.begin(), output_.end(), 
                  [](const fft_complex& item1, const fft_complex& item2 )
                  {
                    return (item1.r < item2.r);
                  })->r;
    
    for(auto& sample : output_)
    {
      //LOG_INFO << sample.r << endl ;//<< " "<< sample.i << endl; 
    }

    // double avg_prob= 0;
    // int k = 0;
    // for(auto& sample : output_)
    // {
    //   double smp = sample.i / output_.size();
    //   if(smp > min_prob)
    //   {
    //     avg_prob+= smp;
    //     k++;
    //   }
    // }
    // if(k!=0){
    //   avg_prob/=k;
    // }
    LOG_INFO << "Resolution: " << resolution_ << endl;

    Stats coarse(formulaMeanAndVariance()), fine(formulaMeanAndVariance(resolution_));
    double ratio = coarse.variance/fine.variance;
    LOG_INFO << "Delta " << delta_ <<endl;
    LOG_INFO << "Coarse mean: " << coarse.mean << ", Coarse variance: " << coarse.variance << endl;
    LOG_INFO << "Fine mean: " << fine.mean << ", fine variance: " << coarse.variance << endl;
    LOG_INFO << "Probability cutoff: " << min_prob << " Ratio: " << ratio <<endl;
    
    UInt k =  0;//-output_.size()/2;
    //unsigned int average_mass = ceil(fine.mean);
    Polynomial pol;
    double p_sum = 0;
    //for(auto& sample : boost::adaptors::reverse(output_))
    for(auto& sample : output_)
    {
      PMember member;
      member.probability = sample.r;
      Int j = k > output_.size()/2 ?  k++ - output_.size() : k++;
      
      if(member.probability > min_prob)
      {
      
        p_sum += member.probability;
        member.power = ratio*
                       ( (j  + average_mass_) * resolution_ - coarse.mean) + fine.mean;
        
        pol.push_back(member);
        LOG_INFO << member.power << " " << member.probability << endl;
      }
      
      
       
    }

    
    LOG_INFO << "Probability sum " << p_sum << endl;

    // // //normalize
    for(auto& point : pol)
    {
       point.probability /= p_sum;
       LOG_INFO << point.power <<" "<<point.probability << endl;
    }

    // // sort(pol.begin(), pol.end(), by_power);
    //merge(pol, resolution_);

  }



  MIDAsFFTID::Stats MIDAsFFTID::formulaMeanAndVariance(double resolution)
  {
    Stats stat = {0,0};

    //throw exception for zero resolution_

    // calculate average
    for(const auto& element : formula_)
    {
      double ave_mw = 0, var_mw = 0;
      for(const auto& iso : element.first->getIsotopeDistribution())
      {
        //round in resolution grid and weight on probability
        ave_mw += round(iso.first / resolution) * resolution * iso.second;
      }

      //calculate variance

      for(const auto& iso : element.first->getIsotopeDistribution())
      {
        //round in resolution grid
        var_mw += iso.second * pow(ave_mw - (round(iso.first / resolution) * resolution), 2);
      }
      
      // find the real variance and mean by scaling with the molecule number in the empirical formula
      stat.variance += element.second * var_mw;
      stat.mean += element.second * ave_mw;
    }

    return stat;

  }

}
