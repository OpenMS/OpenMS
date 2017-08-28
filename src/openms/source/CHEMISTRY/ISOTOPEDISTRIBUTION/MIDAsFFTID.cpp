
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/MIDAsFFTID.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/MATH/MISC/MIDAsFFT.h>

// #define DEBUG

using namespace std;

namespace OpenMS
{
  MIDAsFFTID::MIDAsFFTID(double resolution):
    MIDAs(resolution, 15),
    cutoff_amplitude_factor_(2)
  {

  }


  void MIDAsFFTID::init(const EmpiricalFormula& formula)
  {

#ifdef DEBUG
    LOG_INFO <<"Average mass "<< average_mass_ <<endl;
    LOG_INFO << "Resolution " << resolution_ << endl;
#endif
    UInt k = 0;//-input_.size()/2;
    for (auto& sample : input_)
    {
      Int j = k >= input_.size() / 2?  k++ - input_.size() : k++;
  
      double phi = 0, angle = 0, radius = 1, freq = j * delta_;
      double phase = (2 * Constants::PI * average_mass_ * freq);

      //LOG_INFO << "Delta:" << freq << endl;
      for (const auto& element : formula)
      {
        //Perform temporary calculations on sample data structure
        auto& atoms = element.second;
        
        sample.r = sample.i = 0;
        for (const auto& iso : element.first->getIsotopeDistribution())
        {
          auto mass = round(iso.getMZ() / resolution_);
          auto prob = iso.getIntensity();
          if (not (prob > 0))
          {
            continue;
          }
          phi = 2 * Constants::PI * mass * freq;
          sample.r += prob * cos(phi);
          sample.i += prob * sin(phi);
          
        }
        radius *= pow(hypot(sample.r, sample.i), atoms);
        angle += atoms * atan2(sample.i, sample.r);
        
      }
      
      //After looping assign the value   
      sample.r = radius * cos(angle - phase);
      sample.i = radius * sin(angle - phase);
    }
    
  }

  void MIDAsFFTID::run(const EmpiricalFormula& formula)
  {

    UInt sample_size, k = 0;
    Stats coarse(formulaMeanAndVariance(formula)), 
      fine(formulaMeanAndVariance(formula, resolution_));
    
    double sigma = coarse.variance;
    double range = N * sqrt(1 + sigma);
    mass_range_ = pow(2, ceil(log2(ceil(range))));


    double initial_resolution = resolution_;
    do
    {
      resolution_ = initial_resolution / pow(2,k);  
      sample_size = pow(2, ceil(log2(int(mass_range_ / resolution_))));
      delta_ = 1.0 / sample_size;
      resolution_ = mass_range_ / sample_size;
      k++;
    }while(resolution_ > initial_resolution);
    
    average_mass_ = round(fine.mean) / resolution_;

    fft_complex s = {0,0};
    input_.resize(sample_size, s);
    output_.resize(sample_size, s);
    
#ifdef DEBUG
    LOG_INFO << "Resolution " << resolution_ << endl;
    LOG_INFO << "Coarse Average mass " << coarse.mean << " Variance: " << coarse.variance << endl;
    LOG_INFO << "Fine Average mass " << fine.mean << " Variance: " << fine.variance << endl;
    LOG_INFO << "Sample size: " << input_.size() << "   " << output_.size() << endl;
#endif

    // Create sample of formula
    init(formula);

    UInt size = output_.size() * 2 + 1;
    vector<double> input(size, 0);
    k = 1;

    //Convert samples FFT compatible input
    for(auto& sample : input_)
    {
      if(2*k >= size)
      {
        continue;
      }
      input[2 * k - 1] = sample.r;
      input[2 * k] = sample.i;
      k++;
    }

    //Calculate the inverse fft
    fft(&input.front(), size / 2);
    
    //Dump back the values in the kissfft format
    k = 1;
    for(auto& sample : output_)
    {
      if(2*k >= size)
      {
        continue;
      }
      sample.r = input[2*k-1] / input_.size();
      sample.i = input[2*k] / input_.size();
      k++;
    }


    // Resume normal operation

#ifdef DEBUG
    LOG_INFO << "Inverse FFT done. " <<" Sample size: "<< output_.size() <<endl;
#endif

    double min_prob = -cutoff_amplitude_factor_ * 
      min_element(output_.begin(), output_.end(), 
                  [](const fft_complex& item1, const fft_complex& item2 )
                  {
                    return (item1.r < item2.r);
                  })->r;
    
    double ratio = coarse.variance / fine.variance;

#ifdef DEBUG
    LOG_INFO << "Resolution: " << resolution_ << endl;
    LOG_INFO << "Delta " << delta_ <<endl;
    LOG_INFO << "Coarse mean: " << coarse.mean << ", Coarse variance: " << coarse.variance << endl;
    LOG_INFO << "Fine mean: " << fine.mean << ", fine variance: " << coarse.variance << endl;
    LOG_INFO << "Probability cutoff: " << min_prob << " Ratio: " << ratio <<endl;
#endif

    k = 0;
    Polynomial pol;
    double p_sum = 0;
    for(auto& sample : output_)
    {
      if(sample.r > min_prob)
      {
        Peak1D member;
        member.setIntensity(sample.r);
        Int j = k > input_.size() / 2?  k++ - input_.size() : k++;
        p_sum += member.getIntensity();
        member.setMZ(ratio * ((j + average_mass_) * resolution_ - fine.mean) + coarse.mean);
        pol.push_back(member);
      }
      k++;
    }

    distribution_.assign(pol.begin(), pol.end());
    renormalize();

    //sort by mass
    sort(distribution_.begin(), distribution_.end(), by_power);
    
    //merge(tmp, resolution_);

  }


  MIDAsFFTID::Stats MIDAsFFTID::formulaMeanAndVariance(const EmpiricalFormula& formula, double resolution)
  {
    Stats stat = {0,0};

    //throw exception for zero resolution_

    // calculate average
    for(const auto& element : formula)
    {
      double ave_mw = 0, var_mw = 0;
      for(const auto& iso : element.first->getIsotopeDistribution())
      {
        //round in resolution grid and weight on probability
        double mass = resolution < 1 ? round(iso.getMZ() / resolution) * resolution : iso.getMZ();  
        ave_mw += mass * iso.getIntensity();
      }

      //calculate variance

      for(const auto& iso : element.first->getIsotopeDistribution())
      {
        double mass = resolution < 1 ? round(iso.getMZ() / resolution) * resolution : iso.getMZ();  
        var_mw += iso.getIntensity() * pow(ave_mw - mass, 2);
      }
      
      // find the real variance and mean by scaling with the molecule number in the empirical formula
      stat.variance += element.second * var_mw;
      stat.mean += element.second * ave_mw;
    }

    return stat;
  }

}
