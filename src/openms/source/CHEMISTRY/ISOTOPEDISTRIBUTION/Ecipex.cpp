#include <cstring>
#include <cassert>
#include <iostream>

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/Ecipex.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/MATH/MISC/KissFFTWrapper.h>

#include <kiss_fft.h>

using namespace std;

namespace OpenMS
{

  

  Ecipex::Ecipex() : IsotopePatternGenerator()
  {
    distribution_.push_back(Peak1D(0, 1));
  }

  Ecipex::Ecipex(double threshold, double fft_threshold):
    fft_threshold_(fft_threshold),
    threshold_(threshold)
  {
    distribution_.push_back(Peak1D(0, 1));
  }

  Ecipex::Ecipex(const IsotopeDistribution& isotope_distribution) :
    IsotopePatternGenerator(isotope_distribution)
  {
  }
  
  void Ecipex::run(const EmpiricalFormula& formula)
  {
    for (auto& element : formula) 
    {
      auto pattern = elementIsotopePattern_(
        element.first->getIsotopeDistribution().getContainer(), 
        element.second, 
        fft_threshold_);
      
      convolve_(pattern, threshold_ * threshold_);
    }
    trimIntensities(threshold_);
    sortByMass();
  }



  void Ecipex::sortByIntensityAndNormalize_()
  {
    sortByIntensity();
    
    double max_intensity = distribution_.front().getIntensity();

    transform_([&max_intensity](MassAbundance& m)
               {
                 m.setIntensity(m.getIntensity() / max_intensity);
               });
    
  }

  Ecipex Ecipex::elementIsotopePattern_(const Spectrum& iso_pattern, 
                                       UInt amount, 
                                       double threshold)
  {
    Int dim = iso_pattern.size() - 1;
    if (dim == 0) 
    {
      Ecipex result;
      result[0].setPosition(iso_pattern.front().getMZ() * amount);
      return result;
    }

    if(dim < 0)
    {
      throw Exception::InvalidValue(__FILE__, 
                                    __LINE__, 
                                    OPENMS_PRETTY_FUNCTION, 
                                    "Invalid dimension size", 
                                    String(dim));
    }

    UInt64 edge_len = kiss_fftr_next_fast_size_real(amount + 1);

    vector<Int> dimensions(dim, edge_len);
    FftArray arr{dimensions};

    // array setup
    const auto& iso = iso_pattern;
    arr.scalar_data()[0] = iso.front().getIntensity();
    for (Int i = 0, k = 1; i < dim; i++, k *= edge_len)
    {
      arr.scalar_data()[k] = iso[dim-i].getIntensity();
    }
  
    // forward FFT
    arr.forwardFFT();

    // exponentiation
    for (UInt64 i = 0; i < arr.size(); ++i)
    {
      arr.data()[i] = pow(arr.data()[i], amount);
    }
    // inverse FFT
    arr.inverseFFT();
    for (UInt64 i = 0; i < arr.size(); ++i)
    {
      arr.scalar_data()[i] /= double(arr.size());  // take care of FFT normalization
     }
  
    Ecipex result;
    //Clear the unit convolutional vector
    result.clear();

    vector<UInt64> indices(dim, 0);
    for (UInt64 i = 0; i < arr.size(); ++i) 
    {
      UInt64 k;
      UInt64 n = 0;
      double mass = 0.0;
      for (k = 0; k < indices.size(); ++k) 
      {
        mass += iso_pattern[k + 1].getMZ() * indices[k];
        n += indices[k];
      }

      if (n <= amount && arr.scalar_data()[i] >= threshold) 
      {
        mass += (amount - n) * iso_pattern[0].getMZ();
        result.insert(mass, arr.scalar_data()[i]);
      }

      if (i == arr.size() - 1) 
      {
        break;
      }

      for (k = indices.size() - 1; indices[k] == edge_len - 1; --k)  
      {
      }
        
      indices[k] += 1;
      while (++k < indices.size())
      {
        indices[k] = 0;
      }
    }
    result.sortByIntensityAndNormalize_();
    return result;
  }

  
  void Ecipex::convolve_(IsotopeDistribution& distribution, double threshold)
  {
      
    auto& p1 = *this;
    auto& p2 = distribution;
    
    if(not p1.isNormalized() || not p2.isNormalized())
    {
      throw Exception::InvalidValue(__FILE__, 
                                    __LINE__, 
                                    OPENMS_PRETTY_FUNCTION, 
                                    "Distributions not normalized ", 
                                    "");
    }


    if(p1.isConvolutionUnit())
    {
      p1.set(p2.getContainer());
      return;
    }

    if(p2.isConvolutionUnit())
    {
      return;
    }
    
    Ecipex result;
    // remove convolution unit vector
    result.clear();
    
    UInt64 n1 = p1.size(), n2 = p2.size();

    for(UInt64 i = 0; i < n1; i++)
    {

      for(UInt64 j = 0; j < n2; j++) 
      {
        auto abundance = p1[i].getIntensity() * p2[j].getIntensity();
        if (abundance > threshold) 
        {
          result.insert(p1[i].getMZ() + p2[j].getMZ(), abundance);
        } 
        else 
        {
          if (j == 0) 
          {
            break;
          }
          n2 = j;
        }
      }
    }
    result.sortByIntensityAndNormalize_();
    p1.set(result.getContainer());
  }

}
