#define MAX_ISOTOPIC_COMBINATIONS 1e7 

#include <complex>
#include <cstring>
#include <cassert>
#include <iostream>

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/Ecipex.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <kiss_fft.h>
#include <kiss_fftndr.h>
#include <kiss_fftr.h>

using namespace std;

namespace OpenMS
{
  void print(std::string t, IsotopeDistribution& result)
  {
    for(auto& sample: result.data())
    {
      cout<<t << sample.first <<" " << sample.second << endl;
    }
  }

  class KissFftState 
  {

   public:
    
    KissFftState(const std::vector<Int>& dimensions, bool inverse = false)
      : dims_(dimensions), inverse_(inverse) 
    {
      if (is1d())
      {
        cfg_ = kiss_fftr_alloc(dims_[0], inverse, nullptr, nullptr);
      }
      else
      {
        cfg_ = kiss_fftndr_alloc(&dims_[0], dims_.size(), inverse, nullptr, nullptr);
      }
      if (cfg_ == nullptr) 
      {
        throw Exception::NullPointer(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }
    }

    KissFftState(const KissFftState&) = delete;
    KissFftState operator=(const KissFftState&) = delete;
    
    ~KissFftState() 
    { 
      free(cfg_); 
    }

    void runFFT(void* in, void* out) 
    {
      if (is1d()) 
      { 
        // for some reason kiss_fftndr needs at least 2 dimensions to work
        if (!inverse_)
        {
          kiss_fftr((kiss_fftr_cfg)cfg_, (kiss_fft_scalar*)in, (kiss_fft_cpx*)out);
        }
        else
        {
          kiss_fftri((kiss_fftr_cfg)cfg_, (kiss_fft_cpx*)in, (kiss_fft_scalar*)out);
        }
      } 
      else 
      {
        if (!inverse_)
        {
          kiss_fftndr((kiss_fftndr_cfg)cfg_, (kiss_fft_scalar*)in, (kiss_fft_cpx*)out);
        }
        else
        {
          kiss_fftndri((kiss_fftndr_cfg)cfg_, (kiss_fft_cpx*)in, (kiss_fft_scalar*)out);
        }
      }
    }

   private:
    std::vector<int> dims_;
    bool inverse_;
    void* cfg_;
    bool is1d() const { return dims_.size() == 1; }
  };


  class FftArray 
  {

   public:
    FftArray(const std::vector<Int>& dimensions)
      : dims_(dimensions), data_(nullptr), n_(1) 
    {
      for (const Int& x : dims_)
      {
        n_ *= x;
      }

      if (n_ > MAX_ISOTOPIC_COMBINATIONS)
      {
        throw Exception::InvalidValue(__FILE__, 
                                      __LINE__, 
                                      OPENMS_PRETTY_FUNCTION, 
                                      "Too many isotopic combinations", 
                                      String(n_));
      }

      UInt64 n_bytes = n_ * sizeof(kiss_fft_cpx);

      data_ = reinterpret_cast<std::complex<kiss_fft_scalar>*>(KISS_FFT_MALLOC(n_bytes));
      
      if (data_ == nullptr)
      {
        throw Exception::NullPointer(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }

      std::memset(data_, 0, n_bytes);
    }

    void forwardFFT() 
    {
      KissFftState state(dims_, false);
      state.runFFT(scalar_data(), complex_data());
    }

    void inverseFFT()
    {
      KissFftState state(dims_, true);
      state.runFFT(complex_data(), scalar_data());
    }

    kiss_fft_cpx* complex_data() const 
    { 
      return reinterpret_cast<kiss_fft_cpx*>(data_); 
    }

    kiss_fft_scalar* scalar_data() const 
    {
      return reinterpret_cast<kiss_fft_scalar*>(data_);
    }

    std::complex<kiss_fft_scalar>* data() const 
    { 
      return data_; 
    }

    UInt64 size() const 
    { 
      return n_; 
    }

    ~FftArray() 
    { 
      kiss_fftr_free(data_); 
    }

   private:
    vector<Int> dims_;
    complex<kiss_fft_scalar>* data_;
    UInt64 n_;

  };

  Ecipex::Ecipex() : IsotopeDistribution()
  {
  }

  Ecipex::Ecipex(EmpiricalFormula& formula, double threshold, double fft_threshold):
    formula_(formula),
    fft_threshold_(fft_threshold),
    threshold_(threshold)
  {
  }

  Ecipex::Ecipex(const IsotopeDistribution& isotope_distribution) :
    IsotopeDistribution(isotope_distribution)
  {
  }
  
  void Ecipex::run()
  {
    computeIsotopePattern(threshold_, fft_threshold_);
  }

  void Ecipex::computeIsotopePattern(double threshold, double fft_threshold) 
  {
    for (auto& element : formula_) 
    {
      auto pattern = elementIsotopePattern(
        element.first->getIsotopeDistribution().getContainer(), 
        element.second, 
        fft_threshold);
      
      convolve(pattern, threshold * threshold);
    }
    trimIntensities(threshold);
  }

  void Ecipex::sortAndNormalize()
  {
    sortByIntensity();
    double max_intensity = distribution_.front().second;
    transform_([&max_intensity](MassAbundance& m){
        m.second /= max_intensity;
      });
    
  }

  Ecipex Ecipex::elementIsotopePattern(const Spectrum& iso_pattern, UInt amount, double threshold)
  {
    Int dim = iso_pattern.size() - 1;
    if (dim == 0) 
    {
      Ecipex result;
      result[0].mw = iso_pattern.front().first * amount;
      return result;
    }

    if(dim < 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid dimension size", String(dim));
    }

    

    UInt64 edge_len = kiss_fftr_next_fast_size_real(amount + 1);

    vector<Int> dimensions(dim, edge_len);
    FftArray arr{dimensions};

    // array setup
    const auto& iso = iso_pattern;
    arr.scalar_data()[0] = iso.front().second;
    for (Int i = 0, k = 1; i < dim; i++, k *= edge_len)
    {
      arr.scalar_data()[k] = iso[dim-i].second;
    }
  
    // forward FFT
    arr.forwardFFT();

    // exponentiation
    for (UInt64 i = 0; i < arr.size(); ++i)
    {
      arr.data()[i] = std::pow(arr.data()[i], amount);
    }
    // inverse FFT
    arr.inverseFFT();
    for (UInt64 i = 0; i < arr.size(); ++i)
    {
      arr.scalar_data()[i] /= double(arr.size());  // take care of FFT normalization
    }
  
    Ecipex result;
    result.data().clear();

    vector<UInt64> indices(dim, 0);
    for (UInt64 i = 0; i < arr.size(); ++i) 
    {
      UInt64 k;
      UInt64 n = 0;
      double mass = 0.0;
      for (k = 0; k < indices.size(); ++k) 
      {
        mass += iso_pattern[k + 1].first * indices[k];
        n += indices[k];
      }

      if (n <= amount && arr.scalar_data()[i] >= threshold) 
      {
        mass += (amount - n) * iso_pattern[0].first;

        result.data().push_back(make_pair(mass,arr.scalar_data()[i]));
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
    result.sortAndNormalize();
    return result;
  }

  
  void Ecipex::convolve(IsotopeDistribution& distribution, double threshold)
  {
      
    auto& p1 = *this;
    auto& p2 = distribution;
      
    


    
    if(!p1.isNormalized() || !p2.isNormalized())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Distributions not normalized ", "");
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
    result.data().clear();
    
    UInt64 n1 = p1.size(), n2 = p2.size();

    for(UInt64 i = 0; i < n1; i++)
      for(UInt64 j = 0; j < n2; j++) 
      {
        auto abundance = p1[i].intensity * p2[j].intensity;
        if (abundance > threshold) 
        {
          result.data().push_back(make_pair(p1[i].mw + p2[j].mw, abundance));
        } 
        else 
        {
          if (j == 0) break;
          n2 = j;
        }
      }
    result.sortAndNormalize();
    p1.set(result.getContainer());
  }


}

