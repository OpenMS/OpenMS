
#include <complex>
#include <cstring>

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/Ecipex.h>
#include <OpenMS/CHEMISTRY/Element.h>

#include <kiss_fft.h>
#include <kiss_fftndr.h>
#include <kiss_fftr.h>

using namespace std;

namespace OpenMS
{

  class KissFftState 
  {
    std::vector<int> dims_;
    bool inverse_;
    void* cfg_;
    bool is1d() const { return dims_.size() == 1; }

   public:
    KissFftState(const std::vector<int>& dimensions, bool inverse = false)
      : dims_(dimensions), inverse_(inverse) {
      if (is1d())
        cfg_ = kiss_fftr_alloc(dims_[0], inverse, nullptr, nullptr);
      else
        cfg_ = kiss_fftndr_alloc(&dims_[0], dims_.size(), inverse, nullptr, nullptr);
      if (cfg_ == nullptr) throw std::bad_alloc();
    }

    KissFftState(const KissFftState&) = delete;
    KissFftState operator=(const KissFftState&) = delete;
    ~KissFftState() { free(cfg_); }

    void runFFT(void* in, void* out) {
      if (is1d()) {  // for some reason kiss_fftndr needs at least 2 dimensions to
        // work
        if (!inverse_)
          kiss_fftr((kiss_fftr_cfg)cfg_, (kiss_fft_scalar*)in, (kiss_fft_cpx*)out);
        else
          kiss_fftri((kiss_fftr_cfg)cfg_, (kiss_fft_cpx*)in, (kiss_fft_scalar*)out);
      } else {
        if (!inverse_)
          kiss_fftndr((kiss_fftndr_cfg)cfg_, (kiss_fft_scalar*)in, (kiss_fft_cpx*)out);
        else
          kiss_fftndri((kiss_fftndr_cfg)cfg_, (kiss_fft_cpx*)in, (kiss_fft_scalar*)out);
      }
    }
  };


  class FftArray 
  {
    vector<Int> dims_;
    complex<kiss_fft_scalar>* data_;
    UInt64 n_;

   public:
    FftArray(const std::vector<Int>& dimensions)
      : dims_(dimensions), data_(nullptr), n_(1) {
      for (Int x : dims_)
        n_ *= x;

      if (n_ > 1e7) throw std::runtime_error("too many isotopic combinations");

      UInt64 n_bytes = n_ * sizeof(kiss_fft_cpx);

      data_ = reinterpret_cast<std::complex<kiss_fft_scalar>*>(KISS_FFT_MALLOC(n_bytes));
      if (data_ == nullptr) throw std::bad_alloc();
      std::memset(data_, 0, n_bytes);
    }

    void forwardFFT() {
      KissFftState state(dims_, false);
      state.runFFT(scalar_data(), complex_data());
    }

    void inverseFFT() {
      KissFftState state(dims_, true);
      state.runFFT(complex_data(), scalar_data());
    }

    kiss_fft_cpx* complex_data() const { return reinterpret_cast<kiss_fft_cpx*>(data_); }

    kiss_fft_scalar* scalar_data() const {
      return reinterpret_cast<kiss_fft_scalar*>(data_);
    }

    std::complex<kiss_fft_scalar>* data() const { return data_; }
    size_t size() const { return n_; }

    ~FftArray() { kiss_fftr_free(data_); }
  };

  Ecipex::Ecipex() : MIDAs()
  {}

  Ecipex::Ecipex(const IsotopeDistribution& isotope_distribution) :
    MIDAs(isotope_distribution)
  {
  }
  
  void Ecipex::run()
  {}

  void Ecipex::computeIsotopePattern(double threshold, double fft_threshold) 
  {
    Ecipex result;
    for (auto& element : formula_) 
    {
      auto pattern = elementIsotopePattern(
        element.first->getIsotopeDistribution().getContainer(), 
        element.second, 
        fft_threshold);
      
      distribution_ = result.convolve(pattern, threshold * threshold);
    }
    result.trimIntensities(threshold);
    return;
  }

  void Ecipex::sortAndNormalize()
  {
    sortByIntensity();
    const auto& max_intensity = distribution_.front().second;
    transform_([&max_intensity](MassAbundance& m){
        m.second /= max_intensity;
      });
  }

  Ecipex::Spectrum  Ecipex::elementIsotopePattern(const Spectrum& iso_pattern, UInt amount, double threshold)
  {
    UInt64 dim = iso_pattern.size() - 1;
    if (dim == 0) 
    {
      return Spectrum(1, make_pair(iso_pattern.front().first * amount, 1));
    }

    UInt64 edge_len = kiss_fftr_next_fast_size_real(amount + 1);

    vector<int> dimensions(int(dim), edge_len);
    FftArray arr{dimensions};

    // array setup
    const auto& iso = iso_pattern;
    arr.scalar_data()[0] = iso.front().second;
    for (UInt64 i = 0, k = 1; i < dim; i++, k *= edge_len)
    {
      arr.scalar_data()[k] = iso[dim-1].second;
    }
  
      // forward FFT
      arr.forwardFFT();

      // exponentiation
      for (size_t i = 0; i < arr.size(); ++i)
      {
        arr.data()[i] = std::pow(arr.data()[i], amount);
      }
      // inverse FFT
      arr.inverseFFT();
      for (size_t i = 0; i < arr.size(); ++i)
      {
        arr.scalar_data()[i] /= double(arr.size());  // take care of FFT normalization
      }
  
      Spectrum isotope_pattern;

      vector<UInt64> indices(dim, 0);
      for (UInt64 i = 0; i < arr.size(); ++i) {
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

          isotope_pattern.push_back(make_pair(mass,arr.scalar_data()[i]));
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
      Ecipex result;
      result.set(iso_pattern);
      result.sortAndNormalize();
      return result.getContainer();
    }

    Ecipex::Spectrum Ecipex::convolve(const Spectrum& distribution, double threshold)
    {
    
    }


  }

