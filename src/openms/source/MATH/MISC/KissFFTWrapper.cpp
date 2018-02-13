#include <OpenMS/MATH/MISC/KissFFTWrapper.h>

using namespace std;

namespace OpenMS
{
  KissFftState::KissFftState(const vector<Int>& dimensions, bool inverse)
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

  KissFftState::~KissFftState() 
  { 
    free(cfg_); 
  }

  void KissFftState::runFFT(void* in, void* out) 
  {
    if (is1d()) 
    { 
      // for some reason kiss_fftndr needs at least 2 dimensions to work
      if (not inverse_)
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

  bool KissFftState::is1d() const 
  { 
    return dims_.size() == 1; 
  }


  FftArray::FftArray(const vector<Int>& dimensions)
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

    data_ = reinterpret_cast<complex<kiss_fft_scalar>*>(KISS_FFT_MALLOC(n_bytes));
      
    if (data_ == nullptr)
    {
      throw Exception::NullPointer(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    std::memset(data_, 0, n_bytes);
  }

  void FftArray::forwardFFT() 
  {
    KissFftState state(dims_, false);
    state.runFFT(scalar_data(), complex_data());
  }

  void FftArray::inverseFFT()
  {
    KissFftState state(dims_, true);
    state.runFFT(complex_data(), scalar_data());
  }

  FftArray::~FftArray() 
  { 
    kiss_fftr_free(data_); 
  }

}
