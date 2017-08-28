#define kiss_fft_scalar double 
#define MAX_ISOTOPIC_COMBINATIONS 1e10

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <kiss_fft.h>
#include <kiss_fftndr.h>
#include <kiss_fftr.h>

#include <complex>
#include <vector>
#include <cstring>

namespace OpenMS
{

  class KissFftState 
  {
 public:
    
    KissFftState(const std::vector<Int>& dimensions, bool inverse = false);

    KissFftState(const KissFftState&) = delete;
    KissFftState operator=(const KissFftState&) = delete;
    
    ~KissFftState();

    void runFFT(void* in, void* out);

 private:
    std::vector<Int> dims_;
    bool inverse_;
    void* cfg_;
    bool is1d() const;
  };


  class FftArray 
  {

 public:
    explicit FftArray(const std::vector<Int>& dimensions);
  
    void forwardFFT(); 
    
    void inverseFFT();

    inline kiss_fft_cpx* complex_data() const 
    { 
      return reinterpret_cast<kiss_fft_cpx*>(data_); 
    }

    inline kiss_fft_scalar* scalar_data() const 
    {
      return reinterpret_cast<kiss_fft_scalar*>(data_);
    }

    inline std::complex<kiss_fft_scalar>* data() const 
    { 
      return data_; 
    }

    inline UInt64 size() const 
    { 
      return n_; 
    }

    ~FftArray();

 private:
    std::vector<Int> dims_;
    std::complex<kiss_fft_scalar>* data_;
    UInt64 n_;

  };

}
