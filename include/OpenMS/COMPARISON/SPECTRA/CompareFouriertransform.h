#ifndef COMPAREFOURIERTRANSFORM_H
#define COMPAREFOURIERTRANSFORM_H

#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>

namespace OpenMS
{
	/**
		@brief Compare Discrete Cosines value from a Fourier transformation, also known as Discrete Cosines Transformation
			
		The Direct Cosines Transformation based on the theory of the Fourier transformation.
		In this class the Fast Fourier Transformation(FFT) algorithm of the gsl library is used. FFT  has a run-time complexity of 
		n (log n). To get the Direct Cosines Transformation from a FFT there is preparation necessary. First the 
		input data have to be mirrored. This is necessary, because FFT needs data which have a periodic nature.
		After the computation from FFT only the cosines values are important and stored in the Meta Data Array. So a inverse from these values 
		to get the original spectrum is not available.
		The comparison is done between two Meta Data Arrays, which contain the stored cosines values of there individual spectrum. 
		The advantage of these method is how the comparison works. There is only a sum which have to be count, no multiplication is needed. 
		
		Attention only use the compare function, if the Spectrum was transformed earlier else a error is going to be appear.
		Only use this method of transformation, if you sure there exist enough free memory. 
	
	*/
	
  class CompareFouriertransform : public PeakSpectrumCompareFunctor
  {
  	public:
	
			// @name Constructors and Destructors
			// @{
	    /// default constructor
		  CompareFouriertransform();
	
	    /// copy constructor
		  CompareFouriertransform(const CompareFouriertransform& source);
	
	    /// destructor
	    virtual ~CompareFouriertransform();
			// @}
	
			// @name Operators
			// @{
	    /// assignment operator
	    CompareFouriertransform & operator = (const CompareFouriertransform & source);
		
		  /**
		  	@brief Dummy function,
				
				This function only return 0 for any given PeakSpectrum, please use the other compare operator function
		
		  	@param PeakSpectrum MSSpectrum 
		  	@see MapAlignmentAlgorithmSpectrumAlignment()
		  */
	    double operator () (const PeakSpectrum& )const;
	    /**
	    	@brief compare two MSSpectrums by their Discrete Cosines Transformation.
				
	  		This function compares two given MSSprectrums on their  Discrete Cosines Transformation
	  		First a transformation had to be calculated. Please use the function transform() in this class previously, befor calling this
				function. The comparison works by summing over all elements of both transformation, by subtracts each other coefficient. sum(_i=1)
				^n x_i-y_i. If the sum is zero, both Spectrum are identical in the real part, if the sum is not zero an diversion is done 1/|sum|, 					to get a similarity score.
	  		
	    	@param PeakSpectrum MSSpectrum 
	    	@see MapAlignmentAlgorithmSpectrumAlignment()
	    */
	    double operator () (const PeakSpectrum& spec1 , const PeakSpectrum& spec2 ) const;
	
	    static PeakSpectrumCompareFunctor* create() { return new CompareFouriertransform(); }

  		///
  		static const String getProductName()
  		{
  			return "CompareFouriertransform";
  		}
	protected:
			/**
			 	@brief Search in the MSSpectrum, if a Discrete Fourier transformation occurs, if not a error is going to be throw, else the index 				of the occurrence is returned.
				
				This function gives the position back, which position the transformation was saved in a MetaDataArray. If there is no entry, it 				throws a error, to indicate that first a transformation have do calculated before calling the comparison operator.
				
			 	@param spec  MSSpectrum 
			 	@see MapAlignmentAlgorithmSpectrumAlignment()
			*/
			UInt searchTransformation_(const PeakSpectrum&  spec) const;

			/**
	    	@brief calculate the Discrete Cosines Fourier Transformation.
	       				
	   		This Function transform a given MSSpectrum to an Discrete Cosines Fourier Transformation. It stores only the part of the cosines 					of the FFT in
	   		the MetaDataArray which is a container from the MSSpectrum. Only call this function, if you sure there is no earlier an another 				transformation done over the same MSSpectrum, because it doesn't check if there already exist a transformation.
	          		
	     	@param spec  MSSpectrum 
	     	@see MapAlignmentAlgorithmSpectrumAlignment()
	    */
      void transform_(PeakSpectrum & spec);
  };

}
#endif /*COMPAREFOURIERTRANSFORM_H*/


