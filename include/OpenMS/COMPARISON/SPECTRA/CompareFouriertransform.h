#ifndef COMPAREFOURIERTRANSFORM_H
#define COMPAREFOURIERTRANSFORM_H

#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>

namespace OpenMS
{

  /**
	  @brief Improvement Similarity score of Stein & Scott 

		The details of the score can be found in:
		Signal Maps for Mass Spectrometry-based
		Comparative Proteomics

		

		@ref SteinScottImproveScore_Parameters are explained on a separate page.
		
		
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
	
		/// 
    double operator () (const PeakSpectrum& ) const{return 0;}
    double operator () (const PeakSpectrum& spec1 , const PeakSpectrum& spec2 ) const{
    
    	const DSpectrum<>::MetaDataArrays& temp1 = spec1.getMetaDataArrays();
			
			if(temp1.size()== 0)
			{
				throw "Input needed to be a fouriertransform try first transform()";
				//transform(spec1); 
			}
			
			UInt i=	soc(spec1);

			const DSpectrum<>::MetaDataArrays& temp2 = spec2.getMetaDataArrays();
			if(temp2.size()== 0)
			{
				throw "Second Input needet be a fouriertransfom try first transform of these class";
			}
			UInt j=	soc(spec2);
			if(temp1[i].size() != temp2[j].size())
			{
				std::cout<< temp1[i].size() << temp2[j].size() << std::endl; 
				return 0.0;
			}
			else
			{	Real sum=0;
				for(UInt k=0; k< temp1[i].size();++k)
				{
					std::cout<< temp1[i][k] << " temp1 "<< temp2[j][k] << " temp2 " << std::endl; 
					sum=sum + temp1[i][k]-temp2[j][k];
				}
				std::cout << sum << " summe " << std::endl;
				if(sum !=0)
				{
					return std::abs(1/sum);
				}
				else
					return 1;
			}
		}
	
    UInt soc(const PeakSpectrum&  spec) const
    {
		const DSpectrum<>::MetaDataArrays& temp = spec.getMetaDataArrays();
		UInt i=0;
		while(i< temp.size())
		{
			if(temp[i].getName()=="Fouriertransformation")
			{
				break;
			}
			else
				++i;
		}
		//lesen oder erstellen
		if(i < temp.size() && temp[i].getName()!= "Fouriertransformation")
		{
			 throw " Input needed be a fouriertransfom try first transform of these class";
		
		}
		else
		{
			return i;
		}
	
    }
    template<typename T>
    void transform(T & spec)
    {
    	//leider muß eine Kopie gemacht werden...
    	double* data=  new double [spec.getContainer().size()<<1];
    	bool aflag = false;
    	bool iflag = true;
    	UInt i =0;
    	//spectrum in das array kopieren und duplizieren damit FFT durchgeführt werden kann(perodische funktion)
    	while(!aflag)
    	{	
    		if(i== (spec.getContainer().size()<<1))	//wann muss ich aufhören
    		{
    			aflag=true;
    			break;
    		}
    		if(iflag)//vorderes ende
    		{
    			if(i< spec.getContainer().size())
    			{
    				data[i]= spec.getContainer()[i].getIntensity();
    				++i;
    			}
    			else//spiegeln
    				{
    					iflag= false;
    				}
    		}
    		else//spiegelung
    			{
    				data[i] =spec.getContainer()[(spec.getContainer().size()<<1)-i].getIntensity();
    				++i;
    			}
    	}
    	
    	gsl_fft_real_wavetable * real;
    	gsl_fft_real_workspace * work;
    	work = gsl_fft_real_workspace_alloc (spec.getContainer().size());
    	real = gsl_fft_real_wavetable_alloc (spec.getContainer().size());
    	gsl_fft_real_transform (data,1,spec.getContainer().size(),real, work);
    	gsl_fft_real_wavetable_free (real);
    	gsl_fft_real_workspace_free (work);
    	
    	DSpectrum<>::MetaDataArrays& temp = spec.getMetaDataArrays();
    	i= temp.size();
    	temp.resize(i+1);
    	temp[i].setName("Fouriertransformation");
    	UInt j=0;
    	while(j < spec.getContainer().size())
    	{
    		temp[i].push_back(data[j]);
    		if(j==0) ++j;
    		else j=j+2;
    	}
    	
    	delete[] data;
    }
    
	
    	// @}

		// @name Accessors
		// @{
		///
  static PeakSpectrumCompareFunctor* create() { return new CompareFouriertransform(); }

		///
		static const String getProductName()
		{
			return "CompareFouriertransform";
		}

		// @}

  };

}



#endif /*COMPAREFOURIERTRANSFORM_H*/

