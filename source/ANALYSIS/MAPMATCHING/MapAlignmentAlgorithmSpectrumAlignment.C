// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Vipul Patel $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>

namespace OpenMS
{

	MapAlignmentAlgorithmSpectrumAlignment::MapAlignmentAlgorithmSpectrumAlignment()
	:MapAlignmentAlgorithm(),ProgressLogger()
	{
		setName("MapAlignmentAlgorithmSpectrumAlignment");
		defaults_.setValue("gapcost",2," This is the cost for a open or close gap in the alignment, the input number has to be postive!.",false);
		defaults_.setValue("affinegapcost", 1,"This is the extension cost for a gap. The extenscion cost ist always smaller than the normal gapcost. and have to"
		"fullfill the equation minimum maxscore > 2*extensioncost",false);
		defaults_.setValue("cutoffScore",0.90,"If a score higher or same found, it will be a candidate for pair of a data point for the alignment intervall. (Stepp preAlign)",true);
		defaults_.setValue("scorefunction","SteinScottImproveScore","Sore function for the alignmnent of MSSpectrums",false);
		defaults_.setValidStrings("scorefunction",Factory<PeakSpectrumCompareFunctor>::registeredProducts());
		setLogType(CMD);
		defaultsToParam_();
	}

	MapAlignmentAlgorithmSpectrumAlignment::~MapAlignmentAlgorithmSpectrumAlignment()
	{
	}
	
	void MapAlignmentAlgorithmSpectrumAlignment::alignPeakMaps(std::vector< MSExperiment<> >& peakmaps, std::vector<TransformationDescription>& transformations)
	{
		try
		{ 	
			std::vector<MSSpectrum<>* >spectrum_pointers;			
			msFilter_(peakmaps[0],spectrum_pointers);
			fourierActivation_(spectrum_pointers);

		
			startProgress(0,(peakmaps.size()-1),"alignment");
			for(UInt i = 1 ; i < peakmaps.size();++i )
			{									
				prepareAlign_(spectrum_pointers,peakmaps[i]);
				setProgress(i);
				peakmaps[i].clearMetaDataArrays();
			}
			peakmaps[0].clearMetaDataArrays();
			endProgress();
		}
		catch (Exception::OutOfRange& e) 
		{
			throw Exception::OutOfRange(__FILE__,__LINE__,__PRETTY_FUNCTION__);	
		}
	}
	/**
	 @brief A function to prepare the sequence for the alignment. It calls intern the main function for the alignment.
     
		This function takes two arguments. These Arguments type are two MSExperiments. 
		The first argument should have been filtered, so such only type of MSLevel 1 exist in the Sequence. 
		The second arguments doesn't have to full fill these restriction. It's going to be filter automatically. 
		With these two arguments a precalculation is done, to find some corresponding data points(maximum 4) for building alignments blocks.
		After the alignment a retransformation is done, the new Retention Times appears in the original data.      
    
    The parameters are MSExperiments.
    @param pattern is the template MSExperiment.
    @param aligned is the sequence which has to be align(also MSExperiment).
    @see MapAlignmentAlgorithmSpecturmAlignment()
	*/	

	void MapAlignmentAlgorithmSpectrumAlignment::prepareAlign_(const std::vector<MSSpectrum<>* >& pattern, MSExperiment<>& aligned)
	{
		//tempalign ->container for holding only MSSpectrums with MS-Level 1
		std::vector<MSSpectrum<>*> tempalign;
		msFilter_(aligned,tempalign);
		fourierActivation_(tempalign);

		//if it's possible, built 4 blocks. These can be individually be aligned.						
		std::vector<UInt> alignpoint;
		//saving the first cordinates
		alignpoint.push_back(0); 
		alignpoint.push_back(0);
		//4 blocks : 0-0.25 ,0.25-50,0.50-0.75,1 The data points must have a high similarity score
		for(Real i = 0.25; i<=0.75;i+=0.25)
		{	
			UInt y= (UInt)(tempalign.size() * i);
			UInt x=0;
			Real maxi=-999.0;
			
			for(UInt k = 0; k<pattern.size(); ++k)
			{
				Real s =	scoring_(*pattern[k],*(tempalign[y]));
				if(s > maxi && s > cutoffScore_)
				{
					x=k;
					maxi=s;
				}			
			}
			if(x >= alignpoint[alignpoint.size()-2]+3 && y >= alignpoint[alignpoint.size()-1]+3)
			{
							alignpoint.push_back(x);
							alignpoint.push_back(y);
			}
			
			UInt xn=(UInt)(pattern.size() * i);
			UInt yn=0;
			for(UInt k = 0; k<tempalign.size(); ++k)
			{
				Real s =	scoring_(*pattern[xn],*(tempalign[k]));
				if(s > maxi && s > cutoffScore_)
				{
					yn=k;
					maxi=s;
				}			
			}
			if(xn >= alignpoint[alignpoint.size()-2]+3 && yn >= alignpoint[alignpoint.size()-1]+3)
						{
										alignpoint.push_back(xn);
										alignpoint.push_back(yn);
						}
			//only save possible data points, if they are not already contained 
			
		}

		//save also the endpoint as a data point
		alignpoint.push_back(pattern.size()-1);
		alignpoint.push_back(tempalign.size()-1);
		//the distance of two data points have to be greater than 3, if not the spline would thrown an Expection
		//do a affine gap alignment of the block of the data points x1,y1,x2,y2

		std::vector<int> xcordinate;
		std::vector<double>ycordinate;
		for(UInt i=0; i < alignpoint.size()-2; i+=2)
		{
			affineGapalign_(alignpoint[i],alignpoint[i+1],alignpoint[i+2],alignpoint[i+3], pattern, tempalign,xcordinate,ycordinate);
		}
	/*
		for(UInt i = 0; i< xcordinate.size(); ++i)
		{
			std::cout<< xcordinate[i] << " " << ycordinate[i] << " x  y  ankerpunkte " << std::endl;
		}*/
		//calculate the spline
		calculateSpline_(xcordinate,ycordinate,tempalign,(UInt)0,tempalign.size()-1);
		
	}
	/**
      @brief affine gap cost Alignment

		This Alignment based on Needleman Wunsch Algorithm. 
		To improve the time complexity a banded version was implemented, known as k - alignment. 
		To save some space, the alignment is going to be calculated by position xbegin to xend of one sequence and ybegin 
		and yend by another given sequence. The result of the alignment is stored in the second argument. 
		The first sequence is used to be a template for the alignment.
		
     	 
   	@param xbegin UInt cordinate for the beginning of the template sequence.
   	@param ybegin UInt cordinate for the beginning of the aligend sequence .
   	@param xend UInt cordinate for the end of the template sequence.
   	@param yend UInt cordinate for the end of the aligend sequence.
   	@param pattern is the template MSExperiment.
   	@param aligned is to be align MSExperiment.
   	@see MapAlignmentAlgorithmSpecturmAlignment()
	 */
		void MapAlignmentAlgorithmSpectrumAlignment::affineGapalign_(UInt xbegin, UInt ybegin, UInt xend,UInt yend, const std::vector<MSSpectrum<>* >& pattern,  std::vector<MSSpectrum<>* >& aligned,std::vector<int>& xcordinate, std::vector<double>&ycordinate)
		{	
			//affine gap alignment needs two matrices
		  std::map<UInt, std::map<UInt,Real> > matchmatrix;
		  std::map<UInt, std::map<UInt,Real> > insertmatrix;	
		  //setting the size of column and row	
			
		  UInt n= std::max((xend-xbegin),(yend-ybegin))+1; //column 
		  UInt m= std::min((xend-xbegin),(yend-ybegin))+1; //row
		  //std::cout<< n << " n " << m << " m " << std::endl;
		  Real score_=-999;
		  Int k_=1;
		  
		  bool column_row_orientation= false; 
		  if(n !=(xend-xbegin)+1)
		  {
		  		column_row_orientation = true;
		  }
			//matrix for holding calculated sorces
			std::map<UInt, std::map<UInt,Real> > buffermatrix;
			//calculate the value of k
			k_=	bestk_(pattern, aligned, buffermatrix, column_row_orientation,xbegin,xend, ybegin, yend);
	
			//flag if we have to calculate again the alignment in step k+1
			bool finish = false;
			try
			{
			  while(!finish && k_ <(Int) n+1)	 //if the current score is smaller than the next step -> calculate again
			  {																//best k+1 =(2(k+1)+m-n)d+(n(k+1))matchscore 
																				//if k >= columnsize -> there is no band necessary(normal global alignment)
			  	std::map<UInt,std::map<UInt,UInt> > traceback;
		      //init both  matrices
			  	matchmatrix[0][0]=0.0;
			  	insertmatrix[0][0]=-999.0;
			  	//init column
			  	for(Int i = 1; i<= (Int)(k_+n-m);++i )
			  	{	
			  		Real w = -gap_- ((i-1) * e_);
			  		matchmatrix[i][0]=-999.0;
			  		insertmatrix[i][0]=w;
			  	}
			  	///init row				
			  	for(Int j = 1; j<=(Int) k_ && j<=(Int)m;++j )
			  	{	
			  		Real w = -gap_-(j-1)* e_;
			  		matchmatrix[0][j]=-999.0;
			  		insertmatrix[0][j]=w;
			  	}
			  	for(UInt i = 1; i<= n; ++i)
			  	{	//->band-size
			  		for(Int h= -k_-(n-m); h <= k_ ; ++h) // from -k to k
				    {
			  			//j->row
			  			Int j=i+h;
			  			if(j >= 1 && (UInt)j<=m )
			  			{	
			  				double s=-999.0;
			  				//score calculation, attention column and row a 1 size smaller
			  				try
			  				{
			  					s	=	scoreCalculation_(i,j,xbegin,ybegin,pattern,aligned,buffermatrix ,column_row_orientation);
			  				}
			  				catch (Exception::OutOfRange &e) 
					  		{
			  					throw Exception::OutOfRange(__FILE__,__LINE__,__PRETTY_FUNCTION__);
					  		}
			  				//M(i,j) = max{M(i-1,j-1)+s,I(i-1,j-1)+s}
			  				if(matchmatrix[i-1][j-1] > insertmatrix[i-1][j-1])
			  				{
			  					matchmatrix[i][j]=matchmatrix[i-1][j-1]+s;
			  				}
			  				else
			  				{
			  					matchmatrix[i][j]=insertmatrix[i-1][j-1]+s;
			  				}
			  				Real tempmaximum = matchmatrix[i][j];
			  				//std::cout <<  i << " i " << j << " j " << s << "  Score" <<std::endl;
			  				//calculation of max{M(i-1)-gap_,I(i-1,j)-e_,M(i,j-1)-gap_,I(i,j-1)-e_}
			  				
			  				Real maxima[]={-999.0,-999.0};
			  				//testing if i-1,j is in inside of the band
			  				//std::cout << i << " i " << j << " j " << s <<  "score" << std::endl;
			  				i=i-1;
			  				if(insideBand_(i,j,n,m,k_))
			  				{ 	
			  					i+=1;
			  					Real tempm=matchmatrix[i-1][j]-gap_;
			  					Real tempi=insertmatrix[i-1][j]-e_;
			  					//std::cout << i-1 <<  " i-1 " << j << " j " <<tempm << " Matchmatrix " << tempi << " Insertmatrix " << std::endl;
			  					if(tempm > tempi)
			  					{
			  						maxima[0]=tempm;
			  					}
			  					else
			  					{
			  						maxima[0]=tempi;
			  					}										
			  				}
			  				else	i=i+1;
			  				//testing if i,j-1 is in inside of the band
			  				j=j-1;
			  				if(insideBand_(i,j,n,m,k_))
			  				{	
			  					j+=1;
			  					Real tempm=matchmatrix[i][j-1]-gap_;
			  					Real tempi=insertmatrix[i][j-1]-e_;
			  					//std::cout <<  i <<  " i " << j-1 << " j-1 " << tempm << " Matchmatrix " << tempi << " Insertmatrix " << std::endl;			
			  					if(tempm > tempi)
					      	{
			  						maxima[1]=tempm;
					      	}
			  					else
			  					{
			  						maxima[1]=tempi;
			  					}
			  				}
			  				else j=j+1;
			  				//std::cout << maxima[0]<< " maxima " << maxima[1] << " maxima1" << std::endl;
			  				//std::cout<< std::endl;
			  				if(maxima[0]!=-999.0 || maxima[1]!=-999.0)
			  				{ 
			  					if(maxima[0]>=maxima[1])
			  					{
			  						insertmatrix[i][j]= maxima[0];
			  						if(maxima[0]> tempmaximum)
			  						{
			  							traceback[i][j]= 1;
			  						}
			  					}
			  					else
			  					{
			  						insertmatrix[i][j]= maxima[1];
			  						if(maxima[1]>tempmaximum)
			  						{
			  							traceback[i][j]=2;
			  						}
			  					}
			  				}
			  			}
				    }
			  	}
			  	//std::cout << score_ << " current score " << matchmatrix[n][m] << " new score " <<k_ <<" bandsize "<<  n<< " n " << m << " m " << n <<" dimension of matrix "<< (m-(k_+1))*0.90 -2* gap_ - 2*(k_)*e_ -((Int)n-(Int)m) << "nextscore by extension" <<std::endl;
			  	// hier abstand berechnen (m-(k_+1))*0.9 -2* gap_ - 2*(k_)*e_ -((Int)n-(Int)m)
			  	//only do the traceback, if the actual score is not higher as the old one (m-(k_+1))*0.9 -2* gap_ - (2*(k_+1) +((Int)n-(Int)m))*e
			  	if(score_ >= matchmatrix[n][m] || k_<<1 > (Int)n || matchmatrix[n][m] >= (m-(k_+1))*cutoffScore_ -2* gap_ - (2*(k_) +((Int)m-(Int)n))*e_) //
			  	{
			  		finish= true;
			  		//traceback
			  		bool endtraceback= false;
			  		int i= n;
			  		int j= m;
			  		Real maximum = -999.0;
			  		std::vector<int> xvar;
			  		std::vector<double> yvar;
			  		while(!endtraceback)
			  		{ 
			  			if(i>0 && j>0)
			  			{
			  				Int t = traceback[i][j];								
			  				if(t==0)
			  				{	//std::cout<< column_row_orientation << std::endl;
			  					try
			  					{
			  						if(!column_row_orientation)
			  						{	
			  							/*if(maximum > -999.0)
											{
												xvar.push_back(x);
												yvar.push_back(y);
											}*/
			  							xvar.push_back(j+ybegin-1);
			  							yvar.push_back((*pattern[i+xbegin-1]).getRT());
			  							//aligned[j-1].setRT(pattern[i-1].getRT());
			  						}
			  						else
										{ 	
			  							/*if(maximum > -999.0)
												{
													xvar.push_back(x);
													yvar.push_back(y);
												}*/
												xvar.push_back(i+ybegin-1);
												yvar.push_back((*pattern[j+xbegin-1]).getRT());	
												//aligned[i-1].setRT(pattern[j-1].getRT());
										}
			  						i=i-1;
			  						j=j-1;
			  						maximum = -999.0;
			  					}
			  					catch (Exception::OutOfRange e) 
			  					{
			  						throw Exception::OutOfRange(__FILE__,__LINE__,__PRETTY_FUNCTION__);
			  					}
			  				}
			  				else if( t==1)
								{/*
									if(column_row_orientation)
									{
										if(maximum < matchmatrix[i][j])
										{
											maximum = matchmatrix[i][j];
											x=i-1;
											y=pattern[j-1].getRT();
										}
									}*/
									i=i-1;
								}
								else
								{/*
									if(!column_row_orientation)
									{
										if(maximum < matchmatrix[i][j])
										{
											maximum = matchmatrix[i][j];
											x=j-1;
											y=pattern[i-1].getRT();
										}
									}*/
									j=j-1;
								}
			  			}
			  			else 
			  			{
			  				endtraceback=true;
			  				/*if(maximum > -999.0)
								{
									xvar.push_back(x);
									yvar.push_back(y);
								}*/
			  			}
			  		}
			  		//call the spline on intervall
			  		//csp_(xvar,yvar,aligned,ybegin,yend);
			  		
			   		for(UInt i=0; i <xvar.size(); ++i)
			  		{
			   			
			   			//std::cout<< xvar[xvar.size()-1-i] << " " << std::endl;
			  			xcordinate.push_back(xvar[xvar.size()-1-i]);
			  			ycordinate.push_back(yvar[yvar.size()-1-i]);
			  		}
			   		//std::cout<< xcordinate.size()<< std::endl;
			  	}
			  	score_ = matchmatrix[n][m];
			  	k_ = k_ <<1;
			  }
			}
			catch (Exception::OutOfRange& e) 
			{
				throw Exception::OutOfRange(__FILE__,__LINE__,__PRETTY_FUNCTION__);
			}
		}
		/**
	    @brief filtered the MSLevel to gain only MSLevel 1 
			
		 	The alignment works only on MSLevel 1 data, so a filter has to be run.
		 
			@param peakmap MSExperiment is the sequence which has to be filtered
   		@param spectrum_pointer_container std::vector<MSSpectrum<>* > output container, where pointers of the MSSpectrum are saved(only with MSLevel 1 type)
			
			@exception Exception::IllegalArgument is thrown if no spectra are contained in @p peakmap
		*/
		void MapAlignmentAlgorithmSpectrumAlignment::msFilter_(MSExperiment<>& peakmap,std::vector<MSSpectrum<>* >& spectrum_pointer_container)
		{
			std::vector<UInt> pattern;
			peakmap.updateRanges(-1);
			pattern=peakmap.getMSLevels();
			
			if(pattern.size()!=0)
			{	
				for(UInt i=0; i< peakmap.size();++i)
				{
					if(peakmap[i].getMSLevel()==1)
					{ 
						spectrum_pointer_container.push_back(&(peakmap[i]));							
					}
				}
			}
			else
			{
				throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"No spectra contained");
			}
		}
		/**
	    @brief does the transformation if the Discrete Cosines Fourier Transformation is selected.
			
		 	Call intern the function transform,only if the comparison score function Fourier is selected.
			@param spectrum_pointer_container is the sequence which has to be transform
	   	
		*/
		void MapAlignmentAlgorithmSpectrumAlignment::fourierActivation_(std::vector<MSSpectrum<>* >& spectrum_pointer_container)
		{
			if(c1_->getName()=="CompareFouriertransform")
			{
				for(UInt i=0; i< spectrum_pointer_container.size();++i)
				{
					transform_(*spectrum_pointer_container[i]);
				}
			}
		}
		/**
    	@brief calculate the Discrete Cosines Fourier Transformation.
       				
   		This Function transform a given MSSpectrum to an Discrete Cosines Fourier Transformation. It stores only the part of the cosines of the FFT in
   		the MetaDataArray which is a container from the MSSpectrum. Only call this function, if you sure there is no earlier an another transformation done over the same MSSpectrum, because it doesn't check if there already exist a transformation.
          		
     	@param spec  MSSpectrum 
     	@see MapAlignmentAlgorithmSpectrumAlignment()
    */
		void MapAlignmentAlgorithmSpectrumAlignment::transform_(MSSpectrum<> & spec)
    {
    	MSSpectrum<>::MetaDataArrays& temp = spec.getMetaDataArrays();
			
    	UInt i=0;
    	if(temp.size()>0)
			{
				while(i< temp.size())
	    	{
	    		if(temp[i].getName()=="Fouriertransformation")
	    		{
	    			break;
	    		}
	    		else
	    			{
	    				++i;
	    			}
    		}
			}
			//read or build
			if( (i < temp.size() && temp[i].getName()!= "Fouriertransformation") || i==0 )
			{
			  	//a copy have to be made
				double* data=  new double [spec.getContainer().size()<<1];
				bool aflag = false;
				bool iflag = true;
				Real sum=0;
  			//normalize first the intensity!!!
  			for(UInt k = 0 ;k<spec.getContainer().size(); ++k)
  			{
  				sum+=spec.getContainer()[i].getIntensity();
  			}
				i =0;
				//copy spectrum to the array, and after that mirrow the data, FFT needs perodic function
				while(!aflag)
				{	
					if(i== (spec.getContainer().size()<<1))	
					{
						aflag=true;
						break;
					}
					if(iflag)
					{
						if(i< spec.getContainer().size())
						{
							data[i]= spec.getContainer()[i].getIntensity()/sum;
							++i;
						}
						else
						{
							iflag= false;
						}
					}
					else
					{
						data[i] =spec.getContainer()[(spec.getContainer().size()<<1)-i].getIntensity()/sum;
						++i;
					}
				}
				//calculate the fft by using gsl
				gsl_fft_real_wavetable * real;
				gsl_fft_real_workspace * work;
				work = gsl_fft_real_workspace_alloc (spec.getContainer().size());
				real = gsl_fft_real_wavetable_alloc (spec.getContainer().size());
				gsl_fft_real_transform (data,1,spec.getContainer().size(),real, work);
				gsl_fft_real_wavetable_free (real);
				gsl_fft_real_workspace_free (work);
				//saving the transformation, but only the real part
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
	
    }
		/**
		 @brief function for the test if cell i,j of the grid is inside the band
	
			The function return true if the cell underly these condition:
			-k<=i-j<=k+n-m
			else retun false.
			@param i Int coordinate i
			@param j Int coordinate j
			@param n UInt  size of column
			@param m UInt  size of row
			@param k_ Int  size of k_
			@see MapAlignmentAlgorithmSpectrumAlignment()
		*/
		inline bool  MapAlignmentAlgorithmSpectrumAlignment::insideBand_(UInt i,Int j,UInt n,UInt m,Int k_) 
	  {
	   	if((Int)(-k_)<=(Int)(i-j) &&(Int) (i-j) <=(Int)(k_+n-m))
	   	{
	  		//std::cout << i << " i " << j << " j " << " innerhalb der Bande " << std::endl;
	   		return true;
	   	}
	   	else 
	  	{
	   		//std::cout << i << " i " << j << " j " << "NICHT innerhalb der Bande " << std::endl;
	    	return false;
	    }
	   }
    /**
     	@brief calculate the size of the band for the alignment for two given Sequence
       		
       		This Function calculates the size of the band for the alignment. It takes three samples from the aligned sequence and tries to 
       		find the highscore pairs(matching against the template sequence). The highscore pair with the worst distance is to be chosen as the size of k. 
       		
       		@param pattern const std::vector<T* > vector of pointers of the template sequence
       		@param aligned std::vector<T* > vector of pointers of the aligned sequence
       		@param buffer std::map<UInt, std::map<UInt,Real> >  holds the calculated score of index i,j.
       		@param flag1 bool flag indicate the order of the matrix   		
       		@param xbegin UInt indicate the beginning of the template sequence
       		@param xend UInt indicate the end of the template sequence
       		@param ybegin UInt indicate the beginning of the aligned sequence
       		@param yend UInt indicate the end of the aligned sequence
     	   	@see MapAlignmentAlgorithmSpecturmAlignment()
    */
  
    inline Int  MapAlignmentAlgorithmSpectrumAlignment::bestk_(const std::vector<MSSpectrum<>* >& pattern, std::vector<MSSpectrum<>* >& aligned,std::map<UInt, std::map<UInt,Real> > & buffer,bool column_row_orientation, UInt xbegin,UInt xend, UInt ybegin, UInt yend)
    {	
    	Int ktemp=1;
      for(Real i = 0.25; i<=0.75;i+=0.25)
    	{	
      	UInt temp= (UInt)((yend-ybegin) * i);
    		Real	maxi=-999.0;
    		Real	s=-999.0;
    		for(UInt k = 0; k<=(xend-xbegin); ++k)
    		{
    			UInt x;
    			Int y;
    			if(column_row_orientation)
    			{ 
    				x= (UInt)temp+1;
    				y= k+1;
    				s= scoreCalculation_(x,y,xbegin,ybegin,pattern,aligned,buffer,column_row_orientation);
    			}
    			else
    			{
    				x= (UInt)k+1;
    				y= (Int)(temp+1);
    				s= scoreCalculation_(x,y,xbegin,ybegin,pattern,aligned,buffer,column_row_orientation);
    			}
    			if(s > maxi && s > cutoffScore_)
    			{
    				maxi = s;
    				if(ktemp < std::abs((Int)x-y)+1)
    				{
    					ktemp =std::abs((Int)x-y)+1;
    				}
    			}
    		}
    	}
    	return ktemp;
    }
   /**
      @brief function to calculate a cubicspline to interpolate the Retention time

   	csp(cubic spline) is needed to interpolate the Retention time for the hole length of the sequence.	
   	The data points are the points, in which a match appeared. To get the rest of the Retention times a spline is necessary.
   	
   	@param x std::vector<int> which contain x cordinates
   	@param y  std::vector<double> which contain the retentiontimes
   	@param aligened MSExperiment the aligned sequence
   	@param begin UInt begin of the alignment in the aligned sequence 
   	@param end UInt end of the alignment in the aligned sequence 
   	@see MapAlignmentAlgorithmSpectrumAlignment()
   */
   inline void  MapAlignmentAlgorithmSpectrumAlignment::calculateSpline_(std::vector<int>& x,std::vector<double>& y, std::vector<MSSpectrum<>* >& aligned,UInt begin, UInt end) 
    {
    	if(x.size() >=3)
    	{
    		double *tempx =  new double [x.size()];
    		double *tempy = new double [x.size()];
    		
    		/*for(Int i= x.size()-1 ; i >=0; --i)
    		{
    			tempx[x.size()-1-i]=x[i];
    			tempy[x.size()-1-i]=y[i];
    		}*/
    		for(UInt i= 0 ; i <x.size(); ++i)
    		{
    			tempx[i]=x[i];
    			tempy[i]=y[i];
    		}
    		    		
    	
    		gsl_interp_accel *acc = gsl_interp_accel_alloc();
    		gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, (UInt)x.size());
    		gsl_spline_init(spline,tempx,tempy,(UInt)x.size());

    		for(UInt i = begin; i <= end; ++i)
    		{
    			(*aligned[i]).setRT(gsl_spline_eval(spline,i,acc));
    		}
    		delete[] tempx;
    		delete[] tempy;
    		gsl_spline_free(spline);
    		gsl_interp_accel_free(acc);
    	}
    	else
    	{/*
    		if(x.size()==0)
    		{
    			std::cerr << " For this part, there is no alignment possible " << std::endl;
    		}
    		else
    		{
    			if(x.size()==2)
    			{
    				if((UInt)x[1]== begin &&(UInt) x[0]==end)
    				{
    					UInt distance = std::abs((Int)(begin-end));
    					UInt difference = std::abs((Int)(y[1]-y[0]));
    					Real k= distance/difference; 
    					for(UInt i=0; i<=distance; ++i)
    					{	
    						Real t=y[1]+i*k;
    						(*aligned[begin+i]).setRT(t);
    					}
    				}
    				else if((UInt)x[1]!=begin &&(UInt) x[0]==end)
    				{
    					if(y[1]>(*aligned[begin]).getRT())
    					{
    						x.push_back(begin);
    						y.push_back((*aligned[begin]).getRT());
    						calculateSpline_(x,y,aligned,begin,end);
    					}
    					else
    					{
    						UInt difference = std::abs((Int)((*aligned[x[1]]).getRT()-y[1]));
    						x.push_back(begin);
    						y.push_back((*aligned[x[1]]).getRT()-difference);
    						calculateSpline_(x,y,aligned,begin,end);
    					}
    				}
    				else if((UInt)x[1]==begin &&(UInt) x[0]!=end)
    				{
    					if(y[0]<(*aligned[end]).getRT())
    					{
    						x.push_back(end);
    						y.push_back((*aligned[end]).getRT());
    						ordering_(x,y);
    						calculateSpline_(x,y,aligned,begin,end);
    					}
    					else
    					{
    						UInt difference = std::abs((Int)((*aligned[x[0]]).getRT()-y[0]));
    						x.push_back(end);
    						y.push_back((*aligned[x[0]]).getRT()+difference);
    						ordering_(x,y);
    						calculateSpline_(x,y,aligned,begin,end);
    					}
    				}
    				else if((UInt)x[1]!= begin &&(UInt) x[0]!= end)
    				{
    					if(y[1]<(*aligned[begin]).getRT())
    					{
    						UInt difference = std::abs((Int)((*aligned[x[1]]).getRT()-y[1]));
    						x.push_back(begin);
    						y.push_back((*aligned[x[1]]).getRT()-difference);
    					}
    					else
    					{	
    						x.push_back(begin);
    						y.push_back((*aligned[begin]).getRT());
    										
    					}
    					if(y[0]>(*aligned[end]).getRT())
    					{
    						UInt difference = std::abs((Int)((*aligned[x[0]]).getRT()-y[0]));
    						x.push_back(end);
    						y.push_back((*aligned[x[0]]).getRT()+difference);
    						
    					}
    					else
    					{
    						x.push_back(end);
    						y.push_back((*aligned[end]).getRT());
    					}
    					ordering_(x,y);
    					calculateSpline_(x,y,aligned,begin,end);
    				}
    			}
    		}*/
    	}
    }
   /*
   inline  void  MapAlignmentAlgorithmSpectrumAlignment::ordering_(std::vector<int>& x,std::vector<double>& y )
    {
     	for(UInt i =0 ; i < x.size()-1;++i)
    	{
    		int temp1=x[i+1];
    		double temp2=y[i+1];
    		x[i+1]=x[0];
    		y[i+1]=y[0];
    		x[0]=temp1;
    		y[0]=temp2;
    	}
    	
    }
    */
   
   	/**
         @brief calculate the score of two given MSSpectrums calls intern scoring_
      		
      		This Function calculates the score from to MSSpectrums. These two MSSpectrums are chosen by the coordinate i,j. 
      		I,j indicate the index in the matrix. To find the right index on the sequence, each beginning is also given to the function.
      		A flag indicates which order the matrix contain. The buffermatrix collects the result of the scoring, if the band expand only a lookup of knowing scores is done
      		
      		@param i UInt is a index from the matrix.
      		@param j Int is a index from the matrix.
      		@param patternbegin UInt indicate the beginning of the template sequence
         	@param aligenbegin UInt indicate the beginning of the aligned sequence
         	@param pattern const std::vector<T* > vector of pointers of the template sequence
         	@param aligned std::vector<T* > vector of pointers of the aligned sequence
         	@param buffer std::map<UInt, std::map<UInt,Real> >  holds the calculated score of index i,j.
         	@param flag1 bool flag indicate the order of the matrix
         	@see MapAlignmentAlgorithmSpecturmAlignment()
   */
   inline	Real MapAlignmentAlgorithmSpectrumAlignment::scoreCalculation_(UInt i,Int j, UInt patternbegin, UInt alignbegin ,const std::vector<MSSpectrum<>* >& pattern,  std::vector<MSSpectrum<>* >& aligned,std::map<UInt, std::map<UInt,Real> > & buffer,bool column_row_orientation)
   {
  	 if(!column_row_orientation)
  	 {
  		 if(buffer[i][j]<=0)
  		 {
  			 buffer[i][j]=scoring_(*pattern[i+patternbegin-1], *aligned[j+alignbegin-1]);
  		 }
  		 return buffer[i][j];
   		}
  	 else
  	 {
  		 if(buffer[j][i]<=0)
  		 {
  			 buffer[j][i]=scoring_(*pattern[j+patternbegin-1],*aligned[i+alignbegin-1]);
  		 }
  		 return buffer[j][i];
  	 }
   }
   
	/**
	    @brief return the score of two given MSSpectrums by calling the scorefunction
	 		
	 		@param a is a MSSpectrum.
	   	@param b is a MSSpectrum
	   	@see MapAlignmentAlgorithmSpecturmAlignment()
	*/
	inline	Real MapAlignmentAlgorithmSpectrumAlignment::scoring_(const MSSpectrum<>& a, MSSpectrum<> & b)
	{			
		return c1_->operator ()(a,b);
	}
	void MapAlignmentAlgorithmSpectrumAlignment::updateMembers_()
	{
		gap_	=(Int)param_.getValue("gapcost");
		e_		=(Int)param_.getValue("affinegapcost");
		c1_ = Factory<PeakSpectrumCompareFunctor>::create((String)param_.getValue("scorefunction"));
		cutoffScore_=(Real)param_.getValue("cutoffScore");
	}

} 


