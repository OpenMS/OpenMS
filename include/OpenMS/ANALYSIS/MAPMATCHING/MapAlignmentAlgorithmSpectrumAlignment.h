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
//  version 2.1 of the License, or (at your option) any later version
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Vipul Patel $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMSPECTRUMALIGNMENT_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMSPECTRUMALIGNMENT_H

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>
#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_real.h>

namespace OpenMS
{
	/**
		@brief A map alignment algorithm based on spectrum similarity (dynamic programming). 		
	*/
	class MapAlignmentAlgorithmSpectrumAlignment
	 : public MapAlignmentAlgorithm, 
	 	public ProgressLogger
	{
		public:
			/// Default constructor
			MapAlignmentAlgorithmSpectrumAlignment();

			/// Destructor
			virtual ~MapAlignmentAlgorithmSpectrumAlignment();

			//Docu in base class
		virtual void alignPeakMaps(std::vector< MSExperiment<> >&);
			
			
			///Creates a new instance of this class (for Factory)
			static MapAlignmentAlgorithm* create()
			{
				return new MapAlignmentAlgorithmSpectrumAlignment();
			}
			
			///Returns the product name (for the Factory)
			static String getProductName()
			{
				return "spectrum_alignment";
			}
			
		private:
			///Copy constructor is not implemented -> private
			MapAlignmentAlgorithmSpectrumAlignment(const MapAlignmentAlgorithmSpectrumAlignment& );
			///Assignment operator is not implemented -> private
			MapAlignmentAlgorithmSpectrumAlignment& operator=(const MapAlignmentAlgorithmSpectrumAlignment& );
			
			
			///gap-cost
			Int gap_;
			///extension cost
			Int e_;
			///
			PeakSpectrumCompareFunctor* c1_;
/**
     @brief A function to prepare the sequence for the alignment. It calls intern the main function for the alignment.
     
	This function takes two arguments. These Arguments type are two MSExperiments. 
	The first argument should have been filtered, so such only type of MSLevel 1 exist in the Sequence. 
	The second arguments doesn't have to full fill these restriction. It's going to be filter automatically. 
	With these two arguments a precalculation is done, to find some corresponding data points(maximum 4) for building alignments blocks.
	After the alignment a retransformation is done, the new Retention Times appears in the original data.      
    
    @improvement multithreading(vipulpatel)		

   	
    The parameters are MSExperiments.
    @param pattern is the template MSExperiment.
    @param aligned is the sequence which has to be align(also MSExperiment).
    @see MapAlignmentAlgorithmSpecturmAlignment()
*/				
template<typename T, typename S >
void preparealign_(const std::vector<S* >& pattern, MSExperiment<T>& aligned)
{
	//tempalign ->container for holding only MSSpectrums with MS-Level 1
	std::vector<S*> *tempalign = new std::vector<S*>();
	//update the Ranges for getting the MS-Levels from the function getMSLevels()
	aligned.updateRanges(-1);
	
	
	//substract the MSLevel 1 to tempalign
	if((aligned.getMSLevels()).size()!=0)
	{
		if((aligned.getMSLevels()).size() >1)
		{
			for(UInt i=0; i< aligned.size();++i)
			{
				if(aligned[i].getMSLevel()==1) tempalign->push_back(&aligned[i]);
				if(c1_->getProductName()=="CompareFouriertransform")
				{
					transform_(aligned[i]);
				}
			}
		}
		else {
			for(UInt i=0; i< aligned.size();++i)
				{
					tempalign->push_back(&aligned[i]);
				}
			}
	}
	else
		{
			throw " This is not a valid MSSpectrum ";
		}
	//if it's possible, built 4 blocks. These can be individually be aligned.						
	std::vector<UInt> alignpoint;
	//saving the first cordinates
	alignpoint.push_back(0); 
	alignpoint.push_back(0);
	UInt temp=0;
	UInt x=0;
	Real maxi=-999.0;
	Real s=-999.0;
	//4 blocks : 0-0.25 ,0.25-50,0.50-0.75,1 The data points must have a high similarity score
	for(Real i = 0.25; i<=0.75;)
	{	
		temp= (UInt)(tempalign->size() * i);
		x=0;
		maxi=-999.0;
		s=-999.0;
		for(UInt k = 0; k<pattern.size(); ++k)
		{
			s =	scoring_(*pattern[k],*((*tempalign)[temp]));
			if(s > maxi && s > 0.90)
			{
				x=k;
				maxi=s;
			}			
		}
		//only save possible data points, if they are not already contained 
		if(x >= alignpoint[alignpoint.size()-2]+3 && temp >= alignpoint[alignpoint.size()-1]+3)
		{
			alignpoint.push_back(x);
			alignpoint.push_back(temp);
		}
		i=i+0.25;
	}
	//save also the endpoint as a data point
	alignpoint.push_back(pattern.size()-1);
	alignpoint.push_back(tempalign->size()-1);
	//the distance of two data points have to be greater than 3, if not the spline would thrown an Expection
	UInt i= 0;
	//do a affine gap alignment of the block of the data points x1,y1,x2,y2
	while(i < alignpoint.size()-2)
	{
	//	std::cout<<alignpoint[i]<< " xstart " << alignpoint[i+1]<<" ystart " << alignpoint[i+2]<<" xend " <<alignpoint[i+3]<< " yend "<<std::endl;
		affineGapalign_(alignpoint[i],alignpoint[i+1],alignpoint[i+2],alignpoint[i+3], pattern, (*tempalign));
		i+=2;
	}
	/*
	//retransform from tempalign to align
	i=0;
	UInt j=0;
	bool loop = false;
	while(i < aligned.size())
	{	
		if(i==1 && !loop)loop =true;
		if(i%(aligned.getMSLevels()).size()==0 && loop) ++j;
		//std::cout<< i << " i " << i%patter.size()<< " mod "<< j << " j " << loop << "nachher"<< std::endl;
		if(j> tempalign->size()) break;
		else
		if(aligned[i].getMSLevel()==1)
		{
			double temp = (*(*tempalign)[j]).getRT();
				
				aligned[i].setRT(temp);				
				
		}
		++i;
	}*/
	delete tempalign;
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
template<typename T>
void affineGapalign_(UInt xbegin, UInt ybegin, UInt xend,UInt yend, const std::vector<T* >& pattern,  std::vector<T* >& aligned)
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
  
  bool flag1= false;
  if(n !=(xend-xbegin)+1)
    {
	flag1 = true;
    }
	//matrix for holding calculated sorces
	std::map<UInt, std::map<UInt,Real> > buffermatrix;
	//calculate the value of k
	k_=	bestk_(pattern, aligned, buffermatrix, flag1,xbegin,xend, ybegin, yend);
	//flag if we have to calculate again the alignment in step k+1
	bool finish = false;
	try
	{
	  while(!finish && k_ <(Int) n+1)	 //if the current score is smaller than the next step -> calculate again
	    {								//best k+1 =(2(k+1)+m-n)d+(n(k+1))matchscore 
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
		{	
		  //j->row
		  Int j=i;
		  //->band-size
		  for(Int h= -k_-(n-m); h <= k_ ; ++h) // from -k to k
		    {
			j=i+h;
			if(j >= 1 && (UInt)j<=m )
			{	
			  double s=-999.0;
			  //score calculation, attention column and row a 1 size smaller
			  try
			  	{
			    s= 	scorecalculation_(i,j,xbegin,ybegin,pattern,aligned,buffermatrix ,flag1);
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
			  Real tempm,tempi;
			  Real maxima[]={-999.0,-999.0};
			  //testing if i-1,j is in inside of the band
			  //std::cout << i << " i " << j << " j " << s <<  "score" << std::endl;
			  i=i-1;
			  if(insideband_(i,j,n,m,k_))
			  { 	
				  i+=1;
				  tempm=matchmatrix[i-1][j]-gap_;
				  tempi=insertmatrix[i-1][j]-e_;
				  //std::cout << i-1 <<  " i-1 " << j << " j " <<tempm << " Matchmatrix " << tempi << " Insertmatrix " << std::endl;
				  if(tempm > tempi)
				  {
					  maxima[0]=tempm;
				  }
				  else
				  	{
					  maxima[0]=tempi;
				  	}										
			  }else	i=i+1;
			  //testing if i,j-1 is in inside of the band
			  j=j-1;
			  if(insideband_(i,j,n,m,k_))
			  	{	
				  j+=1;
				  tempm=matchmatrix[i][j-1]-gap_;
				  tempi=insertmatrix[i][j-1]-e_;
				  //std::cout <<  i <<  " i " << j-1 << " j-1 " << tempm << " Matchmatrix " << tempi << " Insertmatrix " << std::endl;			
				  if(tempm > tempi)
			      	{
					  maxima[1]=tempm;
			      	}
				  else
				  	{
					  maxima[1]=tempi;
				  	}
				  }else j=j+1;
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
	if(score_ >= matchmatrix[n][m] || k_<<1 > (Int)n || matchmatrix[n][m] >= (m-(k_+1))*0.9 -2* gap_ - (2*(k_) +((Int)m-(Int)n))*e_)
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
				{	//std::cout<< flag1 << std::endl;
					try
						{
							if(!flag1)
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
					else						
						if( t==1)
						{/*
							if(flag1)
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
								if(!flag1)
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
		//call the spline
		csp_(xvar,yvar,aligned,ybegin,yend);
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
      @brief return the score of two given MSSpectrums by calling the scorefunction
   		
   		@param a is a MSSpectrum.
      	@param b is a MSSpectrum
      	@see MapAlignmentAlgorithmSpecturmAlignment()
*/
template<typename T>
inline	Real scoring_(const T & a, T & b)
	{			
	return c1_->operator ()(a,b);
	}
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
template<typename T>
inline	Real scorecalculation_(UInt& i,Int& j, UInt& patternbegin, UInt& alignbegin ,const std::vector<T* >& pattern,  std::vector<T* >& aligned,std::map<UInt, std::map<UInt,Real> > & buffer,bool& flag1)
			{
			if(!flag1)
				{
				
				if(buffer[i][j]<=0)
				{
					buffer[i][j]=scoring_(*pattern[i+patternbegin-1], *aligned[j+alignbegin-1]);
					
				}
				
				return buffer[i][j];
				}else
				{
					if(buffer[j][i]<=0)
					{
						buffer[j][i]=scoring_(*pattern[j+patternbegin-1],*aligned[i+alignbegin-1]);
					}
					return buffer[j][i];
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
template<typename T >
Int bestk_(const T & pattern, T& aligned,std::map<UInt, std::map<UInt,Real> > & buffer,bool& flag1, UInt& xbegin,UInt& xend, UInt& ybegin, UInt yend)
{	Int ktemp=1;
	Real maxi;
	UInt temp=0;
	Real s;
	for(Real i = 0.25; i<=0.75;)
				{	
							temp= (UInt)((yend-ybegin) * i);
							maxi=-999.0;
							s=-999.0;
							for(UInt k = 0; k<=(xend-xbegin); ++k)
							{
								UInt x;
								Int y;
								if(flag1){ 
									x= (UInt)i+1;//!!!grrrrrrrrrr
									y  = k+1;
									s= scorecalculation_(x,y,xbegin,ybegin,pattern,aligned,buffer,flag1);
								}
								else
								{
									x  = (UInt)k+1;//!!!grrrrrrrrrr
									y= (Int)(i+1);
									s= scorecalculation_(x,y,xbegin,ybegin,pattern,aligned,buffer,flag1);
								}
								
								if(s > maxi && s > 0.90)
								{
									maxi = s;
									if(ktemp < std::abs((Int)x-y)+1)
									{
										ktemp =std::abs((Int)x-y)+1;
									}
								}
								
							}
									
							i=i+0.25;
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
template<typename T>
void csp_(std::vector<int>& x,std::vector<double>& y, std::vector<T* >& aligned,UInt& begin, UInt& end)
{
	if(x.size() >=3)
	{
		double *tempx =  new double [x.size()];
		double *tempy = new double [x.size()];
		
		for(Int i= x.size()-1 ; i >=0; --i)
		{
			tempx[x.size()-1-i]=x[i];
			tempy[x.size()-1-i]=y[i];
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
	{
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
						csp_(x,y,aligned,begin,end);
					}
					else
					{
						UInt difference = std::abs((Int)((*aligned[x[1]]).getRT()-y[1]));
						x.push_back(begin);
						y.push_back((*aligned[x[1]]).getRT()-difference);
						csp_(x,y,aligned,begin,end);
					}
				}
				else if((UInt)x[1]==begin &&(UInt) x[0]!=end)
				{
					if(y[0]<(*aligned[end]).getRT())
					{
						x.push_back(end);
						y.push_back((*aligned[end]).getRT());
						ordering_(x,y);
						csp_(x,y,aligned,begin,end);
					}
					else
					{
						UInt difference = std::abs((Int)((*aligned[x[0]]).getRT()-y[0]));
						x.push_back(end);
						y.push_back((*aligned[x[0]]).getRT()+difference);
						ordering_(x,y);
						csp_(x,y,aligned,begin,end);
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
					csp_(x,y,aligned,begin,end);
				}
			}
			
			
		}
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
bool inline insideband_(UInt & i,Int& j,UInt& n,UInt& m,Int k_) 
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

template<typename T>
    void transform_(T & spec)
    {
			std::vector<DSpectrum<>::MetaDataArray > &temp = spec.getMetaDataArrays();
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
			  	//leider muß eine Kopie gemacht werden...
				double* data=  new double [spec.getContainer().size()<<1];
				bool aflag = false;
				bool iflag = true;
				i =0;
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

void ordering_(std::vector<int>& x,std::vector<double>& y )
{
	int temp1 = 0;
	double temp2 =0.0;
	for(UInt i =0 ; i < x.size()-1;++i)
	{
		temp1=x[i+1];
		temp2=y[i+1];
		x[i+1]=x[0];
		y[i+1]=y[0];
		x[0]=temp1;
		y[0]=temp2;
	}
	
}
	void updateMembers_();


	};
	 
	

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMSPECTRUMALIGNMENT_H


