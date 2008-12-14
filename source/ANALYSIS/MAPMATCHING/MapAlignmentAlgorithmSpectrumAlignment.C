// -*- mode: C++; tab-width: 2; -*-
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
		: MapAlignmentAlgorithm(),
			c1_(0)
	{
		setName("MapAlignmentAlgorithmSpectrumAlignment");
		defaults_.setValue("gapcost",1.0," This Parameter stands for the cost of opining a gap in the Alignment. A Gap means that one Spectrum can not be aligned directly to another Spectrum in the Map. This happens, when the similarity of both spectra a too low or even not present. Imagen as a insert or delete of the spectrum in the map. The gap is necessary for aligning, if we open a gap there is a possibility that an another spectrum can be correct aligned with a higher score as before without gap. But to open a gap is a negative event and has to be punished a bit, so such only in case  it 's a good choice to open a gap, if the score is bad enough. The Parameter is to giving as a positive number, the implementation convert it to a negative number.");
		defaults_.setValue("affinegapcost", 0.5," This Parameter controls the cost of extension a already open gap. The idea behind the affine gapcost lies under the assumption, that it is better to get a long distance of connected gaps than to have a structure gap match gap match.  There for the punishment for the extension of a gap has to be lower than the normal gapcost. If the the result of the aligmnet show high compression, it is a good idea to lower the affine gapcost or the normal gapcost.");
		defaults_.setValue("cutoffScore",0.70,"The Parameter defines the threshold which filtered Spectra, these Spectra are high potential candidate for deciding the interval of a sub-alignment.  Only those pair of Spectra are selected, which has a score higher or same of the threshold.",StringList::create("advanced"));
		defaults_.setValue("bucketsize",10,"Defines the numbers of buckets. It is a quantize of the interval of those points, which defines the main alignment(match points). These points have to filtered, to reduce the amount of points for the calculating a smoother spline curve.",StringList::create("advanced"));
		defaults_.setValue("anchorpoints",10,"Defines the percent of numbers of match points which a selected from one bucket. The high score pairs are previously selected. The reduction of match points helps to get a smoother spline curve.",StringList::create("advanced"));
		defaults_.setValue("debug","false","active the debug mode, there a files written starting with debug prefix.",StringList::create("advanced"));
		defaults_.setValidStrings("debug",StringList::create("true,false"));
		defaults_.setValue("mismatchscore",-1.0,"Defines the score of two Spectra if they have no similarity to each other. ",StringList::create("advanced"));
		defaults_.setValue("scorefunction","SteinScottImproveScore"," The score function is the core of an alignment. The success of an alignment depends mostly of the elected score function. The score function return the similarity of two Spectrum back. The score influence defines later the way of possible traceback. There exist many way of algorithm to calculate the score.");
		defaults_.setValidStrings("scorefunction",Factory<PeakSpectrumCompareFunctor>::registeredProducts());
		defaultsToParam_();

		setLogType(CMD);
	}

	MapAlignmentAlgorithmSpectrumAlignment::~MapAlignmentAlgorithmSpectrumAlignment()
	{
		delete c1_;
	}
	
	void MapAlignmentAlgorithmSpectrumAlignment::alignPeakMaps(std::vector< MSExperiment<> >& peakmaps, std::vector<TransformationDescription>& transformation)
	{
		try
		{ 	
			std::vector<MSSpectrum<>* >spectrum_pointers;
			msFilter_(peakmaps[0],spectrum_pointers);
			fourierActivation_(spectrum_pointers);

		
			startProgress(0,(peakmaps.size()-1),"Alignment");
			for(UInt i = 1 ; i < peakmaps.size();++i )
			{									
				prepareAlign_(spectrum_pointers,peakmaps[i],transformation);
				setProgress(i);
			}
			eraseMetaDataArrayEntry_(spectrum_pointers);
			endProgress();
		}
		catch (Exception::OutOfRange& /*e*/) 
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
		@see MapAlignmentAlgorithmSpectrumAlignment();
	*/	
	void MapAlignmentAlgorithmSpectrumAlignment::prepareAlign_(const std::vector<MSSpectrum<>* >& pattern, MSExperiment<>& aligned,std::vector<TransformationDescription>& transformation)
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

		std::vector<Int> xcoordinate;
		std::vector<Int> xcoordinatepattern;
		std::vector<Real>ycoordinate;
		debugmatrix_.clear();
				
		for(UInt i=0; i < alignpoint.size()-2; i+=2)
		{
			affineGapalign_(alignpoint[i],alignpoint[i+1],alignpoint[i+2],alignpoint[i+3], pattern, tempalign,xcoordinate,ycoordinate,xcoordinatepattern);
		}
		
		//affineGapalign_(0,0,(pattern.size()-1),(tempalign.size()-1), pattern, tempalign,xcoordinate,ycoordinate,xcoordinatepattern);
		if(debug_)
		{
			debugFileCreator_(pattern,tempalign);
		}
		
	/*
		for(UInt i = 0; i< xcoordinate.size(); ++i)
		{
			std::cout<< xcoordinate[i] << " " << ycoordinate[i] << " x  y  anchorpunkte " << std::endl;
		}
		std::cout << std::endl;
		*/
		
		bucketFilter_(pattern, tempalign,xcoordinate,ycoordinate,xcoordinatepattern);
		/*std::cout << xcoordinate.size()<< std::endl;
				for(UInt i = 0; i< xcoordinate.size(); ++i)
				{
					std::cout<< xcoordinate[i] << " " << ycoordinate[i] << " x  y  anchorpunkte " << std::endl;
				}*/
		//calculate the spline
		calculateSpline_(xcoordinate,ycoordinate,tempalign,(UInt)0,tempalign.size()-1,transformation);
		if(xcoordinate[0]!=0)
		{
			for(UInt i =0; i <=(UInt)xcoordinate[0];++i)
			{
				Real rt = (i+1)*(ycoordinate[0]/(xcoordinate[0]+1));
				(*tempalign[i]).setRT(rt);
			}
		}
		if((UInt)xcoordinate[xcoordinate.size()-1]!=tempalign.size()-1)
		{
			Real slope= (ycoordinate[ycoordinate.size()-2]-ycoordinate[ycoordinate.size()-1])/(xcoordinate[xcoordinate.size()-2]-xcoordinate[xcoordinate.size()-1]);
		
			for(UInt i = xcoordinate[xcoordinate.size()-1]; i < tempalign.size();++i)
			{
				Real rt = (i)*slope;
				(*tempalign[i]).setRT(rt);
			}
		}
		//todo plot the result!!!
	/*	
		for(UInt i=0; i< tempalign.size();++i)
		{
			std::cout << (*tempalign[i]).getRT()<< " " << std::endl;
		}
		*/
		eraseMetaDataArrayEntry_(tempalign);
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
		@param xcoordinate std::vector<int> save the postion of anchorpoints
		@param ycoordinate std::vector<double> save the retentiontimes of an anchorpoints
		@param xcoordinatepattern std::vector<int> save the reference position of the anchorpoints from the pattern
		@see MapAlignmentAlgorithmSpectrumAlignment();
	 */
		void MapAlignmentAlgorithmSpectrumAlignment::affineGapalign_(UInt xbegin, UInt ybegin, UInt xend,UInt yend, const std::vector<MSSpectrum<>* >& pattern,  std::vector<MSSpectrum<>* >& aligned,std::vector<int>& xcoordinate, std::vector<Real>&ycoordinate, std::vector<int>& xcoordinatepattern)
		{	
			
			
			
			
			//affine gap alignment needs two matrices
		  std::map<UInt, std::map<UInt,Real> > matchmatrix;
		  std::map<UInt, std::map<UInt,Real> > insertmatrix;	
		  //setting the size of column and row	
			
		  UInt n= std::max((xend-xbegin),(yend-ybegin))+1; //column 
		  UInt m= std::min((xend-xbegin),(yend-ybegin))+1; //row
		 // std::cout<< n << " n " << m << " m " <<  xbegin << " " <<xend<< " " << ybegin << " " << yend <<std::endl;
		  //log the Progress of the subaligmnet
			String temp = "sub-alignment of interval: template sequence " + String(xbegin) + " " + String(xend) + " interval: alignsequence " + String(ybegin) + " " + String(yend);
			startProgress(0,n,temp);
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
		  	std::vector<std::vector<Real> >tempdebugbuffer;//save the cell i, j , matchscore, insertscore, traceback
		  	
			while(!finish && k_ <(Int) n+1)	 //if the current score is smaller than the next step -> calculate again
			  {																//best k+1 =(2(k+1)+m-n)d+(n(k+1))matchscore 
																				//if k >= columnsize -> there is no band necessary(normal global alignment)
			  	tempdebugbuffer.clear();//because a next k run is done so clear the old one
				std::map<UInt,std::map<UInt,UInt> > traceback;
		      //init both  matrices
			  	
			  	for(UInt i = 1; i<= n; ++i)
			  	{	//->band-size
						setProgress(i);
			  		std::vector<Real>debugtmp;
			  		for(Int h= -k_-(n-m); h <= k_ ; ++h) // from -k to k
				    {
			  			//j->row
			  			
			  			Int j=i+h;
			  			if(j >= 1 && (UInt)j<=m )
			  			{	
			  				DoubleReal s=-999.0;
			  				//score calculation, attention column and row a 1 size smaller
			  				try
			  				{
			  					s	=	scoreCalculation_(i,j,xbegin,ybegin,pattern,aligned,buffermatrix ,column_row_orientation);
									//s=s*para_;
									if(debug_)
									{
										std::vector<Real> temp;
										if(!column_row_orientation)
			  						{	
			  							temp.push_back((Real)i+xbegin-1);
											temp.push_back((Real)j+ybegin-1);
											temp.push_back(s);
											temp.push_back(0);
											debugscorematrix_.push_back(temp);
										}
										else
										{
											temp.push_back((Real)j+xbegin-1);
											temp.push_back((Real)i+ybegin-1);
											temp.push_back(s);
											temp.push_back(0);
											debugscorematrix_.push_back(temp);
										}
									}
									//std::cout << s << " s " << std::endl;
			  				}
								catch (Exception::OutOfRange & /*e*/) 
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
			  				//only necessary for the debug mode
			  				if(debug_)
			  				{
			  					if(!column_row_orientation)
			  					{
			  						debugtmp.push_back(i+xbegin-1);
			  						debugtmp.push_back(j+ybegin-1);
			  					}
			  					else
			  					{
			  						debugtmp.push_back(j+xbegin);
			  						debugtmp.push_back(i+ybegin);
			  					}
			  					debugtmp.push_back(matchmatrix[i][j]);
			  				}
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
			  					if(maxima[0]!=-99999.0 || maxima[1]!=-99999.0)
			  				{ 
			  					if(maxima[0]>=maxima[1])
			  					{
			  						insertmatrix[i][j]= maxima[0];
			  						if(debug_)
			  						{
			  							debugtmp.push_back(insertmatrix[i][j]);
			  						}
			  						Real tempcontainer = std::max(matchmatrix[i-1][j],insertmatrix[i-1][j]);
			  						if(tempcontainer >std::max(matchmatrix[i-1][j-1],insertmatrix[i-1][j-1]))//diagonale
			  						{
			  							traceback[i][j]=1;
			  							if(debug_)
			  							{
			  								debugtmp.push_back(1);
			  							}
			  						}
			  						else
			  						{
			  							if(debug_)
			  							{
			  								debugtmp.push_back(0);
			  							}
			  						}
			  					}
			  					else
			  					{
			  						insertmatrix[i][j]= maxima[1];
			  						if(debug_)
			  						{
			  							debugtmp.push_back(insertmatrix[i][j]);
			  						}
			  						Real tempcontainer=std::max(matchmatrix[i][j-1],insertmatrix[i][j-1]);
			  						if(tempcontainer >std::max(matchmatrix[i-1][j-1],insertmatrix[i-1][j-1]))
			  						{
			  							traceback[i][j]=2;
			  							if(debug_)
			  							{
			  								debugtmp.push_back(2);
			  							}
			  						}
			  						else
			  						{
			  							if(debug_)
			  							{
			  								debugtmp.push_back(0);
			  							}
			  						}
			  					}
			  				}
			  			}
				    }
			  		if(debug_)
			  		{
			  			tempdebugbuffer.push_back(debugtmp);
			  			debugtmp.clear();
			  		}
			  	}
/*
					for(UInt i=0; i <=n;++i)
					{
						for(UInt j=0; j<=m;++j)
						{
							std::cout << i << " "<< j << " i j " << matchmatrix[i][j] << " " << insertmatrix[i][j] << " " << traceback[i][j] << " Match Insert Trace" << std::endl;
						}
						std::cout<< std::endl;
					}
					*/
			  	//std::cout << score_ << " current score " << matchmatrix[n][m] << " new score " <<k_ <<" bandsize "<<  n<< " n " << m << " m " << n <<" dimension of matrix "<< (m-(k_+1))*0.90 -2* gap_ - 2*(k_)*e_ -((Int)n-(Int)m) << "nextscore by extension" <<std::endl;
			  	// hier abstand berechnen (m-(k_+1))*0.9 -2* gap_ - 2*(k_)*e_ -((Int)n-(Int)m)
			  	//only do the traceback, if the actual score is not higher as the old one (m-(k_+1))*0.9 -2* gap_ - (2*(k_+1) +((Int)n-(Int)m))*e
			  	if(score_ >= matchmatrix[n][m] || k_<<1 > (Int)n || matchmatrix[n][m] >= (m-(k_+1))*cutoffScore_ -2* gap_ - (2*(k_) +((Int)m-(Int)n))*e_) //
			  	{
			  		finish= true;
			  		if(debug_)
			  		{
			  			for(UInt i =0; i <tempdebugbuffer.size();++i )
			  			{
			  				debugmatrix_.push_back(tempdebugbuffer[i]);
			  			}
			  		}
						//traceback there is now no use of matchmatrix insertmatrix and buffermatrix  so clear them
						buffermatrix.clear();
						matchmatrix.clear();
						insertmatrix.clear();
			  		bool endtraceback= false;
			  		int i= n;
			  		int j= m;
			  		Real maximum = -999.0;
						//container necessary for collecting the positions of both sequence to gain later the correct datapoints for the spline
			  		std::vector<int> xvar;
			  		std::vector<int> xxvar;
			  		std::vector<Real> yvar;
						
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
											//showing the traceback for the debug mode
											if(debug_)
											{
												debugtraceback_.push_back(std::make_pair(i+xbegin-1,j+ybegin-1));
											}
			  							xvar.push_back(j+ybegin-1);
			  							yvar.push_back((*pattern[i+xbegin-1]).getRT());
			  							xxvar.push_back(i+xbegin-1);
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
												xxvar.push_back(j+xbegin-1);
											//showing the traceback for the debug mode
											if(debug_)
											{
												debugtraceback_.push_back(std::make_pair(j+xbegin-1,i+ybegin-1));
											}
			  		
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
			  							//showing the traceback for the debug mode
									if(debug_)
									{
										if(!column_row_orientation)
										{
											debugtraceback_.push_back(std::make_pair(i+xbegin-1,j+ybegin-1));	
										}
										else
										{
										debugtraceback_.push_back(std::make_pair(j+xbegin-1,i+ybegin-1));	
										}
									}
										
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
									if(debug_)
									{
										if(!column_row_orientation)
										{
											debugtraceback_.push_back(std::make_pair(i+xbegin-1,j+ybegin-1));	
										}
										else
										{
											debugtraceback_.push_back(std::make_pair(j+xbegin-1,i+ybegin-1));	
										}
									}
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
			  		for(UInt i=0; i <xvar.size(); ++i)
						{								
							//std::cout<< xvar[xvar.size()-1-i] << " " << std::endl;
							if(xcoordinate.size()>0)
							{
								if(xvar[xvar.size()-1-i]!=xcoordinate[xcoordinate.size()-1])
								{
									xcoordinate.push_back(xvar[xvar.size()-1-i]);
									ycoordinate.push_back(yvar[yvar.size()-1-i]);
									xcoordinatepattern.push_back(xxvar[xxvar.size()-1-i]);
								}
							}
							else
							{
								xcoordinate.push_back(xvar[xvar.size()-1-i]);
								ycoordinate.push_back(yvar[yvar.size()-1-i]);
								xcoordinatepattern.push_back(xxvar[xxvar.size()-1-i]);
							}
						}
			   		//std::cout<< xcoordinate.size()<< std::endl;
						endProgress();
			  	}
			  	score_ = matchmatrix[n][m];
			  	k_ *=2;
			  }
			}
			catch (Exception::OutOfRange& /*e*/) 
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
			@see MapAlignmentAlgorithmSpectrumAlignment();
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
			@see MapAlignmentAlgorithmSpectrumAlignment();
		*/
		inline void MapAlignmentAlgorithmSpectrumAlignment::fourierActivation_(std::vector<MSSpectrum<>* >& spectrum_pointer_container)
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
			@brief Delete entries of the MetaDataArray with was made from CompareFouriertransform 

			This function erase the entries with was done by the CompareFouriertransform function.

			@param spectrum_pointer_container std::vector<MSSpectrum<>* > contains MSSpectrums
			@see MapAlignmentAlgorithmSpecturmAlignment()
	 */
		inline void MapAlignmentAlgorithmSpectrumAlignment::eraseMetaDataArrayEntry_(std::vector<MSSpectrum<>* >& spectrum_pointer_container)
		{
			
			if(c1_->getName()=="CompareFouriertransform")
			{
				for(UInt i=0; i <spectrum_pointer_container.size();++i)
				{
					MSSpectrum<>::MetaDataArrays& temp = (*spectrum_pointer_container[i]).getMetaDataArrays();
					
					if(temp.size()>0)
					{
						MSSpectrum<>::MetaDataArrays::iterator iter;
						iter = temp.begin();
						while(iter!= temp.end())
						{
							if((*iter).getName()=="Fouriertransformation")
							{
								temp.erase(iter);
								break;
							}
							else ++iter;
						}
					}
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
    
			if(i == temp.size() || i==0||temp[i].getName()!= "Fouriertransformation" )
			{
			  	//a copy have to be made
				double* data=  new double [spec.size()<<1];
				bool aflag = false;
				bool iflag = true;
				Real sum=0;
  			//normalize first the intensity!!!
  			for(UInt k = 0 ;k<spec.size(); ++k)
  			{
  				sum+=spec[k].getIntensity();
  			}
				i =0;
				//copy spectrum to the array, and after that mirrow the data, FFT needs perodic function
				while(!aflag)
				{	
					if(i== (spec.size()<<1))	
					{
						aflag=true;
						break;
					}
					if(iflag)
					{
						if(i< spec.size())
						{
							data[i]= spec[i].getIntensity()/sum;
							++i;
						}
						else
						{
							iflag= false;
						}
					}
					else
					{
						data[i] =spec[(spec.size()<<1)-i].getIntensity()/sum;
						++i;
					}
				}
				//calculate the fft by using gsl
				gsl_fft_real_wavetable * real;
				gsl_fft_real_workspace * work;
				work = gsl_fft_real_workspace_alloc (spec.size());
				real = gsl_fft_real_wavetable_alloc (spec.size());
				gsl_fft_real_transform (data,1,spec.size(),real, work);
				gsl_fft_real_wavetable_free (real);
				gsl_fft_real_workspace_free (work);
				//saving the transformation, but only the real part
				i= temp.size();
			  temp.resize(i+1);
			  temp[i].setName("Fouriertransformation");
			  UInt j=0;
			  while(j < spec.size())
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
    	Int ktemp=2;
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
   inline void  MapAlignmentAlgorithmSpectrumAlignment::calculateSpline_(std::vector<int>& x,std::vector<Real>& y, std::vector<MSSpectrum<>* >& aligned,UInt begin, UInt end,std::vector<TransformationDescription>& transformation) 
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
				std::vector<std::pair<Real,Real> >PairVector;
    		for(UInt i= 0 ; i <x.size(); ++i)
    		{
    			tempx[i]=x[i];
    			tempy[i]=y[i];
					//fill pairs with the anchorpoints
					PairVector.push_back(std::make_pair((Real)x[i],(Real)y[i]));
    		}
				//setting transformation by replacing the PairVector
    		TransformationDescription trans;
				trans.setName("spline");
				trans.setPairs(PairVector);
				transformation.push_back(trans);
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
  		 if(buffer[i][j]==0)
  		 {	
				 Real score = scoring_(*pattern[i+patternbegin-1], *aligned[j+alignbegin-1]);
				 if(score==0) score=mismatchscore_; 
				 buffer[i][j]=score;
  		 }
  		 return buffer[i][j];
   		}
  	 else
  	 {
  		 if(buffer[j][i]==0)
  		 {
				 Real score= scoring_(*pattern[j+patternbegin-1],*aligned[i+alignbegin-1]);
				 if(score==0)score=mismatchscore_;
				 buffer[j][i]=score;
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
	
	
	/**
		@brief Help function of preparation function of data points to construct a spline function.

		This function reduced the amount of data values for the next step. The reduction is done by using a number of buckets, where the data points a selected. 
		Within the buckets only defined number a selected to be written back as a data point. 
		The selection within the buckets is done by the scoring function.
		
		@param pattern is the template MSExperiment.
		@param aligned is to be align MSExperiment.
		@param xcoordinate std::vector<int> save the postion of anchorpoints
		@param ycoordinate std::vector<double> save the retentiontimes of an anchorpoints
		@param xcoordinatepattern std::vector<int> save the reference position of the anchorpoints from the pattern
		@see MapAlignmentAlgorithmSpecturmAlignment()
	*/
	inline void MapAlignmentAlgorithmSpectrumAlignment::bucketFilter_(const std::vector<MSSpectrum<>* >& pattern,  std::vector<MSSpectrum<>* >& aligned, std::vector<Int> & xcoordinate, std::vector<Real> & ycoordinate, std::vector<Int>&xcoordinatepattern)
	{
		std::vector<Int>tempxcoordinate;//cordinate aligned
		std::vector<Real>tempycoordinate;//rt aligigned

		std::vector<std::pair<std::pair<Int, Real>, Real> > tempxy;
		UInt size= xcoordinate.size();
		//std::cout <<bucketsize_  << " bucketsize "<<size << " anchorpoints " <<xcoordinate.size() << " xsize()" << std::endl;
		if(bucketsize_ >= xcoordinate.size()) size=xcoordinate.size();
		else size=size/bucketsize_;
		//std::cout << size << " size "<< xcoordinate.size() << " size2 " << std::endl;
		for(UInt i = 0; i<=xcoordinate.size()-size;i+=size)
		{
			std::vector<std::pair<std::pair<Int, Real>, Real> > temp;
			for(UInt j=0; j<size;++j)//inside a bucket
			{
				Real score = scoring_(*pattern[xcoordinatepattern[i+j]],*aligned[xcoordinate[i+j]]);
				
				temp.push_back(std::make_pair(std::make_pair(xcoordinate[i+j],ycoordinate[i+j]),score));		
			}
			/*for(UInt i=0; i < temp.size();++i)
			{ 
				std::cout<< (temp[i].first).first << " " << (temp[i].first).second << " in temp"<< std::endl; 
			}
			std::cout << std::endl;
			*/
			std::sort(temp.begin(),temp.end(),Compare(false));
			Int anchor=(Int)(size*anchorPoints_/100);
			if(anchor<=0) anchor =1;
		
			//std::cout << anchor << " anchorpoints "  << anchorPoints_<< std::endl;
			/*for(UInt i=0; i< temp.size();++i)
			{
				std::cout<< (temp[i].first).first << "first" << (temp[i].first).second << " second" <<std::endl; 
			}*/
			for(UInt k =0; k< (UInt)anchor;++k)
			{
				tempxy.push_back(temp[k]);
			}
		}
		std::sort(tempxy.begin(),tempxy.end(),Compare(true));
		
		xcoordinate.clear();
		ycoordinate.clear();
		for(UInt i =0; i <tempxy.size();++i)
		{
			xcoordinate.push_back((tempxy[i].first).first);
			ycoordinate.push_back((tempxy[i].first).second);
		}
	}
	/**
		@brief Creates files for the debugging

		This function is only active if the debugflag ist true. The debugfileCreator creates following files debugtraceback.txt(gnuplotScript), debugscoreheatmap.r and debugRscript. Debugscoreheatmap.r contains the scores of the Spectra to each other from the alignment and also the traceback. DebugRscript is the R script which reads those data. So both files are only working under R. Start R and type main(location of debugscoreheatmap.r). The output will be a heatmap of each sub-alignment. Debugtraceback.txt shows the way of the Traceback by using gnuplot.

		@param pattern MSExperiment template MSExperiment.
		@param aligned MSExperiment sequence to be align MSExperiment.
		@see MapAlignmentAlgorithmSpecturmAlignment()
	*/
	inline void MapAlignmentAlgorithmSpectrumAlignment::debugFileCreator_(const std::vector<MSSpectrum<>* >& pattern,  std::vector<MSSpectrum<>* >& aligned)
	{
		
		//plotting scores of the alignment
		/*std::ofstream tempfile3;
		Real maximimum=2.0;
		tempfile3.open("debugscore.txt",std::ios::trunc);
		tempfile3 << "set xrange[0:"<< pattern.size()-1<<  "]" << "\n set yrange[0:"<< aligned.size()-1 << "]" << "\n set zrange[0:" 
				<< maximimum << "] \n set view 45,20,1.0,2.5 \n"<< "splot \'-\'" <<std::endl;
		for(UInt i =0; i< debugscorematrix_.size();++i)
		{
			tempfile3<< debugscorematrix_[i][0] << " " << debugscorematrix_[i][1] << " " << debugscorematrix_[i][2] << std::endl;
		}
		tempfile3 << "e" << std::endl;
		tempfile3.close();*/
		
		//gnuplot of the traceback
		std::ofstream myfile;
		myfile.open("debugtraceback.txt",std::ios::trunc);
		myfile << "set xrange[0:"<< pattern.size()-1<<  "]" << "\n set yrange[0:"<< aligned.size()-1 << "] \n plot \'-\' with lines " << std::endl;
		std::sort(debugtraceback_.begin(),debugtraceback_.end(),Compare(false));
		
		for(UInt i =0; i< debugtraceback_.size(); ++i)
		{
			myfile<< debugtraceback_[i].first << " " << debugtraceback_[i].second << std::endl;
			for(UInt p=0;p<debugscorematrix_.size();++p)
			{
				if(debugscorematrix_[p][0] == debugtraceback_[i].first && debugscorematrix_[p][1] == debugtraceback_[i].second)
				{
					debugscorematrix_[p][3]=1;
					break;
				}
				
			}
		}
		myfile <<"e" <<std::endl;
		myfile.close();
		
		
		//R heatplot score of both sequence
		std::map<UInt,std::map<UInt,Real> >debugbuffermatrix;
		
		Real scoremaximum=-2;
		//precalculation for the heatmap
		//getting the maximum score
		for(UInt i=0; i < debugscorematrix_.size();++i)
		{
			if(scoremaximum < debugscorematrix_[i][2]+2) scoremaximum=debugscorematrix_[i][2]+2;
			//shift all score about 2 (to get 0, the default score in the debugbuffermatrix is -2 )
			debugscorematrix_[i][2]+=2;
		}
		//to get the intvall [0,1] divide all score to the global maximum
		for(UInt i=0; i < debugscorematrix_.size();++i)
		{
			if(debugscorematrix_[i][2]!=0) debugscorematrix_[i][2]/=scoremaximum;
		}
	
		
		//write the score in a file
		/*
		for(UInt i=0; i < debugscorematrix_.size();++i)
		{
			debugbuffermatrix[(UInt)debugscorematrix_[i][0]][(UInt)debugscorematrix_[i][1]]=debugscorematrix_[i][2];
		}*/
		std::ofstream scorefile;
		scorefile.open("debugscoreheatmap.r",std::ios::trunc);
		/*
		for(UInt i=0; i < debugbuffermatrix.size();++i)
		{
			for(UInt j=0; j< debugbuffermatrix[i].size(); ++j)
			{
				scorefile<< i << " "<< j << " "<<  debugbuffermatrix[i][j] << std::endl; 	
			}
		}*/
		for(UInt i=0; i < debugscorematrix_.size();++i)
		{
			scorefile<< debugscorematrix_[i][0] << " " << debugscorematrix_[i][1] << " " << debugscorematrix_[i][2] << " "<< debugscorematrix_[i][3] << std::endl;
			
		}
		scorefile.close();
		
		std::ofstream rscript;
		rscript.open("debugRscript.r",std::ios::trunc);
		
		rscript<<"#Name: LoadFile \n #transfer data from file into a matrix \n #Input: Filename \n #Output Matrix \n LoadFile<-function(fname){\n temp<-read.table(fname); \n temp<-as.matrix(temp); \n return(temp); \n } "<< std::endl;
		rscript<<"#Name: ScoreHeatmapPlot \n #plot the score in a way of a heatmap \n #Input: Scorematrix \n #Output Heatmap \n ScoreHeatmapPlot<-function(matrix) { \n xcord<-as.vector(matrix[,1]); \n ycord<-as.vector(matrix[,2]); \n color<-rgb(as.vector(matrix[,4]),as.vector(matrix[,3]),0);\n  plot(xcord,ycord,col=color, main =\"Heatplot of scores included the traceback\" , xlab= \" Template-sequence \", ylab=\" Aligned-sequence \", type=\"p\" ,phc=22)\n } \n main<-function(filenamea) { \n a<-Loadfile(filenamea) \n X11() \n ScoreHeatmapPlot(a) \n  " <<std::endl;
		rscript.close(); 
		/*
		Real matchmaximum=-999.0;
		Real insertmaximum=-999.0;
		for(UInt i =0; i< debugmatrix_.size();++i)
		{
		debugbuffermatrix[debugmatrix_[i][0]][debugmatrix_[i][1]]=(Real)debugmatrix_[i][2];
		if(matchmaximum <debugmatrix_[i][2]) matchmaximum = (Real)debugmatrix_[i][2];
				//debuginsertmatrix[debugmatrix_[i][0]][debugmatrix_[i][1]]=debugmatrix_[i][3];
		if(insertmaximum <debugmatrix_[i][3]) insertmaximum = (Real)debugmatrix_[i][3];
				//debugtracebackmatrix[debugmatrix_[i][0]][debugmatrix_[i][1]]=debugmatrix_[i][4];
	}
			
			//todo option for gnuplot
		myfile<< "#plotting matchmatrix" <<std::endl;
		myfile<< "set multiplot layout 1,1 columnsfirst" << std::endl;
		myfile<< "set title \"Heat Map of the matchmatrix\"  \n "
		"unset key \n"
		"set tic scale 0 \n set palette rgbformula -7,2,-7 \n set cbrange [-999.0:" << matchmaximum <<"] \n set cblabel \"Score\" \n unset cbtics"<< std::endl;
			myfile<< "p \'-\' using 1:2:3 with image" <<std::endl;
			//output matchmatrix
		for(UInt i=0; i< debugbuffermatrix.size();++i)
		{
		if(i!=0) myfile << std::endl;
		for(UInt j=0; j < debugbuffermatrix[0].size();++j)
		{
		myfile<< i << " " << j << " " << debugbuffermatrix[i][j] << std::endl;
	}
	}
		myfile <<"e"<< std::endl;
		myfile << std::endl;
		std::cout << " complete matchmatrix " << std::endl;
		for(UInt i =0; i< debugmatrix_.size();++i)
		{
		debugbuffermatrix[debugmatrix_[i][0]][debugmatrix_[i][1]]=debugmatrix_[i][3];
	}
			
		myfile<< "set title \"Heat Map of the insertmatrix\"  \n "
		"unset key \n"
		"set tic scale 0 \n set palette rgbformula -7,2,-7 \n	set cbrange [-999.0:"<< insertmaximum<<"] \n set cblabel \"Score\" \n unset cbtics"<< std::endl;
		myfile<< "p \'-\' using 1:2:3 with image" <<std::endl;
		for(UInt i=0; i< debugbuffermatrix.size();++i)
		{
		if(i!=0) myfile << std::endl;
		for(UInt j=0; j < debugbuffermatrix[0].size();++j)
		{
		myfile<< i << " " << j << " " << debugbuffermatrix[i][j] << std::endl;
	}
	}
		myfile <<"e"<< std::endl;
		std::cout << "complete insertmatrix " << std::endl;			
		myfile << std::endl;
			
		myfile<< "set multiplot layout 1,1 columnsfirst \n set title \"Heat Map of the tracebackmatrix\" \n unset key \n set tic scale 0 \n set palette rgbformula -7,2,-7 \n set cbrange [0:5] \n set cblabel \"Score\" \n unset cbtics \n set xrange [0: " << pattern.size()-1<< "] \n set yrange[0:" << tempalign.size()<< "] \n set view map \n"<< std::endl;
		myfile<< "p \'-\' using 1:2:3 with image" <<std::endl;
		for(UInt i=0;i <pattern.size();++i)
		{
			for(UInt j=0; j<tempalign.size();++j)
			{
				debugbuffermatrix[i][j]=0;
			}
		}
		for(UInt i =0; i< debugmatrix_.size();++i)
		{
			Real score = 0;
			if(debugmatrix_[i][4] == 0)
			{
				score = 5;				
			}
			else if(debugmatrix_[i][4] == 1)
			{
				score = 3;
			}
			else if(debugmatrix_[i][4] == 2)
			{
				score = 1;
			}
			debugbuffermatrix[debugmatrix_[i][0]][debugmatrix_[i][1]]=score;
		}
		for(UInt i=0; i< debugbuffermatrix.size();++i)
		{
			if(i!=0) myfile << std::endl;
			for(UInt j=0; j < debugbuffermatrix[0].size();++j)
			{
				myfile<< i << " " << j << " " << debugbuffermatrix[i][j] << std::endl;
			}
		}
		myfile <<"e"<< std::endl;
		std::cout << " complete tracematrix " << std::endl;
		myfile << std::endl;
						
		myfile <<"unset multiplot"<<std::endl;
			
			for(UInt i=0; i<debugmatrix_.size();++i)
		{
		myfile<<debugmatrix_[i][0] << " " << debugmatrix_[i][1] << " " << debugmatrix_[i][4] << std::endl; 
	}
		myfile.close();
			*/
		debugmatrix_.clear();
		debugtraceback_.clear();
		debugscorematrix_.clear();	
		
			
	}
	void MapAlignmentAlgorithmSpectrumAlignment::updateMembers_()
	{
		gap_	=(Real)param_.getValue("gapcost");
		e_		=(Real)param_.getValue("affinegapcost");
		if(c1_ ==NULL ||c1_->getName()!=(String)param_.getValue("scorefunction"))
		{
			c1_ = Factory<PeakSpectrumCompareFunctor>::create((String)param_.getValue("scorefunction"));
		}
		
		cutoffScore_=(Real)param_.getValue("cutoffScore");
		bucketsize_=(Int)param_.getValue("bucketsize");
		mismatchscore_=(Real)param_.getValue("mismatchscore");
		anchorPoints_=(Int)param_.getValue("anchorpoints");
		if(anchorPoints_ > 100) anchorPoints_ =100;
		String tmp=(String)param_.getValue("debug");
		if(tmp=="true")
		{
			debug_=true;
		}
		else
		{
			debug_=false;
		}
	}

} 



