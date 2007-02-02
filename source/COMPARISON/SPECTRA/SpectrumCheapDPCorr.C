// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/SPECTRA/SpectrumCheapDPCorr.h>

#include <cmath>

#define SPECTRUMCHEAPDPCORR_DEBUG
#undef  SPECTRUMCHEAPDPCORR_DEBUG

#ifdef SPECTRUMCHEAPDPCORR_DEBUG
	#include <iostream>
#endif

using namespace std;

namespace OpenMS
{
  SpectrumCheapDPCorr::SpectrumCheapDPCorr()
    : PeakSpectrumCompareFunctor(),
			lastconsensus_()
  {
		setName(SpectrumCheapDPCorr::getProductName());
    defaults_.setValue("variation", 0.001);
    defaults_.setValue("int_cnt", 0); 
    defaults_.setValue("keeppeaks", 0);
    factor_ = 0.5;
		defaultsToParam_();
  }

  SpectrumCheapDPCorr::SpectrumCheapDPCorr(const SpectrumCheapDPCorr& source)
    : PeakSpectrumCompareFunctor(source),
			lastconsensus_(source.lastconsensus_),
			factor_(source.factor_)
  {
  }

  SpectrumCheapDPCorr::~SpectrumCheapDPCorr()
  {
  }

  SpectrumCheapDPCorr& SpectrumCheapDPCorr::operator = (const SpectrumCheapDPCorr& source)
  {
		if (this != &source)
		{
    	PeakSpectrumCompareFunctor::operator = (source);
    	lastconsensus_ = source.lastconsensus_;
			factor_ = source.factor_;
		}
    return *this;
  }

  void SpectrumCheapDPCorr::setFactor(double f)
  { 
    if ( f < 1 && f > 0 )
    {
      factor_ = f;
    }
    else
    {
			// TODO exception
      //cerr << "factor should be between 0 and 1, ignored\n";
    }
  }
  
 	double SpectrumCheapDPCorr::operator () (const PeakSpectrum& csa) const
	{
		return operator()(csa, csa);
	}
	
  /**
  looks for peak pairs where there is just one or none possibility for alignment
  and aligns them (if possible). The rest is aligned using dynprog_
  */
  double SpectrumCheapDPCorr::operator () (const PeakSpectrum& x, const PeakSpectrum& y) const
  { 
    double var = (double)param_.getValue("variation");
    double score(0);
    bool keeppeaks_ = (int)param_.getValue("keeppeaks");
    
    lastconsensus_ = PeakSpectrum();
    lastconsensus_.getPrecursorPeak().setPosition((x.getPrecursorPeak().getPosition()[0] + y.getPrecursorPeak().getPosition()[0]) / 2);
    lastconsensus_.getPrecursorPeak().setCharge(x.getPrecursorPeak().getCharge());
		peak_map_.clear();
    
    int xpos = 0;
    int ypos = 0;
		PeakSpectrum::ConstIterator xit = x.getContainer().begin();
		PeakSpectrum::ConstIterator yit = y.getContainer().begin();
    while (xit != x.getContainer().end() && yit != y.getContainer().end())
    {
      double variation = (xit->getPosition()[0] + yit->getPosition()[0]) / 2 * var; 

      //ignore pairs that cannot be paired
      if (fabs(xit->getPosition()[0] - yit->getPosition()[0]) > variation)
      {
        if (xit->getPosition()[0] < yit->getPosition()[0]) // while effizienter?
        {
          DPeak<1> consensuspeak;
          consensuspeak.getPosition()[0] = xit->getPosition()[0];
          consensuspeak.getIntensity() = (xit->getIntensity()) * (1 - factor_);
          if (keeppeaks_) lastconsensus_.getContainer().push_back(consensuspeak);
          ++xit;
          ++xpos;
        }
        else 
        {
          DPeak<1> consensuspeak;
          consensuspeak.getPosition()[0] = yit->getPosition()[0] ;
          consensuspeak.getIntensity()  = (yit->getIntensity()) * (factor_);
          if (keeppeaks_) lastconsensus_.getContainer().push_back(consensuspeak);
          ++yit;
          ++ypos;
        }
      }
      else 
      {
        //x/yrun represents the number of peaks in both spectra that could be paired
        int xrun = 1;
        int yrun = 1;
        while (xit + xrun != x.getContainer().end() && yit + yrun != y.getContainer().end() && 
            (!((xit + xrun - 1)->getPosition()[0] + variation < (yit + yrun)->getPosition()[0]) || 
             !((yit + yrun - 1)->getPosition()[0] + variation < (xit + xrun)->getPosition()[0])))
        {
          if ((yit + yrun - 1)->getPosition()[0] + variation > (xit + xrun)->getPosition()[0])
          {
            xrun++;
          }
          else if ((xit + xrun - 1)->getPosition()[0] + variation > (yit + yrun)->getPosition()[0])
          {
            yrun++;
          }
          else 
          {
            xrun++;
            yrun++;
          }
          if (xit + xrun == x.getContainer().end()) break;
          if (yit + yrun == y.getContainer().end()) break;
        }
        
        //dynamic programming necessary to calculate optimal pairing
        if (xrun > 1 && yrun > 1)
        {
          score += dynprog_(x, y, xpos, xpos + xrun - 1, ypos, ypos + yrun - 1);
          xit = xit + xrun;
          yit = yit + yrun;
          xpos += xrun;
          ypos += yrun;
        }
        //the optimal pairing of 2 peaks is easy...
        else 
        {
          // calculate consensus peak 
          DPeak<1> consensuspeak;
          consensuspeak.getPosition()[0] = (xit->getPosition()[0] * (1 - factor_) + yit->getPosition()[0] * (factor_));
          consensuspeak.getIntensity()  = (xit->getIntensity() * (1 - factor_) + yit->getIntensity() * factor_);
          lastconsensus_.getContainer().push_back(consensuspeak);
					
					if (!peak_map_.has(xit - x.getContainer().begin()))
					{
						peak_map_[xit - x.getContainer().begin()] = yit - y.getContainer().begin();
					}
					else
					{
						peak_map_[xit - x.getContainer().begin()] = yit - y.getContainer().begin() > xit - x.getContainer().begin() ? xit - x.getContainer().begin() : yit - y.getContainer().begin();
					}
          
          variation = (xit->getPosition()[0] + yit->getPosition()[0]) / 2 * var;
          score += comparepeaks_(xit->getPosition()[0], yit->getPosition()[0], xit->getIntensity(), yit->getIntensity());
          ++xit;
          ++yit;
          ++xpos;
          ++ypos;
        }
      }
    }
    factor_ = 0.5;
    return score;
  }

  double SpectrumCheapDPCorr::dynprog_(const PeakSpectrum& x, const PeakSpectrum& y, int xstart, int xend, int ystart, int yend) const
  {
		#ifdef SPECTRUMCHEAPDPCORR_DEBUG
		cerr << "SpectrumCheapDPCorr::dynprog_(const DDiscreteSpectrum<1>& x, const DDiscreteSpectrum<1>& y, " << xstart << ", "<< xend << ", " <<  ystart << ", " << yend << ")" <<  endl;
		#endif
    double var = (double)param_.getValue("variation");
    vector<vector<double> > dparray(xend - xstart + 2, vector<double>(yend - ystart + 2));
    vector<vector<int> > trace(xend - xstart + 2, vector<int>(yend - ystart + 2));
    double align;
    for (int i = 1; i < xend - xstart+2 ; ++i)
    {
      for (int j = 1; j < yend - ystart+2; ++j)
      {
        double variation = (y.getContainer()[ystart + j - 1].getPosition()[0] + x.getContainer()[xstart + i - 1].getPosition()[0]) / 2 * var;
        //positions too different
        if (fabs( x.getContainer()[xstart + i - 1].getPosition()[0] - y.getContainer()[ystart + j - 1].getPosition()[0]) > variation ) align = 0;

        //calculate score of alignment
        else  align = comparepeaks_(x.getContainer()[xstart + i - 1].getPosition()[0], 
																		y.getContainer()[ystart + j - 1].getPosition()[0], 
																		x.getContainer()[xstart + i - 1].getIntensity(), 
																		y.getContainer()[ystart + j - 1].getIntensity());
        //dynaminc programming step
        if ((((dparray[i][j-1]) > (dparray[i-1][j-1]+align)) ? (dparray[i][j-1]) : (dparray[i-1][j-1]+align)) /*== max*/ > dparray[i-1][j])
        {
          if (dparray[i-1][j-1] + align > dparray[i][j-1])
          {
            dparray[i][j] = dparray[i-1][j-1] + align;
            trace[i][j] = 5;
          }
          else
          {
            dparray[i][j] = dparray[i][j-1];
            trace[i][j] = -1;
          }
        }
        else 
        {
          dparray[i][j] = dparray[i-1][j];
          trace[i][j] = 1;
        }
      }
    }

    unsigned int i = xend-xstart+1;
    unsigned int j = yend-ystart+1;
    for (;;) 
    {
      if ( trace[i][j] == 5 ) 
      {
        DPeak<1> consensuspeak;
        consensuspeak.getPosition()[0] = ( y.getContainer()[ystart + j-1].getPosition()[0]*(1-factor_)+x.getContainer()[xstart + i-1].getPosition()[0]*factor_);
        consensuspeak.getIntensity() = ( y.getContainer()[ystart + j-1].getIntensity()*(1-factor_) + x.getContainer()[xstart + i-1].getIntensity()*factor_ );
        lastconsensus_.getContainer().push_back(consensuspeak);
				if (!peak_map_.has(xstart + i-1))
				{
					peak_map_[xstart + i-1] = ystart + j-1;
				}
				else
				{
					peak_map_[xstart + i-1] = ystart + j-1 > peak_map_[xstart + i-1] ? peak_map_[xstart + i-1] : ystart + j-1;
				}
        i--;
        j--;   
      }
      else if ( trace[i][j] == 1 )
      {
        DPeak<1> consensuspeak;
        consensuspeak.getPosition()[0] = x.getContainer()[xstart + i-1].getPosition()[0] ;
        consensuspeak.getIntensity()  = (x.getContainer()[xstart + i-1].getIntensity()) * ( 1-factor_ );
        if ( keeppeaks_ ) lastconsensus_.getContainer().push_back(consensuspeak);
        i--;
      }
      else if ( trace[i][j] == -1 )
      {
        DPeak<1> consensuspeak;
        consensuspeak.getPosition()[0] = y.getContainer()[ystart + j-1].getPosition()[0] ;
        consensuspeak.getIntensity()  = (y.getContainer()[ystart + j-1].getIntensity())*factor_;
        if ( keeppeaks_ ) lastconsensus_.getContainer().push_back(consensuspeak);
        j--;
      }
      if ( !i || !j ) break;
    }

    return dparray[xend-xstart+1][yend-ystart+1]; 
  }
 
  const PeakSpectrum& SpectrumCheapDPCorr::lastconsensus() const 
	{
 		return lastconsensus_;
	}

	HashMap<Size, Size> SpectrumCheapDPCorr::getPeakMap() const
	{
		return peak_map_;
	}

  /**
   @param posa position of peak a
   @param posb position of peak b
   @param inta intensity of peak a
   @param intb intensity of peak b
   @return score
   */
  double SpectrumCheapDPCorr::comparepeaks_(double posa, double posb, double inta, double intb) const
  {
    double variation = (posa + posb) / 2 * (double)param_.getValue("variation");
		unsigned int int_cnt = (unsigned int)param_.getValue("int_cnt");
    if (int_cnt == 0)
   	{
      return gsl_ran_gaussian_pdf(posa - posb, variation) * inta * intb;
    }
    else
		{
			if (int_cnt == 1)
    	{
      	return gsl_ran_gaussian_pdf(posa - posb, variation) * sqrt(inta * intb);
   		}
    	else 
			{
				if (int_cnt == 2)
    		{
      		return gsl_ran_gaussian_pdf(posa - posb, variation) * (inta + intb);
    		}
    		else 
				{
					if (int_cnt == 3)
    			{
      			return max(0.0, gsl_ran_gaussian_pdf(posa - posb, variation) * ((inta + intb) / 2 - fabs(inta - intb)));
    			}
    			else
    			{
						// TODO exception
      			// cerr << "int_cnt is not in [0,1,2,3]\n";
      			return -1;
    			}
				}
			}
		}
  }
}
