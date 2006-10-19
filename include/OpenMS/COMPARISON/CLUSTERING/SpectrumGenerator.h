// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer:  $
// --------------------------------------------------------------------------

#ifndef OPENMS_COMPARISON_CLUSTERING_SPECTRUMGENERATOR_H
#define OPENMS_COMPARISON_CLUSTERING_SPECTRUMGENERATOR_H

#include <OpenMS/KERNEL/MSSpectrum.h>

#include <vector>
#include <map>

namespace OpenMS
{
  /**
  	@brief SpectrumGenerator creates theoretical tandem-MS spectra
  */
  class SpectrumGenerator
  {
  public:
    /// Access to instance 
    static SpectrumGenerator* instance(); 
    
    /**
    	@brief 
    	
    	The modifications list is possibly not complete!
  		
  		Here is only what was encountered during minimal usage.
    	
	    @param csit iterator to amno acid character (possibly with modification)
	    @param seq amino acid sequence where @p csit points to
	    @return residuemass (considers modifications)
    */
    double residuemass(String::const_iterator& csit, const String& seq) const ;
    
    /**
    	@brief 
    	
			@param seq peptide sequence
    	@return mass of peptide
    */
    double getPeptidemass(const String& seq) const;
    
    /** 
    	@brief create simple theoretical spectrum
    	
    	@param seq peptide sequence
    	@return theoretical stick spectrum
    */
    MSSpectrum< DPeak<1> >* getSpectrum(const String& seq) const;

    /**
    	@brief 
    	
    	@param mz position of peak
    	@param intensity of peak
    	@return theoretical isotope peaks, based on the probability of atom isotopes and a the medium atom composition of an amino acid
    */
    std::vector<std::pair<double,double> > isotopepeaks(double mz, double intensity) const;
  
  private:
    ///@brief use instance
    SpectrumGenerator();
    
    
    static SpectrumGenerator* instance_;

    /// create isotope peaks?
    bool isotopes_;
    /// create peaks for neutral losses?
    bool neutrallosses_;

    /// monoisiotopic masses
    std::map<char,double> residuemasses_;

    /// frequency of the residues
    std::map<char,double> residuefrequency_;

    /// mass differences in modified residues
    std::map<char,std::map<char,double> > modifications_;

    /// element count for all residues
    std::map<char,std::map<char,double> > elements_per_residue_;

    /// probablilty for isotopes from elements
    std::map<char,std::vector<double> > isotopeprob_;

    /// dont create peaks below this mz threshold ( quadrupol )
    double low_range_cutoff_ ;
  };
}
#endif //OPENMS_COMPARISON_CLUSTERING_SPECTRUMGENERATOR_H
