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
//
#include <OpenMS/COMPARISON/CLUSTERING/SpectrumGenerator.h>

#include <OpenMS/METADATA/Identification.h>

#include <cmath>

using namespace std;

namespace OpenMS
{
  SpectrumGenerator* SpectrumGenerator::instance_ = 0;
  
  SpectrumGenerator::SpectrumGenerator()
  {
    isotopes_ = 1;
    neutrallosses_ = 0; 
    low_range_cutoff_ = 200;
    
    residuemasses_['G'] = 57.02;
    residuemasses_['A'] = 71.04;
    residuemasses_['S'] = 87.03;
    residuemasses_['P'] = 97.05;
    residuemasses_['V'] = 99.07;
    residuemasses_['T'] = 101.05;
    residuemasses_['C'] = 103.01;
    residuemasses_['L'] = 113.08;
    residuemasses_['I'] = 113.08;
    residuemasses_['X'] = 113.08; // L or I ( sequest uses this sometimes )
    residuemasses_['N'] = 114.04;
    residuemasses_['B'] = 114.54; // N and D (avg)
    residuemasses_['D'] = 115.03;
    residuemasses_['Q'] = 128.06;
    residuemasses_['K'] = 128.09;
    residuemasses_['Z'] = 128.55; // Q and E (avg)
    residuemasses_['E'] = 129.04;
    residuemasses_['M'] = 131.04;
    residuemasses_['H'] = 137.06;
    residuemasses_['U'] = 149.91; // selenocysteine
    residuemasses_['F'] = 147.07;
    residuemasses_['R'] = 156.1;
    residuemasses_['Y'] = 163.06;
    residuemasses_['W'] = 186.08;
    modifications_['*'].insert(make_pair('M',16)); //oxidized methionine
    modifications_['#'].insert(make_pair('C',71.03)); // acrylocysteine
   
    // values from Gilis et al 2001
    residuefrequency_.insert(make_pair('X',0));
    residuefrequency_.insert(make_pair('B',0));
    residuefrequency_.insert(make_pair('Z',0));
    residuefrequency_.insert(make_pair('A',0.0780));
    residuefrequency_.insert(make_pair('R',0.0523));
    residuefrequency_.insert(make_pair('D',0.0519));
    residuefrequency_.insert(make_pair('N',0.0437));
    residuefrequency_.insert(make_pair('C',0.0110));
    residuefrequency_.insert(make_pair('E',0.0672));
    residuefrequency_.insert(make_pair('Q',0.0345));
    residuefrequency_.insert(make_pair('G',0.0677));
    residuefrequency_.insert(make_pair('H',0.0203));
    residuefrequency_.insert(make_pair('I',0.0695));
    residuefrequency_.insert(make_pair('L',0.1015));
    residuefrequency_.insert(make_pair('K',0.0632));
    residuefrequency_.insert(make_pair('M',0.0229));
    residuefrequency_.insert(make_pair('F',0.0439));
    residuefrequency_.insert(make_pair('P',0.0426));
    residuefrequency_.insert(make_pair('S',0.0646));
    residuefrequency_.insert(make_pair('T',0.0512));
    residuefrequency_.insert(make_pair('W',0.0109));
    residuefrequency_.insert(make_pair('Y',0.0330));
    residuefrequency_.insert(make_pair('V',0.0701));
    // median amino acid
    residuefrequency_.insert(make_pair('u',0));
    
    // values from www.ionsource.com/tutorial/isotopes/comp.htm
    elements_per_residue_['A'].insert(make_pair('C',3));
    elements_per_residue_['A'].insert(make_pair('H',5));
    elements_per_residue_['A'].insert(make_pair('N',1));
    elements_per_residue_['A'].insert(make_pair('O',1));
    elements_per_residue_['A'].insert(make_pair('S',0));
    elements_per_residue_['R'].insert(make_pair('C',6));
    elements_per_residue_['R'].insert(make_pair('H',12));
    elements_per_residue_['R'].insert(make_pair('N',4));
    elements_per_residue_['R'].insert(make_pair('O',1));
    elements_per_residue_['R'].insert(make_pair('S',0));
    elements_per_residue_['D'].insert(make_pair('C',4));
    elements_per_residue_['D'].insert(make_pair('H',6));
    elements_per_residue_['D'].insert(make_pair('N',2));
    elements_per_residue_['D'].insert(make_pair('O',2));
    elements_per_residue_['D'].insert(make_pair('S',0));
    elements_per_residue_['N'].insert(make_pair('C',4));
    elements_per_residue_['N'].insert(make_pair('H',5));
    elements_per_residue_['N'].insert(make_pair('N',1));
    elements_per_residue_['N'].insert(make_pair('O',3));
    elements_per_residue_['N'].insert(make_pair('S',0));
    elements_per_residue_['C'].insert(make_pair('C',3));
    elements_per_residue_['C'].insert(make_pair('H',4));
    elements_per_residue_['C'].insert(make_pair('N',1));
    elements_per_residue_['C'].insert(make_pair('O',1));
    elements_per_residue_['C'].insert(make_pair('S',1));
    elements_per_residue_['E'].insert(make_pair('C',5));
    elements_per_residue_['E'].insert(make_pair('H',7));
    elements_per_residue_['E'].insert(make_pair('N',1));
    elements_per_residue_['E'].insert(make_pair('O',3));
    elements_per_residue_['E'].insert(make_pair('S',0));
    elements_per_residue_['Q'].insert(make_pair('C',5));
    elements_per_residue_['Q'].insert(make_pair('H',8));
    elements_per_residue_['Q'].insert(make_pair('N',2));
    elements_per_residue_['Q'].insert(make_pair('O',2));
    elements_per_residue_['Q'].insert(make_pair('S',0));
    elements_per_residue_['G'].insert(make_pair('C',2));
    elements_per_residue_['G'].insert(make_pair('H',3));
    elements_per_residue_['G'].insert(make_pair('N',1));
    elements_per_residue_['G'].insert(make_pair('O',1));
    elements_per_residue_['G'].insert(make_pair('S',0));
    elements_per_residue_['H'].insert(make_pair('C',6));
    elements_per_residue_['H'].insert(make_pair('H',7));
    elements_per_residue_['H'].insert(make_pair('N',3));
    elements_per_residue_['H'].insert(make_pair('O',1));
    elements_per_residue_['H'].insert(make_pair('S',0));
    elements_per_residue_['I'].insert(make_pair('C',6));
    elements_per_residue_['I'].insert(make_pair('H',11));
    elements_per_residue_['I'].insert(make_pair('N',1));
    elements_per_residue_['I'].insert(make_pair('O',1));
    elements_per_residue_['I'].insert(make_pair('S',0));
    elements_per_residue_['L'].insert(make_pair('C',6));
    elements_per_residue_['L'].insert(make_pair('H',11));
    elements_per_residue_['L'].insert(make_pair('N',1));
    elements_per_residue_['L'].insert(make_pair('O',1));
    elements_per_residue_['I'].insert(make_pair('S',0));
    elements_per_residue_['L'].insert(make_pair('C',6));
    elements_per_residue_['L'].insert(make_pair('H',11));
    elements_per_residue_['L'].insert(make_pair('N',1));
    elements_per_residue_['L'].insert(make_pair('O',1));
    elements_per_residue_['L'].insert(make_pair('S',0));
    elements_per_residue_['K'].insert(make_pair('C',6));
    elements_per_residue_['K'].insert(make_pair('H',12));
    elements_per_residue_['K'].insert(make_pair('N',2));
    elements_per_residue_['K'].insert(make_pair('O',1));
    elements_per_residue_['K'].insert(make_pair('S',0));
    elements_per_residue_['M'].insert(make_pair('C',5));
    elements_per_residue_['M'].insert(make_pair('H',9));
    elements_per_residue_['M'].insert(make_pair('N',1));
    elements_per_residue_['M'].insert(make_pair('O',1));
    elements_per_residue_['M'].insert(make_pair('S',1));
    elements_per_residue_['F'].insert(make_pair('C',5));
    elements_per_residue_['F'].insert(make_pair('H',9));
    elements_per_residue_['F'].insert(make_pair('N',1));
    elements_per_residue_['F'].insert(make_pair('O',1));
    elements_per_residue_['F'].insert(make_pair('S',0));
    elements_per_residue_['P'].insert(make_pair('C',5));
    elements_per_residue_['P'].insert(make_pair('H',7));
    elements_per_residue_['P'].insert(make_pair('N',1));
    elements_per_residue_['P'].insert(make_pair('O',1));
    elements_per_residue_['P'].insert(make_pair('S',0));
    elements_per_residue_['S'].insert(make_pair('C',3));
    elements_per_residue_['S'].insert(make_pair('H',5));
    elements_per_residue_['S'].insert(make_pair('N',1));
    elements_per_residue_['S'].insert(make_pair('O',2));
    elements_per_residue_['S'].insert(make_pair('S',0));
    elements_per_residue_['T'].insert(make_pair('C',4));
    elements_per_residue_['T'].insert(make_pair('H',7));
    elements_per_residue_['T'].insert(make_pair('N',1));
    elements_per_residue_['T'].insert(make_pair('O',2));
    elements_per_residue_['T'].insert(make_pair('S',0));
    elements_per_residue_['W'].insert(make_pair('C',11));
    elements_per_residue_['W'].insert(make_pair('H',10));
    elements_per_residue_['W'].insert(make_pair('N',2));
    elements_per_residue_['W'].insert(make_pair('O',1));
    elements_per_residue_['W'].insert(make_pair('S',0));
    elements_per_residue_['Y'].insert(make_pair('C',9));
    elements_per_residue_['Y'].insert(make_pair('H',9));
    elements_per_residue_['Y'].insert(make_pair('N',1));
    elements_per_residue_['Y'].insert(make_pair('O',2));
    elements_per_residue_['Y'].insert(make_pair('S',0));
    elements_per_residue_['V'].insert(make_pair('C',5));
    elements_per_residue_['V'].insert(make_pair('H',9));
    elements_per_residue_['V'].insert(make_pair('N',1));
    elements_per_residue_['V'].insert(make_pair('O',1));
    elements_per_residue_['V'].insert(make_pair('S',0));
    
    // calculate monoisotopic residue mass for median amino acid 'u'
    for ( map<char,double>::const_iterator cmit = residuemasses_.begin(); cmit != residuemasses_.end(); ++cmit )
    {
      residuemasses_['u'] += residuefrequency_[cmit->first] * cmit->second;
    }

    // calculate elements per residue for the median amino acid 'u'
    for ( map<char,map<char,double> >::const_iterator cmit = elements_per_residue_.begin(); cmit != elements_per_residue_.end(); ++cmit )
    {
      for ( map<char,double>::const_iterator cmit2 = cmit->second.begin(); cmit2 != cmit->second.end(); ++cmit2 )
      {
        elements_per_residue_['u'][cmit2->first] += cmit2->second * residuefrequency_[cmit->first];
      }
    }
    
    // isotopeprob_[AA][0] = light isotope, isotopeprob_[AA][1] = isotope with +1 neutron ... to +3 neutron
    // values from http:://physics.nist.gov/PhysRefData/Compositions/index.html
    isotopeprob_['H'].push_back(0.999885);
    isotopeprob_['H'].push_back(0.000115);
    isotopeprob_['H'].push_back(0);
    isotopeprob_['H'].push_back(0);

    isotopeprob_['C'].push_back(0.9893);
    isotopeprob_['C'].push_back(0.0107);
    isotopeprob_['C'].push_back(0);
    isotopeprob_['C'].push_back(0);
    
    isotopeprob_['N'].push_back(0.99632);
    isotopeprob_['N'].push_back(0.00368);
    isotopeprob_['N'].push_back(0);
    isotopeprob_['N'].push_back(0);
    
    isotopeprob_['O'].push_back(0.99757);
    isotopeprob_['O'].push_back(0.00038);
    isotopeprob_['O'].push_back(0.00205);
    isotopeprob_['O'].push_back(0);
    
    isotopeprob_['S'].push_back(0.9493);
    isotopeprob_['S'].push_back(0.0076);
    isotopeprob_['S'].push_back(0.0429);
    isotopeprob_['S'].push_back(0.0002);
  }

  SpectrumGenerator* SpectrumGenerator::instance()
  {
    if (!instance_)
    {
      instance_ = new SpectrumGenerator();
    }
    return instance_;
  }

  /**
  \param seq peptide sequence
  \return peptidemass
  */
  double SpectrumGenerator::getPeptidemass(const String& seq) const
  {
    double mass = 18; // C and N- Terminal H and OH 
    for (String::const_iterator csit = seq.begin(); csit != seq.end(); ++csit)
    {
      // modification
      mass += SpectrumGenerator::instance()->residuemass(csit,seq);
    }
    return mass;
  }
          
  /**
  the modifications list is possibly not complete! <br>
  here is only what was encountered during minimal usage <br>
  \param sit iterator where amino acid lies
  \param seq String where sit points to
  \return mass ( considers modifications which are )
  */
  double SpectrumGenerator::residuemass(String::const_iterator& sit, const String& seq) const
  {

    map<char,double>::const_iterator cmit = residuemasses_.find(*sit);
    double mass = 0;
    if ( cmit == residuemasses_.end() ) 
    {
      cerr << "unknown amino acid:  " << *sit << " in " << seq << endl;
    }
    else
    {
      mass = cmit->second;
    }
    if ((sit+1) != seq.end() && modifications_.find(*(sit+1)) != modifications_.end() )
    {
      //debug
      if (modifications_.find(*(sit+1))->second.find(*sit) == modifications_.find(*(sit+1))->second.end())
      {
        char message[2];
        message[0] = *sit;
        message[1] = 0;
        throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"unknown modification",message);
      }
      mass += modifications_.find(*(sit+1))->second.find(*sit)->second;
      ++sit;
    }
    return mass;
  }

  /**
  \param seq peptide sequence
  \return theoretical stick spectrum
  */
  MSSpectrum< DPeak<1> >* SpectrumGenerator::getspectrum(const String& seq) const
  {
    // current method only for tryptic digest, creates spectra from peptides
    // with charge state 2, 1 H+ C-Terminal, 1 mobile
    // resulting fragments have charge +1
    MSSpectrum< DPeak<1> >* specp = new MSSpectrum< DPeak<1> >();
    // sequest dta have mass of charge state 1
    double mass = getPeptidemass(seq)+2;
    map<double,pair<double,String> > fragments;
    
    specp->getPrecursorPeak().getPosition()[0] = mass;
    specp->getPrecursorPeak().getCharge() = 2;  
    
    double bion = 1;   
    double yion = mass -1 ; 
    
    // Values from Tabb anal chem 2003
   
    for(String::const_iterator csit = seq.begin(); csit != seq.end(); ++csit)
    {
      double aamass = this->residuemass(csit,seq);
      bion += aamass;
      yion -= aamass;

      // b-ion
      fragments[bion].first += 46;
      fragments[bion].second += "b_ion";
      // y-ion
      fragments[yion].first += 100;
      fragments[bion].second += "y_ion";

      if ( neutrallosses_ )
      {
        // neutral losses
        // ammonia
        fragments[bion-17].first += 46/2;
        fragments[bion-17].second += "ammonium_loss_from_b_ion";
        fragments[yion-17].first += 100/5;
        fragments[yion-17].second += "ammonium_loss_from_y_ion";
        // water
        fragments[bion-18].first += 46/3;
        fragments[bion-18].second += "water_loss_from_b_ion";
        fragments[yion-18].first += 100/10;
        fragments[yion-18].second += "water_loss_from_y_ion";
      }
    }

    for (map<double,pair<double,String> >::const_iterator cmit = fragments.begin(); cmit != fragments.end(); ++cmit)
    {
      double mz = cmit->first;
      double intensity = cmit->second.first;
      if ( mz < low_range_cutoff_ ) continue;
      if ( isotopes_ ) 
      {
        vector<pair<double,double> > peaks = isotopepeaks(mz,intensity);
        uint i = 0;
        for ( vector<pair<double,double> >::const_iterator cvit = peaks.begin(); cvit != peaks.end(); ++cvit )
        {
          DPeak<1> peak;
          stringstream ss;
          ss << cmit->second.second << "_iso_" << i++;
          peak.setMetaValue("peaktype",ss.str());
          peak.getPosition()[0] = cvit->first;
          peak.getIntensity() = cvit->second;
          specp->getContainer().push_back(peak);
        }
      }
      else
      {
        DPeak<1> peak;
        peak.setMetaValue("peaktype",cmit->second.second);
        peak.getPosition()[0] = mz;
        peak.getIntensity() = intensity;
        specp->getContainer().push_back(peak);
      }
    }
  
    Identification dbs;
    dbs.insertPeptideHit(PeptideHit(0.5,"probability", 1, seq));
    specp->getIdentifications().push_back(dbs);
    
    return specp;
  }
  
  /**
  \param mz position of the peak
  \param intensity of the peak
  \return theoretical isotope peaks, based on the probability of atom isotopes and a the medium atom composition of an amino acid
  */
  std::vector<std::pair<double,double> > SpectrumGenerator::isotopepeaks(double mz, double intensity) const
  {
    // how many median AAs 'u' fit into the fragment?
    double nru = mz/residuemasses_.find('u')->second;

    map<char,double> nrelements;
    nrelements['H'] = elements_per_residue_.find('u')->second.find('H')->second * nru;
    nrelements['C'] = elements_per_residue_.find('u')->second.find('C')->second * nru;
    nrelements['S'] = elements_per_residue_.find('u')->second.find('S')->second * nru;
    nrelements['N'] = elements_per_residue_.find('u')->second.find('N')->second * nru;
    nrelements['N'] = elements_per_residue_.find('u')->second.find('N')->second * nru;

    // calculate frequencies of isotopes 0...3
    vector<double> isotopefreqs;
    // calculate frequency of isotope 0
    double iso0 = 1;
    for ( map<char,double>::const_iterator cmit = elements_per_residue_.find('u')->second.begin(); cmit != elements_per_residue_.find('u')->second.end(); ++cmit )
    {
      iso0 *= pow(isotopeprob_.find(cmit->first)->second[0],nrelements[cmit->first]);
    }
    isotopefreqs.push_back(iso0);
   
    // calculate frequency of isotope 1
    double iso1 = 0;
    for ( map<char,double>::const_iterator cmit = elements_per_residue_.find('u')->second.begin(); cmit != elements_per_residue_.find('u')->second.end(); ++cmit )
    {
      double isotemp = isotopeprob_.find(cmit->first)->second[1];
      isotemp *= pow(isotopeprob_.find(cmit->first)->second[0],nrelements[cmit->first] -1 );
      iso1 += isotemp * nrelements[cmit->first];
    }
    isotopefreqs.push_back(iso1);
    
    // calculate frequency of isotope 2
    double iso2 = 0;
    for ( map<char,double>::const_iterator cmit = elements_per_residue_.find('u')->second.begin(); cmit != elements_per_residue_.find('u')->second.end(); ++cmit )
    {
      double isotemp = pow(isotopeprob_.find(cmit->first)->second[1],2);
      isotemp *= pow(isotopeprob_.find(cmit->first)->second[0],nrelements[cmit->first] -2 );
      iso2 += isotemp * ((nrelements[cmit->first]*(nrelements[cmit->first]-1))/2);
      isotemp = isotopeprob_.find(cmit->first)->second[2];
      isotemp *= pow(isotopeprob_.find(cmit->first)->second[0],nrelements[cmit->first] -1 );
      iso2 += isotemp * nrelements[cmit->first];
    }
    isotopefreqs.push_back(iso2);

    vector<pair<double,double> > result;
    double isosum = 0;
    for ( uint i = 0; i < isotopefreqs.size(); ++i )
    {
      isosum += isotopefreqs[i];
    }
    // isotope intensities are not scaled to look like in real spectra
    for ( uint i = 0; i < isotopefreqs.size(); ++i )
    {
      result.push_back(make_pair(mz+i,(isotopefreqs[i]/isosum)*intensity));
    }
    return result;
  }
                  
}

