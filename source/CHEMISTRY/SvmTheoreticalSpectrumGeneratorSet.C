// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Sandro Andreotti $
// $Authors: Sandro Andreotti $
// --------------------------------------------------------------------------

#include<OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGeneratorSet.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>



namespace OpenMS
{

  // Default constructor
  SvmTheoreticalSpectrumGeneratorSet::SvmTheoreticalSpectrumGeneratorSet(){}

  // Copy constructor
  SvmTheoreticalSpectrumGeneratorSet::SvmTheoreticalSpectrumGeneratorSet(const SvmTheoreticalSpectrumGeneratorSet& source):
      simulators_(source.simulators_)
  {}

  //Destructor
  SvmTheoreticalSpectrumGeneratorSet::~SvmTheoreticalSpectrumGeneratorSet()
  {}

  // Assignment operator
  SvmTheoreticalSpectrumGeneratorSet& SvmTheoreticalSpectrumGeneratorSet::operator =(const SvmTheoreticalSpectrumGeneratorSet& rhs)
  {
    if (this != &rhs)
    {
      simulators_=rhs.simulators_;
    }
    return *this;
  }

  // Generate the MS/MS according to the given probabilistic model
  void SvmTheoreticalSpectrumGeneratorSet::simulate(RichPeakSpectrum &spectrum, const AASequence &peptide, const gsl_rng *rng, Size precursor_charge)
  {
    std::map<Size, SvmTheoreticalSpectrumGenerator>::iterator it=simulators_.find(precursor_charge);
    if(it!=simulators_.end())
    {
      it->second.simulate(spectrum, peptide, rng, precursor_charge);
    }
    else
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Invalid Precursor charge, no Model available", String(precursor_charge));
    }
  }

  //Load a trained Svm and Prob. models
  void SvmTheoreticalSpectrumGeneratorSet::load(String filename)
  {
    if (! File::readable( filename ) )
    { // look in OPENMS_DATA_PATH
      filename = File::find( filename );
    }
    TextFile file(filename);

    Param sim_param = SvmTheoreticalSpectrumGenerator().getDefaults();
    for(Size line_num=1; line_num<file.size(); ++line_num)
    {
      String line(file[line_num]);
      std::vector<String>spl;
      line.split(":",spl);
      Int precursor_charge=spl[0].toInt();

      if(spl.size()!=2 || precursor_charge<1)
      {
        OpenMS::Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, file[line_num]," Invalid entry in SVM model File");
      }

      //load the model into the map
      sim_param.setValue("model_file_name", File::path(filename)+"/"+spl[1]);
      simulators_[precursor_charge].setParameters(sim_param);
      simulators_[precursor_charge].load();
    }
  }

  //Return precursor charges for which a model is contained in the set
  void SvmTheoreticalSpectrumGeneratorSet::getSupportedCharges(std::set<Size>&charges)
  {
    charges.clear();
    std::map<Size, SvmTheoreticalSpectrumGenerator>::const_iterator it;
    for(it=simulators_.begin(); it!=simulators_.end(); ++it)
    {
      charges.insert(it->first);
    }
  }

  //return a modifiable reference to the SVM model with given charge. If charge is not supported throw exception
  SvmTheoreticalSpectrumGenerator & SvmTheoreticalSpectrumGeneratorSet::getSvmModel(Size prec_charge)
  {
    std::map<Size, SvmTheoreticalSpectrumGenerator>::iterator it = simulators_.find(prec_charge);
    if(it==simulators_.end())
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Invalid Precursor charge, no Model available", String(prec_charge));
    }
    return it->second;
  }
}


