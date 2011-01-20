// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/Fitter1D.h>

// include derived classes here
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedIsotopeFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaIsotopeFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaGaussFitter1D.h>

namespace OpenMS
{
    Fitter1D::Fitter1D()
    : DefaultParamHandler("Fitter1D")
    {
      defaults_.setValue("interpolation_step",0.2,"Sampling rate for the interpolation of the model function.", StringList::create("advanced"));
      defaults_.setValue("statistics:mean",1.0,"Centroid position of the model.", StringList::create("advanced"));
      defaults_.setValue("statistics:variance",1.0,"The variance of the model.", StringList::create("advanced"));
      defaults_.setValue("tolerance_stdev_bounding_box",3.0,"Bounding box has range [minimim of data, maximum of data] enlarged by tolerance_stdev_bounding_box times the standard deviation of the data.", StringList::create("advanced"));
                        
      defaultsToParam_();
    }
    
    Fitter1D::Fitter1D(const Fitter1D& source)
    : DefaultParamHandler(source)
    {
      setParameters( source.getParameters() );
      updateMembers_();
    }

    Fitter1D& Fitter1D::operator = (const Fitter1D& source)
    {
      if (&source == this) return *this;
                
      DefaultParamHandler::operator = (source);
      setParameters( source.getParameters() );
      updateMembers_();
    
      return *this;
    }
    
    void Fitter1D::registerChildren()
    {
      Factory< Fitter1D >::registerProduct(GaussFitter1D::getProductName(), &GaussFitter1D::create);
      Factory< Fitter1D >::registerProduct(BiGaussFitter1D::getProductName(), &BiGaussFitter1D::create);
      Factory< Fitter1D >::registerProduct(IsotopeFitter1D::getProductName(), &IsotopeFitter1D::create);
      Factory< Fitter1D >::registerProduct(LmaIsotopeFitter1D::getProductName(),&LmaIsotopeFitter1D::create);
      Factory< Fitter1D >::registerProduct(ExtendedIsotopeFitter1D::getProductName(), &ExtendedIsotopeFitter1D::create);
		  Factory< Fitter1D >::registerProduct(EmgFitter1D::getProductName(), &EmgFitter1D::create);
      Factory< Fitter1D >::registerProduct(LmaGaussFitter1D::getProductName(), &LmaGaussFitter1D::create);
    }
    
    void Fitter1D::updateMembers_()
    {
      tolerance_stdev_box_ = param_.getValue( "tolerance_stdev_bounding_box" );
      interpolation_step_ = param_.getValue("interpolation_step");
      statistics_.setMean(param_.getValue("statistics:mean"));
      statistics_.setVariance(param_.getValue("statistics:variance"));
    }
          
} // namespace OpenMS

