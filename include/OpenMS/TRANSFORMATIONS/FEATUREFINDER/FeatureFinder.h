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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

namespace OpenMS
{
  class BaseModelFitter;

	/**
		@defgroup FeatureFinder FeatureFinder
		
		@brief The feature detection algorithm in OpenMS.
		
		@ingroup Transformations
	*/
	
	/**
		@brief The base class of the Feature Finding algorithm.
		
		By feature finding, we understand the search for 
		peptides or other chemical compounds different from noise
		in a MS map. 
		
		Our algorithm consists of three different stages: seeding,
		extension and model fitting. In the seeding stage, we mark
		data points with high intensity as seeds. During the extension
		phase, these seeds are extended by collecting other data points
		that lie within a certain distance from the seed and have an intensity
		significantly larger than the background noise in the data set.
		The seed together with the data points collected during its
		extension comprise the region of a potential feature.
		
		Finally, we fit a hypothetical model of the compound we want
		to detect to the feature region. If the correlation between the
		distribution of the data points and the model is too low, the
		feature is discarded. Furthermore, we can remove datapoints
		having a very low probability given the model from the feature
		region.
		
		This class reads the parameters from the Param class,
		initialises the traits class and the seeder, extender and
		model fitting classes.
		
		@note The input data should be sorted and not contain any tandem spectra.
	  	
		@ingroup FeatureFinder
	*/
	class FeatureFinder
	{
		public:
	    typedef std::vector<BaseSeeder*> SeederVector;
	    typedef std::vector<BaseExtender*> ExtenderVector;
	    typedef std::vector<BaseModelFitter*> FitterVector;
			typedef FeaFiTraits::MapType MapType;
	
	    /// Default constructor.
	    FeatureFinder();
	    /// copy constructor
	    FeatureFinder(const FeatureFinder& source);
	    /// destructor
	    virtual ~FeatureFinder();
	    /// assignment operator
	    FeatureFinder& operator = (const FeatureFinder& source);
	    /// test for equality
	    bool operator == (const FeatureFinder& rhs) const;
	
	    /**
				@brief add seeder to FeatureFinder  
				
				@param name unique name (see BaseSeeder) for Seeder class derived from BaseSeeder 
				@param param settings for the seeder
	    */
	    void addSeeder(const String& name, const Param* param=0);
	    /**
		 		@brief add extender to FeatureFinder  
		 		
				@param name unique name (see BaseExtender) for Extender class derived from BaseExtender
				@param param settings for the extender
	    */
	    void addExtender(const String& name, const Param* param=0);
	
	    /**
	  		@brief add modelfitter to FeatureFinder  
	  		
	  		@param name unique name (see BaseFitter) for Fitter class derived from BaseFitter
	  		@param param settings for the fitter
	    */
	    void addFitter(const String& name, const Param* param=0);
	
	    ///	 remove first seeder with name @p name
	    void removeSeeder(const String& name);
	    /// remove first extender with name @p name
	    void removeExtender(const String& name);
	    /// remove first modelfitter with name @p name
	    void removeFitter(const String& name);
	
	    /**
	    	@brief set FeatureFinder param file  
	     	@return false if param file is not valid 
	    */
	    bool setParam(const Param& param);
	
	
	    /// Sets data using instance of MSExperimentExtern.
	    template <class SpectrumIteratorType>
	    void setData(const SpectrumIteratorType& begin, const SpectrumIteratorType& end, UnsignedInt buffer_size)
		  { 
	    	traits_->setData(begin, end, buffer_size); 
	    }
		
	    /// start feature finding
	    const DFeatureMap<2>& run();
	
	    friend std::ostream& operator << (std::ostream& os, const FeatureFinder& finder);
	
		protected:
	
	    inline bool setModule(const String& module)
	    {
	      Param mod_param = param_.copy(module+":",true);
	      if (mod_param.empty()) return false;
	      Param::ConstIterator it;
	      for (it=mod_param.begin(); it!=mod_param.end(); it++)
	      {
	        if (it->first.substr(it->first.size()-2,2) == "ID")
	        {
	          String impl_name = it->second;
	          std::vector<String> substrings;
	          String(it->first).split(':',substrings);
	          String node = substrings[substrings.size()-2];
	          Param sub_param;
	          sub_param = mod_param.copy(node+":",true);
	          sub_param.remove("ID");
	
	          if (module=="Seeders")
	          {
	          	addSeeder(impl_name, &sub_param);
	          }
	          else if (module=="Extenders")
	          {
	          	addExtender(impl_name, &sub_param);
	          }
	          else if (module=="ModelFitters")
	          {
	          	addFitter(impl_name, &sub_param);
	        	}
	        }
	      }
	      return true;
	    }
	
	    Param param_;
	    FeaFiTraits* traits_;
	    SeederVector seeders_;
	    ExtenderVector extenders_;
	    FitterVector fitters_;
	
	};

	/// Print the module names of a FeatureFinder object to a stream.
	std::ostream& operator << (std::ostream& os, const FeatureFinder& finder);

}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_H
