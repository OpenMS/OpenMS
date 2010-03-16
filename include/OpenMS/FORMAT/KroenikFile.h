// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_KROENIKFILE_H
#define OPENMS_FORMAT_KROENIKFILE_H

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <fstream>
#include <vector>

namespace OpenMS
{
 	/**
 		@brief File adapter for Kroenik (HardKloer sibling) files.
 		
  	The first line is the header and contains the column names:<br>
    File	First Scan	Last Scan	Num of Scans	Charge	Monoisotopic Mass	Base Isotope Peak	Best Intensity	Summed Intensity	First RTime	Last RTime	Best RTime	Best Correlation	Modifications

    Every subsequent line is a feature.
  	
  	@ingroup FileIO
  */
  class OPENMS_DLLAPI KroenikFile
  {
    public:
      /// Default constructor
      KroenikFile();
			/// Destructor
      virtual ~KroenikFile();
      
      /**
 				@brief Loads a Kroenik file into a featureXML.
 				
 				The content of the file is stored in @p features.

				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
      */
      template <typename FeatureMapType>
      void load(const String& filename, FeatureMapType& feature_map)
      {
        // load input
				TextFile input(filename);
		
				// reset map
        FeatureMapType fmap;
				feature_map = fmap;
			
				for (Size i=1; i<input.size(); ++i)
				{
					String line = input[i];
				
					//split lines: File,	First Scan,	Last Scan,	Num of Scans,	Charge,	Monoisotopic Mass, Base Isotope Peak,	Best Intensity,	Summed Intensity,	First RTime,	Last RTime,	Best RTime,	Best Correlation,	Modifications
					std::vector< String > parts;
					line.split('\t', parts);
					
					if (parts.size() != 14)
					{
						throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", 
              String("Failed parsing in line ") + String(i+1) + ": missing 14 tab-separated entries (got " + String(parts.size()) + ")\nLine was: " + line);
					}
					//create feature
					Feature f;
					f.setCharge(parts[4].toInt());
					f.setMZ(parts[5].toDouble()/f.getCharge() + Constants::PROTON_MASS_U);
					f.setRT(parts[11].toDouble());
					f.setOverallQuality(parts[12].toDouble());
					f.setIntensity(parts[8].toDouble());
					ConvexHull2D hull;
					ConvexHull2D::PointType point;
					
					point.setX(parts[9].toDouble());
					point.setY(f.getMZ());
					hull.addPoint(point);

					point.setX(parts[9].toDouble());
					point.setY(f.getMZ()+3.0/(DoubleReal)f.getCharge());
					hull.addPoint(point);

					point.setX(parts[10].toDouble());
					point.setY(f.getMZ()+3.0/(DoubleReal)f.getCharge());
					hull.addPoint(point);

					point.setX(parts[10].toDouble());
					point.setY(f.getMZ());
					hull.addPoint(point);
					
					point.setX(parts[9].toDouble());
					point.setY(f.getMZ());
					hull.addPoint(point);
					
					std::vector< ConvexHull2D > hulls;
					hulls.push_back(hull);
					f.setConvexHulls(hulls);
					f.setMetaValue("Mass",parts[5].toDouble());
					f.setMetaValue("FirstScan",parts[1].toDouble());
					f.setMetaValue("LastScan",parts[2].toInt());
					f.setMetaValue("NumOfScans",parts[3].toDouble());
					f.setMetaValue("AveragineModifications",parts[13]);
					feature_map.push_back(f);
				}

				LOG_INFO << "Hint: The convex hulls are approximated in m/z dimension (Kroenik lacks this information)!\n";
      }

      /**
      	@brief Stores a featureXML as a Kroenik file.
      	
        NOT IMPLEMENTED

				@exception Exception::UnableToCreateFile is thrown if the file could not be created
      */
      template <typename SpectrumType>
      void store(const String& filename, const SpectrumType& spectrum) const
      {
        throw Exception::NotImplemented (__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }
  };
} // namespace OpenMS

#endif // OPENMS_FORMAT_KROENIKFILE_H

