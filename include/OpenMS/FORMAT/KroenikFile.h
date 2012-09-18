// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    File,  First Scan,  Last Scan,  Num of Scans,  Charge,  Monoisotopic Mass,  Base Isotope Peak,  Best Intensity,  Summed Intensity,  First RTime,  Last RTime,  Best RTime,  Best Correlation,  Modifications

    Every subsequent line is a feature.
    
    All properties in the file are converted to Feature properties, whereas "First Scan", "Last Scan", "Num of Scans" and "Modifications" are stored as 
    metavalues with the following names "FirstScan", "LastScan", "NumOfScans" and "AveragineModifications".

    The width in m/z of the overall convex hull of each feature is set to 3 Th in lack of a value provided by the Kroenik file.

    @note Kroenik files are Tab (\\t) separated files.
  	
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
				
					//split lines: File,  First Scan,  Last Scan,  Num of Scans,  Charge,  Monoisotopic Mass,  Base Isotope Peak,  Best Intensity,  Summed Intensity,  First RTime,  Last RTime,  Best RTime,  Best Correlation,  Modifications
					std::vector< String > parts;
					line.split('\t', parts);
					
					if (parts.size() != 14)
					{
						throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", 
              String("Failed parsing in line ") + String(i+1) + ": missing 14 tab-separated entries (got " + String(parts.size()) + ")\nLine was: '" + line + "'");
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
        std::cerr << "Store() for KroenikFile not implemented. Filename was: " << filename << ", spec of size " << spectrum.size() << "\n";
        throw Exception::NotImplemented (__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }
  };
} // namespace OpenMS

#endif // OPENMS_FORMAT_KROENIKFILE_H

