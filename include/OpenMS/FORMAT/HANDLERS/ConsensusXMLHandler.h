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
// $Maintainer: Eva Lange$
// --------------------------------------------------------------------------


#ifndef OPENMS_FORMAT_HANDLERS_CONSENSUSXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_CONSENSUSXMLHANDLER_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/StarAlignment.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLSchemes.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/ConsensusPeak.h>

#include <fstream>

// STL includes
#include <iostream>


namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief XML Handler for a consensusXML.
     
    */
    template < typename AlignmentT >
    class ConsensusXMLHandler
    	: public XMLHandler
    {
      public:
        typedef typename AlignmentT::ElementContainerType ElementContainerType;
        typedef typename AlignmentT::ElementType ElementType;
        typedef typename AlignmentT::ConsensusElementType ConsensusElementType;

      
        typedef typename GridCell::MappingVector MappingVector;
        typedef DPosition<2> PositionType;

      	/// Constructor
        ConsensusXMLHandler(ConsensusMap<ConsensusElementType>& consensus_map , const String& filename, bool load_elements_maps = true)
            : XMLHandler(filename),
            consensus_map_(&consensus_map),
            act_cons_element_(),
            calignment_(0),
            consensus_element_range_(false),
            feature_map_flag_(true),
            load_elements_maps_(load_elements_maps)
        {
        }

        /// Copy constructor
        ConsensusXMLHandler(const StarAlignment< ConsensusElementType >& alignment, const String& filename)
            : XMLHandler(filename),
            consensus_map_(0),
            act_cons_element_(),
            calignment_(&alignment),
            consensus_element_range_(false),
            feature_map_flag_(true)
        {
        }

        /// Destructor
        virtual ~ConsensusXMLHandler()
        {
        }

        // Docu in base class
        virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);

        // Docu in base class
        virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);

        // Docu in base class
        virtual void characters(const XMLCh* const chars, unsigned int length);

        ///Writes the contents to a stream
        void writeTo(std::ostream& os);



      protected:
        ConsensusMap<ConsensusElementType>* consensus_map_;
        ConsensusElementType act_cons_element_;
        // StartAlignment pointer for writing
        const StarAlignment< ConsensusElementType >* calignment_;
        bool consensus_element_range_;
        bool feature_map_flag_;
        bool consensus_map_flag_;
        bool load_elements_maps_;
        typename ConsensusElementType::PositionType pos_;
        typename ConsensusElementType::IntensityType it_;
        typename ConsensusElementType::PositionBoundingBoxType pos_range_;
        typename ConsensusElementType::IntensityBoundingBoxType it_range_;

        /// This function fills the members of a picked peak of type OutputPeakType.
        template <typename ConsensusElementT >
        void loadFile_(const String& /* file_name */, UInt /* id */, const ConsensusElementT& /* c */) throw (Exception::FileNotFound, Exception::ParseError)
        {
          std::cout << "Read no file " << std::endl;
        }

        void writeCellList_(std::ostream& os, const Grid& grid)
        {
          // write features with their attributes
          for (UInt s=0; s < grid.size(); s++)
          {
            const GridCell& cell = grid[s];

            os << "\t\t<cell nr=\"" << s << "\">" << std::endl;
            os << "\t\t\t\t<range rtMin=\"" << cell.min()[0] << "\" rtMax=\"" << cell.max()[0]
            << "\" mzMin=\"" << cell.min()[1] << "\" mzMax=\"" << cell.max()[1] << "\"/>\n";

            os << "\t\t\t\t<mappinglist>\n";
            MappingVector mappings = cell.getMappings();

            typename MappingVector::const_iterator citer = mappings.begin();
            UInt i = 0;
            while (citer != mappings.end() )
            {

              if (i==0)
              {
                os << "\t\t\t\t\t<rtMapping name=\"" << (*citer)->getName() << "\">\n";
              }
              else
                if (i == 1)
                {
                  os << "\t\t\t\t\t<mzMapping name=\"" << (*citer)->getName() << "\">\n";
                }

              Param map_param = (*citer)->getParameters();
              Param::ParamIterator piter = map_param.begin();
              while (piter != map_param.end())
              {
                os << "\t\t\t\t\t\t<param name=\"" << piter.getName() << "\" value=\"" << piter->value << "\"/>\n";
                piter++;
              }

              if (i==0)
              {
                os << "\t\t\t\t\t</rtMapping>\n";
              }
              else
                if (i == 1)
                {
                  os << "\t\t\t\t\t</mzMapping>\n";
                }

              citer++;
              ++i;
            }

            os << "\t\t\t\t</mappinglist>\n";
            os << "\t\t\t</cell>\n";

          } // end for ( features )
        }
    };

    //--------------------------------------------------------------------------------

    template < typename  AlignmentT >
    void ConsensusXMLHandler<AlignmentT>::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
    	//std::cout << "END: " << sm_.convert(qname) << std::endl;
    	static XMLCh* s_consensuselement = xercesc::XMLString::transcode("consensusElement");
      
      if (equal(qname,s_consensuselement))
      {
				consensus_map_->push_back(act_cons_element_);
      }
    }

    template < typename AlignmentT >
    void ConsensusXMLHandler<AlignmentT>::characters(const XMLCh* const /*chars*/, unsigned int /*length*/)
  	{
  	}

    template < typename  AlignmentT >
    void ConsensusXMLHandler<AlignmentT>::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
    {
      //std::cout << "BEGIN: " << sm_.convert(qname) << std::endl;
      static XMLCh* s_consensuselement = xercesc::XMLString::transcode("consensusElement");
      static XMLCh* s_count = xercesc::XMLString::transcode("count");
      static XMLCh* s_name = xercesc::XMLString::transcode("name");
      static XMLCh* s_maplist = xercesc::XMLString::transcode("mapList");
      static XMLCh* s_maptype = xercesc::XMLString::transcode("mapType");
      static XMLCh* s_map = xercesc::XMLString::transcode("map");
      static XMLCh* s_element = xercesc::XMLString::transcode("element");
      static XMLCh* s_range = xercesc::XMLString::transcode("range");
      static XMLCh* s_centroid = xercesc::XMLString::transcode("centroid");
      static XMLCh* s_rt = xercesc::XMLString::transcode("rt");
      static XMLCh* s_mz = xercesc::XMLString::transcode("mz");
      static XMLCh* s_it = xercesc::XMLString::transcode("it");
      static XMLCh* s_id = xercesc::XMLString::transcode("id");
      static XMLCh* s_rtmin = xercesc::XMLString::transcode("rtMin");
      static XMLCh* s_rtmax = xercesc::XMLString::transcode("rtMax");
      static XMLCh* s_mzmin = xercesc::XMLString::transcode("mzMin");
      static XMLCh* s_mzmax = xercesc::XMLString::transcode("mzMax");
      static XMLCh* s_itmin = xercesc::XMLString::transcode("itMin");
      static XMLCh* s_itmax = xercesc::XMLString::transcode("itMax");
      
      String tmp_str;
      if (equal(qname,s_maplist))
      {
        tmp_str = attributeAsString_(attributes,s_count);
        if (tmp_str != "")
        {
          UInt count = asUInt_(tmp_str);
          consensus_map_->getMapVector().resize(count);
          consensus_map_->getFilenames().resize(count);
        }      	
      }
      else if (equal(qname,s_maptype))
    	{
				tmp_str = attributeAsString_(attributes,s_name);
				if (tmp_str != "")
				{
					if (tmp_str == "feature_map")
					{
						feature_map_flag_ = true;
						consensus_map_flag_ = false;
					}
					else if (tmp_str == "consensus_map")
					{
						consensus_map_flag_ = true;
						feature_map_flag_ = false;
					}		
				}
			}
      else if (equal(qname,s_map))
    	{
        tmp_str = attributeAsString_(attributes,s_id);
        if (tmp_str != "")
        {
          UInt id = asUInt_(tmp_str);
          tmp_str = attributeAsString_(attributes,s_name);

          // load FeatureMapXML
          if (feature_map_flag_ && load_elements_maps_)
          {
            loadFile_(tmp_str,id,act_cons_element_);
          }
          // load MzData
          else
          {
            if (consensus_map_flag_ && load_elements_maps_)
            {
              loadFile_(tmp_str,id,act_cons_element_);
            }
            else if (load_elements_maps_)
            {
              loadFile_(tmp_str,id,act_cons_element_);
            }
          }
          consensus_map_->getFilenames()[id]=tmp_str;
        }
    	}
      else if (equal(qname,s_consensuselement))
    	{
        act_cons_element_ = ConsensusElementType();
        consensus_element_range_ = true;
    	}
      else if (equal(qname,s_centroid))
    	{
          tmp_str = attributeAsString_(attributes,s_rt);
          if (tmp_str != "")
          {
            pos_[Peak2D::RT] = asDouble_(tmp_str);
          }

          tmp_str = attributeAsString_(attributes,s_mz);
          if (tmp_str != "")
          {
            pos_[Peak2D::MZ] = asDouble_(tmp_str);
          }

          tmp_str = attributeAsString_(attributes,s_it);
          if (tmp_str != "")
          {
            it_ = asDouble_(tmp_str);
          }

    	}
      else if (equal(qname,s_range))
    	{
        if (consensus_element_range_)
        {
          tmp_str = attributeAsString_(attributes, s_rtmin);
          if (tmp_str != "")
          {
            pos_range_.setMinX(asDouble_(tmp_str));

            tmp_str = attributeAsString_(attributes, s_rtmax);
            if (tmp_str != "")
            {
              pos_range_.setMaxX(asDouble_(tmp_str));

              tmp_str = attributeAsString_(attributes, s_mzmin);
              if (tmp_str != "")
              {
                pos_range_.setMinY(asDouble_(tmp_str));

                tmp_str = attributeAsString_(attributes, s_mzmax);
                if (tmp_str != "")
                {
                  pos_range_.setMaxY(asDouble_(tmp_str));

                  tmp_str = attributeAsString_(attributes, s_itmin);
                  if (tmp_str != "")
                  {
                    it_range_.setMin(asDouble_(tmp_str));

                    tmp_str = attributeAsString_(attributes, s_itmax);
                    if (tmp_str != "")
                    {
                      it_range_.setMax(asDouble_(tmp_str));

                      consensus_element_range_ = false;
                    }
                  }
                }
              }
            }
          }
        }
    	}
      else if (equal(qname,s_element))
    	{
        if (load_elements_maps_)
        {
          IndexTuple< ElementContainerType > act_index_tuple;
          tmp_str = attributeAsString_(attributes, s_map);
          if (tmp_str != "")
          {
            UInt map_index = asUInt_(tmp_str);
            tmp_str = attributeAsString_(attributes, s_id);

            if (tmp_str != "")
            {
              UInt element_index = asUInt_(tmp_str);

              act_index_tuple.setMapIndex(map_index);
              act_index_tuple.setElementIndex(element_index);

              tmp_str = attributeAsString_(attributes, s_rt);
              PositionType pos;
              pos[0] = asDouble_(tmp_str);
              tmp_str = attributeAsString_(attributes, s_mz);
              pos[1] = asDouble_(tmp_str);

              act_index_tuple.setTransformedPosition(pos);
              act_index_tuple.setElement((*((consensus_map_->getMapVector())[map_index]))[element_index]);
              act_cons_element_.insert(act_index_tuple);
            }
          }
        }
        act_cons_element_.getPosition() = pos_;
        act_cons_element_.getPositionRange() = pos_range_;
        act_cons_element_.setIntensity(it_);
        act_cons_element_.getIntensityRange() = it_range_;
    	}
    }



    template < typename  AlignmentT  >
    void ConsensusXMLHandler<AlignmentT>::writeTo(std::ostream& os)
    {
      os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
      << "<consensusXML xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/ConsensusXML_1_0.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";

      const std::vector< ElementContainerType* >& map_vector = calignment_->getElementMapVector();
      const std::vector< String >& name_vector = calignment_->getFileNames();
      os << "\t<mapList count=\"" << map_vector.size() << "\">\n";
      os << "\t<mapType name=\"" << calignment_->getMapType() << "\"/>\n";

      // write aligned maps (mapList)
      UInt n = map_vector.size();
      for (UInt s=0; s < n; s++)
      {
        os << "\t\t<map id=\"" << s << "\" name =\"" << name_vector[s] << "\" count=\"" << map_vector[s]->size() << " \"/>\n";
      } // end for ( mapList)
      os << "\t</mapList>\n";

      os << "\t<alignmentMethod name=\"StarAlignmemt\">\n";
      os << "\t\t<matchingAlgorithm name=\"" << calignment_->getParameters().getValue("matching_algorithm:type") << "\"/>\n";
      os << "\t\t<consensusAlgorithm name=\"delaunay\"/>\n";
      os << "\t</alignmentMethod>\n";

      UInt ref_index = calignment_->getReferenceMapIndex();
      os << "\t<alignmentNewickTree> " << calignment_->getAlignmentTree() << "</alignmentNewickTree>\n";

      os << "\t<transformationList>\n";
      os << "\t\t<transformation id=\"0\" name =\"IdentityTransformation\"/>\n";
      const std::vector< Grid >& transformation_vector = calignment_->getTransformationVector();
      n = transformation_vector.size();
      UInt j=0;
      for (UInt s=0; s < n; ++s)
      {
        if (s != ref_index)
        {
          os << "\t\t<transformation id=\"" << (j+1) << "\" name=\"AffineTransformation\">\n";
          os << "\t";
          writeCellList_(os,transformation_vector[s]);
          os << "\t\t</transformation>\n";
          ++j;
        }
      }
      os << "\t</transformationList>\n";

      const std::vector < ConsensusElementType >& final_consensus_map = calignment_->getFinalConsensusMap();
      n=final_consensus_map.size();
      os << "\t<consensusElementList>\n";
      for (UInt i = 0; i < n; ++i)
      {
        os << "\t\t<consensusElement id=\""<< i << "\">\n";
        os << "\t\t\t<centroid rt=\"" << final_consensus_map[i].getRT()
        << "\" mz=\"" << final_consensus_map[i].getMZ()
        << "\" it=\"" << final_consensus_map[i].getIntensity() <<"\"/>\n";
        os << "\t\t\t<range rtMin=\"" << final_consensus_map[i].getPositionRange().min()[0]
        << "\" rtMax=\"" << final_consensus_map[i].getPositionRange().max()[0]
        << "\" mzMin=\"" << final_consensus_map[i].getPositionRange().min()[1]
        << "\" mzMax=\"" << final_consensus_map[i].getPositionRange().max()[1]
        << "\" itMin=\"" << final_consensus_map[i].getIntensityRange().min()
        << "\" itMax=\"" << final_consensus_map[i].getIntensityRange().max() <<"\"/>\n";

        os << "\t\t\t<groupedElementList>\n";
        for (typename ConsensusElementType::Group::const_iterator it = final_consensus_map[i].begin(); it != final_consensus_map[i].end(); ++it)
        {
          os  << "\t\t\t\t<element id=\"" << it->getElementIndex()
          << "\" map=\"" << it->getMapIndex()
          << "\" rt=\"" << it->getTransformedPosition()[0]
          << "\" mz=\"" << it->getTransformedPosition()[1]
          << "\" it=\"" << it->getElement().getIntensity() << "\"/>\n";
        }
        os << "\t\t\t</groupedElementList>\n";
        os << "\t\t</consensusElement>\n";
      }
      os << "\t</consensusElementList>\n";
      os << "</consensusXML>"<< std::endl;
    }

    /// Load the peaks
    template <>
    template <>
    void ConsensusXMLHandler< StarAlignment< ConsensusFeature< FeatureMap< > > > >::loadFile_< ConsensusFeature< FeatureMap< > > >(const String& file_name, UInt id, const ConsensusFeature< FeatureMap< > >& /* c */ ) throw (Exception::FileNotFound, Exception::ParseError);

    // load MzData
    template <>
    template <>
    void ConsensusXMLHandler< StarAlignment< ConsensusPeak< DPeakArray<Peak2D> > > >::loadFile_< ConsensusPeak< DPeakArray<Peak2D> > >( const String& file_name, UInt id, const ConsensusPeak< DPeakArray<Peak2D> >& /* c */) throw (Exception::FileNotFound, Exception::ParseError);

    // load consensusXML
    template <>
    template <>
    void ConsensusXMLHandler< StarAlignment< ConsensusFeature< ConsensusMap< ConsensusFeature< FeatureMap< > > > > > >::loadFile_<ConsensusFeature< ConsensusMap< ConsensusFeature< FeatureMap< > > > > >(const String& file_name, UInt id, const ConsensusFeature< ConsensusMap< ConsensusFeature< FeatureMap< > > > >& /* c */) throw (Exception::FileNotFound, Exception::ParseError);

  } // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_CONSENSUSXMLHANDLER_H
