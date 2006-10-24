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
// $Maintainer: Eva Lange$
// --------------------------------------------------------------------------


#ifndef OPENMS_FORMAT_HANDLERS_CONSENSUSXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_CONSENSUSXMLHANDLER_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/StarAlignment.h>
#include <OpenMS/FORMAT/DFeatureMapFile.h>
#include <OpenMS/FORMAT/UniqueIdGenerator.h>
#include <OpenMS/FORMAT/HANDLERS/SchemaHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLSchemes.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/util/TransService.hpp>

// STL includes
#include <iostream>


namespace OpenMS
{
  namespace Internal
  {

    /**
      
      @brief XML Handler for a StarAlignment
     
      This class can be used to save the content of a
      StarAlignment object into an XML file. 
    */
    template < typename AlignmentT >
    class ConsensusXMLHandler
          : public SchemaHandler
    {
      public:

        /// Defines the coordinates of peaks / features.
        enum DimensionId
        {
          RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
          MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
      };

        /**
        @name Type definitions
        */
        //@{
        typedef typename AlignmentT::ElementContainerType ElementContainerType;
        typedef typename AlignmentT::ElementType ElementType;
        typedef typename AlignmentT::ConsensusMapType ConsensusMapType;
        typedef typename AlignmentT::ConsensusElementType ConsensusElementType;
        typedef typename ConsensusElementType::TraitsType TraitsType;

        typedef DGrid<2> GridType;
        typedef DGridCell<2> GridCellType;
        typedef typename GridCellType::MappingVector MappingVector;
        typedef DPosition<2,TraitsType> PositionType;

        /**@name Constructors and destructor */
        //@{
        ///
        ConsensusXMLHandler(ConsensusMap<ConsensusElementType>& consensus_map , const String& filename)
            : SchemaHandler(TAG_NUM,MAP_NUM,filename),
            consensus_map_(&consensus_map),
            calignment_(0),
            act_cons_element_(0),
            act_index_tuple_(0),
            consensus_element_range_(false)
        {
          fillMaps_(Schemes::ConsensusXML[schema_]); // fill maps with current schema
        }

        ///
        ConsensusXMLHandler(const StarAlignment< ConsensusElementType >& alignment, const String& filename)
            : SchemaHandler(TAG_NUM,MAP_NUM,filename),
            calignment_(&alignment),
            act_cons_element_(0),
            act_index_tuple_(0),
            consensus_element_range_(false)
        {
          fillMaps_(Schemes::ConsensusXML[schema_]); // fill maps with current schema
        }

        ///
        virtual ~ConsensusXMLHandler()
        { }
        //@}

        // Docu in base class
        virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);

        // Docu in base class
        virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);

        // Docu in base class
        virtual void characters(const XMLCh* const chars, const unsigned int length);

        ///Writes the contents to a stream
        void writeTo(std::ostream& os);

      protected:
        ConsensusElementType* act_cons_element_;
        ConsensusMap<ConsensusElementType>* consensus_map_;
        IndexTuple< ElementContainerType >* act_index_tuple_;
        // StartAlignment pointer for writing
        const StarAlignment< ConsensusElementType >* calignment_;
        bool consensus_element_range_;
        typename ConsensusElementType::PositionType pos_;
        typename ConsensusElementType::IntensityType it_;
				typename ConsensusElementType::PositionBoundingBoxType pos_range_;
				typename ConsensusElementType::IntensityBoundingBoxType it_range_;

        /** @brief indices for tags used by StarAlignment

          Used to access is_parser_in_tag_.
          If you add tags, also add them to XMLSchemes.h.
          Add no elements to the enum after TAG_NUM.
        */
        enum Tags { TAGNULL, CONSENSUSXML, MAPLIST, MAP, ALIGNMENT, ALIGNMENTMETHOD,
                    MATCHINGALGORITHM, CONSENSUSALGORITHM, ALIGNMENTNEWICKTREE, TRANSFORMATIONLIST, TRANSFORMATION,
                    CELL, RANGE, PARAMETERS, CONSENSUSELEMENTLIST, CONSENSUSELEMENT, CENTROID, GROUPEDELEMENTLIST, ELEMENT, TAG_NUM};

        /** @brief indices for enum2str-maps used by DFeatureMapFile

          Used to access enum2str_().
          If you add maps, also add them to XMLSchemes.h.
          Add no elements to the enum after MAP_NUM.
          Each map corresponds to a string in XMLSchemes.h.
        */
        enum MapTypes { TAGMAP, MAP_NUM };


        void writeCellList_(std::ostream& os, const GridType& grid)
        {
          // write features with their attributes
          for (UnsignedInt s=0; s < grid.size(); s++)
          {
            const GridCellType& cell = grid[s];

            os << "\t\t<cell nr=\"" << s << "\">" << std::endl;
            os << "\t\t\t\t<range rtMin=\"" << cell.min()[0] << "\" rtMax=\"" << cell.max()[0]
            << "\" mzMin=\"" << cell.min()[1] << "\" mzMax=\"" << cell.max()[1] << "\"/>\n";

            os << "\t\t\t\t<mappinglist>\n";
            MappingVector mappings = cell.getMappings();

            typename MappingVector::const_iterator citer = mappings.begin();
            UnsignedInt i = 0;
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

              Param map_param = (*citer)->getParam();
              Param::ConstIterator piter = map_param.begin();
              while (piter != map_param.end())
              {
                os << "\t\t\t\t\t\t<param name=\"" << piter->first << "\" value=\"" << piter->second << "\"/>\n";
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

    template < typename  ConsensusElementT >
    void ConsensusXMLHandler<ConsensusElementT>::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
      int tag = str2enum_(TAGMAP,xercesc::XMLString::transcode(qname),"closing tag");  // index of current tag
      is_parser_in_tag_[tag] = false;

      // Do something depending on the tag
      switch(tag)
      {
          case CONSENSUSELEMENT:
          consensus_map_->push_back(*act_cons_element_);
          delete act_cons_element_;
          break;     	
      }
    }


    template < typename  ConsensusElementT >
    void ConsensusXMLHandler<ConsensusElementT>::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
    {
      int tag = str2enum_(TAGMAP,xercesc::XMLString::transcode(qname),"opening tag"); // index of current tag
      is_parser_in_tag_[tag] = true;

      // Do something depending on the tag
      switch(tag)
      {
          case MAPLIST:    
          if (attributes.getIndex(xercesc::XMLString::transcode("count")) != -1)
          	{
          		UnsignedInt count = asUnsignedInt_(xercesc::XMLString::transcode(attributes.getValue(xercesc::XMLString::transcode("count"))));
            	consensus_map_->getMapVector().resize(count);
            	consensus_map_->getFilenames().resize(count);
            }
          break;
          case MAP:
          if (attributes.getIndex(xercesc::XMLString::transcode("id")) != -1)
          	{
          		UnsignedInt id = asUnsignedInt_(xercesc::XMLString::transcode(attributes.getValue(xercesc::XMLString::transcode("id"))));
          		if (attributes.getIndex(xercesc::XMLString::transcode("name")) != -1)
          			{
          				String act_filename = xercesc::XMLString::transcode(attributes.getValue(xercesc::XMLString::transcode("name")));
         					DFeatureMapFile feature_file;
          				feature_file.load(act_filename,(consensus_map_->getMapVector())[id]);
          				consensus_map_->getFilenames()[id]=act_filename;
          			}
        		}
          break;
          case CONSENSUSELEMENT:
          act_cons_element_ = new ConsensusElementType;
          consensus_element_range_ = true;
          break;
          case CENTROID:
          if (attributes.getIndex(xercesc::XMLString::transcode("rt"))!=-1)
          {
            pos_[RT] = asDouble_(xercesc::XMLString::transcode(attributes.getValue(xercesc::XMLString::transcode("rt"))));
          }

          if (attributes.getIndex(xercesc::XMLString::transcode("mz"))!=-1)
          {
            pos_[MZ] = asDouble_(xercesc::XMLString::transcode(attributes.getValue(xercesc::XMLString::transcode("mz"))));
          }

          if (attributes.getIndex(xercesc::XMLString::transcode("it")) != -1)
          {
            it_ = asDouble_(xercesc::XMLString::transcode(attributes.getValue(xercesc::XMLString::transcode("it"))));
          }
          break;
          case RANGE:
          if (consensus_element_range_)
          {
            if (attributes.getIndex(xercesc::XMLString::transcode("rtMin")) != -1)
            {
              pos_range_.setMinX(asDouble_(xercesc::XMLString::transcode(attributes.getValue(xercesc::XMLString::transcode("rtMin")))));

              if (attributes.getIndex(xercesc::XMLString::transcode("rtMax")) != -1)
              {
                pos_range_.setMaxX(asDouble_(xercesc::XMLString::transcode(attributes.getValue(xercesc::XMLString::transcode("rtMax")))));

                if (attributes.getIndex(xercesc::XMLString::transcode("mzMin")) != -1)
                {
                  pos_range_.setMinY(asDouble_(xercesc::XMLString::transcode(attributes.getValue(xercesc::XMLString::transcode("mzMin")))));

                  if (attributes.getIndex(xercesc::XMLString::transcode("mzMax")) != -1)
                  {
                    pos_range_.setMaxY(asDouble_(xercesc::XMLString::transcode(attributes.getValue(xercesc::XMLString::transcode("mzMax")))));

                    if (attributes.getIndex(xercesc::XMLString::transcode("itMin")) != -1)
                    {
                      it_range_.setMin(asDouble_(xercesc::XMLString::transcode(attributes.getValue(xercesc::XMLString::transcode("itMin")))));

                      if (attributes.getIndex(xercesc::XMLString::transcode("itMax")) != -1)
                      {
                        it_range_.setMax(asDouble_(xercesc::XMLString::transcode(attributes.getValue(xercesc::XMLString::transcode("itMax")))));

                        consensus_element_range_ = false;
                      }
                    }
                  }
                }
              }
            }
          }
          break;
          case ELEMENT:
          act_index_tuple_ = new IndexTuple< ElementContainerType >;
          if (attributes.getIndex(xercesc::XMLString::transcode("map")) != -1)
          {
            UnsignedInt map_index = asUnsignedInt_(xercesc::XMLString::transcode(attributes.getValue(xercesc::XMLString::transcode("map"))));
            if (attributes.getIndex(xercesc::XMLString::transcode("id")) != -1)
            {
              UnsignedInt element_index = asUnsignedInt_(xercesc::XMLString::transcode(attributes.getValue(xercesc::XMLString::transcode("id"))));

              act_index_tuple_->getMapIndex() = map_index;
              act_index_tuple_->getElementIndex() = element_index;
              act_index_tuple_->setElement(((consensus_map_->getMapVector())[map_index])[element_index]);
              act_cons_element_->insert(*act_index_tuple_);
              act_cons_element_->getPosition() = pos_;
              act_cons_element_->getPositionRange() = pos_range_;
              act_cons_element_->getIntensity() = it_;
              act_cons_element_->getIntensityRange() = it_range_;
            }
          }
          break;
      }
    }

    template < typename ConsensusElementT >
    void ConsensusXMLHandler<ConsensusElementT>::characters(const XMLCh* const chars, const unsigned int /*length*/)
    {
      // find the tag that the parser is in right now
      for (Size i=0; i<is_parser_in_tag_.size(); i++)
      {
        if (is_parser_in_tag_[i])
        {
          switch(i)
          {
              
          }
        }
      }
    }


    template < typename  ConsensusElementT  >
    void ConsensusXMLHandler<ConsensusElementT>::writeTo(std::ostream& os)
    {
      os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
      << "<consensusXML>\n";

      const std::vector< ElementContainerType* >& map_vector = calignment_->getElementMapVector();
      const std::vector< String >& name_vector = calignment_->getFileNames();
      os << "\t<mapList count=\"" << map_vector.size() << "\">\n";

      // write aligned maps (mapList)
      UnsignedInt n = map_vector.size();
      for (UnsignedInt s=0; s < n; s++)
      {
        os << "\t\t<map id=\"" << s << "\" name =\"" << name_vector[s] << "\" count=\"" << map_vector[s]->size() << " \"/>\n";
      } // end for ( mapList)
      os << "\t</mapList>\n";

      os << "\t<alignmentMethod name=\"StarAlignmemt\">\n";
      os << "\t\t<matchingAlgorithm name=\"" << calignment_->getParam().getValue("matchingAlgorithm") << "\"/>\n";
      os << "\t\t<consensusAlgorithm name=\"" << calignment_->getParam().getValue("consensusAlgorithm") << "\"/>\n";
      os << "\t</alignmentMethod>\n";

      UnsignedInt ref_index = calignment_->getReferenceMapIndex();
      os << "\t<alignmentNewickTree> " << calignment_->getAlignmentTree() << "</alignmentNewickTree>\n";

      os << "\t<transformationList>\n";
      os << "\t\t<transformation id=\"0\" name =\"IdentityTransformation\"/>\n";
      const std::vector< GridType >& transformation_vector = calignment_->getTransformationVector();
      n = transformation_vector.size();
      UnsignedInt j=0;
      for (UnsignedInt s=0; s < n; ++s)
      {
        if (s != ref_index)
        {

          os << "\t\t<transformation id=\"" << (j+1) << "\" name=\"AffineTransformation\"/>\n";
          os << "\t";
          writeCellList_(os,transformation_vector[s]);
          ++j;
        }
      }
      os << "\t</transformationList>\n";

      const std::vector < ConsensusElementType >& final_consensus_map = calignment_->getFinalConsensusMap();
      n=final_consensus_map.size();
      os << "\t<consensusElementList>\n";
      for (UnsignedInt i = 0; i < n; ++i)
      {
        os << "\t\t<consensusElement id=\""<< i << "\">\n";
        os << "\t\t\t<centroid rt=\"" << final_consensus_map[i].getPosition()[0]
        << "\" mz=\"" << final_consensus_map[i].getPosition()[1]
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
          << "\" rt=\"" << it->getElement().getPosition()[0]
          << "\" mz=\"" << it->getElement().getPosition()[1]
          << "\" it=\"" << it->getElement().getIntensity() << "\"/>\n";
        }
        os << "\t\t\t</groupedElementList>\n";
        os << "\t\t</consensusElement>\n";
      }
      os << "\t</consensusElementList>\n";
      os << "</consensusXML>"<< std::endl;
    }


  } // namespace Internal
} // namespace OpenMS

#endif
