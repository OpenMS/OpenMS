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
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/UniqueIdGenerator.h>
#include <OpenMS/FORMAT/HANDLERS/SchemaHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLSchemes.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/ConsensusPeak.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/SYSTEM/File.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/util/TransService.hpp>

#include <fstream>

// STL includes
#include <iostream>


namespace OpenMS
{
  namespace Internal
  {
    /**
      
      @brief XML Handler for a consensusXML.
     
      This class can be used to load the content of a consensusXML file into a consensusMap, as well as to
      save the content of a StarAlignment object into an XML file. 
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
            act_cons_element_(),
            calignment_(0),
            consensus_element_range_(false),
            feature_map_flag_(true)
        {
          fillMaps_(Schemes::ConsensusXML[schema_]); // fill maps with current schema
          setMaps_(TAGMAP, ATTMAP);
        }

        ///
        ConsensusXMLHandler(const StarAlignment< ConsensusElementType >& alignment, const String& filename)
            : SchemaHandler(TAG_NUM,MAP_NUM,filename),
            consensus_map_(0),
            act_cons_element_(),
            calignment_(&alignment),
            consensus_element_range_(false),
            feature_map_flag_(true)
        {
          fillMaps_(Schemes::ConsensusXML[schema_]); // fill maps with current schema
          setMaps_(TAGMAP, ATTMAP);
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
        ConsensusMap<ConsensusElementType>* consensus_map_;
        ConsensusElementType act_cons_element_;
        // StartAlignment pointer for writing
        const StarAlignment< ConsensusElementType >* calignment_;
        bool consensus_element_range_;
        bool feature_map_flag_;
        bool consensus_map_flag_;
        typename ConsensusElementType::PositionType pos_;
        typename ConsensusElementType::IntensityType it_;
        typename ConsensusElementType::PositionBoundingBoxType pos_range_;
        typename ConsensusElementType::IntensityBoundingBoxType it_range_;

        /** @brief indices for tags used by StarAlignment

          Used to access is_parser_in_tag_.
          If you add tags, also add them to XMLSchemes.h.
          Add no elements to the enum after TAG_NUM.
        */
        enum Tags { TAGNULL, CONSENSUSXML, MAPLIST, MAPTYPE, MAP, ALIGNMENT, ALIGNMENTMETHOD,
                    MATCHINGALGORITHM, CONSENSUSALGORITHM, ALIGNMENTNEWICKTREE, TRANSFORMATIONLIST, TRANSFORMATION,
                    CELL, RANGE, PARAMETERS, CONSENSUSELEMENTLIST, CONSENSUSELEMENT, CENTROID, GROUPEDELEMENTLIST, ELEMENT, TAG_NUM};

        /** @brief indices for attributes used by StarAlignment

          If you add attributes, also add them to XMLSchemes.h.
          Add no elements to the enum after ATT_NUM.
        */
        enum Attributes { ATTNULL, COUNT, NAME, ID, RT_ATT, MZ_ATT, IT, RTMIN, RTMAX, MZMIN, MZMAX, ITMIN, ITMAX, MAP_ATT, ATT_NUM};

        /** @brief indices for enum2str-maps used by DFeatureMapFile

          Used to access enum2str_().
          If you add maps, also add them to XMLSchemes.h.
          Add no elements to the enum after MAP_NUM.
          Each map corresponds to a string in XMLSchemes.h.
        */
        enum MapTypes { TAGMAP, ATTMAP, MAP_NUM };

        /// This function fills the members of a picked peak of type OutputPeakType.
        template <typename ConsensusElementT >
        void loadFile_(const String& /* file_name */, UnsignedInt /* id */, const ConsensusElementT& /* c */) throw (Exception::FileNotFound, Exception::ParseError)
        {}

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

    template < typename  AlignmentT >
    void ConsensusXMLHandler<AlignmentT>::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
      int tag = leaveTag(qname);

      // Do something depending on the tag
      switch(tag)
      {
          case CONSENSUSELEMENT:
          consensus_map_->push_back(act_cons_element_);
          break;
      }
    }
    
    template < typename AlignmentT >
    void ConsensusXMLHandler<AlignmentT>::characters(const XMLCh* const /*chars*/, const unsigned int /*length*/)
  {}

    template < typename  AlignmentT >
    void ConsensusXMLHandler<AlignmentT>::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
    {
      int tag = enterTag(qname, attributes);

      String tmp_str;
      // Do something depending on the tag
      switch(tag)
      {
          case MAPLIST:
          tmp_str = getAttributeAsString(COUNT);
          if (tmp_str != "")
          {
            UnsignedInt count = asUnsignedInt_(tmp_str);
            consensus_map_->getMapVector().resize(count);
            consensus_map_->getFilenames().resize(count);
          }
          break;
          case MAPTYPE:
          tmp_str = getAttributeAsString(NAME);
          if (tmp_str != "")
          {
            if (getAttributeAsString(ID) == "feature_map")
            {
              feature_map_flag_ = true;
              consensus_map_flag_ = false;
            }
            else
              if (getAttributeAsString(ID) == "consensus_map")
              {
                consensus_map_flag_ = true;
                feature_map_flag_ = false;
              }

          }
          break;
          case MAP:
          tmp_str = getAttributeAsString(ID);
          if (tmp_str != "")
          {
            UnsignedInt id = asUnsignedInt_(tmp_str);
            tmp_str = getAttributeAsString(NAME);
            if (tmp_str != "")
            {
              String act_filename = tmp_str;
              // load FeatureMapXML
              if (feature_map_flag_)
              {
                loadFile_(act_filename,id,act_cons_element_);
              }
              // load MzData
              else
              {
                if (consensus_map_flag_)
                {
                  loadFile_(act_filename,id,act_cons_element_);
                }
                else
                {
                  loadFile_(act_filename,id,act_cons_element_);
                }
              }
              consensus_map_->getFilenames()[id]=act_filename;
            }
          }
          break;
          case CONSENSUSELEMENT:
          act_cons_element_ = ConsensusElementType();
          consensus_element_range_ = true;
          break;
          case CENTROID:
          tmp_str = getAttributeAsString(RT_ATT);
          if (tmp_str != "")
          {
            pos_[RT] = asDouble_(tmp_str);
          }

          tmp_str = getAttributeAsString(MZ_ATT);
          if (tmp_str != "")
          {
            pos_[MZ] = asDouble_(tmp_str);
          }

          tmp_str = getAttributeAsString(IT);
          if (tmp_str != "")
          {
            it_ = asDouble_(tmp_str);
          }
          break;
          case RANGE:
          if (consensus_element_range_)
          {
            tmp_str = getAttributeAsString(RTMIN);
            if (tmp_str != "")
            {
              pos_range_.setMinX(asDouble_(tmp_str));

              tmp_str = getAttributeAsString(RTMAX);
              if (tmp_str != "")
              {
                pos_range_.setMaxX(asDouble_(tmp_str));

                tmp_str = getAttributeAsString(MZMIN);
                if (tmp_str != "")
                {
                  pos_range_.setMinY(asDouble_(tmp_str));

                  tmp_str = getAttributeAsString(MZMAX);
                  if (tmp_str != "")
                  {
                    pos_range_.setMaxY(asDouble_(tmp_str));

                    tmp_str = getAttributeAsString(ITMIN);
                    if (tmp_str != "")
                    {
                      it_range_.setMin(asDouble_(tmp_str));

                      tmp_str = getAttributeAsString(ITMAX);
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
          break;
          case ELEMENT:
          IndexTuple< ElementContainerType > act_index_tuple;
          tmp_str = getAttributeAsString(MAP_ATT);
          if (tmp_str != "")
          {
            UnsignedInt map_index = asUnsignedInt_(tmp_str);
            tmp_str = getAttributeAsString(ID);
            if (tmp_str != "")
            {
              UnsignedInt element_index = asUnsignedInt_(tmp_str);

              act_index_tuple.getMapIndex() = map_index;
              act_index_tuple.getElementIndex() = element_index;
              act_index_tuple.setElement(((consensus_map_->getMapVector())[map_index])[element_index]);
              act_cons_element_.insert(act_index_tuple);
              act_cons_element_.getPosition() = pos_;
              act_cons_element_.getPositionRange() = pos_range_;
              act_cons_element_.getIntensity() = it_;
              act_cons_element_.getIntensityRange() = it_range_;
            }
          }
          break;
      }
    }



    template < typename  AlignmentT  >
    void ConsensusXMLHandler<AlignmentT>::writeTo(std::ostream& os)
    {
      os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
      << "<consensusXML>\n";

      const std::vector< ElementContainerType* >& map_vector = calignment_->getElementMapVector();
      const std::vector< String >& name_vector = calignment_->getFileNames();
      os << "\t<mapList count=\"" << map_vector.size() << "\">\n";
      os << "\t<map_type name=\"" << calignment_->getMapType() << "\"/>\n";

      // write aligned maps (mapList)
      UnsignedInt n = map_vector.size();
      for (UnsignedInt s=0; s < n; s++)
      {
        os << "\t\t<map id=\"" << s << "\" name =\"" << name_vector[s] << "\" count=\"" << map_vector[s]->size() << " \"/>\n";
      } // end for ( mapList)
      os << "\t</mapList>\n";

      os << "\t<alignmentMethod name=\"StarAlignmemt\">\n";
      os << "\t\t<matching_algorithm name=\"" << calignment_->getParam().getValue("matching_algorithm") << "\"/>\n";
      os << "\t\t<consensus_algorithm name=\"" << calignment_->getParam().getValue("consensus_algorithm") << "\"/>\n";
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
