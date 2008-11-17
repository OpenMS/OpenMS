// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//									 OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//	Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
//
//	This library is free software; you can redistribute it and/or
//	modify it under the terms of the GNU Lesser General Public
//	License as published by the Free Software Foundation; either
//	version 2.1 of the License, or (at your option) any later version.
//
//	This library is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
//	Lesser General Public License for more details.
//
//	You should have received a copy of the GNU Lesser General Public
//	License along with this library; if not, write to the Free Software
//	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA	02111-1307	USA
//
// --------------------------------------------------------------------------
// $Maintainer: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/FeatureXMLHandler.h>

using namespace std;

namespace OpenMS
{
	namespace Internal
	{
		void FeatureXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
		{
			static const XMLCh* s_dim = xercesc::XMLString::transcode("dim");
			static const XMLCh* s_name = xercesc::XMLString::transcode("name");
			static const XMLCh* s_version = xercesc::XMLString::transcode("version");
			static const XMLCh* s_value = xercesc::XMLString::transcode("value");
			static const XMLCh* s_type = xercesc::XMLString::transcode("type");
			static const XMLCh* s_completion_time = xercesc::XMLString::transcode("completion_time");
			static const XMLCh* s_id = xercesc::XMLString::transcode("id");

			String tag = sm_.convert(qname);
			String parent_tag;
			if (open_tags_.size()!=0) parent_tag = open_tags_.back();
			open_tags_.push_back(tag);
			
			//for downward compatibility, all tags in the old description must be ignored
			if (in_description_) return;
			
			if (tag=="description")
			{
				in_description_ = true;
			}
			else if (tag=="feature")
			{
				// create new feature at apropriate level
				updateCurrentFeature_(true);
			}
			else if (tag=="subordinate")
			{ // this is not safe towards malformed xml!
				++subordinate_feature_level_;
			}
			else if (tag=="featureList")
			{
				if (options_.getMetadataOnly()) throw EndParsingSoftly(__FILE__,__LINE__,__PRETTY_FUNCTION__);
			}
			else if (tag=="quality" || tag=="hposition" || tag=="position")
			{
				dim_ = attributeAsInt_(attributes,s_dim);
			}
			else if (tag=="convexhull")
			{
				current_chull_ = ConvexHull2D();
			}
			else if (tag=="hullpoint")
			{
				hull_position_ = DPosition<2>::zero;
			}
			else if (tag=="model")
			{
				model_desc_ = new ModelDescription<2>();
				param_.clear();
				model_desc_->setName(attributeAsString_(attributes,s_name));
			}
			else if (tag=="param")
			{
				String name = attributeAsString_(attributes,s_name);
				String value = attributeAsString_(attributes,s_value);
				if (name != "" && value != "") param_.setValue(name, value);
			}
			else if (tag == "userParam")
			{
				String name = attributeAsString_(attributes,s_name);
				String type = attributeAsString_(attributes,s_type);

				if (parent_tag=="feature")
				{
					if(type=="int")
					{
						last_meta_->setMetaValue(name, attributeAsInt_(attributes,s_value));
					}
					else if (type=="float")
					{
						last_meta_->setMetaValue(name, attributeAsDouble_(attributes,s_value));
					}
					else if (type=="string")
					{
						last_meta_->setMetaValue(name, (String)attributeAsString_(attributes,s_value));
					}
					else
					{
						throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Invalid userParam type '") + type + "'" );
					}
				}
			}
			else if (tag=="featureMap")
			{
				//check file version against schema version
				String file_version="1.0"; // default schema is 1.0
				optionalAttributeAsString_(file_version,attributes,s_version);
				if (file_version.toDouble()>version_.toDouble())
				{
					warning(String("The XML file (") + file_version +") is newer than the parser (" + version_ + "). This might lead to undefinded program behaviour.");
				}
				//handle file id
				String id;
				if (optionalAttributeAsString_(id, attributes, s_id))
				{
					map_->setIdentifier(id);
				}
			}
			else if (tag=="dataProcessing")
			{
				DataProcessing tmp;
				tmp.setCompletionTime(asDateTime_(attributeAsString_(attributes, s_completion_time)));
				map_->getDataProcessing().push_back(tmp);
				last_meta_ = &(map_->getDataProcessing().back());
			}
			else if (tag=="software" && parent_tag=="dataProcessing")
			{
				map_->getDataProcessing().back().getSoftware().setName(attributeAsString_(attributes, s_name));
				map_->getDataProcessing().back().getSoftware().setVersion(attributeAsString_(attributes, s_version));
			}
			else if (tag=="processingAction" && parent_tag=="dataProcessing")
			{
				String name = attributeAsString_(attributes, s_name);
				for (UInt i=0; i< DataProcessing::SIZE_OF_PROCESSINGACTION; ++i)
				{
					if (name == DataProcessing::NamesOfProcessingAction[i])
					{
						map_->getDataProcessing().back().getProcessingActions().insert((DataProcessing::ProcessingAction)i);
					}
				}
			}
		}

		void FeatureXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
		{
			static const XMLCh* s_feature = xercesc::XMLString::transcode("feature");
			static const XMLCh* s_model = xercesc::XMLString::transcode("model");
			static const XMLCh* s_description = xercesc::XMLString::transcode("description");
			static const XMLCh* s_hullpoint = xercesc::XMLString::transcode("hullpoint");
			static const XMLCh* s_convexhull = xercesc::XMLString::transcode("convexhull");
			static const XMLCh* s_subordinate = xercesc::XMLString::transcode("subordinate");

			open_tags_.pop_back();
			
			//for downward compatibility, all tags in the old description must be ignored
			if (equal_(qname,s_description))
			{
				in_description_ = false;
			}
			if (in_description_) return;

			if (equal_(qname,s_feature))
			{
				if ((!options_.hasRTRange() || options_.getRTRange().encloses(current_feature_->getRT()))
						&&	(!options_.hasMZRange() || options_.getMZRange().encloses(current_feature_->getMZ()))
						&&	(!options_.hasIntensityRange() || options_.getIntensityRange().encloses(current_feature_->getIntensity())))
				{
				}
				else
				{
					// this feature does not pass the restrictions --> remove it
					if (subordinate_feature_level_==0)
					{
						map_->pop_back();
					}
					else
					{
						Feature* f1;
						if (!map_->empty())
						{
							f1 = &(map_->back());
						}
						else
						{
							throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Feature", "<Feature> with unexpected location." );
						}

						for (Int level = 1; level<subordinate_feature_level_;++level)
						{
							f1 = &(f1->getSubordinates().back());
						}
						// delete the offending feature
						f1->getSubordinates().pop_back();
					}
				}
				updateCurrentFeature_(false);
			}
			else if (equal_(qname,s_model))
			{
				model_desc_->setParam(param_);
				current_feature_->setModelDescription(*model_desc_);
				delete model_desc_;
			}
			else if (equal_(qname,s_hullpoint))
			{
				current_chull_.addPoint(hull_position_);
			}
			else if (equal_(qname,s_convexhull))
			{
				current_feature_->getConvexHulls().push_back(current_chull_);
			}
			else if (equal_(qname,s_subordinate))
			{ // this is not safe towards malformed xml!
				--subordinate_feature_level_;
				if (subordinate_feature_level_ < 0) throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "subordinate", "Too many closing tags for </subordinate>." );
				// reset current_feature
				updateCurrentFeature_(false);
			}
		}

		void FeatureXMLHandler::characters(const XMLCh* const chars, unsigned int /*length*/)
		{
			//for downward compatibility, all tags in the old description must be ignored
			if (in_description_) return;
			
			String& current_tag = open_tags_.back();
			if (current_tag == "intensity")
			{
				current_feature_->setIntensity(asDouble_(sm_.convert(chars)));
			}
			else if (current_tag == "position")
			{
				current_feature_->getPosition()[dim_] = asDouble_(sm_.convert(chars));
			}
			else if (current_tag == "quality")
			{
				current_feature_->setQuality(dim_, asDouble_(sm_.convert(chars)));
			}
			else if (current_tag == "overallquality")
			{
				current_feature_->setOverallQuality(asDouble_(sm_.convert(chars)));
			}
			else if (current_tag == "charge")
			{
				current_feature_->setCharge(asInt_(chars));
			}
			else if (current_tag == "hposition")
			{
				hull_position_[dim_] = asDouble_(sm_.convert(chars));
			}
		}


		void FeatureXMLHandler::writeTo(ostream& os)
		{
			os.precision(writtenDigits<DoubleReal>());

			os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
				 << "<featureMap version=\"" << version_ << "\"";
			if (cmap_->getIdentifier()!="")
			{
				os << " id=\"" << cmap_->getIdentifier() << "\"";
			}
			os << " xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/FeatureXML_1_3.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";
			
			//write data processing
			for (UInt i=0; i< cmap_->getDataProcessing().size(); ++i)
			{
				const DataProcessing& processing = cmap_->getDataProcessing()[i];
				os << "\t<dataProcessing completion_time=\"" << processing.getCompletionTime().get() << "\">\n";
				os << "\t\t<software name=\"" << processing.getSoftware().getName() << "\" version=\"" << processing.getSoftware().getVersion() << "\" />\n";
				for (std::set<DataProcessing::ProcessingAction>::const_iterator it = processing.getProcessingActions().begin(); it!=processing.getProcessingActions().end(); ++it)
				{
					os << "\t\t<processingAction name=\"" << DataProcessing::NamesOfProcessingAction[*it] << "\" />\n";
				}
				writeUserParam_ ("userParam", os, processing, 2);
				os << "\t</dataProcessing>\n";
			}
			
			// write features with their corresponding attributes
			os << "\t<featureList count=\"" << cmap_->size() << "\">\n";
			for (UInt s=0; s<cmap_->size(); s++)
			{
				const Feature& feat = (*cmap_)[s];
				writeFeature_(os, feat, "f_", s, 0);
			}

			os << "\t</featureList>\n";
			os << "</featureMap>\n";
		}

		void FeatureXMLHandler::setOptions(const PeakFileOptions& options)
		{
			options_ = options;
		}

		void FeatureXMLHandler::writeFeature_(ostream& os, const Feature& feat, const String& identifier_prefix, const UInt identifier, const UInt& indentation_level)
		{
			String indent = String(indentation_level,'\t');

			os << indent << "\t\t<feature id=\"" << identifier_prefix << identifier << "\">\n";
			for (UInt i=0; i<2;i++)
			{
				os << indent <<	"\t\t\t<position dim=\"" << i << "\">" << precisionWrapper(feat.getPosition()[i]) << "</position>\n";
			}
			os << indent << "\t\t\t<intensity>" << precisionWrapper(feat.getIntensity()) << "</intensity>\n";
			for (UInt i=0; i<2;i++)
			{
				os << indent << "\t\t\t<quality dim=\"" << i << "\">" << precisionWrapper(feat.getQuality(i)) << "</quality>\n";
			}
			os << indent << "\t\t\t<overallquality>" << precisionWrapper(feat.getOverallQuality()) << "</overallquality>\n";
			os << indent << "\t\t\t<charge>" << feat.getCharge() << "</charge>\n";

			// write model description
			ModelDescription<2> desc = feat.getModelDescription();
			if (!desc.getName().empty() || !desc.getParam().empty())
			{
				os << indent << "\t\t\t<model name=\"" << desc.getName() << "\">\n";
				Param modelp = desc.getParam();
				Param::ParamIterator piter = modelp.begin();
				while (piter != modelp.end())
				{
					os << indent << "\t\t\t\t<param name=\"" << piter.getName() << "\" value=\"" << piter->value << "\"/>\n";
					piter++;
				}
				os << indent << "\t\t\t</model>\n";
			}

			// write convex hull
			vector<ConvexHull2D> hulls = feat.getConvexHulls();
			vector<ConvexHull2D>::iterator citer = hulls.begin();

			UInt hulls_count = hulls.size();

			for (UInt i=0;i<hulls_count; i++)
			{
				os << indent << "\t\t\t<convexhull nr=\"" << i << "\">\n";

				ConvexHull2D current_hull = hulls[i];
				UInt hull_size	= current_hull.getPoints().size();

				for (UInt j=0;j<hull_size;j++)
				{
					os << indent << "\t\t\t\t<hullpoint>\n";

					DPosition<2> pos = current_hull.getPoints()[j];
					UInt pos_size = pos.size();
					for (UInt k=0; k<pos_size; k++)
					{
						os << indent << "\t\t\t\t\t<hposition dim=\"" << k << "\">" << precisionWrapper(pos[k]) << "</hposition>\n";
					}

					os << indent << "\t\t\t\t</hullpoint>\n";
				} // end for (..hull_size..)

				os << indent << "\t\t\t</convexhull>\n";
			} // end	for ( ... hull_count..)

			if (!feat.getSubordinates().empty())
			{
				os << indent << "\t\t\t<subordinate>\n";
				UInt identifier_subordinate = 0;
				for (std::size_t i=0;i<feat.getSubordinates().size();++i)
				{
					writeFeature_(os, feat.getSubordinates()[i], identifier_prefix+identifier+"_", identifier_subordinate, indentation_level+2);
					++identifier_subordinate;
				}
				os << indent << "\t\t\t</subordinate>\n";
			}


			writeUserParam_("userParam", os, feat, 3+indentation_level);

			os << indent << "\t\t</feature>\n";
		}

		void FeatureXMLHandler::updateCurrentFeature_(bool create)
		{
			if (subordinate_feature_level_==0)
			{
				if (create)
				{
					map_->push_back(Feature());
					current_feature_ = &map_->back();
					last_meta_ =  &map_->back();
				}
				else
				{
					if (map_->empty())
					{
						current_feature_ = 0;
					last_meta_ =  0;
					}
					else
					{
						current_feature_ = &map_->back();
						last_meta_ =  &map_->back();
					}
				}
				return;
			}

			Feature* f1 = 0;
			if (map_->empty())
			{
				// do NOT throw an exception here. this is a valid case!	e.g. the
				// only one feature in a map was discarded during endElement(), thus
				// the map_ is empty() now and we cannot assign a current_feature,
				// because there is none!
				current_feature_ = 0;
				last_meta_ = 0;
				return;
			}
			else
			{
				f1 = &map_->back();
			}

			for (Int level = 1; level<subordinate_feature_level_;++level)
			{
				// if all features of the current level are discarded (due to
				// range-restrictions etc), then the current feature is the one which
				// is one level up
				if (f1->getSubordinates().empty())
				{
					current_feature_ = f1;
					last_meta_ = f1;
					return;
				}
				f1 = &f1->getSubordinates().back();
			}
			if (create)
			{
				f1->getSubordinates().push_back(Feature());
				current_feature_ = &f1->getSubordinates().back();
				last_meta_ = &f1->getSubordinates().back();
				return;
			}
			else
			{
				if (f1->getSubordinates().empty())
				{
					current_feature_ = 0;
					last_meta_ = 0;
					return;
				}
				else
				{
					current_feature_ = &f1->getSubordinates().back();
					last_meta_ = &f1->getSubordinates().back();
					return;
				}
			}
		}

	} //namespace Internal
} // namespace OpenMS
