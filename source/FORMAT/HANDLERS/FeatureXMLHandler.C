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
// $Maintainer: Ole Schulz-Trieglaff $
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
			static const XMLCh* s_value = xercesc::XMLString::transcode("value");
			static const XMLCh* s_type = xercesc::XMLString::transcode("type");
			static const XMLCh* s_int = xercesc::XMLString::transcode("int");
			static const XMLCh* s_float = xercesc::XMLString::transcode("float");
			static const XMLCh* s_string = xercesc::XMLString::transcode("string");
			static const XMLCh* s_featuremap = xercesc::XMLString::transcode("featureMap");
			static const XMLCh* s_id = xercesc::XMLString::transcode("id");

			String tag = sm_.convert(qname);
			open_tags_.push_back(tag);

			// collect Experimental Settings
			if (in_description_)
			{
				exp_sett_ << '<' << sm_.convert(qname);
				UInt n=attributes.getLength();
				for (UInt i=0; i<n; ++i)
				{
					exp_sett_ << ' ' << sm_.convert(attributes.getQName(i)) << "=\""	<< sm_.convert(attributes.getValue(i)) << '\"';
				}
				exp_sett_ << '>';
				return;
			}

			if (tag=="feature")
			{
				// create new feature at apropriate level
				updateCurrentFeature_(true);
			}
			else if (tag=="subordinate")
			{ // this is not safe towards malformed xml!
				++subordinate_feature_level_;
			}
			else if (tag=="description")
			{
				exp_sett_ << "<description>";
				in_description_ = true;
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
				const XMLCh* value = attributes.getValue(s_value);
				const XMLCh* type = attributes.getValue(s_type);
				String name = sm_.convert(attributes.getValue(s_name));

				if(*type==*s_int)
				{
					current_feature_->setMetaValue(name, xercesc::XMLString::parseInt(value));
				}
				else if (*type==*s_float)
				{
					current_feature_->setMetaValue(name, atof(sm_.convert(value)) );
				}
				else if (*type==*s_string)
				{
					current_feature_->setMetaValue(name, (String)sm_.convert(value));
				}
				else
				{
					throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Invalid userParam type '") + sm_.convert(type) + "'" );
				}
			}
			else if (equal_(qname,s_featuremap))
			{
				//check file version against schema version
				String file_version="1.0"; // default schema is 1.0
				optionalAttributeAsString_(file_version,attributes,"version");
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
		}

		void FeatureXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
		{
			static const XMLCh* s_description = xercesc::XMLString::transcode("description");
			static const XMLCh* s_feature = xercesc::XMLString::transcode("feature");
			static const XMLCh* s_model = xercesc::XMLString::transcode("model");
			static const XMLCh* s_hullpoint = xercesc::XMLString::transcode("hullpoint");
			static const XMLCh* s_convexhull = xercesc::XMLString::transcode("convexhull");
			static const XMLCh* s_subordinate = xercesc::XMLString::transcode("subordinate");
			// std::cout << "end tag: '" << sm_.convert(qname) <<"' - '"<< sm_.convert(s_description) << "'" << std::endl;

			open_tags_.pop_back();

			// collect Experimental Settings
			if (in_description_)
			{
				exp_sett_ << "</" << sm_.convert(qname) << ">\n";
				if (!equal_(qname,s_description)) return;
			}

			if (equal_(qname,s_description))
			{
				in_description_ = false;
				// call MzDataExpSett parser
				xercesc::XMLPlatformUtils::Initialize();
				xercesc::SAX2XMLReader* parser = xercesc::XMLReaderFactory::createXMLReader();
				parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces,false);
				parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpacePrefixes,false);

				MzDataExpSettHandler handler( *((ExperimentalSettings*)map_),file_);
				handler.resetErrors();
				parser->setContentHandler(&handler);
				parser->setErrorHandler(&handler);

				String tmp(exp_sett_.str().c_str());

				xercesc::MemBufInputSource source((const XMLByte*)(tmp.c_str()), tmp.size(), "dummy");
				parser->parse(source);
				delete(parser);
			}
			else if (equal_(qname,s_feature))
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
			// std::cout << "characters: "	 << sm_.convert(chars) << std::endl;
			// collect Experimental Settings
			if (in_description_)
			{
				exp_sett_ << sm_.convert(chars);
				return;
			}

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
			os.precision(written_digits_doublereal);

			os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
				 << "<featureMap version=\"" << version_ << "\"";
			if (cmap_->getIdentifier()!="")
			{
				os << " id=\"" << cmap_->getIdentifier() << "\"";
			}
			os << " xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/FeatureXML_1_3.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";

			// delegate control to ExperimentalSettings handler
			Internal::MzDataExpSettHandler handler(*((const ExperimentalSettings*)cmap_),"");
			handler.writeTo(os);

			os << "\t<featureList count=\"" << cmap_->size() << "\">\n";

			// write features with their corresponding attributes
			for (UInt s=0; s<cmap_->size(); s++)
			{
				const Feature& feat = (*cmap_)[s];
				writeFeature_(os, feat, "f_", s, 0);

			} // end for ( features )

			os << "\t</featureList>\n</featureMap>\n";
		}

		void FeatureXMLHandler::setOptions(const PeakFileOptions& options)
		{
			options_ = options;
		}

		void FeatureXMLHandler::writeFeature_(ostream& os, const Feature& feat, const String& identifier_prefix, const UInt identifier, const UInt& intendation_level)
		{
			String intend = String(intendation_level,'\t');

			os << intend << "\t\t<feature id=\"" << identifier_prefix << identifier << "\">\n";
			for (UInt i=0; i<2;i++)
			{
				os << intend <<	"\t\t\t<position dim=\"" << i << "\">" << feat.getPosition()[i] << "</position>\n";
			}
			os << intend << "\t\t\t<intensity>" << feat.getIntensity() << "</intensity>\n";
			for (UInt i=0; i<2;i++)
			{
				os << intend << "\t\t\t<quality dim=\"" << i << "\">" << feat.getQuality(i) << "</quality>\n";
			}
			os << intend << "\t\t\t<overallquality>" << feat.getOverallQuality() << "</overallquality>\n";
			os << intend << "\t\t\t<charge>" << feat.getCharge() << "</charge>\n";

			// write model description
			ModelDescription<2> desc = feat.getModelDescription();
			if (!desc.getName().empty() || !desc.getParam().empty())
			{
				os << intend << "\t\t\t<model name=\"" << desc.getName() << "\">\n";
				Param modelp = desc.getParam();
				Param::ParamIterator piter = modelp.begin();
				while (piter != modelp.end())
				{
					os << intend << "\t\t\t\t<param name=\"" << piter.getName() << "\" value=\"" << piter->value << "\"/>\n";
					piter++;
				}
				os << intend << "\t\t\t</model>\n";
			}

			// write convex hull
			vector<ConvexHull2D> hulls = feat.getConvexHulls();
			vector<ConvexHull2D>::iterator citer = hulls.begin();

			UInt hulls_count = hulls.size();

			for (UInt i=0;i<hulls_count; i++)
			{
				os << intend << "\t\t\t<convexhull nr=\"" << i << "\">\n";

				ConvexHull2D current_hull = hulls[i];
				UInt hull_size	= current_hull.getPoints().size();

				for (UInt j=0;j<hull_size;j++)
				{
					os << intend << "\t\t\t\t<hullpoint>\n";

					DPosition<2> pos = current_hull.getPoints()[j];
					UInt pos_size = pos.size();
					for (UInt k=0; k<pos_size; k++)
					{
						os << intend << "\t\t\t\t\t<hposition dim=\"" << k << "\">" << pos[k] << "</hposition>\n";
					}

					os << intend << "\t\t\t\t</hullpoint>\n";
				} // end for (..hull_size..)

				os << intend << "\t\t\t</convexhull>\n";
			} // end	for ( ... hull_count..)

			if (!feat.getSubordinates().empty())
			{
				os << intend << "\t\t\t<subordinate>\n";
				UInt identifier_subordinate = 0;
				for (std::size_t i=0;i<feat.getSubordinates().size();++i)
				{
					writeFeature_(os, feat.getSubordinates()[i], identifier_prefix+identifier+"_", identifier_subordinate, intendation_level+2);
					++identifier_subordinate;
				}
				os << intend << "\t\t\t</subordinate>\n";
			}

			writeUserParam_("userParam", os, feat, 3+intendation_level);

			os << intend << "\t\t</feature>\n";
		}

		void FeatureXMLHandler::updateCurrentFeature_(const bool create)
		{
			if (subordinate_feature_level_==0)
			{
				if (create)
				{
					map_->push_back(Feature());
					current_feature_ = &map_->back();
				}
				else
				{
					if (map_->empty())
					{
						current_feature_ = 0;
					}
					else
					{
						current_feature_ = &map_->back();
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
					return;
				}
				f1 = &f1->getSubordinates().back();
			}
			if (create)
			{
				f1->getSubordinates().push_back(Feature());
				current_feature_ = &f1->getSubordinates().back();
				return;
			}
			else
			{
				if (f1->getSubordinates().empty())
				{
					current_feature_ = 0;
					return;
				}
				else
				{
					current_feature_ = &f1->getSubordinates().back();
					return;
				}
			}
		}

	} //namespace Internal
} // namespace OpenMS
