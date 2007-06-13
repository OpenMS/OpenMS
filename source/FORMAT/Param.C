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
// $Maintainer: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/Param.h>

#include <iostream>
#include <fstream>

using namespace std;

#include <OpenMS/FORMAT/HANDLERS/ParamXMLHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

namespace OpenMS
{

	using namespace Exception;

	Param::Param(): values_(), descriptions_(), inheritance_steps_max(15)
	{
	}

	Param::Param(const Param& rhs): values_(rhs.values_), descriptions_(rhs.descriptions_)
	{
	}

	Param::~Param()
	{
		
	}

	Param& Param::operator = (const Param& rhs)
	{
		values_ = rhs.values_;
		descriptions_ = rhs.descriptions_;
		return *this;
	}

	bool Param::operator == (const Param& rhs) const
	{
		return (values_ == rhs.values_);
	}
	
	
	// get/set the values
	void Param::setValue(const String& key, Int value, const String& description)
	{
		values_[key] = value;
		if(!description.empty())
		{
			descriptions_[key]=description;
		}
	}
	
	void Param::setValue(const String& key, float value, const String& description)
	{
		values_[key] = value;
		if(!description.empty())
		{
			descriptions_[key]=description;
		}
	}

	void Param::setValue(const String& key, double value, const String& description)
	{
		values_[key] = value;
		if(!description.empty())
		{
			descriptions_[key]=description;
		}
	}

	void Param::setValue(const String& key, const String& value, const String& description)
	{
		values_[key] = value;
		if(!description.empty())
		{
			descriptions_[key]=description;
		}
	}
	
	const DataValue& Param::getValue(const String& key) const
	{
		map<String, DataValue>::const_iterator it=values_.find(key);
		if (it!=values_.end())
		{
			return it->second;
		}
		return DataValue::EMPTY;
			
	}
	
	const String& Param::getDescription(const String& key) const
	{
		map<String, String>::const_iterator it=descriptions_.find(key);
		if (it!=descriptions_.end())
		{
			return it->second;
		}
		return String::EMPTY;
			
	}
	
	void Param::insert(String prefix, const Param& para)
	{
		if (prefix.empty() )
		{
			for(map<String,DataValue>::const_iterator it = para.values_.begin(); it != para.values_.end();++it)
			{
				values_[it->first]=it->second;
			}
			for(map<String, String>::const_iterator it = para.descriptions_.begin(); it != para.descriptions_.end();++it)
			{
				descriptions_[it->first]=it->second;
			}
		}
		else
		{
			prefix.ensureLastChar(':');
			for(map<String,DataValue>::const_iterator it = para.values_.begin(); it != para.values_.end();++it)
			{
				values_[prefix+it->first]=it->second;
			}
			for(map<String,String>::const_iterator it = para.descriptions_.begin(); it != para.descriptions_.end();++it)
			{
				descriptions_[prefix+it->first]=it->second;
			}
		}
	}

	void Param::setDefaults(const Param& defaults, String prefix, bool showMessage)
	{
		if ( !prefix.empty() )
		{
			prefix.ensureLastChar(':');
		}	
		
		for(map<String,DataValue>::const_iterator it = defaults.values_.begin(); it != defaults.values_.end();++it)
		{
			if (values_.find(prefix+it->first)==values_.end())
			{
				if (showMessage) cout << "Setting " << prefix+it->first << " to " << it->second << endl;
				values_[prefix+it->first]=it->second;
			}
		}
		for(map<String,String>::const_iterator it = defaults.descriptions_.begin(); it != defaults.descriptions_.end();++it)
		{
			if (descriptions_.find(prefix+it->first)==descriptions_.end())
			{
				descriptions_[prefix+it->first]=it->second;
			}
		}
	}
	
	void Param::remove(const String& prefix)
	{
		map<String,DataValue>::iterator it = values_.lower_bound(prefix);
		while (it!=values_.end())
		{
			if (it->first.substr(0,prefix.size())==prefix)
			{
				values_.erase(it);
			}
			else
			{
				break;
			}
			it = values_.lower_bound(prefix);
		}

		//delete descriptions
		map<String,String>::iterator it2 = descriptions_.lower_bound(prefix);
		while (it2!=descriptions_.end())
		{
			if (it2->first.substr(0,prefix.size())==prefix)
			{
				descriptions_.erase(it2);
			}
			else
			{
				break;
			}
			it2 = descriptions_.lower_bound(prefix);
		}

	}

	Param Param::copy(const String& prefix, bool remove_prefix, String new_prefix) const
	{
		if (!new_prefix.empty())
		{
			new_prefix.ensureLastChar(':');
		}
		
		Param out;
		String key;
		for ( map<String,DataValue>::const_iterator it = values_.lower_bound(prefix);
					(it != values_.end()) && (it->first.size() >= prefix.size()) && (it->first.substr(0,prefix.size())==prefix);
					++it
				)
		{
			//remove old prefix
			if (remove_prefix)
			{
				key = it->first.substr(prefix.size(),it->first.size() - prefix.size());
			}
			else
			{
				key = it->first;
			}
			
			// add new prefix
			if (!new_prefix.empty())
			{
				key = new_prefix + key;
			}
			
			out.values_[key]=it->second;
		}
		for ( map<String,String>::const_iterator it = descriptions_.lower_bound(prefix);
					(it != descriptions_.end()) && (it->first.size() >= prefix.size()) && (it->first.substr(0,prefix.size())==prefix);
					++it
				)
		{
			//remove old prefix
			if (remove_prefix)
			{
				key = it->first.substr(prefix.size(),it->first.size() - prefix.size());
			}
			else
			{
				key = it->first;
			}
			
			// add new prefix
			if (!new_prefix.empty())
			{
				key = new_prefix + key;
			}
			
			out.descriptions_[key]=it->second;
		}
		return out;
	}

	Param Param::copyWithInherit(const String& old_prefix, const String& new_prefix) const
	{
		if ( *old_prefix.rbegin() != ':' )
		{
			return copy( old_prefix, true, new_prefix );
		}
		else
		{
			Param result = copy( old_prefix, true, "" );
			int inheritance_steps = 0;
			DataValue const * inherit_path_value = & result.getValue("inherit");
			while ( ! inherit_path_value->isEmpty() )
			{
				String inherit_path = inherit_path_value->toString();
				if ( ++inheritance_steps > inheritance_steps_max )
				{
					throw Exception::ParseError
						( __FILE__, __LINE__, __PRETTY_FUNCTION__,
							"<ITEM name=\"inherit\" ... />",
							String("Too many inheritance steps (")+inheritance_steps_max+" allowed).  Perhaps there is a cycle?"
							+" old_prefix="+old_prefix+ " new_prefix="+new_prefix+" inherit_path="+inherit_path
						);
				}
				result.remove("inherit");
				result.setDefaults( copy(inherit_path+':',true,"" ), "", false);
				inherit_path_value = & result.getValue("inherit");
			}
			return result.copy("",true,new_prefix);
		}
	}


	void Param::store(const String& filename) const throw (Exception::UnableToCreateFile)
	{
		String up, down ,key, key_without_prefix, new_prefix ,type, prefix = "";
		UInt common, level=1;
		
		ofstream os;
		os.open (filename.c_str(), ofstream::out);
		
		if(!os)
		{
			 throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}
		
  	os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
  	os << "<PARAMETERS>\n";
		for(map<String,DataValue>::const_iterator it = values_.begin(); it != values_.end();++it)
		{
			//init variables
			key = it->first;
			if (key.find(":")==string::npos)
			{
				key_without_prefix = key;
				new_prefix = "";
			}
			else
			{
				key_without_prefix = key.substr(key.rfind(":")+1,key.size());
				new_prefix = key.substr(0,key.rfind(":")+1);
			} 
			//common prefix
			common=0;
			for (UInt i=0;i<min(key.size(),prefix.size());++i)
			{
				if (prefix[i]!=key[i])
				{
					break;
				}
				if (prefix[i]==':')
				{
					common = i+1;
				}
			}
			//cout << "key_wo: "<<key_without_prefix<<endl;
			//cout << "key   : "<<key<<endl;
			//cout << "prefix: "<<prefix<<endl;
			//cout << "|||   : "<<key.substr(0, common)<<endl;
			//write down
			down = prefix.substr(common, prefix.size());
			if (down!="")
			{
 				//cout << "  <-  : "<<down<<endl;
				for (UInt i = 0; i < down.size();++i)
				{
					if (down[i]==':')
					{
						--level;
						String tmp = string (2*level,' ');
						os << tmp.c_str();
						os << "</NODE>\n";						
					}
				}	
			}
			
			//write up
			up = key.substr(common, key.size()-common-key_without_prefix.size());
			String nodepath = key.substr(0,common);
			if (!nodepath.empty())
			{
				nodepath = nodepath.substr(0,-1);
			}
			if (up!="")
			{
				//cout << "  ->  : "<<up<<endl;
				while (up != "")
				{
					//name
					String nodename = up.substr(0,up.find(":"));
					os << String (2*level,' ') << "<NODE name=\"" << nodename << "\"";
					
					//description
					if (nodepath.empty())
					{
						nodepath = nodename;
					}
					else
					{
						nodepath = nodepath + ":" + nodename;
					}
					//cout << "NODE: '" << nodename << "' Path: '" << nodepath << "' Common: '" << key.substr(0,common) << "'" << endl;
					map<String, String>::const_iterator iter=descriptions_.find(nodepath);
					if(iter!=descriptions_.end())
					{
						//replace double quotes and newlines in descriptions
						String tmp = iter->second;
						tmp.substitute('"','\'');
						tmp.substitute("\n","#br#");
						//cout << "DESCRIPTION: " << iter->second << endl;
						os << " description=\"" << tmp <<"\"";
					}
					
					//closing tag and loop updates
					os << ">\n";
					++level;
					up = up.substr(nodename.size()+1,up.size());
				}				
			}
			
			//write item
			String tmp = string (2*level,' ');
			os << tmp.c_str();
			if (it->second.valueType()==DataValue::INTVALUE || it->second.valueType()==DataValue::LONVALUE || it->second.valueType()==DataValue::SHOVALUE  )
			{
				type = "int";
			}
			if (it->second.valueType()==DataValue::FLOVALUE || it->second.valueType()==DataValue::DOUVALUE )
			{
				type = "float";
			}
			if (it->second.valueType()==DataValue::STRVALUE )
			{
				type = "string";
			}
			
			if(it->second.valueType()!=DataValue::EMPTYVALUE)
			{
				map<String, String>::const_iterator iter=descriptions_.find(key);
				if(iter!=descriptions_.end())
				{
					String d = iter->second;
					d.substitute('"','\'');
					d.substitute("\n","#br#");
					tmp = "<ITEM name=\""+key_without_prefix+"\" value=\""+it->second.toString()+"\" type=\""+type+"\" description=\""+d+"\" />\n";
				}
				else
				{
					tmp = "<ITEM name=\""+key_without_prefix+"\" value=\""+it->second.toString()+"\" type=\""+type+"\" />\n";
				}
				
				os << tmp.c_str();					
			}
				
			//set new prefix
			prefix = new_prefix;
			//cout << "newprf: "<<new_prefix<<endl;
			//cout << endl;				
		}
		
		//close remaining prefix tags
		down = prefix;
		if (down!="")
		{
//				cout << "  <-  : "<<down<<endl;
			for (UInt i = 0; i < down.size();++i)
			{
				if (down[i]==':')
				{
					--level;
					String tmp = string (2*level,' ');
					os << tmp.c_str();
					os << "</NODE>\n";						
				}
			}	
		}
		
		os << "</PARAMETERS>\n";
    os.close();
	}
	
	void Param::load(const String& filename) throw (FileNotFound,ParseError)
	{
   	//try to open file
		if (!File::exists(filename))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
		
		// initialize parser
		try 
		{
			xercesc::XMLPlatformUtils::Initialize();
		}
		catch (const xercesc::XMLException& toCatch) 
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Error during initialization: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
	  }

		xercesc::SAX2XMLReader* parser = xercesc::XMLReaderFactory::createXMLReader();
		parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces,false);
		parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpacePrefixes,false);
		Internal::ParamXMLHandler handler(values_, descriptions_, filename);
		parser->setContentHandler(&handler);
		parser->setErrorHandler(&handler);
		
		xercesc::LocalFileInputSource source( xercesc::XMLString::transcode(filename.c_str()) );
		try 
    {
    	parser->parse(source);
    	delete(parser);
    }
    catch (const xercesc::XMLException& toCatch) 
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("XMLException: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
    }
    catch (const xercesc::SAXException& toCatch) 
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("SAXException: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
    }
	}

	void Param::parseCommandLine(const int argc , char**argv, String prefix)
	{
		//determine prefix
		if (prefix!="")
		{
			prefix.ensureLastChar(':');
		}
		
		//parse arguments
    String arg,arg1;
    for(int i = 1; i < argc; ++i )
    { 
      //load the current and next argument:  arg and arg1 ("" at the last argument)
      arg = argv[i];
      if (i+1<argc)
      {
      	arg1 = argv[i+1];
      }
    	else
    	{
    		arg1 = "";
    	}
    	
    	//cout << "Parse: '"<< arg << "' '" << arg1 << "'" << endl;
    	
      //flag (option without text argument)
      if(arg[0]  == '-' && (arg1[0] == '-' || arg1 ==""))
      {
	      	values_[prefix+arg] = "";
      }
      //option with argument
      else if(arg[0]  == '-')
      {
	      	values_[prefix+arg] = arg1;
	      	++i;
      }      
      //just text arguments (not preceded by an option)
      else
      {
      	if (values_[prefix+"misc"].isEmpty())
      	{
      		values_[prefix+"misc"] = arg;
      	}
      	else
      	{
      		values_[prefix+"misc"] = string(values_[prefix+"misc"])+" "+arg;
      	}
      }
    }
	}

	void Param::parseCommandLine(const int argc , char **argv, const map<String, String>& options_with_argument, const std::map<String, String>& options_without_argument, const String& misc, const String& unknown)
	{
		//determine misc key
    String misc_key = misc;

		//determine unknown key
    String unknown_key = unknown;

		//parse arguments
    String arg,arg1;
    for(int i = 1; i < argc; ++i )
    { 
      //load the current and next argument:  arg and arg1 ("" at the last argument)
      arg = argv[i];
      if (i+1<argc)
      {
      	arg1 = argv[i+1];
      }
    	else
    	{
    		arg1 = "";
    	}

			//without argument
			if (options_without_argument.find(arg)!=options_without_argument.end())
			{
				values_[options_without_argument.find(arg)->second] = "true";
			}
			//with argument
			else if (options_with_argument.find(arg)!=options_with_argument.end())
			{
				//next argument is not a option
				if (arg1[0]!='-')
				{
					values_[options_with_argument.find(arg)->second] = arg1;
					++i;
				}
				//next argument is a option
				else
				{
					values_[options_with_argument.find(arg)->second] = "";
				}
			}
			//unknown option
			else if (arg[0]=='-')
			{
      	if (values_[unknown_key].isEmpty())
      	{
      		values_[unknown_key] = arg;
      	}
      	else
      	{
      		values_[unknown_key] = string(values_[unknown_key])+" "+arg;
      	}				
			}
			//just text argument
			else
			{
	    	if (values_[misc_key].isEmpty())
	    	{
	    		values_[misc_key] = arg;
	    	}
	    	else
	    	{
	    		values_[misc_key] = string(values_[misc_key])+" "+arg;
	    	}				
			}
    }
	}

	ostream& operator << (ostream& os, const Param& param)
 	{
		for (map<String,DataValue>::const_iterator it = param.values_.begin(); it != param.values_.end();++it)
		{
			map<String, String>::const_iterator iter=param.descriptions_.find(it->first);
			String description;
			if (iter!=param.descriptions_.end())
			{
				description=" :"+iter->second;
			}
		
			os << "\""<<it->first<< "\"  ->  \""<< it->second.toString()<< "\"" <<description<<endl;
		}
		return os;
	}

	UInt Param::size() const
	{
		return values_.size();
	}

	bool Param::empty() const
	{
		return values_.empty();
	}

	void Param::clear()
	{
		values_.clear();
		descriptions_.clear();
	}

	void Param::checkDefaults(const String& name, const Param& defaults, String prefix, std::ostream& os) const
	{
		//Extract right parameters
		map<String,DataValue> check_values;
		map<String,String> check_descriptions;
		if ( prefix.empty() )
		{
			check_values = values_;
			check_descriptions=descriptions_;
		}	
		else
		{
			prefix.ensureLastChar(':');
			check_values = copy(prefix,true).values_;
			check_descriptions = copy(prefix,true).descriptions_;
		}
		//check
		for(map<String,DataValue>::const_iterator it = check_values.begin(); it != check_values.end();++it)
		{
			if (defaults.values_.find(it->first)==defaults.values_.end())
			{
				os << "Warning: " << name << " received the unknown parameter '" << it->first << "'";
				if (!prefix.empty()) os << " in '" << prefix << "'";
				os << "!" << endl;
			}
		}
	}

	void Param::setDescription(const String& location, const String& description)
	{
		map<String,DataValue>::iterator it = values_.lower_bound(location);
		if (it!=values_.end() && description!="")
		{
			descriptions_[location] = description;
		}
	}

}	//namespace
