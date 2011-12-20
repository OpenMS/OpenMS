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
// $Maintainer: Erhan Kenar $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cctype> // for "isalpha"

#include <OpenMS/FORMAT/HANDLERS/ParamXMLHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace std;
using namespace OpenMS::Exception;
using namespace OpenMS::Internal;

namespace OpenMS
{

	//********************************* ParamEntry **************************************
	Param::ParamEntry::ParamEntry()
		: name(),
			description(),
			value(),
			tags(),
			min_float(-std::numeric_limits<DoubleReal>::max()),
			max_float(std::numeric_limits<DoubleReal>::max()),
			min_int(-std::numeric_limits<Int>::max()),
			max_int(std::numeric_limits<Int>::max()),
			valid_strings()
	{
	}

	Param::ParamEntry::ParamEntry(const String& n, const DataValue& v, const String& d, const StringList& t)
		: name(n),
			description(d),
			value(v),
			tags(),
			min_float(-std::numeric_limits<DoubleReal>::max()),
			max_float(std::numeric_limits<DoubleReal>::max()),
			min_int(-std::numeric_limits<Int>::max()),
			max_int(std::numeric_limits<Int>::max()),
			valid_strings()
	{
		//add tags
		for (Size i=0; i<t.size(); ++i)
		{
			tags.insert(t[i]);
		}
		//check name
		if (name.has(':'))
		{
			cerr << "Error ParamEntry name must not contain ':' characters!" << endl;
		}
	}
	
	Param::ParamEntry::~ParamEntry()
	{
	}
	
	bool Param::ParamEntry::isValid(String& message) const
	{
			if (value.valueType()==DataValue::STRING_VALUE)
			{
        if (valid_strings.size()!=0)
				{
          bool ok = false;
          if (std::find(valid_strings.begin(),valid_strings.end(), value) != valid_strings.end())
          {
            ok = true;
          }
          else if (std::find(tags.begin(), tags.end(), "input file") != tags.end() || std::find(tags.begin(), tags.end(), "output file") != tags.end())
          {
            //do not check restrictions on file names for now
            ok = true;
          }

          if (!ok)
          {
            String valid;
            valid.concatenate(valid_strings.begin(),valid_strings.end(),",");
            message = "Invalid string parameter value '"+(String)value+"' for parameter '"+name+"' given! Valid values are: '"+valid+"'.";
            return false;
          }
				}
			}
			else if(value.valueType()==DataValue::STRING_LIST)
			{
				String str_value;
				StringList ls_value = (StringList) value;
				for (Size i = 0; i < ls_value.size(); ++i)
				{
					str_value = ls_value[i];
					
          if (valid_strings.size()!=0)
					{
            bool ok = false;
            if (std::find(valid_strings.begin(),valid_strings.end(), str_value) != valid_strings.end())
            {
              ok = true;
            }
            else if (std::find(tags.begin(), tags.end(), "input file") != tags.end() || std::find(tags.begin(), tags.end(), "output file") != tags.end())
            {
              //do not check restrictions on file names for now
              ok = true;
            }

            if (!ok)
            {
              String valid;
              valid.concatenate(valid_strings.begin(),valid_strings.end(),",");
              message = "Invalid string parameter value '"+str_value+"' for parameter '"+name+"' given! Valid values are: '"+valid+"'.";
              return false;
            }
					}
				}	
			}
			else if (value.valueType()==DataValue::INT_VALUE)
			{
				Int tmp = value;
				if ((min_int != -std::numeric_limits<Int>::max() && tmp < min_int) || (max_int!=std::numeric_limits<Int>::max() && tmp > max_int))
				{
					message = String("Invalid integer parameter value '")+String(tmp)+"' for parameter '"+name+"' given! The valid range is: ["+min_int+":"+max_int+"].";
					return false;
				}
			}
			else if(value.valueType()==DataValue::INT_LIST)
			{
				Int int_value;
				IntList ls_value =(IntList) value;
				for (Size i = 0; i < ls_value.size(); ++i)
				{
					int_value = ls_value[i];
					if ((min_int != -std::numeric_limits<Int>::max() && int_value < min_int) || (max_int!=std::numeric_limits<Int>::max() && int_value > max_int))
					{
						message = String("Invalid integer parameter value '")+int_value+"' for parameter '"+name+"' given! The valid range is: ["+min_int+":"+max_int+"].";
						return false;
					}
				}
			}
			else if (value.valueType()==DataValue::DOUBLE_VALUE)
			{
				DoubleReal tmp = value;
				if ((min_float!=-std::numeric_limits<DoubleReal>::max() && tmp < min_float) || (max_float!=std::numeric_limits<DoubleReal>::max() && tmp > max_float))
				{
					message = String("Invalid double parameter value '")+tmp+"' for parameter '"+name+"' given! The valid range is: ["+min_float+":"+max_float+"].";
					return false;
			}
			}
			else if(value.valueType()==DataValue::DOUBLE_LIST)
			{
				DoubleReal dou_value; 
				DoubleList ls_value = (DoubleList)value;
				for (Size i = 0; i < ls_value.size(); ++i)
				{
					dou_value = ls_value[i];
					if ((min_float!=-std::numeric_limits<DoubleReal>::max() && dou_value < min_float) || (max_float!=std::numeric_limits<DoubleReal>::max() && dou_value > max_float))
					{
						message = String("Invalid double parameter value '")+dou_value+"' for parameter '"+name+"' given! The valid range is: ["+min_float+":"+max_float+"].";
						return false;
					}	
				}
			}	
			return true;
	}
	
	bool Param::ParamEntry::operator==(const ParamEntry& rhs) const
	{
		return name==rhs.name && value==rhs.value;
	}
	
	//********************************* ParamNode **************************************
	Param::ParamNode::ParamNode()
		: name(),
			description(),
			entries(),
			nodes()
	{
	}

	Param::ParamNode::ParamNode(const String& n, const String& d)
		: name(n),
			description(d),
			entries(),
			nodes()
	{
		if (name.has(':'))
		{
			cerr << "Error ParamNode name must not contain ':' characters!" << endl;
		}
	}

	Param::ParamNode::~ParamNode()
	{
	}

	bool Param::ParamNode::operator==(const ParamNode& rhs) const
	{
		if (name!=rhs.name || entries.size()!=rhs.entries.size() || nodes.size()!=rhs.nodes.size()) return false;
		
		//order of sections / entries should not matter
		for (Size i=0; i< entries.size(); ++i)
		{
			if (find(rhs.entries.begin(),rhs.entries.end(),entries[i])==rhs.entries.end())
			{
				return false;
			}
		}
		for (Size i=0; i< nodes.size(); ++i)
		{
			if (find(rhs.nodes.begin(),rhs.nodes.end(),nodes[i])==rhs.nodes.end())
			{
				return false;
			}
		}
		
		return true;
	}

	Param::ParamNode::EntryIterator Param::ParamNode::findEntry(const String& name)
	{
		for(EntryIterator it = entries.begin(); it!=entries.end(); ++it)
		{
			if(it->name==name)
			{
				return it;
			}
		}
		return entries.end();
	}
	
	Param::ParamNode::NodeIterator Param::ParamNode::findNode(const String& name)
	{
		for(NodeIterator it = nodes.begin(); it!=nodes.end(); ++it)
		{
			if(it->name==name)
			{
				return it;
			}
		}
		return nodes.end();
	}
	
	Param::ParamNode* Param::ParamNode::findParentOf(const String& name)
	{
		//cout << "findParentOf nodename: " << this->name << " - nodes: " << this->nodes.size() << " - find: "<< name << endl;
		if (!name.has(':')) // we are in the right child
		{
			//check if a node or entry prefix match
			for (Size i=0; i<nodes.size();++i)
			{
				if (nodes[i].name.hasPrefix(name)) return this;
			}
			for (Size i=0; i<entries.size();++i)
			{
				if (entries[i].name.hasPrefix(name)) return this;
			}
			return 0;
		}
		else //several subnodes to browse through
		{
			String prefix = name.prefix(':');
			//cout << " - Prefix: '" << prefix << "'" << endl;
			NodeIterator it = findNode(prefix);
			if (it==nodes.end()) //subnode not found
			{
				return 0;
			}
			//recursively call findNode for the rest of the path
			String new_name = name.substr(it->name.size()+1);
			//cout << " - Next name: '" << new_name << "'" << endl;
			return it->findParentOf(new_name);
		}
	}

	Param::ParamEntry* Param::ParamNode::findEntryRecursive(const String& name)
	{
		ParamNode* parent = findParentOf(name);
		if (parent==0) return 0;

		EntryIterator it = parent->findEntry(suffix(name));
		if (it==parent->entries.end()) return 0;
		
		return &(*it);
	}
	
	void Param::ParamNode::insert(const ParamNode& node, const String& prefix)
	{
		//cerr << "INSERT NODE  " << node.name << " (" << prefix << ")" << endl;
		String prefix2 = prefix + node.name;
		
		ParamNode* insert_node = this;
		while (prefix2.has(':'))
		{
			String name = prefix2.prefix(':');
			//check if the node already exists
			NodeIterator it = insert_node->findNode(name);
			if (it!=insert_node->nodes.end()) //exists
			{
				insert_node = &(*it);
			}
			else //create it
			{
				insert_node->nodes.push_back(ParamNode(name,""));
				insert_node = &(insert_node->nodes.back());
				//cerr << " - Created new node: " << insert_node->name << endl;
			}
			//remove prefix
			prefix2 = prefix2.substr(name.size()+1);
		}
		
		//check if the node already exists
		NodeIterator it = insert_node->findNode(prefix2);
		if (it!=insert_node->nodes.end()) //append nodes and entries
		{
			for (ConstNodeIterator it2=node.nodes.begin(); it2!=node.nodes.end(); ++it2)
			{
				it->insert(*it2);
			}
			for (ConstEntryIterator it2=node.entries.begin(); it2!=node.entries.end(); ++it2)
			{
				it->insert(*it2);
			}
			if (it->description=="" || node.description!="") //replace description if not empty in new node
			{
				it->description = node.description;
			}
		}
		else //insert it
		{
			Param::ParamNode tmp(node);
			tmp.name = prefix2;
			insert_node->nodes.push_back(tmp);
		}
	}

	void Param::ParamNode::insert(const ParamEntry& entry, const String& prefix)
	{
		//cerr << "INSERT ENTRY " << entry.name << " (" << prefix << ")" << endl;
		String prefix2 = prefix + entry.name;
		//cerr << " - inserting: " << prefix2 << endl;
		
		ParamNode* insert_node = this;
		while (prefix2.has(':'))
		{
			String name = prefix2.prefix(':');
			//cerr << " - looking for node: " << name << endl;
			//look up if the node already exists
			NodeIterator it = insert_node->findNode(name);
			if (it!=insert_node->nodes.end()) //exists
			{
				insert_node = &(*it);
			}
			else //create it
			{
				insert_node->nodes.push_back(ParamNode(name,""));
				insert_node = &(insert_node->nodes.back());
				//cerr << " - Created new node: " << insert_node->name << endl;
			}
			//remove prefix
			prefix2 = prefix2.substr(name.size()+1);
			//cerr << " - new prefix: " << prefix2 << endl;
		}
		
		//check if the entry already exists
		//cerr << " - final entry name: " << prefix2 << endl;
		EntryIterator it = insert_node->findEntry(prefix2);
		if (it!=insert_node->entries.end()) //overwrite entry
		{
			it->value = entry.value;
			it->tags = entry.tags;
			if (it->description=="" || entry.description!="") //replace description if not empty in new entry
			{
				it->description = entry.description;
			}
		}
		else //insert it
		{
			Param::ParamEntry tmp(entry);
			tmp.name = prefix2;
			insert_node->entries.push_back(tmp);
		}
	}
	
	Size Param::ParamNode::size() const
	{
		Size subnode_size = 0;
		for (vector<ParamNode>::const_iterator it=nodes.begin(); it!=nodes.end(); ++it)
		{
			subnode_size += it->size();
		}
		return entries.size() + subnode_size;
	}
	
	String Param::ParamNode::suffix(const String& key) const
	{
		if (key.has(':'))
		{
			return key.suffix(':');
		}
		return key;
	}
	
	//********************************* Param **************************************
	
	Param::Param()
		: XMLFile("/SCHEMAS/Param_1_3.xsd","1.3"),
			root_("ROOT","") 
	{
	}

	Param::Param(const Param& rhs)
		: Internal::XMLFile(rhs),
			root_(rhs.root_)
	{
	}

	Param::~Param()
	{	
	}

	Param& Param::operator=(const Param& rhs)
	{
		root_ = rhs.root_;
		return *this;
	}

	Param::Param(const ParamNode& node)
		: root_(node) 
	{
		root_.name="ROOT";
		root_.description="";
	}

	bool Param::operator==(const Param& rhs) const
	{
		return root_ == rhs.root_;
	}	
	
	void Param::setValue(const String& key, const DataValue& value, const String& description, const StringList& tags)
	{
		root_.insert(ParamEntry("",value,description,tags),key);
	}

	void Param::setValidStrings(const String& key, const std::vector<String>& strings)
	{
		ParamEntry& entry = getEntry_(key);
		//check if correct parameter type
		if (entry.value.valueType()!=DataValue::STRING_VALUE && entry.value.valueType() != DataValue::STRING_LIST) 
		{
			throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,key);
		}
		//check for commas
		for (Size i=0; i<strings.size(); ++i)
		{
			if (strings[i].has(','))
			{
				throw InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Comma characters in Param string restrictions are not allowed!");
			}
		}
		entry.valid_strings = strings;
	}

	void Param::setMinInt(const String& key, Int min)
	{
		ParamEntry& entry = getEntry_(key);
		if (entry.value.valueType()!=DataValue::INT_VALUE && entry.value.valueType()!=DataValue::INT_LIST) 
		{
			throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,key);
		}
		entry.min_int = min;
	}

	void Param::setMaxInt(const String& key, Int max)
	{
		ParamEntry& entry = getEntry_(key);
		if (entry.value.valueType()!=DataValue::INT_VALUE && entry.value.valueType()!=DataValue::INT_LIST) 
		{
			throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,key);
		}
		entry.max_int = max;
	}

	void Param::setMinFloat(const String& key, DoubleReal min)
	{
		ParamEntry& entry = getEntry_(key);
		if (entry.value.valueType()!=DataValue::DOUBLE_VALUE && entry.value.valueType() !=DataValue::DOUBLE_LIST) 
		{
			throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,key);
		}
		entry.min_float = min;
	}

	void Param::setMaxFloat(const String& key, DoubleReal max)
	{
		ParamEntry& entry = getEntry_(key);
		if (entry.value.valueType()!=DataValue::DOUBLE_VALUE && entry.value.valueType() !=DataValue::DOUBLE_LIST) 
		{
			throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,key);
		}
		entry.max_float = max;
	}

	const DataValue& Param::getValue(const String& key) const
	{
		return getEntry_(key).value;
	}
	
	const String& Param::getSectionDescription(const String& key) const
	{
		//This variable is used instead of String::EMPTY as the method is used in
		//statis initialization und thus cannot rely on String::EMPTY been initialized.
		static String empty;

		ParamNode* node = root_.findParentOf(key);
		if (node==0)
		{
			return empty;
		}
		
		Param::ParamNode::NodeIterator it = node->findNode(node->suffix(key));
		if (it==node->nodes.end())
		{
			return empty;
		}
		
		return it->description;
	}
	
	void Param::insert(const String& prefix, const Param& param)
	{
		//cerr << "INSERT PARAM (" << prefix << ")" << endl;
		for (Param::ParamNode::NodeIterator it=param.root_.nodes.begin(); it!=param.root_.nodes.end(); ++it)
		{
			root_.insert(*it,prefix);
		}
		for (Param::ParamNode::EntryIterator it=param.root_.entries.begin(); it!=param.root_.entries.end(); ++it)
		{
			root_.insert(*it,prefix);
		}
	}

	void Param::setDefaults(const Param& defaults, const String& prefix, bool showMessage)
	{
		String prefix2 = prefix;
		if (prefix2!="")
		{
			prefix2.ensureLastChar(':');
		}
		
		String pathname;
		for(Param::ParamIterator it = defaults.begin(); it != defaults.end();++it)
		{
			if (!exists(prefix2 + it.getName()))
			{
				if (showMessage) cerr << "Setting " << prefix2+it.getName() << " to " << it->value << endl;
				String name = prefix2+it.getName();
				root_.insert(ParamEntry("", it->value, it->description), name);
				//copy tags
				for (set<String>::const_iterator tag_it=it->tags.begin(); tag_it!=it->tags.end(); ++tag_it)
				{
					addTag(name,*tag_it);
				}
				//copy restrictions
				if (it->value.valueType()==DataValue::STRING_VALUE || it->value.valueType()==DataValue::STRING_LIST)
				{
					setValidStrings(name,it->valid_strings);
				}
				else if (it->value.valueType()==DataValue::INT_VALUE || it->value.valueType() == DataValue::INT_LIST)
				{
					setMinInt(name,it->min_int);
					setMaxInt(name,it->max_int);
				}
				else if (it->value.valueType()==DataValue::DOUBLE_VALUE || it->value.valueType()==DataValue::DOUBLE_LIST)
				{
					setMinFloat(name,it->min_float);
					setMaxFloat(name,it->max_float);
				}
			}
			
			//copy section descriptions
			const std::vector< ParamIterator::TraceInfo >& trace = it.getTrace();
			for(std::vector< ParamIterator::TraceInfo >::const_iterator it2 = trace.begin(); it2!=trace.end(); ++it2)
			{
				if (it2->opened)
				{
					pathname += it2->name + ":";
				}
				else
				{
					pathname.resize(pathname.size() - it2->name.size() -1);
				}
				String real_pathname = pathname.chop(1); //remove ':' at the end
				if (real_pathname != "")
				{
					String description_old = "";
					String description_new = "";
					
					description_old = getSectionDescription(prefix+real_pathname);
					description_new = defaults.getSectionDescription(real_pathname);
					if (description_old=="")
					{
						//cerr << "## Setting description of " << prefix+real_pathname << " to"<< endl;
						//cerr << "## " << description_new << endl;
						setSectionDescription(prefix2+real_pathname, description_new);
					}
				}
			}
		}
	}

	void Param::remove(const String& key)
	{
    String keyname = key;
		if (key.hasSuffix(':')) // delete section
		{
      keyname = key.chop(1);

      ParamNode* node_parent = root_.findParentOf(keyname);
			if (node_parent!=0)
			{
				Param::ParamNode::NodeIterator it = node_parent->findNode(node_parent->suffix(keyname));
				if (it!=node_parent->nodes.end())
				{
          String name = it->name;
          node_parent->nodes.erase(it); // will automatically delete subnodes
          if (node_parent->nodes.empty()  && node_parent->entries.empty())
          {
            // delete last section name (could be partial)
            remove( keyname.chop(name.size()) ); // keep last ':' to indicate deletion of a section
          }
				}
			}
    }
    else
    {
		  ParamNode* node = root_.findParentOf(keyname);
		  if (node!=0)
		  {
        String entryname = node->suffix(keyname); // get everything beyond last ':'
        Param::ParamNode::EntryIterator it = node->findEntry(entryname);
			  if (it!=node->entries.end())
			  {
          node->entries.erase(it); // delete entry
          if (node->nodes.empty()  && node->entries.empty())
          {
            // delete if section is now empty
            remove( keyname.chop(entryname.size()) ); // keep last ':' to indicate deletion of a section
          }
			  }
		  }
    }
	}

	void Param::removeAll(const String& prefix)
	{
		if (prefix.hasSuffix(':')) //we have to delete one node only (and its subnodes)
		{
			ParamNode* node = root_.findParentOf(prefix.chop(1));
			if (node!=0)
			{
				Param::ParamNode::NodeIterator it = node->findNode(node->suffix(prefix.chop(1)));
				if (it!=node->nodes.end())
				{
          String name = it->name;
          node->nodes.erase(it); // will automatically delete subnodes
          if (node->nodes.empty()  && node->entries.empty())
          {
            // delete last section name (could be partial)
            removeAll(prefix.chop(name.size()+1)); // '+1' for the tailing ':'
          }
				}
			}
		}
		else //we have to delete all entries and nodes starting with the prefix
		{
			ParamNode* node = root_.findParentOf(prefix);
			if (node!=0)
			{			
				String suffix = node->suffix(prefix); // name behind last ":"
				
				for (Param::ParamNode::NodeIterator it = node->nodes.begin(); it!=node->nodes.end();/*do nothing*/)
				{
					if (it->name.hasPrefix(suffix))
					{
						it = node->nodes.erase(it);
					}
					else if (it!=node->nodes.end())
					{
						++it;
					}
				}
				for (Param::ParamNode::EntryIterator it = node->entries.begin(); it!=node->entries.end();/*do nothing*/)
				{
					if (it->name.hasPrefix(suffix))
					{
						it = node->entries.erase(it);
					}
					else if (it!=node->entries.end())
					{
						++it;
					}
				}
        // the parent node might now be empty (delete it as well - otherwise the trace might be broken)
        if (node->nodes.empty() && node->entries.empty())
        {
          // delete last section name (could be partial)
          removeAll(prefix.chop(suffix.size()));
        }
			}
		}
	}

	Param Param::copy(const String& prefix, bool remove_prefix) const
	{
		ParamNode out("ROOT","");

		ParamNode* node = root_.findParentOf(prefix);
		if (node==0)
		{
			return Param();
		}
		//we have to copy this node only
		if (prefix.hasSuffix(':'))
		{
			if (remove_prefix)
			{
				out = *node;
			}
			else
			{
				out.insert(*node,prefix.chop(node->name.size()+1));
			}
		}
		else //we have to copy all entries and nodes starting with the right suffix
		{
			String suffix = node->suffix(prefix);
			for (Param::ParamNode::NodeIterator it = node->nodes.begin(); it!=node->nodes.end(); ++it)
			{
				if (it->name.hasPrefix(suffix))
				{
					if (remove_prefix)
					{
						ParamNode tmp = *it;
						tmp.name = tmp.name.substr(suffix.size());
						out.insert(tmp);
					}
					else
					{
						out.insert(*it,prefix.chop(suffix.size()));
					}
				}
			}
			for (Param::ParamNode::EntryIterator it = node->entries.begin(); it!=node->entries.end(); ++it)
			{
				if (it->name.hasPrefix(suffix))
				{
					if (remove_prefix)
					{
						ParamEntry tmp = *it;
						tmp.name = tmp.name.substr(suffix.size());
						out.insert(tmp);
					}
					else
					{
						out.insert(*it,prefix.chop(suffix.size()));
					}
				}
			}
		}

		return Param(out);
	}

  void Param::store(const String& filename) const
  {
    //open file
    ofstream os_;
    ostream* os_ptr;
    if ( filename != "-" )
    {
      os_.open (filename.c_str(), ofstream::out);
      if(!os_)
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
      }
      os_ptr = &os_;
    }
    else
    {
      os_ptr = &std::cout;
    }

    //write to file stream
    writeXMLToStream(os_ptr);

    os_.close();
  }

  void Param::writeXMLToStream(ostream* os_ptr) const
	{
    // hint: the handling of 'getTrace()' is vulnerable to an unpruned tree (a path of nodes, but no entries in them), i.e.
    //       too many closing tags are written to the INI file, but no openening ones.
    //       This currently cannot happen, as removeAll() was fixed to prune the tree, just keep it in mind.

    ostream& os = *os_ptr;

		os.precision(writtenDigits<DoubleReal>());
		
  	os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
  	os << "<PARAMETERS version=\"" << getVersion() << "\" xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/Param_1_3.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";
		String indentation = "  ";
		ParamIterator it = begin();
		while(it != end())
		{
			//init variables
			string key = it.getName();
			
			//write opened/closed nodes
			const std::vector< ParamIterator::TraceInfo >& trace = it.getTrace();
			for(std::vector< ParamIterator::TraceInfo >::const_iterator it2 = trace.begin(); it2!=trace.end(); ++it2)
			{
				if (it2->opened) //opened node
				{
					String d = it2->description;
					//d.substitute('"','\'');
					d.substitute("\n","#br#");
					//d.substitute("<","&lt;");
					//d.substitute(">","&gt;");
					os << indentation  << "<NODE name=\"" << writeXMLEscape(it2->name) << "\" description=\"" << writeXMLEscape(d) << "\">" << "\n";
					indentation += "  ";
				}
				else //closed node
				{
					indentation.resize(indentation.size()-2);
					os << indentation << "</NODE>" << "\n";
				}
			}
			
			//write item
			if(it->value.valueType()!=DataValue::EMPTY_VALUE)
			{
				DataValue::DataType value_type = it->value.valueType();
				//write opening tag
				switch(value_type)
				{
					case DataValue::INT_VALUE:
						os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << it->value.toString() << "\" type=\"int\"";
						break;
					case DataValue::DOUBLE_VALUE:
						os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << it->value.toString() << "\" type=\"float\"";
						break;
					case DataValue::STRING_VALUE:
						os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << writeXMLEscape(it->value.toString()) << "\" type=\"string\"";
						break;
					case DataValue::STRING_LIST:
						os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << "\" type=\"string\"";
						break;
					case DataValue::INT_LIST:
						os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << "\" type=\"int\"";
						break;
					case DataValue::DOUBLE_LIST:
						os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << "\" type=\"float\"";
						break;
					default:
						break;
				};
				
				//replace all critical characters in description
				String d = it->description;
				//d.substitute("\"","'");
				d.substitute("\n","#br#");
				//d.substitute("<","&lt;");
				//d.substitute(">","&gt;");
				os << " description=\"" << writeXMLEscape(d) << "\"";
				
				//tags
				if (!it->tags.empty())
				{
					String list;
					for (set<String>::const_iterator tag_it=it->tags.begin(); tag_it!=it->tags.end(); ++tag_it)
					{
						if (!list.empty()) list += ",";
						list += *tag_it;
					}
					os << " tags=\"" << writeXMLEscape(list) << "\"";
				}

				//restrictions
				String restrictions = "";
				switch(value_type)
				{
					case DataValue::INT_VALUE:
					case DataValue::INT_LIST:	
						{
							bool min_set = (it->min_int!=-numeric_limits<Int>::max());
							bool max_set = (it->max_int!=numeric_limits<Int>::max());
							if (max_set || min_set)
							{
								if (min_set)
								{
									restrictions += String(it->min_int);
								}
								restrictions += ':';
								if (max_set)
								{
									restrictions += String(it->max_int);
								}
							}
						}
						break;
					case DataValue::DOUBLE_VALUE:
					case DataValue::DOUBLE_LIST:
						{
							bool min_set = (it->min_float!=-numeric_limits<DoubleReal>::max());
							bool max_set = (it->max_float!=numeric_limits<DoubleReal>::max());
							if (max_set || min_set)
							{
								if (min_set)
								{
									restrictions += String(it->min_float);
								}
								restrictions += ':';
								if (max_set)
								{
									restrictions += String(it->max_float);
								}
							}
						}
						break;
					case DataValue::STRING_VALUE:
					case DataValue::STRING_LIST:
						if (it->valid_strings.size()!=0)
						{
							restrictions.concatenate(it->valid_strings.begin(),it->valid_strings.end(),",");
						}
						break;
					default:
						break;
				};
				if (restrictions!="")
				{
					os << " restrictions=\"" << writeXMLEscape(restrictions) << "\"";
				}
				
				//finish opening tag
				switch(value_type)
				{
					case DataValue::INT_VALUE:
					case DataValue::DOUBLE_VALUE:
					case DataValue::STRING_VALUE:
						os << " />" <<  "\n";	
						break;
					case DataValue::STRING_LIST:
						{
							os << ">" <<  "\n";
							const StringList& list = it->value;
							for (Size i=0; i<list.size();++i)
							{
								os << indentation << "  <LISTITEM value=\"" << writeXMLEscape(list[i]) << "\"/>" << "\n";	
							}
							os << indentation << "</ITEMLIST>" << "\n";	
						}
						break;
					case DataValue::INT_LIST:
						{
							os << ">" <<  "\n";
							const IntList& list = it->value;
							for (Size i=0; i<list.size();++i)
							{
								os << indentation << "  <LISTITEM value=\"" << list[i] << "\"/>" << "\n";	
							}
							os << indentation << "</ITEMLIST>" << "\n";	
						}
						break;
					case DataValue::DOUBLE_LIST:
						{
							os << ">" <<  "\n";
							const DoubleList& list = it->value;
							for (Size i=0; i<list.size();++i)
							{
								os << indentation << "  <LISTITEM value=\"" << list[i] << "\"/>" << "\n";	
							}
							os << indentation << "</ITEMLIST>" << "\n";	
						}
						break;
					default:
						break;
				};
			}
			++it;
		}
		
    // if we had tags ...
    if (begin()!=end())
    {
		  //close remaining tags
		  const std::vector< ParamIterator::TraceInfo >& trace = it.getTrace();
		  for(std::vector< ParamIterator::TraceInfo >::const_iterator it2 = trace.begin(); it2!=trace.end(); ++it2)
		  {
        Size ss = indentation.size();
        indentation.resize(ss-2);
			  os << indentation << "</NODE>" << "\n";	
		  }
    }
		
    os << "</PARAMETERS>" << endl; // forces a flush
	}
	
	void Param::load(const String& filename)
	{
		Internal::ParamXMLHandler handler(*this, filename, schema_version_);
		parse_(filename, &handler);
	}


	void Param::parseCommandLine(const int argc , const char** argv, const String& prefix)
	{
		//determine prefix
		String prefix2 = prefix;
		if (prefix2!="")
		{
			prefix2.ensureLastChar(':');
		}
		
		//parse arguments
    String arg,arg1 ;
    for(int i = 1; i < argc; ++i )
    {
      //load the current and next argument:  arg and arg1 ("" at the last argument)
      arg = argv[i];
      arg1 = "";
      if (i+1<argc)
      {
      	arg1 = argv[i+1];
      }
    	
    	//it is a option when it starts with a '-' and the second character is not a number
    	bool arg_is_option = false;
    	if (arg.size()>=2 && arg[0]=='-' && arg[1]!='0' && arg[1]!='1' && arg[1]!='2' && arg[1]!='3' && arg[1]!='4' && arg[1]!='5' && arg[1]!='6' && arg[1]!='7' && arg[1]!='8' && arg[1]!='9') arg_is_option = true;
    	bool arg1_is_option = false;
    	if (arg1.size()>=2 && arg1[0]=='-' && arg1[1]!='0' && arg1[1]!='1' && arg1[1]!='2' && arg1[1]!='3' && arg1[1]!='4' && arg1[1]!='5' && arg1[1]!='6' && arg1[1]!='7' && arg1[1]!='8' && arg1[1]!='9') arg1_is_option = true;
    	
    	//cout << "Parse: '"<< arg << "' '" << arg1 << "'" << endl;
    	
      //flag (option without text argument)
      if(arg_is_option && arg1_is_option)
      {
	    	root_.insert(ParamEntry(arg,String(),""),prefix2);
      }
      //option with argument
      else if(arg_is_option && !arg1_is_option)
      {
      	root_.insert(ParamEntry(arg,arg1,""),prefix2);
      	++i;
      }      
      //just text arguments (not preceded by an option)
      else
      {

      	ParamEntry* misc_entry = root_.findEntryRecursive(prefix2+"misc");
      	if (misc_entry==0)
      	{
      		StringList sl;
      		sl << arg;
					// create "misc"-Node: 
      		root_.insert(ParamEntry("misc",sl,""),prefix2);
      	}
      	else
      	{
      		StringList sl = (StringList)misc_entry->value;
      		sl << arg;
      		misc_entry->value = sl;
      	}				
      }
    }
	}


	void Param::parseCommandLine(const int argc , const char** argv,const Map<String, String>& options_with_one_argument,const Map<String, String>& options_without_argument,const Map<String,String>& options_with_multiple_argument, const String& misc, const String& unknown)
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
      arg1 = "";
      if (i+1<argc)
      {
      	arg1 = argv[i+1];
      }
      
    	//it is a option when it starts with a '-' and the second character is not a number
    	bool arg_is_option = false;
    	if (arg.size()>=2 && arg[0]=='-' && arg[1]!='0' && arg[1]!='1' && arg[1]!='2' && arg[1]!='3' && arg[1]!='4' && arg[1]!='5' && arg[1]!='6' && arg[1]!='7' && arg[1]!='8' && arg[1]!='9') arg_is_option = true;
    	bool arg1_is_option = false;
    	if (arg1.size()>=2 && arg1[0]=='-' && arg1[1]!='0' && arg1[1]!='1' && arg1[1]!='2' && arg1[1]!='3' && arg1[1]!='4' && arg1[1]!='5' && arg1[1]!='6' && arg1[1]!='7' && arg1[1]!='8' && arg1[1]!='9') arg1_is_option = true;
    	

			//with multpile argument
			if(options_with_multiple_argument.has(arg))
			{
				//next argument is an option
				if(arg1_is_option)
				{
					root_.insert(ParamEntry("",StringList(),""),options_with_multiple_argument.find(arg)->second);
				}		
				//next argument is not an option	
				else
				{
					StringList sl;
					int j=(i+1);
					while(j< argc && !(arg1.size()>=2 && arg1[0]=='-' && arg1[1]!='0' && arg1[1]!='1' && arg1[1]!='2' && arg1[1]!='3' && arg1[1]!='4' && arg1[1]!='5' && arg1[1]!='6' && arg1[1]!='7' && arg1[1]!='8' && arg1[1]!='9'))
					{
						sl << arg1;
						++j;
						if (j< argc) arg1 = argv[j];
					}
					
					root_.insert(ParamEntry("",sl,""),options_with_multiple_argument.find(arg)->second);
					i = j-1;
				}
			}
			//without argument
			else if (options_without_argument.has(arg))
			{
				root_.insert(ParamEntry("",String("true"),""),options_without_argument.find(arg)->second);
			}
			//with one argument
			else if (options_with_one_argument.has(arg))

			{
				//next argument is not an option
				if (!arg1_is_option)
				{
					root_.insert(ParamEntry("",arg1,""),options_with_one_argument.find(arg)->second);
					++i;
				}
				//next argument is an option
				else
				{

					root_.insert(ParamEntry("",String(),""),options_with_one_argument.find(arg)->second);
				}
			}
			//unknown option
			else if (arg_is_option)
			{
      	ParamEntry* unknown_entry = root_.findEntryRecursive(unknown);
      	if (unknown_entry==0)
      	{
      		StringList sl;
      		sl << arg;
      		root_.insert(ParamEntry("",sl,""),unknown);
      	}
      	else
      	{
      		StringList sl = (StringList)unknown_entry->value;
      		sl << arg;
      		unknown_entry->value = sl;
      	}		
			}
			//just text argument
			else
			{
      	ParamEntry* misc_entry = root_.findEntryRecursive(misc);
      	if (misc_entry==0)
      	{
      		StringList sl;
      		sl << arg;
					// create "misc"-Node: 
      		root_.insert(ParamEntry("",sl,""),misc);
      	}
      	else
      	{
      		StringList sl = (StringList)misc_entry->value;
      		sl << arg;
      		misc_entry->value = sl;
      	}				
			}
    }
	}


	void Param::parseCommandLine(const int argc, const char** argv, const vector<ParameterInformation>& parameters, const String& misc, const String& unknown)
	{
		// prepare map of parameters:
		typedef map<String, vector<ParameterInformation>::const_iterator> ParamMap;
		ParamMap param_map;
		for (vector<ParameterInformation>::const_iterator it = parameters.begin(); it != parameters.end(); ++it)
		{
			param_map["-" + it->name] = it;
		}

		// list to store "misc"/"unknown" items:
		map<String, StringList> misc_unknown;

		list<String> queue; // queue for arguments
		// we parse the arguments in reverse order, so that we have arguments already when we encounter the option that uses them!
		for (int i = argc - 1; i > 0; --i)
		{
			String arg = argv[i];
			// options start with "-" or "--" followed by a letter:
			bool is_option = (arg.size() >= 2) && (arg[0] == '-') && (isalpha(arg[1]) || ((arg[1] == '-') && (arg.size() >= 3) &&  isalpha(arg[2])));
			if (is_option) // process content of the queue
			{
				ParamMap::iterator pos = param_map.find(arg);
				if (pos != param_map.end()) // parameter is defined
				{
					DataValue value;
					if (pos->second->type == ParameterInformation::FLAG) // flag
					{
						value = String("true");
					}
					else // option with argument(s)
					{
						switch (pos->second->default_value.valueType())
						{
						  case DataValue::STRING_VALUE:
								if (queue.empty()) value = String();
								else value = queue.front();
								break;
						  case DataValue::INT_VALUE:
								if (!queue.empty()) value = queue.front().toInt();
								break;
						  case DataValue::DOUBLE_VALUE:
								if (!queue.empty()) value = queue.front().toDouble();
								break;
						  case DataValue::STRING_LIST:
							{
								vector<String> arg_list(queue.begin(), queue.end());
								value = StringList(arg_list);
								queue.clear();
								break;
							}
						  case DataValue::INT_LIST:
							{
								IntList arg_list;
								for (list<String>::iterator it = queue.begin(); it != queue.end(); ++it)
								{
									arg_list << it->toInt();
								}
								value = arg_list;
								queue.clear();
								break;
							}
						  case DataValue::DOUBLE_LIST:
							{
								DoubleList arg_list;
								for (list<String>::iterator it = queue.begin(); it != queue.end(); ++it)
								{
									arg_list << it->toDouble();
								}
								value = arg_list;
								queue.clear();
								break;
							}
						  case DataValue::EMPTY_VALUE:
								break;
						}
						if (!queue.empty()) queue.pop_front(); // argument was already used
					}
					root_.insert(ParamEntry("", value, ""), pos->second->name);
				}
				else // unknown argument -> append to "unknown" list
				{
					misc_unknown[unknown] << arg;
				}
				// rest of the queue is just text -> insert into "misc" list:
				StringList& misc_list = misc_unknown[misc];
				misc_list.insert(misc_list.begin(), queue.begin(), queue.end());
				queue.clear();
			}
			else // more arguments
			{
				queue.push_front(arg); // order in the queue is not reversed!
			}
		}
		// remaining items in the queue are leading text arguments:
		StringList& misc_list = misc_unknown[misc];
		misc_list.insert(misc_list.begin(), queue.begin(), queue.end());
		
		// store "misc"/"unknown" items, if there were any:
		for (map<String, StringList>::iterator it = misc_unknown.begin();
				 it != misc_unknown.end(); ++it)
		{	
			if (it->second.empty()) continue;
			ParamEntry* entry = root_.findEntryRecursive(it->first);
			if (entry == 0) // create new node
			{
				root_.insert(ParamEntry("", it->second, ""), it->first);
			}
			else
			{
				StringList items = entry->value;
				items.insert(items.end(), it->second.begin(), it->second.end());
				entry->value = items;
			}
		}
	}


	ostream& operator << (ostream& os, const Param& param)
 	{
		for (Param::ParamIterator it = param.begin(); it!=param.end(); ++it)
		{
			String prefix = it.getName().chop(it->name.size()+1);
			if (prefix!="")
			{
				prefix += "|";
			}
			os << '"' << prefix << it->name << "\" -> \"" << it->value << '"';
			if (it->description!="")
			{
				os << " (" << it->description << ")";
			}
			os << endl;
		}
		return os;
	}

	Size Param::size() const
	{
		return root_.size();
	}

	bool Param::empty() const
	{
		return size()==0;
	}

	void Param::clear()
	{
		root_ = ParamNode("ROOT","");
	}

	void Param::checkDefaults(const String& name, const Param& defaults, const String& prefix, std::ostream& os) const
	{
		//Extract right parameters
		String prefix2 = prefix;
		if (prefix2!="")
		{
			prefix2.ensureLastChar(':');
		}
		Param check_values = copy(prefix2,true);
		
		//check
		for(ParamIterator it = check_values.begin(); it != check_values.end();++it)
		{
			//unknown parameter
			if (!defaults.exists(it.getName()))
			{
				os << "Warning: " << name << " received the unknown parameter '" << it.getName() << "'";
				if (!prefix2.empty()) os << " in '" << prefix2 << "'";
				os << "!" << endl;
			}

			//different types
			ParamEntry* default_value = defaults.root_.findEntryRecursive(prefix2+it.getName());
			if (default_value==0) continue;
			if (default_value->value.valueType()!=it->value.valueType())
			{
				String d_type;
				if (default_value->value.valueType()==DataValue::STRING_VALUE) d_type = "string";
        if (default_value->value.valueType()==DataValue::STRING_LIST) d_type = "string list";
				if (default_value->value.valueType()==DataValue::EMPTY_VALUE) d_type = "empty";
				if (default_value->value.valueType()==DataValue::INT_VALUE) d_type = "integer";
        if (default_value->value.valueType()==DataValue::INT_LIST) d_type = "integer list";
				if (default_value->value.valueType()==DataValue::DOUBLE_VALUE) d_type = "float";
        if (default_value->value.valueType()==DataValue::DOUBLE_LIST) d_type = "float list";
				String p_type;
				if (it->value.valueType()==DataValue::STRING_VALUE) p_type = "string";
        if (it->value.valueType()==DataValue::STRING_LIST) p_type = "string list";
				if (it->value.valueType()==DataValue::EMPTY_VALUE) p_type = "empty";
				if (it->value.valueType()==DataValue::INT_VALUE) p_type = "integer";
        if (it->value.valueType()==DataValue::INT_LIST) p_type = "integer list";
				if (it->value.valueType()==DataValue::DOUBLE_VALUE) p_type = "float";
        if (it->value.valueType()==DataValue::DOUBLE_LIST) p_type = "float list";

				throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,name+": Wrong parameter type '"+p_type+"' for "+d_type+" parameter '"+it.getName()+"' given!");
			}
			//parameter restrictions
			ParamEntry pe= *default_value;
			pe.value = it->value;
			String s;
			if (!	pe.isValid(s))	throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,name+": "+s);
		}
	}

  void Param::update(const Param& old_version, const bool report_new_params, const bool only_update_old, Logger::LogStream& stream)
  {
		// augment
		for(Param::ParamIterator it = old_version.begin(); it != old_version.end();++it)
		{
			if (this->exists(it.getName()))
			{
				// param 'version': do not override!
				if (it.getName().hasSuffix(":version"))
				{	
					if (this->getValue(it.getName()) != it->value)
					{
						stream << "Warning: for ':version' entry, augmented and Default Ini-File differ in value. Default value will not be altered!\n";
					}
					continue;
				}
				// param 'type': do not override!
        else if (it.getName().hasSuffix(":type") && 
                 it.getName().toQString().count(':')==2) // only for TOPP type (e.g. PeakPicker:1:type), any other 'type' param is ok
				{	
					if (this->getValue(it.getName()) != it->value)
					{
						stream << "Warning: for ':type' entry, augmented and Default Ini-File differ in value. Default value will not be altered!\n";
					}
					continue;
				}
				
				// all other parameters:
				Param::ParamEntry entry = this->getEntry (it.getName());
				if (entry.value.valueType() == it->value.valueType())
				{
					if (entry.value != it->value)
					{
						// check entry for consistency (in case restrictions have changed)							
						entry.value = it->value;
						String s;
						if (entry.isValid(s))
						{
							// overwrite default value
							stream << "Overriding Default-Parameter '" << it.getName() << "' with new value " << it->value << "\n"; 
							this->setValue(it.getName(),it->value, entry.description, this->getTags(it.getName()));
						}
						else
						{
							stream << "Parameter '" << it.getName() << "' does not fit into new restriction settings! Ignoring..."; 
						}
					}
					else
					{
						// value stayed the same .. nothing to be done
					}
				}
				else
				{
					stream << "Parameter '" << it.getName() << "' has changed value type! Ignoring...\n"; 
				}
			}
			else
			{
        if (!only_update_old)
        {
				  stream << "Deprecated Parameter '" << it.getName() << "' given in old parameter file! Ignoring...\n"; 
        }
        else
        {
          stream << "Deprecated Parameter '" << it.getName() << "' given in old parameter file! But leaving it there..." << std::endl; 
          Param::ParamEntry entry = old_version.getEntry (it.getName());
          String prefix = "";
          if (it.getName().has(':')) prefix = it.getName().substr(0, 1+it.getName().find_last_of(':'));
          this->root_.insert(entry, prefix);//->setValue(it.getName(), entry.value, entry.description, entry.tags);
        }
			}
		}

    // print new parameters (unique to this Param, but not in old one)
    // list new parameters not known to old version (just nice to know)
	  for(Param::ParamIterator it = this->begin(); it != this->end();/* do nothing */)
	  {
		  if (!old_version.exists(it.getName()))
      {
        if (report_new_params)
        {
          stream << "Information: New Parameter '" << it.getName() << "' not contained in old parameter file.\n"; 
        }
        if (only_update_old) 
        {
          this->removeAll(it.getName());
          it = this->begin(); // reset it, as it just became invalid
          continue;
        }
      }
      ++it;
    }

  }

	void Param::setSectionDescription(const String& key, const String& description)
	{
		ParamNode* node = root_.findParentOf(key);
		if (node==0)
		{
			throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,key);
		}
		
		Param::ParamNode::NodeIterator it = node->findNode(node->suffix(key));
		if (it==node->nodes.end())
		{
			throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,key);
		}
		it->description = description;
	}

	Param::ParamIterator Param::begin() const
	{
		return ParamIterator(root_);
	}

	Param::ParamIterator Param::end() const
	{
		return ParamIterator();
	}
		
	Param::ParamIterator::ParamIterator()
		: root_(0),
			current_(0),
      stack_(),
      trace_()
	{	
	}
	
	Param::ParamIterator::ParamIterator(const Param::ParamNode& root)
		: root_(&root),
			current_(-1),
      stack_(),
      trace_()
	{
		//Empty Param => begin == end iterator
		if (root_->entries.empty() && root_->nodes.empty())
		{
			root_ = 0;
			return;
		}
		
		//find first entry
		stack_.push_back(root_);
		operator++();
	}
	
	Param::ParamIterator::~ParamIterator()
	{
	}
	
	const Param::ParamEntry& Param::ParamIterator::operator*()
  {
  	return stack_.back()->entries[current_];
  }
  
  const Param::ParamEntry* Param::ParamIterator::operator->()
  {
  	return &(stack_.back()->entries[current_]);
  }

  Param::ParamIterator Param::ParamIterator::operator++(Int)
  {
		ParamIterator tmp(*this);
		++(*this);
		return tmp;
  }
  
  Param::ParamIterator& Param::ParamIterator::operator++()
  {
  	if (root_==0) return *this;
  	
  	trace_.clear();
  	while(true)
  	{
    	const Param::ParamNode* node = stack_.back();
    		
    	//cout << "############ operator++ #### " << node->name << " ## " << current_ <<endl;
    		
    	//check if there is a next entry in the current node
    	if (current_+1<(Int)node->entries.size())
    	{
    		//cout << " - next entry" <<endl;
    		++current_;
    		return *this;
    	}
    	//visit subnodes after entries
    	else if (!node->nodes.empty())
    	{
    		current_ = -1;
    		stack_.push_back(&(node->nodes[0]));
    		//cout << " - entering into: " << node->nodes[0].name <<endl;
    		//track changes (enter a node)
    		trace_.push_back(TraceInfo(node->nodes[0].name,node->nodes[0].description,true));
    		
    		continue;
    	}
    	//go back in tree until the node we came from is not the last subnode
    	//of the current node. Go into the next subnode.
    	else 
    	{
    		while(true)
	  		{
	    		const Param::ParamNode* last = node;
    			stack_.pop_back();
    		  //cout << " - stack size: " << stack_.size() << endl;
    		  //we have reached the end
		  		if (stack_.empty())
		  		{
		  			//cout << " - reached the end" << endl;
		  			root_ = 0;
		  			return *this;
		  		}
	    		node = stack_.back();
			
					//cout << " - last was: " << last->name << endl;
					//cout << " - descended to: " << node->name << endl;
							
					//track changes (leave a node)
					trace_.push_back(TraceInfo(last->name,last->description,false));
					
					//check of new subtree is accessable
					UInt next_index = (last - &(node->nodes[0])) +1;
					if (next_index < node->nodes.size()) 
					{
						current_=-1;
						stack_.push_back(&(node->nodes[next_index]));
						//cout << " - entering into: " << node->nodes[next_index].name  << endl;
						//track changes (enter a node)
						trace_.push_back(TraceInfo(node->nodes[next_index].name,node->nodes[next_index].description,true));
						break;
					}
	  		}
	  	}
  	}

  	return *this;
  }
	
	bool Param::ParamIterator::operator==(const ParamIterator& rhs) const
  {
  	return (root_==0 && rhs.root_==0) || (stack_==rhs.stack_ && current_==rhs.current_);
  }
	
	bool Param::ParamIterator::operator!=(const ParamIterator& rhs) const
  {
  	return !operator==(rhs);
  }
  
  String Param::ParamIterator::getName() const
  {
  	String tmp;
  	for (std::vector<const Param::ParamNode*>::const_iterator it = stack_.begin()+1; it!=stack_.end();++it)
  	{
  		tmp += (*it)->name + ':';
  	}
  	return tmp + stack_.back()->entries[current_].name;
  }
	
	const std::vector< Param::ParamIterator::TraceInfo >& Param::ParamIterator::getTrace() const
  {
  	return trace_;
  }

	const Param::ParamEntry& Param::getEntry(const String& key) const
	{
		return getEntry_(key);
	}
	
	const String& Param::getDescription(const String& key) const
	{
		return getEntry_(key).description;
	}

	void Param::addTag(const String& key, const String& tag)
	{
		if (tag.has(','))
		{
			throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Param tags may not contain comma characters",tag); 
		}
		getEntry_(key).tags.insert(tag);
	}
	
	void Param::addTags(const String& key, const StringList& tags)
	{
		ParamEntry& entry = getEntry_(key);
		for (Size i=0; i!=tags.size(); ++i)
		{
			if (tags[i].has(','))
			{
				throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Param tags may not contain comma characters",tags[i]); 
			}
			entry.tags.insert(tags[i]);
		}
	}

	StringList Param::getTags(const String& key) const
	{
		ParamEntry& entry = getEntry_(key);
		StringList list;
		for(set<String>::const_iterator it=entry.tags.begin(); it!=entry.tags.end(); ++it)
		{
			list.push_back(*it);
		}
		return list;
	}
	
	void Param::clearTags(const String& key)
	{
		getEntry_(key).tags.clear();
	}

	bool Param::hasTag(const String& key, const String& tag) const
	{
		return getEntry_(key).tags.count(tag);
	}
	
	bool Param::exists(const String& key) const
	{
		return root_.findEntryRecursive(key);
	}

	Param::ParamEntry& Param::getEntry_(const String& key) const
	{
		ParamEntry* entry = root_.findEntryRecursive(key);
		if (entry==0)
		{
			throw ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,key);
		}
		
		return *entry;
	}

}	//namespace
