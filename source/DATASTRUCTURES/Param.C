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

#include <OpenMS/DATASTRUCTURES/Param.h>

#include <iostream>
#include <fstream>

#include <OpenMS/FORMAT/HANDLERS/ParamXMLHandler.h>

using namespace std;
using namespace OpenMS::Exception;

namespace OpenMS
{

	//********************************* ParamEntry **************************************
	Param::ParamEntry::ParamEntry()
		: name(),
			description(),
			value(),
			advanced(true)
	{
	}

	Param::ParamEntry::ParamEntry(const String& n, const DataValue& v, const String& d, bool a)
		: name(n),
			description(d),
			value(v),
			advanced(a)
	{
		if (name.has(':'))
		{
			cerr << "Error ParamEntry name must not contain ':' characters!" << endl;
		}
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

	bool Param::ParamNode::operator==(const ParamNode& rhs) const
	{
		if (name!=rhs.name || entries.size()!=rhs.entries.size() || nodes.size()!=rhs.nodes.size()) return false;
		
		//order of sections / entries should not matter
		for(UInt i=0; i< entries.size(); ++i)
		{
			if (find(rhs.entries.begin(),rhs.entries.end(),entries[i])==rhs.entries.end())
			{
				return false;
			}
		}
		for(UInt i=0; i< nodes.size(); ++i)
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
			for(UInt i=0; i<nodes.size();++i)
			{
				if (nodes[i].name.hasPrefix(name)) return this;
			}
			for(UInt i=0; i<entries.size();++i)
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
			it->nodes.insert(it->nodes.end(), node.nodes.begin(), node.nodes.end());
			it->entries.insert(it->entries.end(), node.entries.begin(), node.entries.end());
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
		//cerr << "ParamEntry::Insert(ParamEntry)" << endl;
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
			it->advanced = entry.advanced;
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
	
	UInt Param::ParamNode::size() const
	{
		UInt subnode_size = 0;
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
		: XMLFile(OPENMS_PATH"/data/SCHEMAS/Param_1_0.xsd"),
			root_("ROOT",""), 
			inheritance_steps_max(15)
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

	Param& Param::operator = (const Param& rhs)
	{
		root_ = rhs.root_;
		return *this;
	}

	Param::Param(const ParamNode& node)
		: root_(node), 
			inheritance_steps_max(15)
	{
		root_.name="ROOT";
		root_.description="";
	}

	bool Param::operator == (const Param& rhs) const
	{
		return root_ == rhs.root_;
	}	
	
	void Param::setValue(const String& key, Int value, const String& description, bool advanced)
	{
		root_.insert(ParamEntry("",DataValue(value),description,advanced),key);
	}
	
	void Param::setValue(const String& key, float value, const String& description, bool advanced)
	{
		root_.insert(ParamEntry("",DataValue(value),description,advanced),key);
	}

	void Param::setValue(const String& key, double value, const String& description, bool advanced)
	{
		root_.insert(ParamEntry("",DataValue(value),description,advanced),key);
	}

	void Param::setValue(const String& key, const String& value, const String& description, bool advanced)
	{
		root_.insert(ParamEntry("",DataValue(value),description,advanced),key);
	}
	
	const DataValue& Param::getValue(const String& key) const throw (ElementNotFound<String>)
	{
		ParamEntry* entry = root_.findEntryRecursive(key);
		if (entry==0)
		{
			throw ElementNotFound<String>(__FILE__,__LINE__,__PRETTY_FUNCTION__,key);
		}
		
		return entry->value;
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
	
	void Param::insert(String prefix, const Param& param)
	{
		for (Param::ParamNode::NodeIterator it=param.root_.nodes.begin(); it!=param.root_.nodes.end(); ++it)
		{
			root_.insert(*it,prefix);
		}
		for (Param::ParamNode::EntryIterator it=param.root_.entries.begin(); it!=param.root_.entries.end(); ++it)
		{
			root_.insert(*it,prefix);
		}
	}

	void Param::setDefaults(const Param& defaults, String prefix, bool showMessage)
	{
		if ( !prefix.empty() )
		{
			prefix.ensureLastChar(':');
		}	
		String pathname;
		for(Param::ParamIterator it = defaults.begin(); it != defaults.end();++it)
		{
			if (!exists(prefix + it.getName()))
			{
				if (showMessage) cerr << "Setting " << prefix+it.getName() << " to " << it->value << endl;
				root_.insert(ParamEntry("", it->value, it->description, it->advanced), prefix+it.getName());
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
				String real_pathname = pathname.substr(0,-1); //remove ':' at the end
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
						setSectionDescription(prefix+real_pathname, description_new);
					}
				}
			}
		}
	}
	
	void Param::remove(const String& prefix)
	{
		if (prefix.hasSuffix(':')) //we have to delete one node only
		{
			ParamNode* node = root_.findParentOf(prefix.substr(0,-1));
			if (node!=0)
			{
				Param::ParamNode::NodeIterator it = node->findNode(node->suffix(prefix.substr(0,-1)));
				if (it!=node->nodes.end())
				{
					node->nodes.erase(it);
				}
			}
		}
		else //we have to delete all entries and nodes starting with the suffix
		{
			ParamNode* node = root_.findParentOf(prefix);
			if (node!=0)
			{			
				String suffix = node->suffix(prefix);
				
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
				out.insert(*node,prefix.substr(0,-(node->name.size()+1)));
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
						out.insert(*it,prefix.substr(0,-suffix.size()));
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
						out.insert(*it,prefix.substr(0,-suffix.size()));
					}
				}
			}
		}

		return Param(out);
	}

	Param Param::copyWithInherit(const String& prefix) const
	{
		return copy(prefix,true);
//		if ( !prefix.hasSuffix(':') )
//		{
//			return copy(prefix, true);
//		}
//		else
//		{
//			Param result = copy(prefix, true);
//			Int inheritance_steps = 0;
//			const DataValue* inherit_path_value = &DataValue::EMPTY;
//			if (result.exists("inherit"))
//			{
//				inherit_path_value = &(result.getValue("inherit")); 
//			}
//			while (!inherit_path_value->isEmpty())
//			{
//				String inherit_path = inherit_path_value->toString();
//				if ( ++inheritance_steps > inheritance_steps_max )
//				{
//					throw Exception::ParseError ( __FILE__, __LINE__, __PRETTY_FUNCTION__, "<ITEM name=\"inherit\" ... />", String("Too many inheritance steps (")+inheritance_steps_max+" allowed).  Perhaps there is a cycle? prefix="+prefix+" inherit_path="+inherit_path);
//				}
//				result.remove("inherit");
//				result.setDefaults(copy(inherit_path+':',true), "", false);
//				
//				const DataValue* inherit_path_value = &DataValue::EMPTY;
//				if (result.exists("inherit"))
//				{
//					inherit_path_value = &(result.getValue("inherit")); 
//				}
//			}
//			return result.copy("",true);
//		}
	}


	void Param::store(const String& filename) const throw (Exception::UnableToCreateFile)
	{
		//open file
		ofstream os;
		os.open (filename.c_str(), ofstream::out);
		if(!os)
		{
			 throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}
		
  	os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
  	os << "<PARAMETERS  xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/Param_1_0.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";
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
					d.substitute('"','\'');
					d.substitute("\n","#br#");
					d.substitute("<","&lt;");
					d.substitute(">","&gt;");
					os << indentation  << "<NODE name=\"" << it2->name << "\" description=\"" << d << "\">" << endl;
					indentation += "  ";
				}
				else //closed node
				{
					indentation.resize(indentation.size()-2);
					os << indentation << "</NODE>" << endl;
				}
			}
			
			//write item
			if(it->value.valueType()!=DataValue::EMPTYVALUE)
			{				
				os << indentation << "<ITEM name=\"" << it->name << "\" value=\"" << it->value.toString() << "\" type=\"";
				
				switch(it->value.valueType())
				{
					case DataValue::INTVALUE:
					case DataValue::LONVALUE:
					case DataValue::SHOVALUE:
						os << "int";
						break;
					case DataValue::FLOVALUE:
					case DataValue::DOUVALUE:
						os << "float";
						break;
					case DataValue::STRVALUE:
						os << "string";
						break;
					default:
						break;
				};
				//replace all critical characters in description
				String d = it->description;
				d.substitute("\"","'");
				d.substitute("\n","#br#");
				d.substitute("<","&lt;");
				d.substitute(">","&gt;");
				os << "\" description=\"" << d << "\"";
				//advanced parameter handling
				if (it->advanced)
				{
					os << " advanced=\"true\"";
				}
				else
				{
					os << " advanced=\"false\"";
				}
				os << " />" <<  endl;	
			}
			++it;
		}
		
		//close remaining tags
		const std::vector< ParamIterator::TraceInfo >& trace = it.getTrace();
		for(std::vector< ParamIterator::TraceInfo >::const_iterator it2 = trace.begin(); it2!=trace.end(); ++it2)
		{
			indentation.resize(indentation.size()-2);
			os << indentation << "</NODE>" << endl;	
		}
		
		os << "</PARAMETERS>\n";
    os.close();
	}
	
	void Param::load(const String& filename) throw (FileNotFound,ParseError)
	{
		Internal::ParamXMLHandler handler(*this, filename);
		parse_(filename, &handler);
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
    	
    	//it is a option when it starts with a '-' and the second character is not a number
    	bool arg1_is_option = false;
    	bool arg_is_option = false;
    	if (arg1[0]=='-' && arg1[1]!='0' && arg1[1]!='1' && arg1[1]!='2' && arg1[1]!='3' && arg1[1]!='4' && arg1[1]!='5' && arg1[1]!='6' && arg1[1]!='7' && arg1[1]!='8' && arg1[1]!='9') arg1_is_option = true;
    	if (arg[0]=='-' && arg[1]!='0' && arg[1]!='1' && arg[1]!='2' && arg[1]!='3' && arg[1]!='4' && arg[1]!='5' && arg[1]!='6' && arg[1]!='7' && arg[1]!='8' && arg[1]!='9') arg_is_option = true;
    	
    	//cout << "Parse: '"<< arg << "' '" << arg1 << "'" << endl;
    	
      //flag (option without text argument)
      if(arg_is_option && arg1_is_option)
      {
	    	root_.insert(ParamEntry(arg,"","",false),prefix);
      }
      //option with argument
      else if(arg_is_option && !arg1_is_option)
      {
      	root_.insert(ParamEntry(arg,arg1,"",false),prefix);
      	++i;
      }      
      //just text arguments (not preceded by an option)
      else
      {
      	ParamEntry* misc_entry = root_.findEntryRecursive(prefix+"misc");
      	if (misc_entry==0)
      	{
      		root_.insert(ParamEntry("misc",arg,"",false),prefix);
      	}
      	else
      	{
      		misc_entry->value = String(misc_entry->value)+" "+arg;
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

    	//it is a option when it starts with a '-' and the second character is not a number
    	bool arg1_is_option = false;
    	bool arg_is_option = false;
    	if (arg1[0]=='-' && arg1[1]!='0' && arg1[1]!='1' && arg1[1]!='2' && arg1[1]!='3' && arg1[1]!='4' && arg1[1]!='5' && arg1[1]!='6' && arg1[1]!='7' && arg1[1]!='8' && arg1[1]!='9') arg1_is_option = true;
    	if (arg[0]=='-' && arg[1]!='0' && arg[1]!='1' && arg[1]!='2' && arg[1]!='3' && arg[1]!='4' && arg[1]!='5' && arg[1]!='6' && arg[1]!='7' && arg[1]!='8' && arg[1]!='9') arg_is_option = true;
    	

			//without argument
			if (options_without_argument.find(arg)!=options_without_argument.end())
			{
				root_.insert(ParamEntry("","true","",false),options_without_argument.find(arg)->second);
			}
			//with argument
			else if (options_with_argument.find(arg)!=options_with_argument.end())
			{
				//next argument is not a option
				if (!arg1_is_option)
				{
					root_.insert(ParamEntry("",arg1,"",false),options_with_argument.find(arg)->second);
					++i;
				}
				//next argument is a option
				else
				{
					root_.insert(ParamEntry("","","",false),options_with_argument.find(arg)->second);
				}
			}
			//unknown option
			else if (arg_is_option)
			{
      	ParamEntry* unknown_entry = root_.findEntryRecursive(unknown);
      	if (unknown_entry==0)
      	{
      		root_.insert(ParamEntry("",arg,"",false),unknown);
      	}
      	else
      	{
      		unknown_entry->name = unknown_entry->value+" "+arg;
      	}		
			}
			//just text argument
			else
			{
      	ParamEntry* misc_entry = root_.findEntryRecursive(misc);
      	if (misc_entry==0)
      	{
      		root_.insert(ParamEntry("",arg,"",false),misc);
      	}
      	else
      	{
      		misc_entry->value = String(misc_entry->value)+" "+arg;
      	}				
			}
    }
	}

	ostream& operator << (ostream& os, const Param& param)
 	{
		for (Param::ParamIterator it = param.begin(); it!=param.end(); ++it)
		{
			String prefix = it.getName().substr(0,-(it->name.size()+1));
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

	UInt Param::size() const
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

	void Param::checkDefaults(const String& name, const Param& defaults, String prefix, std::ostream& os) const
	{
		//Extract right parameters
		if ( !prefix.empty() )
		{
			prefix.ensureLastChar(':');
		}	
		Param check_values = copy(prefix,true);
		
		//check
		for(ParamIterator it = check_values.begin(); it != check_values.end();++it)
		{
			if (!defaults.exists(it.getName()))
			{
				os << "Warning: " << name << " received the unknown parameter '" << it.getName() << "'";
				if (!prefix.empty()) os << " in '" << prefix << "'";
				os << "!" << endl;
			}
		}
	}

	void Param::setSectionDescription(const String& key, const String& description)  throw (ElementNotFound<String>)
	{
		ParamNode* node = root_.findParentOf(key);
		if (node==0)
		{
			throw ElementNotFound<String>(__FILE__,__LINE__,__PRETTY_FUNCTION__,key);
		}
		
		Param::ParamNode::NodeIterator it = node->findNode(node->suffix(key));
		if (it==node->nodes.end())
		{
			throw ElementNotFound<String>(__FILE__,__LINE__,__PRETTY_FUNCTION__,key);
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
			current_(0)
	{	
	}
	
	Param::ParamIterator::ParamIterator(const Param::ParamNode& root)
		: root_(&root),
			current_(-1)
	{
		//Empty Param => begin == end iterator
		if (root_->entries.size()==0 && root_->nodes.size()==0)
		{
			root_ = 0;
			return;
		}
		
		//find first entry
		stack_.push_back(root_);
		operator++();
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
    	if (current_+1<node->entries.size())
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

	const Param::ParamEntry& Param::getEntry(const String& key) const throw (ElementNotFound<String>)
	{
		ParamEntry* entry = root_.findEntryRecursive(key);
		if (entry==0)
		{
			throw ElementNotFound<String>(__FILE__,__LINE__,__PRETTY_FUNCTION__,key);
		}
		
		return *entry;
	}
	
	const String& Param::getDescription(const String& key) const throw (ElementNotFound<String>)
	{
		ParamEntry* entry = root_.findEntryRecursive(key);
		if (entry==0)
		{
			throw ElementNotFound<String>(__FILE__,__LINE__,__PRETTY_FUNCTION__,key);
		}
		
		return entry->description;
	}
	
	bool Param::isAdvancedParameter(const String& key) const throw (ElementNotFound<String>)
	{
		ParamEntry* entry = root_.findEntryRecursive(key);
		if (entry==0)
		{
			throw ElementNotFound<String>(__FILE__,__LINE__,__PRETTY_FUNCTION__,key);
		}
		
		return entry->advanced;
	}
	
	bool Param::exists(const String& key) const
	{
		return root_.findEntryRecursive(key);
	}

}	//namespace
