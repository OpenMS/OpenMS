// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/Param.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <algorithm>

namespace OpenMS
{

  //********************************* ParamEntry **************************************
  Param::ParamEntry::ParamEntry() :
    name(),
    description(),
    value(),
    tags(),
    min_float(-std::numeric_limits<double>::max()),
    max_float(std::numeric_limits<double>::max()),
    min_int(-std::numeric_limits<int>::max()),
    max_int(std::numeric_limits<int>::max()),
    valid_strings()
  {
  }

  Param::ParamEntry::ParamEntry(const std::string& n, const ParamValue& v, const std::string& d, const std::vector<std::string>& t) :
    name(n),
    description(d),
    value(v),
    tags(),
    min_float(-std::numeric_limits<double>::max()),
    max_float(std::numeric_limits<double>::max()),
    min_int(-std::numeric_limits<int>::max()),
    max_int(std::numeric_limits<int>::max()),
    valid_strings()
  {
    //add tags
    for (size_t i = 0; i < t.size(); ++i)
    {
      tags.insert(t[i]);
    }
    //check name
    if (name.find(':') != std::string::npos)
    {
      std::cerr << "Error ParamEntry name must not contain ':' characters!" << std::endl;
    }
  }

  Param::ParamEntry::~ParamEntry() = default;

  bool Param::ParamEntry::isValid(std::string& message) const
  {
    if (value.valueType() == ParamValue::STRING_VALUE)
    {
      if (!valid_strings.empty())
      {
        bool ok = false;
        if (std::find(valid_strings.begin(), valid_strings.end(), value) != valid_strings.end())
        {
          ok = true;
        }
        else if (std::find(tags.begin(), tags.end(), "input file") != tags.end() 
          || std::find(tags.begin(), tags.end(), "output file") != tags.end()
          || std::find(tags.begin(), tags.end(), "output prefix") != tags.end())
        {
          //do not check restrictions on file names for now
          ok = true;
        }

        if (!ok)
        {
          std::string valid = valid_strings.front();
          for (auto it = valid_strings.begin() + 1, end = valid_strings.end(); it != end; ++it) {
              valid += "," + *it;
          }
          message = "Invalid string parameter value '" + value.toString() + "' for parameter '" + name + "' given! Valid values are: '" + valid + "'.";
          return false;
        }
      }
    }
    else if (value.valueType() == ParamValue::STRING_LIST)
    {
      std::string str_value;
      std::vector<std::string> ls_value = value;
      for (size_t i = 0; i < ls_value.size(); ++i)
      {
        str_value = ls_value[i];

        if (!valid_strings.empty())
        {
          bool ok = false;
          if (std::find(valid_strings.begin(), valid_strings.end(), str_value) != valid_strings.end())
          {
            ok = true;
          }
          else if (std::find(tags.begin(), tags.end(), "input file") != tags.end() 
            || std::find(tags.begin(), tags.end(), "output file") != tags.end())
          {
            //do not check restrictions on file names for now
            ok = true;
          }

          if (!ok)
          {
            std::string valid = valid_strings.front();
            for (auto it = valid_strings.begin() + 1, end = valid_strings.end(); it != end; ++it) {
                valid += "," + *it;
            }
            message = "Invalid string parameter value '" + str_value + "' for parameter '" + name + "' given! Valid values are: '" + valid + "'.";
            return false;
          }
        }
      }
    }
    else if (value.valueType() == ParamValue::INT_VALUE)
    {
      int tmp = value;
      if ((min_int != -std::numeric_limits<int>::max() && tmp < min_int) || (max_int != std::numeric_limits<int>::max() && tmp > max_int))
      {
        message = "Invalid integer parameter value '" + std::to_string(tmp) + "' for parameter '" + name + "' given! The valid range is: [" + std::to_string(min_int) + ":" + std::to_string(max_int) + "].";
        return false;
      }
    }
    else if (value.valueType() == ParamValue::INT_LIST)
    {
      int int_value;
      std::vector<int> ls_value = value;
      for (size_t i = 0; i < ls_value.size(); ++i)
      {
        int_value = ls_value[i];
        if ((min_int != -std::numeric_limits<int>::max() && int_value < min_int) || (max_int != std::numeric_limits<int>::max() && int_value > max_int))
        {
          message = "Invalid integer parameter value '" + std::to_string(int_value) + "' for parameter '" + name + "' given! The valid range is: [" + std::to_string(min_int) + ":" + std::to_string(max_int) + "].";
          return false;
        }
      }
    }
    else if (value.valueType() == ParamValue::DOUBLE_VALUE)
    {
      double tmp = value;
      if ((min_float != -std::numeric_limits<double>::max() && tmp < min_float) || (max_float != std::numeric_limits<double>::max() && tmp > max_float))
      {
        message = "Invalid double parameter value '" + std::to_string(tmp) + "' for parameter '" + name + "' given! The valid range is: [" + std::to_string(min_float)+ ":" + std::to_string(max_float) + "].";
        return false;
      }
    }
    else if (value.valueType() == ParamValue::DOUBLE_LIST)
    {
      std::vector<double> ls_value = value;
      for (size_t i = 0; i < ls_value.size(); ++i)
      {
        double dou_value = ls_value[i];
        if ((min_float != -std::numeric_limits<double>::max() && dou_value < min_float) || (max_float != std::numeric_limits<double>::max() && dou_value > max_float))
        {
          message = "Invalid double parameter value '" + std::to_string(dou_value) + "' for parameter '" + name + "' given! The valid range is: [" + std::to_string(min_float) + ":" + std::to_string(max_float) + "].";
          return false;
        }
      }
    }
    return true;
  }

  bool Param::ParamEntry::operator==(const ParamEntry& rhs) const
  {
    return name == rhs.name && value == rhs.value;
  }

  //********************************* ParamNode **************************************
  Param::ParamNode::ParamNode() :
    name(),
    description(),
    entries(),
    nodes()
  {
  }

  Param::ParamNode::ParamNode(const std::string& n, const std::string& d) :
    name(n),
    description(d),
    entries(),
    nodes()
  {
      if (name.find(':') != std::string::npos) 
      {
        std::cerr << "Error ParamNode name must not contain ':' characters!" << std::endl;
      }
  }

  Param::ParamNode::~ParamNode() = default;

  bool Param::ParamNode::operator==(const ParamNode& rhs) const
  {
    if (name != rhs.name || entries.size() != rhs.entries.size() || nodes.size() != rhs.nodes.size())
    {
      return false;
    }
    //order of sections / entries should not matter
    for (size_t i = 0; i < entries.size(); ++i)
    {
      if (find(rhs.entries.begin(), rhs.entries.end(), entries[i]) == rhs.entries.end())
      {
        return false;
      }
    }
    for (size_t i = 0; i < nodes.size(); ++i)
    {
      if (find(rhs.nodes.begin(), rhs.nodes.end(), nodes[i]) == rhs.nodes.end())
      {
        return false;
      }
    }

    return true;
  }

  Param::ParamNode::EntryIterator Param::ParamNode::findEntry(const std::string& local_name)
  {
    for (EntryIterator it = entries.begin(); it != entries.end(); ++it)
    {
      if (it->name == local_name)
      {
        return it;
      }
    }
    return entries.end();
  }

  Param::ParamNode::NodeIterator Param::ParamNode::findNode(const std::string& local_name)
  {
    for (NodeIterator it = nodes.begin(); it != nodes.end(); ++it)
    {
      if (it->name == local_name)
      {
        return it;
      }
    }
    return nodes.end();
  }

  Param::ParamNode* Param::ParamNode::findParentOf(const std::string& local_name)
  {
    //cout << "findParentOf nodename: " << this->name << " - nodes: " << this->nodes.size() << " - find: "<< name << std::endl;
    if (local_name.find(':') != std::string::npos) //several subnodes to browse through
    {
        size_t pos = local_name.find(':');
        std::string prefix = local_name.substr(0, pos);

        //cout << " - Prefix: '" << prefix << "'" << std::endl;
        NodeIterator it = findNode(prefix);
        if (it == nodes.end()) //subnode not found
        {
          return nullptr;
        }
        //recursively call findNode for the rest of the path
        std::string new_name = local_name.substr(it->name.size() + 1);
        //cout << " - Next name: '" << new_name << "'" << std::endl;
        return it->findParentOf(new_name);
    }
    else // we are in the right child
    {
        //check if a node or entry prefix match
        for (size_t i = 0; i < nodes.size(); ++i)
        {
            if (nodes[i].name.compare(0, local_name.size(), local_name) == 0)
            {
              return this;
            }
        }
        for (size_t i = 0; i < entries.size(); ++i)
        {
            if (entries[i].name.compare(0, local_name.size(), local_name) == 0)
            {
            return this;
            }
        }
        return nullptr;
    }
  }

  Param::ParamEntry* Param::ParamNode::findEntryRecursive(const std::string& local_name)
  {
    ParamNode* parent = findParentOf(local_name);
    if (parent == nullptr)
    {
      return nullptr;
    }

    EntryIterator it = parent->findEntry(suffix(local_name));
    if (it == parent->entries.end())
    {
      return nullptr;
    }

    return &(*it);
  }

  void Param::ParamNode::insert(const ParamNode& node, const std::string& prefix)
  {
    //std::cerr << "INSERT NODE  " << node.name << " (" << prefix << ")" << std::endl;
    std::string prefix2 = prefix + node.name;

    ParamNode* insert_node = this;
    while (prefix2.find(':') != std::string::npos)
    {
      size_t pos = prefix2.find(':');
      std::string local_name = prefix2.substr(0, pos);
      //check if the node already exists
      NodeIterator it = insert_node->findNode(local_name);
      if (it != insert_node->nodes.end()) //exists
      {
        insert_node = &(*it);
      }
      else //create it
      {
        insert_node->nodes.emplace_back(local_name, "");
        insert_node = &(insert_node->nodes.back());
        //std::cerr << " - Created new node: " << insert_node->name << std::endl;
      }
      //remove prefix
      prefix2 = prefix2.substr(local_name.size() + 1);
    }

    // check if the node exists as ParamEntry
    EntryIterator entry_it = insert_node->findEntry(prefix2);
    if (entry_it != insert_node->entries.end()) {
      std::string message = "Duplicate option \""
                            + prefix
                            + "\" into \""
                            + name
                            + "\", should not be added as ParamNode and ParamEntry at the same time (1).";
      throw Exception::InternalToolError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
    }

    //check if the node already exists
    NodeIterator it = insert_node->findNode(prefix2);
    if (it != insert_node->nodes.end()) //append nodes and entries
    {
      for (ConstNodeIterator it2 = node.nodes.begin(); it2 != node.nodes.end(); ++it2)
      {
        it->insert(*it2);
      }
      for (ConstEntryIterator it2 = node.entries.begin(); it2 != node.entries.end(); ++it2)
      {
        it->insert(*it2);
      }
      if (it->description.empty() || !node.description.empty()) //replace description if not empty in new node
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

  void Param::ParamNode::insert(const ParamEntry& entry, const std::string& prefix)
  {
    //std::cerr << "INSERT ENTRY " << entry.name << " (" << prefix << ")" << std::endl;
    std::string prefix2 = prefix + entry.name;
    //std::cerr << " - inserting: " << prefix2 << std::endl;

    ParamNode* insert_node = this;
    while (prefix2.find(':') != std::string::npos)
    {
      size_t pos = prefix2.find(':');
      std::string local_name = prefix2.substr(0, pos);
      //std::cerr << " - looking for node: " << name << std::endl;
      //look up if the node already exists
      NodeIterator it = insert_node->findNode(local_name);
      if (it != insert_node->nodes.end()) //exists
      {
        insert_node = &(*it);
      }
      else //create it
      {
        insert_node->nodes.emplace_back(local_name, "");
        insert_node = &(insert_node->nodes.back());
        //std::cerr << " - Created new node: " << insert_node->name << std::endl;
      }
      //remove prefix
      prefix2 = prefix2.substr(local_name.size() + 1);
      //std::cerr << " - new prefix: " << prefix2 << std::endl;
    }

    // check if the entry exists as ParamNode
    NodeIterator node_it = insert_node->findNode(prefix2);
    if (node_it != insert_node->nodes.end())
    {
      std::string message = "Duplicate option \""
                            + prefix
                            + "\" into \""
                            + name
                            + "\", should not be added as ParamNode and ParamEntry at the same time (2).";
      throw Exception::InternalToolError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
    }

    //check if the entry already exists
    //std::cerr << " - final entry name: " << prefix2 << std::endl;
    EntryIterator it = insert_node->findEntry(prefix2);
    if (it != insert_node->entries.end()) //overwrite entry
    {
      it->value = entry.value;
      it->tags = entry.tags;
      if (it->description.empty() || !entry.description.empty()) //replace description if not empty in new entry
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

  size_t Param::ParamNode::size() const
  {
    size_t subnode_size = 0;
    for (std::vector<ParamNode>::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
      subnode_size += it->size();
    }
    return entries.size() + subnode_size;
  }

  std::string Param::ParamNode::suffix(const std::string& key) const
  {
    size_t pos = key.rfind(':');
    if (pos != std::string::npos)
    {
      return key.substr(++pos);
    }
    return key;
  }

  //********************************* Param **************************************

  Param::Param() :
    root_("ROOT", "")
  {
  }

  Param::~Param() = default;

  Param::Param(const ParamNode& node) :
    root_(node)
  {
    root_.name = "ROOT";
    root_.description = "";
  }

  bool Param::operator==(const Param& rhs) const
  {
    return root_ == rhs.root_;
  }

  void Param::setValue(const std::string& key, const ParamValue& value, const std::string& description, const std::vector<std::string>& tags)
  {
    root_.insert(ParamEntry("", value, description, tags), key);
  }

  void Param::setValidStrings(const std::string& key, const std::vector<std::string>& strings)
  {
    ParamEntry& entry = getEntry_(key);
    //check if correct parameter type
    if (entry.value.valueType() != ParamValue::STRING_VALUE && entry.value.valueType() != ParamValue::STRING_LIST)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, key);
    }
    //check for commas
    for (size_t i = 0; i < strings.size(); ++i)
    {
      if (strings[i].find(',') != std::string::npos)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Comma characters in Param string restrictions are not allowed!");
      }
    }
    entry.valid_strings = strings;
  }

  const std::vector<std::string>& Param::getValidStrings(const std::string& key) const
  {
    ParamEntry& entry = getEntry_(key);
    // check if correct parameter type
    if (entry.value.valueType() != ParamValue::STRING_VALUE && entry.value.valueType() != ParamValue::STRING_LIST)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, key);
    }
    return entry.valid_strings;
  }

  void Param::setMinInt(const std::string& key, int min)
  {
    ParamEntry& entry = getEntry_(key);
    if (entry.value.valueType() != ParamValue::INT_VALUE && entry.value.valueType() != ParamValue::INT_LIST)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, key);
    }
    entry.min_int = min;
  }

  void Param::setMaxInt(const std::string& key, int max)
  {
    ParamEntry& entry = getEntry_(key);
    if (entry.value.valueType() != ParamValue::INT_VALUE && entry.value.valueType() != ParamValue::INT_LIST)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, key);
    }
    entry.max_int = max;
  }

  void Param::setMinFloat(const std::string& key, double min)
  {
    ParamEntry& entry = getEntry_(key);
    if (entry.value.valueType() != ParamValue::DOUBLE_VALUE && entry.value.valueType() != ParamValue::DOUBLE_LIST)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, key);
    }
    entry.min_float = min;
  }

  void Param::setMaxFloat(const std::string& key, double max)
  {
    ParamEntry& entry = getEntry_(key);
    if (entry.value.valueType() != ParamValue::DOUBLE_VALUE && entry.value.valueType() != ParamValue::DOUBLE_LIST)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, key);
    }
    entry.max_float = max;
  }

  const ParamValue& Param::getValue(const std::string& key) const
  {
    return getEntry_(key).value;
  }

  const std::string& Param::getSectionDescription(const std::string& key) const
  {
    //This variable is used instead of String::EMPTY as the method is used in
    //static initialization and thus cannot rely on String::EMPTY been initialized.
    static std::string empty;

    ParamNode* node = root_.findParentOf(key);
    if (node == nullptr)
    {
      return empty;
    }

    Param::ParamNode::NodeIterator it = node->findNode(node->suffix(key));
    if (it == node->nodes.end())
    {
      return empty;
    }

    return it->description;
  }

  void Param::insert(const std::string& prefix, const Param& param)
  {
    //std::cerr << "INSERT PARAM (" << prefix << ")" << std::endl;
    for (Param::ParamNode::NodeIterator it = param.root_.nodes.begin(); it != param.root_.nodes.end(); ++it)
    {
      root_.insert(*it, prefix);
    }
    for (Param::ParamNode::EntryIterator it = param.root_.entries.begin(); it != param.root_.entries.end(); ++it)
    {
      root_.insert(*it, prefix);
    }
  }

  void Param::setDefaults(const Param& defaults, const std::string& prefix, bool showMessage)
  {
    std::string prefix2 = prefix;
    if (!prefix2.empty())
    {
      if (prefix2.back() != ':')
      {
        prefix2 += ':';
      }
    }

    std::string pathname;
    for (Param::ParamIterator it = defaults.begin(); it != defaults.end(); ++it)
    {
      if (!exists(prefix2 + it.getName()))
      {
        if (showMessage)
          std::cerr << "Setting " << prefix2 + it.getName() << " to " << it->value << std::endl;
        std::string name = prefix2 + it.getName();
        root_.insert(ParamEntry("", it->value, it->description), name);
        //copy tags
        for (std::set<std::string>::const_iterator tag_it = it->tags.begin(); tag_it != it->tags.end(); ++tag_it)
        {
          addTag(name, *tag_it);
        }
        //copy restrictions
        if (it->value.valueType() == ParamValue::STRING_VALUE || it->value.valueType() == ParamValue::STRING_LIST)
        {
          setValidStrings(name, it->valid_strings);
        }
        else if (it->value.valueType() == ParamValue::INT_VALUE || it->value.valueType() == ParamValue::INT_LIST)
        {
          setMinInt(name, it->min_int);
          setMaxInt(name, it->max_int);
        }
        else if (it->value.valueType() == ParamValue::DOUBLE_VALUE || it->value.valueType() == ParamValue::DOUBLE_LIST)
        {
          setMinFloat(name, it->min_float);
          setMaxFloat(name, it->max_float);
        }
      }

      //copy section descriptions
      const std::vector<ParamIterator::TraceInfo>& trace = it.getTrace();
      for (std::vector<ParamIterator::TraceInfo>::const_iterator it2 = trace.begin(); it2 != trace.end(); ++it2)
      {
        if (it2->opened)
        {
          pathname += it2->name + ":";
        }
        else
        {
          pathname.resize(pathname.size() - it2->name.size() - 1);
        }
        std::string real_pathname = pathname.substr(0, pathname.length() - 1); //remove ':' at the end
        if (!real_pathname.empty())
        {
          std::string description_old = getSectionDescription(prefix + real_pathname);
          const std::string& description_new = defaults.getSectionDescription(real_pathname);
          if (description_old.empty())
          {
            //std::cerr << "## Setting description of " << prefix+real_pathname << " to"<< std::endl;
            //std::cerr << "## " << description_new << std::endl;
            setSectionDescription(prefix2 + real_pathname, description_new);
          }
        }
      }
    }
  }

  void Param::remove(const std::string& key)
  {
    std::string keyname = key;
    if (!key.empty() && key.back() == ':') // delete section
    {
      keyname = key.substr(0, key.length() - 1);

      ParamNode* node_parent = root_.findParentOf(keyname);
      if (node_parent != nullptr)
      {
        Param::ParamNode::NodeIterator it = node_parent->findNode(node_parent->suffix(keyname));
        if (it != node_parent->nodes.end())
        {
          std::string name = it->name;
          node_parent->nodes.erase(it); // will automatically delete subnodes
          if (node_parent->nodes.empty()  && node_parent->entries.empty())
          {
            // delete last section name (could be partial)
            remove(keyname.substr(0, keyname.size() - name.size())); // keep last ':' to indicate deletion of a section
          }
        }
      }
    }
    else
    {
      ParamNode* node = root_.findParentOf(keyname);
      if (node != nullptr)
      {
        std::string entryname = node->suffix(keyname); // get everything beyond last ':'
        Param::ParamNode::EntryIterator it = node->findEntry(entryname);
        if (it != node->entries.end())
        {
          node->entries.erase(it); // delete entry
          if (node->nodes.empty()  && node->entries.empty())
          {
            // delete if section is now empty
            remove(keyname.substr(0, keyname.length() - entryname.length())); // keep last ':' to indicate deletion of a section
          }
        }
      }
    }
  }

  void Param::removeAll(const std::string& prefix)
  {
    if (!prefix.empty() && prefix.back() == ':')//we have to delete one node only (and its subnodes)
    {
      ParamNode* node = root_.findParentOf(prefix.substr(0, prefix.length() - 1));
      if (node != nullptr)
      {
        Param::ParamNode::NodeIterator it = node->findNode(node->suffix(prefix.substr(0, prefix.length() - 1)));
        if (it != node->nodes.end())
        {
          std::string name = it->name;
          node->nodes.erase(it); // will automatically delete subnodes
          if (node->nodes.empty()  && node->entries.empty())
          {
            // delete last section name (could be partial)
            removeAll(prefix.substr(0, prefix.length() - name.length() - 1));// '-1' for the tailing ':'
          }
        }
      }
    }
    else //we have to delete all entries and nodes starting with the prefix
    {
      ParamNode* node = root_.findParentOf(prefix);
      if (node != nullptr)
      {
        std::string suffix = node->suffix(prefix); // name behind last ":"

        for (Param::ParamNode::NodeIterator it = node->nodes.begin(); it != node->nodes.end(); /*do nothing*/)
        {
          if (it->name.compare(0, suffix.length(), suffix) == 0)
          {
            it = node->nodes.erase(it);
          }
          else if (it != node->nodes.end())
          {
            ++it;
          }
        }
        for (Param::ParamNode::EntryIterator it = node->entries.begin(); it != node->entries.end(); /*do nothing*/)
        {
          if (it->name.compare(0, suffix.size(), suffix) == 0)
          {
            it = node->entries.erase(it);
          }
          else if (it != node->entries.end())
          {
            ++it;
          }
        }
        // the parent node might now be empty (delete it as well - otherwise the trace might be broken)
        if (node->nodes.empty() && node->entries.empty())
        {
          // delete last section name (could be partial)
          removeAll(prefix.substr(0, prefix.size() - suffix.size()));
        }
      }
    }
  }

  Param Param::copySubset(const Param& subset) const
  {
    ParamNode out("ROOT", "");

    for (const auto& entry : subset.root_.entries)
    {
      const auto& n = root_.findEntry(entry.name);
      if (n == root_.entries.end())
      {
        OPENMS_LOG_WARN << "Warning: Trying to copy non-existent parameter entry " << entry.name << std::endl;
      }
      else
      {
        out.insert(*n);
      }
    }

    for (const auto& node : subset.root_.nodes)
    {
      const auto& n = root_.findNode(node.name);
      if (n == root_.nodes.end())
      {
        OPENMS_LOG_WARN << "Warning: Trying to copy non-existent parameter node " << node.name << std::endl;
      }
      else
      {
        out.insert(*n);
      }
    }
    return Param(out);
  }

  Param Param::copy(const std::string& prefix, bool remove_prefix) const
  {
    ParamNode out("ROOT", "");

    ParamNode* node = root_.findParentOf(prefix);
    if (node == nullptr)
    {
      return Param();
    }
    //we have to copy this node only
    if (!prefix.empty() && prefix.back() == ':')
    {
      if (remove_prefix)
      {
        out = *node;
      }
      else
      {
        out.insert(*node, prefix.substr(0, prefix.size() - node->name.size() - 1));
      }
    }
    else //we have to copy all entries and nodes starting with the right suffix
    {
      std::string suffix = node->suffix(prefix);
      for (Param::ParamNode::NodeIterator it = node->nodes.begin(); it != node->nodes.end(); ++it)
      {
        if (it->name.compare(0, suffix.size(), suffix) == 0)
        {
          if (remove_prefix)
          {
            ParamNode tmp = *it;
            tmp.name = tmp.name.substr(suffix.size());
            out.insert(tmp);
          }
          else
          {
            out.insert(*it, prefix.substr(0, prefix.size() - suffix.size()));
          }
        }
      }
      for (Param::ParamNode::EntryIterator it = node->entries.begin(); it != node->entries.end(); ++it)
      {
        if (it->name.compare(0, suffix.size(), suffix) == 0)
        {
          if (remove_prefix)
          {
            ParamEntry tmp = *it;
            tmp.name = tmp.name.substr(suffix.size());
            out.insert(tmp);
          }
          else
          {
            out.insert(*it, prefix.substr(0, prefix.size() - suffix.size()));
          }
        }
      }
    }

    return Param(out);
  }

  void Param::parseCommandLine(const int argc, const char** argv, const std::string& prefix)
  {
    //determine prefix
    std::string prefix2 = prefix;
    if (!prefix2.empty())
    {
      //prefix2.ensureLastChar(':');
      if (prefix2.back() != ':')
      {
        prefix2.append(1, ':');
      }
    }

    //parse arguments
    std::string arg, arg1;
    for (int i = 1; i < argc; ++i)
    {
      //load the current and next argument:  arg and arg1 ("" at the last argument)
      arg = argv[i];
      arg1 = "";
      if (i + 1 < argc)
      {
        arg1 = argv[i + 1];
      }

      //it is a option when it starts with a '-' and the second character is not a number
      bool arg_is_option = false;
      if (arg.size() >= 2 && arg[0] == '-' && arg[1] != '0' && arg[1] != '1' && arg[1] != '2' && arg[1] != '3' && arg[1] != '4' && arg[1] != '5' && arg[1] != '6' && arg[1] != '7' && arg[1] != '8' && arg[1] != '9')
      {
        arg_is_option = true;
      }
      bool arg1_is_option = false;
      if (arg1.size() >= 2 && arg1[0] == '-' && arg1[1] != '0' && arg1[1] != '1' && arg1[1] != '2' && arg1[1] != '3' && arg1[1] != '4' && arg1[1] != '5' && arg1[1] != '6' && arg1[1] != '7' && arg1[1] != '8' && arg1[1] != '9')
      {
        arg1_is_option = true;
      }
      //cout << "Parse: '"<< arg << "' '" << arg1 << "'" << std::endl;

      //flag (option without text argument)
      if (arg_is_option && arg1_is_option)
      {
        root_.insert(ParamEntry(arg, std::string(), ""), prefix2);
      }
      //option with argument
      else if (arg_is_option && !arg1_is_option)
      {
        root_.insert(ParamEntry(arg, arg1, ""), prefix2);
        ++i;
      }
      //just text arguments (not preceded by an option)
      else
      {

        ParamEntry* misc_entry = root_.findEntryRecursive(prefix2 + "misc");
        if (misc_entry == nullptr)
        {
          std::vector<std::string> sl;
          sl.push_back(arg);
          // create "misc"-Node:
          root_.insert(ParamEntry("misc", sl, ""), prefix2);
        }
        else
        {
          std::vector<std::string> sl = misc_entry->value;
          sl.push_back(arg);
          misc_entry->value = sl;
        }
      }
    }
  }

  void Param::parseCommandLine(const int argc, const char** argv, const std::map<std::string, std::string>& options_with_one_argument, const std::map<std::string, std::string>& options_without_argument, const std::map<std::string, std::string>& options_with_multiple_argument, const std::string& misc, const std::string& unknown)
  {
    //determine misc key
    
    //determine unknown key
    
    //parse arguments
    std::string arg, arg1;
    for (int i = 1; i < argc; ++i)
    {
      //load the current and next argument:  arg and arg1 ("" at the last argument)
      arg = argv[i];
      arg1 = "";
      if (i + 1 < argc)
      {
        arg1 = argv[i + 1];
      }

      //it is a option when it starts with a '-' and the second character is not a number
      bool arg_is_option = false;
      if (arg.size() >= 2 && arg[0] == '-' && arg[1] != '0' && arg[1] != '1' && arg[1] != '2' && arg[1] != '3' && arg[1] != '4' && arg[1] != '5' && arg[1] != '6' && arg[1] != '7' && arg[1] != '8' && arg[1] != '9')
      {
        arg_is_option = true;
      }
      bool arg1_is_option = false;
      if (arg1.size() >= 2 && arg1[0] == '-' && arg1[1] != '0' && arg1[1] != '1' && arg1[1] != '2' && arg1[1] != '3' && arg1[1] != '4' && arg1[1] != '5' && arg1[1] != '6' && arg1[1] != '7' && arg1[1] != '8' && arg1[1] != '9')
      {
        arg1_is_option = true;
      }

      //with multiple argument
      if (options_with_multiple_argument.find(arg) != options_with_multiple_argument.end())
      {
        //next argument is an option
        if (arg1_is_option)
        {
          root_.insert(ParamEntry("", std::vector<std::string>(), ""), options_with_multiple_argument.find(arg)->second);
        }
        //next argument is not an option
        else
        {
          std::vector<std::string> sl;
          int j = (i + 1);
          while (j < argc && !(arg1.size() >= 2 && arg1[0] == '-' && arg1[1] != '0' && arg1[1] != '1' && arg1[1] != '2' && arg1[1] != '3' && arg1[1] != '4' && arg1[1] != '5' && arg1[1] != '6' && arg1[1] != '7' && arg1[1] != '8' && arg1[1] != '9'))
          {
            sl.push_back(arg1);
            ++j;
            if (j < argc)
            {
              arg1 = argv[j];
            }
          }

          root_.insert(ParamEntry("", sl, ""), options_with_multiple_argument.find(arg)->second);
          i = j - 1;
        }
      }
      //without argument
      else if (options_without_argument.find(arg) != options_without_argument.end())
      {
        root_.insert(ParamEntry("", "true", ""), options_without_argument.find(arg)->second);
      }
      //with one argument
      else if (options_with_one_argument.find(arg) != options_with_one_argument.end())
      {
        //next argument is not an option
        if (!arg1_is_option)
        {
          root_.insert(ParamEntry("", arg1, ""), options_with_one_argument.find(arg)->second);
          ++i;
        }
        //next argument is an option
        else
        {

          root_.insert(ParamEntry("", std::string(), ""), options_with_one_argument.find(arg)->second);
        }
      }
      //unknown option
      else if (arg_is_option)
      {
        ParamEntry* unknown_entry = root_.findEntryRecursive(unknown);
        if (unknown_entry == nullptr)
        {
          std::vector<std::string> sl;
          sl.push_back(arg);
          root_.insert(ParamEntry("", sl, ""), unknown);
        }
        else
        {
          std::vector<std::string> sl = unknown_entry->value;
          sl.push_back(arg);
          unknown_entry->value = sl;
        }
      }
      //just text argument
      else
      {
        ParamEntry* misc_entry = root_.findEntryRecursive(misc);
        if (misc_entry == nullptr)
        {
          std::vector<std::string> sl;
          sl.push_back(arg);
          // create "misc"-Node:
          root_.insert(ParamEntry("", sl, ""), misc);
        }
        else
        {
          std::vector<std::string> sl = misc_entry->value;
          sl.push_back(arg);
          misc_entry->value = sl;
        }
      }
    }
  }

  std::ostream& operator<<(std::ostream& os, const Param& param)
  {
    for (Param::ParamIterator it = param.begin(); it != param.end(); ++it)
    {
      os << '"';
      if (it.getName().length() > it->name.length() + 1)
      {
        os << it.getName().substr(0, it.getName().length() - it->name.length() - 1) << "|";
      }
      os  << it->name << "\" -> \"" << it->value << '"';
      if (!it->description.empty())
      {
        os << " (" << it->description << ")";
      }
      os << std::endl;
    }
    return os;
  }

  size_t Param::size() const
  {
    return root_.size();
  }

  bool Param::empty() const
  {
    return size() == 0;
  }

  void Param::clear()
  {
    root_ = ParamNode("ROOT", "");
  }

  void Param::checkDefaults(const std::string& name, const Param& defaults, const std::string& prefix) const
  {
    //Extract right parameters
    std::string prefix2 = prefix;
    if (!prefix2.empty())
    {
      if (prefix2.back() != ':') 
      {
          prefix2 += ':';
      }
    }
    Param check_values = copy(prefix2, true);

    //check
    for (ParamIterator it = check_values.begin(); it != check_values.end(); ++it)
    {
      //unknown parameter
      if (!defaults.exists(it.getName()))
      {
        OPENMS_LOG_WARN << "Warning: " << name << " received the unknown parameter '" << it.getName() << "'";
        if (!prefix2.empty())
        {
          OPENMS_LOG_WARN << " in '" << prefix2 << "'";
        }
        OPENMS_LOG_WARN << "!" << std::endl;
      }

      //different types
      ParamEntry* default_value = defaults.root_.findEntryRecursive(prefix2 + it.getName());
      if (default_value == nullptr)
      {
        continue;
      }
      if (default_value->value.valueType() != it->value.valueType())
      {
        std::string d_type;
        if (default_value->value.valueType() == ParamValue::STRING_VALUE)
        {
          d_type = "string";
        }
        if (default_value->value.valueType() == ParamValue::STRING_LIST)
        {
          d_type = "string list";
        }
        if (default_value->value.valueType() == ParamValue::EMPTY_VALUE)
        {
          d_type = "empty";
        }
        if (default_value->value.valueType() == ParamValue::INT_VALUE)
        {
          d_type = "integer";
        }
        if (default_value->value.valueType() == ParamValue::INT_LIST)
        {
          d_type = "integer list";
        }
        if (default_value->value.valueType() == ParamValue::DOUBLE_VALUE)
        {
          d_type = "float";
        }
        if (default_value->value.valueType() == ParamValue::DOUBLE_LIST)
        {
          d_type = "float list";
        }
        std::string p_type;
        if (it->value.valueType() == ParamValue::STRING_VALUE)
        {
          p_type = "string";
        }
        if (it->value.valueType() == ParamValue::STRING_LIST)
        {
          p_type = "string list";
        }
        if (it->value.valueType() == ParamValue::EMPTY_VALUE)
        {
          p_type = "empty";
        }
        if (it->value.valueType() == ParamValue::INT_VALUE)
        {
          p_type = "integer";
        } 
        if (it->value.valueType() == ParamValue::INT_LIST)
        {
          p_type = "integer list";
        }
        if (it->value.valueType() == ParamValue::DOUBLE_VALUE)
        {
          p_type = "float";
        }
        if (it->value.valueType() == ParamValue::DOUBLE_LIST)
        {
          p_type = "float list";
        }

        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name + ": Wrong parameter type '" + p_type + "' for " + d_type + " parameter '" + it.getName() + "' given!");
      }
      //parameter restrictions
      ParamEntry pe = *default_value;
      pe.value = it->value;
      std::string s;
      if (!pe.isValid(s))
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name + ": " + s);
    }
  }

  Param::ParamIterator Param::findFirst(const std::string& leaf) const
  {
    for (Param::ParamIterator it = this->begin(); it != this->end(); ++it)
    {
      std::string suffix = ":" + leaf;
      if (!(suffix.length() > it.getName().length()) &&
          it.getName().compare(it.getName().length() - suffix.length(), suffix.length(), suffix) == 0)
      {
        return it;
      }
    }
    return this->end();
  }

  Param::ParamIterator Param::findNext(const std::string& leaf, const ParamIterator& start_leaf) const
  {
    // start at NEXT entry
    Param::ParamIterator it = start_leaf;
    if (it != this->end()) ++it;

    for (; it != this->end(); ++it)
    {
      std::string suffix = ":" + leaf;
      if (!(suffix.length() > it.getName().length()) &&
          it.getName().compare(it.getName().length() - suffix.length(), suffix.length(), suffix) == 0)
      {
        return it;
      }
    }
    return this->end();
  }

  bool Param::update(const Param& p_outdated, const bool add_unknown)
  {
    return update(p_outdated, add_unknown, OpenMS_Log_warn);
  }

  bool Param::update(const Param& p_outdated, const bool add_unknown, Logger::LogStream& stream)
  {
    bool fail_on_invalid_values = false;
    bool fail_on_unknown_parameters = false;
    return update(p_outdated, true, add_unknown, fail_on_invalid_values, fail_on_unknown_parameters, stream);
  }

  
  bool Param::update(const Param& p_outdated, bool verbose, const bool add_unknown, bool fail_on_invalid_values, bool fail_on_unknown_parameters, Logger::LogStream& stream)
  {
    bool is_update_success(true);
    // augment
    for (Param::ParamIterator it = p_outdated.begin(); it != p_outdated.end(); ++it)
    {
      Param::ParamEntry new_entry; // entry of new location (used to retain new description)
      std::string target_name; // fully qualified name in new param

      if (this->exists(it.getName()))
      {
        // param 'version': do not override!
        std::string suffix = ":version";
        if (!(suffix.length() > it.getName().length()) &&
            it.getName().compare(it.getName().length() - suffix.length(), suffix.length(), suffix) == 0)
        {
          if (this->getValue(it.getName()) != it->value)
          {
OPENMS_THREAD_CRITICAL(LOGSTREAM)
            stream << "Warning: for ':version' entry, augmented and Default Ini-File differ in value. Default value will not be altered!\n";
          }
          continue;
        }
        // param 'type': do not override!
        else if (suffix = ":type",
                !(suffix.length() > it.getName().length()) &&
                it.getName().compare(it.getName().length() - suffix.length(), suffix.length(), suffix) == 0) // only for TOPP type (e.g. PeakPicker:1:type), any other 'type' param is ok
        {
          size_t first = it.getName().find(':');
          if (first != std::string::npos &&
              it.getName().find(':', first+1) != std::string::npos)
          {
            if (this->getValue(it.getName()) != it->value) 
            {
                OPENMS_THREAD_CRITICAL(LOGSTREAM)
                stream << "Warning: for ':type' entry, augmented and Default Ini-File differ in value. Default value will not be altered!\n";
            }
            continue;
          }
        }

        // all other parameters:
        new_entry = this->getEntry(it.getName());
        target_name = it.getName();

      }
      else // outdated param non-existent in new param
      {
        // search by suffix in new param. Only match complete names, e.g. myname will match newsection:myname, but not newsection:othermyname
        Param::ParamEntry l1_entry = p_outdated.getEntry(it.getName());
        // since the outdated param with full path does not exist within new param,
        // we will never find the new entry by using exists() as above, thus
        // its safe to modify it here

        ParamIterator it_match = this->findFirst(l1_entry.name);
        if (it_match != this->end())
        {
          // make sure the same leaf name does not exist at any other position
          if (this->findNext(l1_entry.name, it_match) == this->end())
          {
OPENMS_THREAD_CRITICAL(LOGSTREAM)
            stream << "Found '" << it.getName() << "' as '" << it_match.getName() << "' in new param." << std::endl;
            new_entry = this->getEntry(it_match.getName());
            target_name = it_match.getName();
          }
        }

        if (target_name.empty()) // no mapping was found
        {
          if (fail_on_unknown_parameters)
          {
OPENMS_THREAD_CRITICAL(LOGSTREAM)
            stream << "Unknown (or deprecated) Parameter '" << it.getName() << "' given in outdated parameter file!" << std::endl;
            is_update_success = false;
          }
          else if (add_unknown)
          {
OPENMS_THREAD_CRITICAL(LOGSTREAM)
            stream << "Unknown (or deprecated) Parameter '" << it.getName() << "' given in outdated parameter file! Adding to current set." << std::endl;
            Param::ParamEntry local_entry = p_outdated.getEntry(it.getName());
            std::string prefix = "";
            if (it.getName().find(':') != std::string::npos)
            {
              prefix = it.getName().substr(0, 1 + it.getName().find_last_of(':'));
            }
            this->root_.insert(local_entry, prefix); //->setValue(it.getName(), local_entry.value, local_entry.description, local_entry.tags);
          }
          else if (verbose)
          {
OPENMS_THREAD_CRITICAL(LOGSTREAM)
            stream << "Unknown (or deprecated) Parameter '" << it.getName() << "' given in outdated parameter file! Ignoring parameter. " << std::endl;
          }
          continue;
        }
      }

      // do the actual updating (we found a matching pair)
      if (new_entry.value.valueType() == it->value.valueType())
      {
        if (new_entry.value != it->value)
        {
          // check entry for consistency (in case restrictions have changed)
          ParamValue default_value = new_entry.value;
          new_entry.value = it->value;
          std::string validation_result;
          if (new_entry.isValid(validation_result))
          {
            // overwrite default value
            if (verbose) 
            {
OPENMS_THREAD_CRITICAL(LOGSTREAM)
                stream << "Default-Parameter '" << target_name << "' overridden: '" << default_value << "' --> '" << it->value << "'!" << std::endl;
            }
            this->setValue(target_name, it->value, new_entry.description, this->getTags(target_name));
          }
          else
          {
OPENMS_THREAD_CRITICAL(LOGSTREAM)
            stream << validation_result;
            if (fail_on_invalid_values)
            {
OPENMS_THREAD_CRITICAL(LOGSTREAM)
              stream << " Updating failed!" << std::endl;
              is_update_success = false;
            }
            else
            {
OPENMS_THREAD_CRITICAL(LOGSTREAM)
              stream << " Ignoring invalid value (using new default '" << default_value << "')!" << std::endl;
              new_entry.value = default_value;
            }
          }
        }
        else
        {
          // value stayed the same .. nothing to be done
        }
      }
      else
      {
OPENMS_THREAD_CRITICAL(LOGSTREAM)
        stream << "Parameter '" << it.getName() << "' has changed value type!\n";
        if (fail_on_invalid_values)
        {
OPENMS_THREAD_CRITICAL(LOGSTREAM)
          stream << " Updating failed!" << std::endl;
          is_update_success = false;
        } 
        else
        {
OPENMS_THREAD_CRITICAL(LOGSTREAM)
          stream << " Ignoring invalid value (using new default)!" << std::endl;
        }
      }

    } // next param in outdated tree

    return is_update_success;
  }

  void Param::merge(const OpenMS::Param& toMerge)
  {
    // keep track of the path inside the param tree
    std::string pathname;

    // augment
    for (Param::ParamIterator it = toMerge.begin(); it != toMerge.end(); ++it)
    {
      std::string prefix = "";
      if (it.getName().find(':') != std::string::npos)
        prefix = it.getName().substr(0, 1 + it.getName().find_last_of(':'));

      // we care only about values that do not exist already
      if (!this->exists(it.getName()))
      {
        Param::ParamEntry entry = *it;
        OPENMS_LOG_DEBUG << "[Param::merge] merging " << it.getName() << std::endl;
        this->root_.insert(entry, prefix);
      }

      //copy section descriptions
      const std::vector<ParamIterator::TraceInfo>& trace = it.getTrace();
      for (std::vector<ParamIterator::TraceInfo>::const_iterator traceIt = trace.begin(); traceIt != trace.end(); ++traceIt)
      {
        if (traceIt->opened)
        {
          OPENMS_LOG_DEBUG << "[Param::merge] extending param trace " << traceIt->name << " (" << pathname << ")" << std::endl;
          pathname += traceIt->name + ":";
        }
        else
        {
          OPENMS_LOG_DEBUG << "[Param::merge] reducing param trace " << traceIt->name << " (" << pathname << ")" << std::endl;
          std::string suffix = traceIt->name + ":";
          if (suffix.size() <= pathname.size() && pathname.compare(pathname.size() - suffix.size(), suffix.size(), suffix) == 0)
            pathname.resize(pathname.size() - traceIt->name.size() - 1);
        }
        std::string real_pathname = pathname.substr(0, pathname.size() - 1);//remove ':' at the end
        if (!real_pathname.empty())
        {
          std::string description_old = getSectionDescription(prefix + real_pathname);
          const std::string& description_new = toMerge.getSectionDescription(real_pathname);
          if (description_old.empty())
          {
            setSectionDescription(real_pathname, description_new);
          }
        }
      }

    }
  }

  void Param::setSectionDescription(const std::string& key, const std::string& description)
  {
    ParamNode* node = root_.findParentOf(key);
    if (node == nullptr)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, key);
    }

    Param::ParamNode::NodeIterator it = node->findNode(node->suffix(key));
    if (it == node->nodes.end())
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, key);
    }
    it->description = description;
  }

  void Param::addSection(const std::string& key, const std::string& description)
  {
    root_.insert(ParamNode("",description),key);
  }

  Param::ParamIterator Param::begin() const
  {
    return ParamIterator(root_);
  }

  Param::ParamIterator Param::end() const
  {
    return ParamIterator();
  }

  Param::ParamIterator::ParamIterator() :
    root_(nullptr),
    current_(0),
    stack_(),
    trace_()
  {
  }

  Param::ParamIterator::ParamIterator(const Param::ParamNode& root) :
    root_(&root),
    current_(-1),
    stack_(),
    trace_()
  {
    //Empty Param => begin == end iterator
    if (root_->entries.empty() && root_->nodes.empty())
    {
      root_ = nullptr;
      return;
    }

    //find first entry
    stack_.push_back(root_);
    operator++();
  }

  Param::ParamIterator::~ParamIterator() = default;

  const Param::ParamEntry& Param::ParamIterator::operator*()
  {
    return stack_.back()->entries[current_];
  }

  const Param::ParamEntry* Param::ParamIterator::operator->()
  {
    return &(stack_.back()->entries[current_]);
  }

  Param::ParamIterator Param::ParamIterator::operator++(int)
  {
    ParamIterator tmp(*this);
    ++(*this);
    return tmp;
  }

  Param::ParamIterator& Param::ParamIterator::operator++()
  {
    if (root_ == nullptr)
    {
      return *this;
    }
    trace_.clear();
    while (true)
    {
      const Param::ParamNode* node = stack_.back();

      //std::cout << "############ operator++ #### " << node->name << " ## " << current_ << std::endl;

      //check if there is a next entry in the current node
      if (current_ + 1 < (int)node->entries.size())
      {
        //std::cout << " - next entry" << std::endl;
        ++current_;
        return *this;
      }
      //visit subnodes after entries
      else if (!node->nodes.empty())
      {
        current_ = -1;
        stack_.push_back(&(node->nodes[0]));
        //std::cout << " - entering into: " << node->nodes[0].name << std::endl;
        //track changes (enter a node)
        trace_.emplace_back(node->nodes[0].name, node->nodes[0].description, true);

        continue;
      }
      //go back in tree until the node we came from is not the last subnode
      //of the current node. Go into the next subnode.
      else
      {
        while (true)
        {
          const Param::ParamNode* last = node;
          stack_.pop_back();
          //std::cout << " - stack size: " << stack_.size() << std::endl;
          //we have reached the end
          if (stack_.empty())
          {
            //std::cout << " - reached the end" << std::endl;
            root_ = nullptr;
            return *this;
          }
          node = stack_.back();

          //std::cout << " - last was: " << last->name << std::endl;
          //std::cout << " - descended to: " << node->name << std::endl;

          //track changes (leave a node)
          if (!trace_.empty() && trace_.back().name == last->name && trace_.back().opened) // was empty subnode
          {
            trace_.pop_back();
          }
          else
          {
            trace_.emplace_back(last->name, last->description, false);
          }

          //check of new subtree is accessible
          unsigned int next_index = (last - &(node->nodes[0])) + 1;
          if (next_index < node->nodes.size())
          {
            current_ = -1;
            stack_.push_back(&(node->nodes[next_index]));
            //cout << " - entering into: " << node->nodes[next_index].name  << endl;
            //track changes (enter a node)
            trace_.emplace_back(node->nodes[next_index].name, node->nodes[next_index].description, true);
            break;
          }
        }
      }
    }
  }

  bool Param::ParamIterator::operator==(const ParamIterator& rhs) const
  {
    return (root_ == nullptr && rhs.root_ == nullptr) || (stack_ == rhs.stack_ && current_ == rhs.current_);
  }

  bool Param::ParamIterator::operator!=(const ParamIterator& rhs) const
  {
    return !operator==(rhs);
  }

  std::string Param::ParamIterator::getName() const
  {
    std::string tmp;
    for (std::vector<const Param::ParamNode*>::const_iterator it = stack_.begin() + 1; it != stack_.end(); ++it)
    {
      tmp += (*it)->name + ':';
    }
    return tmp + stack_.back()->entries[current_].name;
  }

  const std::vector<Param::ParamIterator::TraceInfo>& Param::ParamIterator::getTrace() const
  {
    return trace_;
  }

  const Param::ParamEntry& Param::getEntry(const std::string& key) const
  {
    return getEntry_(key);
  }

  ParamValue::ValueType Param::getValueType(const std::string& key) const
  {
    return getEntry_(key).value.valueType();
  }

  const std::string& Param::getDescription(const std::string& key) const
  {
    return getEntry_(key).description;
  }

  void Param::addTag(const std::string& key, const std::string& tag)
  {
    if (tag.find(',') != std::string::npos)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Param tags may not contain comma characters", tag);
    }
    getEntry_(key).tags.insert(tag);
  }

  void Param::addTags(const std::string& key, const std::vector<std::string>& tags)
  {
    ParamEntry& entry = getEntry_(key);
    for (size_t i = 0; i != tags.size(); ++i)
    {
      if (tags[i].find(',') != std::string::npos)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Param tags may not contain comma characters", tags[i]);
      }
      entry.tags.insert(tags[i]);
    }
  }

  std::vector<std::string> Param::getTags(const std::string& key) const
  {
    ParamEntry& entry = getEntry_(key);
    std::vector<std::string> list;
    for (std::set<std::string>::const_iterator it = entry.tags.begin(); it != entry.tags.end(); ++it)
    {
      list.push_back(*it);
    }
    return list;
  }

  void Param::clearTags(const std::string& key)
  {
    getEntry_(key).tags.clear();
  }

  bool Param::hasTag(const std::string& key, const std::string& tag) const
  {
    return getEntry_(key).tags.count(tag);
  }

  bool Param::exists(const std::string& key) const
  {
    return root_.findEntryRecursive(key);
  }

  bool Param::hasSection(const std::string &key) const
  {
    if (key.back() == ':')
    {
      // Remove trailing colon from key
      return root_.findParentOf(key.substr(0, key.size() - 1)) != nullptr;
    }
    else
    {
      return root_.findParentOf(key) != nullptr;
    }
  }

  Param::ParamEntry& Param::getEntry_(const std::string& key) const
  {
    ParamEntry* entry = root_.findEntryRecursive(key);
    if (entry == nullptr)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, key);
    }

    return *entry;
  }

} //namespace
