// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/Param.h>

#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/DATASTRUCTURES/Map.h>

#include <QtCore/QString>
#include <fstream>

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
    min_int(-std::numeric_limits<Int>::max()),
    max_int(std::numeric_limits<Int>::max()),
    valid_strings()
  {
  }

  Param::ParamEntry::ParamEntry(const String& n, const DataValue& v, const String& d, const StringList& t) :
    name(n),
    description(d),
    value(v),
    tags(),
    min_float(-std::numeric_limits<double>::max()),
    max_float(std::numeric_limits<double>::max()),
    min_int(-std::numeric_limits<Int>::max()),
    max_int(std::numeric_limits<Int>::max()),
    valid_strings()
  {
    //add tags
    for (Size i = 0; i < t.size(); ++i)
    {
      tags.insert(t[i]);
    }
    //check name
    if (name.has(':'))
    {
      std::cerr << "Error ParamEntry name must not contain ':' characters!" << std::endl;
    }
  }

  Param::ParamEntry::~ParamEntry()
  {
  }

  bool Param::ParamEntry::isValid(String& message) const
  {
    if (value.valueType() == DataValue::STRING_VALUE)
    {
      if (valid_strings.size() != 0)
      {
        bool ok = false;
        if (std::find(valid_strings.begin(), valid_strings.end(), value) != valid_strings.end())
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
          valid.concatenate(valid_strings.begin(), valid_strings.end(), ",");
          message = "Invalid string parameter value '" + static_cast<String>(value) + "' for parameter '" + name + "' given! Valid values are: '" + valid + "'.";
          return false;
        }
      }
    }
    else if (value.valueType() == DataValue::STRING_LIST)
    {
      String str_value;
      StringList ls_value = value;
      for (Size i = 0; i < ls_value.size(); ++i)
      {
        str_value = ls_value[i];

        if (valid_strings.size() != 0)
        {
          bool ok = false;
          if (std::find(valid_strings.begin(), valid_strings.end(), str_value) != valid_strings.end())
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
            valid.concatenate(valid_strings.begin(), valid_strings.end(), ",");
            message = "Invalid string parameter value '" + str_value + "' for parameter '" + name + "' given! Valid values are: '" + valid + "'.";
            return false;
          }
        }
      }
    }
    else if (value.valueType() == DataValue::INT_VALUE)
    {
      Int tmp = value;
      if ((min_int != -std::numeric_limits<Int>::max() && tmp < min_int) || (max_int != std::numeric_limits<Int>::max() && tmp > max_int))
      {
        message = String("Invalid integer parameter value '") + String(tmp) + "' for parameter '" + name + "' given! The valid range is: [" + min_int + ":" + max_int + "].";
        return false;
      }
    }
    else if (value.valueType() == DataValue::INT_LIST)
    {
      Int int_value;
      IntList ls_value = value;
      for (Size i = 0; i < ls_value.size(); ++i)
      {
        int_value = ls_value[i];
        if ((min_int != -std::numeric_limits<Int>::max() && int_value < min_int) || (max_int != std::numeric_limits<Int>::max() && int_value > max_int))
        {
          message = String("Invalid integer parameter value '") + int_value + "' for parameter '" + name + "' given! The valid range is: [" + min_int + ":" + max_int + "].";
          return false;
        }
      }
    }
    else if (value.valueType() == DataValue::DOUBLE_VALUE)
    {
      double tmp = value;
      if ((min_float != -std::numeric_limits<double>::max() && tmp < min_float) || (max_float != std::numeric_limits<double>::max() && tmp > max_float))
      {
        message = String("Invalid double parameter value '") + tmp + "' for parameter '" + name + "' given! The valid range is: [" + min_float + ":" + max_float + "].";
        return false;
      }
    }
    else if (value.valueType() == DataValue::DOUBLE_LIST)
    {
      DoubleList ls_value = value;
      for (Size i = 0; i < ls_value.size(); ++i)
      {
        double dou_value = ls_value[i];
        if ((min_float != -std::numeric_limits<double>::max() && dou_value < min_float) || (max_float != std::numeric_limits<double>::max() && dou_value > max_float))
        {
          message = String("Invalid double parameter value '") + dou_value + "' for parameter '" + name + "' given! The valid range is: [" + min_float + ":" + max_float + "].";
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

  Param::ParamNode::ParamNode(const String& n, const String& d) :
    name(n),
    description(d),
    entries(),
    nodes()
  {
    if (name.has(':'))
    {
      std::cerr << "Error ParamNode name must not contain ':' characters!" << std::endl;
    }
  }

  Param::ParamNode::~ParamNode()
  {
  }

  bool Param::ParamNode::operator==(const ParamNode& rhs) const
  {
    if (name != rhs.name || entries.size() != rhs.entries.size() || nodes.size() != rhs.nodes.size())
      return false;

    //order of sections / entries should not matter
    for (Size i = 0; i < entries.size(); ++i)
    {
      if (find(rhs.entries.begin(), rhs.entries.end(), entries[i]) == rhs.entries.end())
      {
        return false;
      }
    }
    for (Size i = 0; i < nodes.size(); ++i)
    {
      if (find(rhs.nodes.begin(), rhs.nodes.end(), nodes[i]) == rhs.nodes.end())
      {
        return false;
      }
    }

    return true;
  }

  Param::ParamNode::EntryIterator Param::ParamNode::findEntry(const String& local_name)
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

  Param::ParamNode::NodeIterator Param::ParamNode::findNode(const String& local_name)
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

  Param::ParamNode* Param::ParamNode::findParentOf(const String& local_name)
  {
    //cout << "findParentOf nodename: " << this->name << " - nodes: " << this->nodes.size() << " - find: "<< name << std::endl;
    if (!local_name.has(':')) // we are in the right child
    {
      //check if a node or entry prefix match
      for (Size i = 0; i < nodes.size(); ++i)
      {
        if (nodes[i].name.hasPrefix(local_name))
          return this;
      }
      for (Size i = 0; i < entries.size(); ++i)
      {
        if (entries[i].name.hasPrefix(local_name))
          return this;
      }
      return nullptr;
    }
    else //several subnodes to browse through
    {
      String prefix = local_name.prefix(':');
      //cout << " - Prefix: '" << prefix << "'" << std::endl;
      NodeIterator it = findNode(prefix);
      if (it == nodes.end()) //subnode not found
      {
        return nullptr;
      }
      //recursively call findNode for the rest of the path
      String new_name = local_name.substr(it->name.size() + 1);
      //cout << " - Next name: '" << new_name << "'" << std::endl;
      return it->findParentOf(new_name);
    }
  }

  Param::ParamEntry* Param::ParamNode::findEntryRecursive(const String& local_name)
  {
    ParamNode* parent = findParentOf(local_name);
    if (parent == nullptr)
      return nullptr;

    EntryIterator it = parent->findEntry(suffix(local_name));
    if (it == parent->entries.end())
      return nullptr;

    return &(*it);
  }

  void Param::ParamNode::insert(const ParamNode& node, const String& prefix)
  {
    //std::cerr << "INSERT NODE  " << node.name << " (" << prefix << ")" << std::endl;
    String prefix2 = prefix + node.name;

    ParamNode* insert_node = this;
    while (prefix2.has(':'))
    {
      String local_name = prefix2.prefix(':');
      //check if the node already exists
      NodeIterator it = insert_node->findNode(local_name);
      if (it != insert_node->nodes.end()) //exists
      {
        insert_node = &(*it);
      }
      else //create it
      {
        insert_node->nodes.push_back(ParamNode(local_name, ""));
        insert_node = &(insert_node->nodes.back());
        //std::cerr << " - Created new node: " << insert_node->name << std::endl;
      }
      //remove prefix
      prefix2 = prefix2.substr(local_name.size() + 1);
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
      if (it->description == "" || node.description != "") //replace description if not empty in new node
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
    //std::cerr << "INSERT ENTRY " << entry.name << " (" << prefix << ")" << std::endl;
    String prefix2 = prefix + entry.name;
    //std::cerr << " - inserting: " << prefix2 << std::endl;

    ParamNode* insert_node = this;
    while (prefix2.has(':'))
    {
      String local_name = prefix2.prefix(':');
      //std::cerr << " - looking for node: " << name << std::endl;
      //look up if the node already exists
      NodeIterator it = insert_node->findNode(local_name);
      if (it != insert_node->nodes.end()) //exists
      {
        insert_node = &(*it);
      }
      else //create it
      {
        insert_node->nodes.push_back(ParamNode(local_name, ""));
        insert_node = &(insert_node->nodes.back());
        //std::cerr << " - Created new node: " << insert_node->name << std::endl;
      }
      //remove prefix
      prefix2 = prefix2.substr(local_name.size() + 1);
      //std::cerr << " - new prefix: " << prefix2 << std::endl;
    }

    //check if the entry already exists
    //std::cerr << " - final entry name: " << prefix2 << std::endl;
    EntryIterator it = insert_node->findEntry(prefix2);
    if (it != insert_node->entries.end()) //overwrite entry
    {
      it->value = entry.value;
      it->tags = entry.tags;
      if (it->description == "" || entry.description != "") //replace description if not empty in new entry
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
    for (std::vector<ParamNode>::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
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

  Param::Param() :
    root_("ROOT", "")
  {
  }

  Param::~Param()
  {
  }

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

  void Param::setValue(const String& key, const DataValue& value, const String& description, const StringList& tags)
  {
    root_.insert(ParamEntry("", value, description, tags), key);
  }

  void Param::setValidStrings(const String& key, const std::vector<String>& strings)
  {
    ParamEntry& entry = getEntry_(key);
    //check if correct parameter type
    if (entry.value.valueType() != DataValue::STRING_VALUE && entry.value.valueType() != DataValue::STRING_LIST)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, key);
    }
    //check for commas
    for (Size i = 0; i < strings.size(); ++i)
    {
      if (strings[i].has(','))
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Comma characters in Param string restrictions are not allowed!");
      }
    }
    entry.valid_strings = strings;
  }

  void Param::setMinInt(const String& key, Int min)
  {
    ParamEntry& entry = getEntry_(key);
    if (entry.value.valueType() != DataValue::INT_VALUE && entry.value.valueType() != DataValue::INT_LIST)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, key);
    }
    entry.min_int = min;
  }

  void Param::setMaxInt(const String& key, Int max)
  {
    ParamEntry& entry = getEntry_(key);
    if (entry.value.valueType() != DataValue::INT_VALUE && entry.value.valueType() != DataValue::INT_LIST)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, key);
    }
    entry.max_int = max;
  }

  void Param::setMinFloat(const String& key, double min)
  {
    ParamEntry& entry = getEntry_(key);
    if (entry.value.valueType() != DataValue::DOUBLE_VALUE && entry.value.valueType() != DataValue::DOUBLE_LIST)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, key);
    }
    entry.min_float = min;
  }

  void Param::setMaxFloat(const String& key, double max)
  {
    ParamEntry& entry = getEntry_(key);
    if (entry.value.valueType() != DataValue::DOUBLE_VALUE && entry.value.valueType() != DataValue::DOUBLE_LIST)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, key);
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
    //static initialization and thus cannot rely on String::EMPTY been initialized.
    static String empty;

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

  void Param::insert(const String& prefix, const Param& param)
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

  void Param::setDefaults(const Param& defaults, const String& prefix, bool showMessage)
  {
    String prefix2 = prefix;
    if (prefix2 != "")
    {
      prefix2.ensureLastChar(':');
    }

    String pathname;
    for (Param::ParamIterator it = defaults.begin(); it != defaults.end(); ++it)
    {
      if (!exists(prefix2 + it.getName()))
      {
        if (showMessage)
          std::cerr << "Setting " << prefix2 + it.getName() << " to " << it->value << std::endl;
        String name = prefix2 + it.getName();
        root_.insert(ParamEntry("", it->value, it->description), name);
        //copy tags
        for (std::set<String>::const_iterator tag_it = it->tags.begin(); tag_it != it->tags.end(); ++tag_it)
        {
          addTag(name, *tag_it);
        }
        //copy restrictions
        if (it->value.valueType() == DataValue::STRING_VALUE || it->value.valueType() == DataValue::STRING_LIST)
        {
          setValidStrings(name, it->valid_strings);
        }
        else if (it->value.valueType() == DataValue::INT_VALUE || it->value.valueType() == DataValue::INT_LIST)
        {
          setMinInt(name, it->min_int);
          setMaxInt(name, it->max_int);
        }
        else if (it->value.valueType() == DataValue::DOUBLE_VALUE || it->value.valueType() == DataValue::DOUBLE_LIST)
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
        String real_pathname = pathname.chop(1); //remove ':' at the end
        if (real_pathname != "")
        {
          String description_old = getSectionDescription(prefix + real_pathname);
          String description_new = defaults.getSectionDescription(real_pathname);
          if (description_old == "")
          {
            //std::cerr << "## Setting description of " << prefix+real_pathname << " to"<< std::endl;
            //std::cerr << "## " << description_new << std::endl;
            setSectionDescription(prefix2 + real_pathname, description_new);
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
      if (node_parent != nullptr)
      {
        Param::ParamNode::NodeIterator it = node_parent->findNode(node_parent->suffix(keyname));
        if (it != node_parent->nodes.end())
        {
          String name = it->name;
          node_parent->nodes.erase(it); // will automatically delete subnodes
          if (node_parent->nodes.empty()  && node_parent->entries.empty())
          {
            // delete last section name (could be partial)
            remove(keyname.chop(name.size())); // keep last ':' to indicate deletion of a section
          }
        }
      }
    }
    else
    {
      ParamNode* node = root_.findParentOf(keyname);
      if (node != nullptr)
      {
        String entryname = node->suffix(keyname); // get everything beyond last ':'
        Param::ParamNode::EntryIterator it = node->findEntry(entryname);
        if (it != node->entries.end())
        {
          node->entries.erase(it); // delete entry
          if (node->nodes.empty()  && node->entries.empty())
          {
            // delete if section is now empty
            remove(keyname.chop(entryname.size())); // keep last ':' to indicate deletion of a section
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
      if (node != nullptr)
      {
        Param::ParamNode::NodeIterator it = node->findNode(node->suffix(prefix.chop(1)));
        if (it != node->nodes.end())
        {
          String name = it->name;
          node->nodes.erase(it); // will automatically delete subnodes
          if (node->nodes.empty()  && node->entries.empty())
          {
            // delete last section name (could be partial)
            removeAll(prefix.chop(name.size() + 1)); // '+1' for the tailing ':'
          }
        }
      }
    }
    else //we have to delete all entries and nodes starting with the prefix
    {
      ParamNode* node = root_.findParentOf(prefix);
      if (node != nullptr)
      {
        String suffix = node->suffix(prefix); // name behind last ":"

        for (Param::ParamNode::NodeIterator it = node->nodes.begin(); it != node->nodes.end(); /*do nothing*/)
        {
          if (it->name.hasPrefix(suffix))
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
          if (it->name.hasPrefix(suffix))
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
          removeAll(prefix.chop(suffix.size()));
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

  Param Param::copy(const String& prefix, bool remove_prefix) const
  {
    ParamNode out("ROOT", "");

    ParamNode* node = root_.findParentOf(prefix);
    if (node == nullptr)
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
        out.insert(*node, prefix.chop(node->name.size() + 1));
      }
    }
    else //we have to copy all entries and nodes starting with the right suffix
    {
      String suffix = node->suffix(prefix);
      for (Param::ParamNode::NodeIterator it = node->nodes.begin(); it != node->nodes.end(); ++it)
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
            out.insert(*it, prefix.chop(suffix.size()));
          }
        }
      }
      for (Param::ParamNode::EntryIterator it = node->entries.begin(); it != node->entries.end(); ++it)
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
            out.insert(*it, prefix.chop(suffix.size()));
          }
        }
      }
    }

    return Param(out);
  }

  void Param::parseCommandLine(const int argc, const char** argv, const String& prefix)
  {
    //determine prefix
    String prefix2 = prefix;
    if (prefix2 != "")
    {
      prefix2.ensureLastChar(':');
    }

    //parse arguments
    String arg, arg1;
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
        arg_is_option = true;
      bool arg1_is_option = false;
      if (arg1.size() >= 2 && arg1[0] == '-' && arg1[1] != '0' && arg1[1] != '1' && arg1[1] != '2' && arg1[1] != '3' && arg1[1] != '4' && arg1[1] != '5' && arg1[1] != '6' && arg1[1] != '7' && arg1[1] != '8' && arg1[1] != '9')
        arg1_is_option = true;

      //cout << "Parse: '"<< arg << "' '" << arg1 << "'" << std::endl;

      //flag (option without text argument)
      if (arg_is_option && arg1_is_option)
      {
        root_.insert(ParamEntry(arg, String(), ""), prefix2);
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
          StringList sl;
          sl.push_back(arg);
          // create "misc"-Node:
          root_.insert(ParamEntry("misc", sl, ""), prefix2);
        }
        else
        {
          StringList sl = misc_entry->value;
          sl.push_back(arg);
          misc_entry->value = sl;
        }
      }
    }
  }

  void Param::parseCommandLine(const int argc, const char** argv, const Map<String, String>& options_with_one_argument, const Map<String, String>& options_without_argument, const Map<String, String>& options_with_multiple_argument, const String& misc, const String& unknown)
  {
    //determine misc key
    String misc_key = misc;

    //determine unknown key
    String unknown_key = unknown;

    //parse arguments
    String arg, arg1;
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
        arg_is_option = true;
      bool arg1_is_option = false;
      if (arg1.size() >= 2 && arg1[0] == '-' && arg1[1] != '0' && arg1[1] != '1' && arg1[1] != '2' && arg1[1] != '3' && arg1[1] != '4' && arg1[1] != '5' && arg1[1] != '6' && arg1[1] != '7' && arg1[1] != '8' && arg1[1] != '9')
        arg1_is_option = true;


      //with multiple argument
      if (options_with_multiple_argument.has(arg))
      {
        //next argument is an option
        if (arg1_is_option)
        {
          root_.insert(ParamEntry("", StringList(), ""), options_with_multiple_argument.find(arg)->second);
        }
        //next argument is not an option
        else
        {
          StringList sl;
          int j = (i + 1);
          while (j < argc && !(arg1.size() >= 2 && arg1[0] == '-' && arg1[1] != '0' && arg1[1] != '1' && arg1[1] != '2' && arg1[1] != '3' && arg1[1] != '4' && arg1[1] != '5' && arg1[1] != '6' && arg1[1] != '7' && arg1[1] != '8' && arg1[1] != '9'))
          {
            sl.push_back(arg1);
            ++j;
            if (j < argc)
              arg1 = argv[j];
          }

          root_.insert(ParamEntry("", sl, ""), options_with_multiple_argument.find(arg)->second);
          i = j - 1;
        }
      }
      //without argument
      else if (options_without_argument.has(arg))
      {
        root_.insert(ParamEntry("", String("true"), ""), options_without_argument.find(arg)->second);
      }
      //with one argument
      else if (options_with_one_argument.has(arg))
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

          root_.insert(ParamEntry("", String(), ""), options_with_one_argument.find(arg)->second);
        }
      }
      //unknown option
      else if (arg_is_option)
      {
        ParamEntry* unknown_entry = root_.findEntryRecursive(unknown);
        if (unknown_entry == nullptr)
        {
          StringList sl;
          sl.push_back(arg);
          root_.insert(ParamEntry("", sl, ""), unknown);
        }
        else
        {
          StringList sl = unknown_entry->value;
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
          StringList sl;
          sl.push_back(arg);
          // create "misc"-Node:
          root_.insert(ParamEntry("", sl, ""), misc);
        }
        else
        {
          StringList sl = misc_entry->value;
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
      String prefix = it.getName().chop(it->name.size() + 1);
      if (prefix != "")
      {
        prefix += "|";
      }
      os << '"' << prefix << it->name << "\" -> \"" << it->value << '"';
      if (it->description != "")
      {
        os << " (" << it->description << ")";
      }
      os << std::endl;
    }
    return os;
  }

  Size Param::size() const
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

  void Param::checkDefaults(const String& name, const Param& defaults, const String& prefix) const
  {
    //Extract right parameters
    String prefix2 = prefix;
    if (prefix2 != "")
    {
      prefix2.ensureLastChar(':');
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
          OPENMS_LOG_WARN << " in '" << prefix2 << "'";
        OPENMS_LOG_WARN << "!" << std::endl;
      }

      //different types
      ParamEntry* default_value = defaults.root_.findEntryRecursive(prefix2 + it.getName());
      if (default_value == nullptr)
        continue;
      if (default_value->value.valueType() != it->value.valueType())
      {
        String d_type;
        if (default_value->value.valueType() == DataValue::STRING_VALUE)
          d_type = "string";
        if (default_value->value.valueType() == DataValue::STRING_LIST)
          d_type = "string list";
        if (default_value->value.valueType() == DataValue::EMPTY_VALUE)
          d_type = "empty";
        if (default_value->value.valueType() == DataValue::INT_VALUE)
          d_type = "integer";
        if (default_value->value.valueType() == DataValue::INT_LIST)
          d_type = "integer list";
        if (default_value->value.valueType() == DataValue::DOUBLE_VALUE)
          d_type = "float";
        if (default_value->value.valueType() == DataValue::DOUBLE_LIST)
          d_type = "float list";
        String p_type;
        if (it->value.valueType() == DataValue::STRING_VALUE)
          p_type = "string";
        if (it->value.valueType() == DataValue::STRING_LIST)
          p_type = "string list";
        if (it->value.valueType() == DataValue::EMPTY_VALUE)
          p_type = "empty";
        if (it->value.valueType() == DataValue::INT_VALUE)
          p_type = "integer";
        if (it->value.valueType() == DataValue::INT_LIST)
          p_type = "integer list";
        if (it->value.valueType() == DataValue::DOUBLE_VALUE)
          p_type = "float";
        if (it->value.valueType() == DataValue::DOUBLE_LIST)
          p_type = "float list";

        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name + ": Wrong parameter type '" + p_type + "' for " + d_type + " parameter '" + it.getName() + "' given!");
      }
      //parameter restrictions
      ParamEntry pe = *default_value;
      pe.value = it->value;
      String s;
      if (!pe.isValid(s))
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name + ": " + s);
    }
  }

  Param::ParamIterator Param::findFirst(const String& leaf) const
  {
    for (Param::ParamIterator it = this->begin(); it != this->end(); ++it)
    {
      if (it.getName().hasSuffix(String(":") + leaf))
      {
        return it;
      }
    }
    return this->end();
  }

  Param::ParamIterator Param::findNext(const String& leaf, const ParamIterator& start_leaf) const
  {
    // start at NEXT entry
    Param::ParamIterator it = start_leaf;
    if (it != this->end()) ++it;

    for (; it != this->end(); ++it)
    {
      if (it.getName().hasSuffix(String(":") + leaf))
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
      String target_name; // fully qualified name in new param

      if (this->exists(it.getName()))
      {
        // param 'version': do not override!
        if (it.getName().hasSuffix(":version"))
        {
          if (this->getValue(it.getName()) != it->value)
          {
OPENMS_THREAD_CRITICAL(oms_log)
            stream << "Warning: for ':version' entry, augmented and Default Ini-File differ in value. Default value will not be altered!\n";
          }
          continue;
        }
        // param 'type': do not override!
        else if (it.getName().hasSuffix(":type") &&
                 it.getName().toQString().count(':') == 2) // only for TOPP type (e.g. PeakPicker:1:type), any other 'type' param is ok
        {
          if (this->getValue(it.getName()) != it->value)
          {
OPENMS_THREAD_CRITICAL(oms_log)
            stream << "Warning: for ':type' entry, augmented and Default Ini-File differ in value. Default value will not be altered!\n";
          }
          continue;
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
OPENMS_THREAD_CRITICAL(oms_log)
            stream << "Found '" << it.getName() << "' as '" << it_match.getName() << "' in new param." << std::endl;
            new_entry = this->getEntry(it_match.getName());
            target_name = it_match.getName();
          }
        }

        if (target_name.empty()) // no mapping was found
        {
          if (fail_on_unknown_parameters)
          {
OPENMS_THREAD_CRITICAL(oms_log)
            stream << "Unknown (or deprecated) Parameter '" << it.getName() << "' given in outdated parameter file!" << std::endl;
            is_update_success = false;
          }
          else if (add_unknown)
          {
OPENMS_THREAD_CRITICAL(oms_log)
            stream << "Unknown (or deprecated) Parameter '" << it.getName() << "' given in outdated parameter file! Adding to current set." << std::endl;
            Param::ParamEntry local_entry = p_outdated.getEntry(it.getName());
            String prefix = "";
            if (it.getName().has(':'))
            {
              prefix = it.getName().substr(0, 1 + it.getName().find_last_of(':'));
            }
            this->root_.insert(local_entry, prefix); //->setValue(it.getName(), local_entry.value, local_entry.description, local_entry.tags);
          }
          else
          {
OPENMS_THREAD_CRITICAL(oms_log)
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
          DataValue default_value = new_entry.value;
          new_entry.value = it->value;
          String validation_result;
          if (new_entry.isValid(validation_result))
          {
            // overwrite default value
            if (verbose) 
            {
OPENMS_THREAD_CRITICAL(oms_log)
                stream << "Default-Parameter '" << target_name << "' overridden: '" << default_value << "' --> '" << it->value << "'!" << std::endl;
            }
            this->setValue(target_name, it->value, new_entry.description, this->getTags(target_name));
          }
          else
          {
OPENMS_THREAD_CRITICAL(oms_log)
            stream << validation_result;
            if (fail_on_invalid_values)
            {
OPENMS_THREAD_CRITICAL(oms_log)
              stream << " Updating failed!" << std::endl;
              is_update_success = false;
            }
            else
            {
OPENMS_THREAD_CRITICAL(oms_log)
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
OPENMS_THREAD_CRITICAL(oms_log)
        stream << "Parameter '" << it.getName() << "' has changed value type!\n";
        if (fail_on_invalid_values)
        {
OPENMS_THREAD_CRITICAL(oms_log)
          stream << " Updating failed!" << std::endl;
          is_update_success = false;
        } 
        else
        {
OPENMS_THREAD_CRITICAL(oms_log)
          stream << " Ignoring invalid value (using new default)!" << std::endl;
        }
      }

    } // next param in outdated tree

    return is_update_success;
  }

  void Param::merge(const OpenMS::Param& toMerge)
  {
    // keep track of the path inside the param tree
    String pathname;

    // augment
    for (Param::ParamIterator it = toMerge.begin(); it != toMerge.end(); ++it)
    {
      String prefix = "";
      if (it.getName().has(':'))
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
          if (pathname.hasSuffix(traceIt->name + ":"))
            pathname.resize(pathname.size() - traceIt->name.size() - 1);
        }
        String real_pathname = pathname.chop(1); //remove ':' at the end
        if (real_pathname != "")
        {
          String description_old = getSectionDescription(prefix + real_pathname);
          String description_new = toMerge.getSectionDescription(real_pathname);
          if (description_old == "")
          {
            setSectionDescription(real_pathname, description_new);
          }
        }
      }

    }
  }

  void Param::setSectionDescription(const String& key, const String& description)
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

  void Param::addSection(const String& key, const String& description)
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
    if (root_ == nullptr)
      return *this;

    trace_.clear();
    while (true)
    {
      const Param::ParamNode* node = stack_.back();

      //cout << "############ operator++ #### " << node->name << " ## " << current_ <<endl;

      //check if there is a next entry in the current node
      if (current_ + 1 < (Int)node->entries.size())
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
        trace_.push_back(TraceInfo(node->nodes[0].name, node->nodes[0].description, true));

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
          //cout << " - stack size: " << stack_.size() << endl;
          //we have reached the end
          if (stack_.empty())
          {
            //cout << " - reached the end" << endl;
            root_ = nullptr;
            return *this;
          }
          node = stack_.back();

          //cout << " - last was: " << last->name << endl;
          //cout << " - descended to: " << node->name << endl;

          //track changes (leave a node)
          trace_.push_back(TraceInfo(last->name, last->description, false));

          //check of new subtree is accessible
          UInt next_index = (last - &(node->nodes[0])) + 1;
          if (next_index < node->nodes.size())
          {
            current_ = -1;
            stack_.push_back(&(node->nodes[next_index]));
            //cout << " - entering into: " << node->nodes[next_index].name  << endl;
            //track changes (enter a node)
            trace_.push_back(TraceInfo(node->nodes[next_index].name, node->nodes[next_index].description, true));
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

  String Param::ParamIterator::getName() const
  {
    String tmp;
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
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Param tags may not contain comma characters", tag);
    }
    getEntry_(key).tags.insert(tag);
  }

  void Param::addTags(const String& key, const StringList& tags)
  {
    ParamEntry& entry = getEntry_(key);
    for (Size i = 0; i != tags.size(); ++i)
    {
      if (tags[i].has(','))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Param tags may not contain comma characters", tags[i]);
      }
      entry.tags.insert(tags[i]);
    }
  }

  StringList Param::getTags(const String& key) const
  {
    ParamEntry& entry = getEntry_(key);
    StringList list;
    for (std::set<String>::const_iterator it = entry.tags.begin(); it != entry.tags.end(); ++it)
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
    if (entry == nullptr)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, key);
    }

    return *entry;
  }

} //namespace
