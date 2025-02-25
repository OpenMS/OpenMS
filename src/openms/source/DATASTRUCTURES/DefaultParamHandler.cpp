// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

using namespace std;

namespace OpenMS
{
  DefaultParamHandler::DefaultParamHandler(const String& name) :
    param_(),
    defaults_(),
    subsections_(),
    error_name_(name),
    check_defaults_(true),
    warn_empty_defaults_(true)
  {

  }

  DefaultParamHandler::DefaultParamHandler(const DefaultParamHandler& rhs) = default;

  DefaultParamHandler& DefaultParamHandler::operator=(const DefaultParamHandler& rhs)
  {
    if (&rhs == this)
    {
      return *this;
    }
    //copy members
    param_ = rhs.param_;
    defaults_ = rhs.defaults_;
    subsections_ = rhs.subsections_;
    error_name_ = rhs.error_name_;
    check_defaults_ = rhs.check_defaults_;
    warn_empty_defaults_ = rhs.warn_empty_defaults_;

    return *this;
  }

  bool DefaultParamHandler::operator==(const DefaultParamHandler& rhs) const
  {
    return param_ == rhs.param_ &&
           defaults_ == rhs.defaults_ &&
           subsections_ == rhs.subsections_ &&
           error_name_ == rhs.error_name_ &&
           check_defaults_ == rhs.check_defaults_ &&
           warn_empty_defaults_ == rhs.warn_empty_defaults_;
  }

  DefaultParamHandler::~DefaultParamHandler() = default;

  void DefaultParamHandler::setParameters(const Param& param)
  {
    //set defaults and apply new parameters
    Param tmp(param);
    tmp.setDefaults(defaults_);
    param_ = tmp;

    if (check_defaults_)
    {
      if (defaults_.empty() && warn_empty_defaults_)
      {
        OPENMS_LOG_WARN << "Warning: No default parameters for DefaultParameterHandler '" << error_name_ << "' specified!" << endl;
      }

      //remove registered subsections
      for (vector<String>::const_iterator it = subsections_.begin(); it != subsections_.end(); ++it)
      {
        tmp.removeAll(*it + ':');
      }

      //check defaults
      tmp.checkDefaults(error_name_, defaults_);
    }

    //do necessary changes to other member variables
    updateMembers_();
  }

  void DefaultParamHandler::defaultsToParam_()
  {
    //check if a description is given for all defaults
    bool description_missing = false;
    String missing_parameters;
    for (Param::ParamIterator it = defaults_.begin(); it != defaults_.end(); ++it)
    {
      //cout << "Name: " << it->getName() << endl;
      if (it->description.empty())
      {
        description_missing = true;
        missing_parameters += it.getName() + ",";
        break;
      }
    }
    if (description_missing)
    {
      cerr << "Warning: no default parameter description for parameters '" << missing_parameters << "' of DefaultParameterHandler '" << error_name_ << "' given!" << endl;
    }
    param_.setDefaults(defaults_);
    updateMembers_();
  }

  void DefaultParamHandler::updateMembers_()
  {

  }

  const Param& DefaultParamHandler::getParameters() const
  {
    return param_;
  }

  const Param& DefaultParamHandler::getDefaults() const
  {
    return defaults_;
  }

  const String& DefaultParamHandler::getName() const
  {
    return error_name_;
  }

  void DefaultParamHandler::setName(const String& name)
  {
    error_name_ = name;
  }

  const std::vector<String>& DefaultParamHandler::getSubsections() const
  {
    return subsections_;
  }

  void DefaultParamHandler::writeParametersToMetaValues(const Param& write_this, MetaInfoInterface& write_here, const String& prefix)
  {
    String prefix_(prefix);
    if (!prefix_.empty())
    {
      if (prefix_.compare(prefix_.size() - 1, 1, ":") != 0) // ends with colon?
      {
        prefix_ += ":";
      }
    }
    for (auto it = write_this.begin(); it != write_this.end(); it++)
    {
      write_here.setMetaValue(prefix_ + (*it).name, DataValue((*it).value));
    }
  }

} // namespace OpenMS
