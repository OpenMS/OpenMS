// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <OpenMS/DATASTRUCTURES/Param.h>

using namespace std;

namespace OpenMS
{
  DefaultParamHandler::DefaultParamHandler(const String & name) :
    param_(),
    defaults_(),
    subsections_(),
    error_name_(name),
    check_defaults_(true),
    warn_empty_defaults_(true)
  {

  }

  DefaultParamHandler::DefaultParamHandler(const DefaultParamHandler & rhs) :
    param_(rhs.param_),
    defaults_(rhs.defaults_),
    subsections_(rhs.subsections_),
    error_name_(rhs.error_name_),
    check_defaults_(rhs.check_defaults_),
    warn_empty_defaults_(rhs.warn_empty_defaults_)
  {
  }

  DefaultParamHandler & DefaultParamHandler::operator=(const DefaultParamHandler & rhs)
  {
    if (&rhs == this)
      return *this;

    //copy members
    param_ = rhs.param_;
    defaults_ = rhs.defaults_;
    subsections_ = rhs.subsections_;
    error_name_ = rhs.error_name_;
    check_defaults_ = rhs.check_defaults_;
    warn_empty_defaults_ = rhs.warn_empty_defaults_;

    return *this;
  }

  bool DefaultParamHandler::operator==(const DefaultParamHandler & rhs) const
  {
    return param_ == rhs.param_ &&
           defaults_ == rhs.defaults_ &&
           subsections_ == rhs.subsections_ &&
           error_name_ == rhs.error_name_ &&
           check_defaults_ == rhs.check_defaults_ &&
           warn_empty_defaults_ == rhs.warn_empty_defaults_;
  }

  DefaultParamHandler::~DefaultParamHandler()
  {
  }

  void DefaultParamHandler::setParameters(const Param & param)
  {
    //set defaults and apply new parameters
    Param tmp(param);
    tmp.setDefaults(defaults_);
    param_ = tmp;

    if (check_defaults_)
    {
      if (defaults_.empty() && warn_empty_defaults_)
      {
        LOG_WARN << "Warning: No default parameters for DefaultParameterHandler '" << error_name_ << "' specified!" << endl;
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
      if (it->description == "")
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

  const Param & DefaultParamHandler::getParameters() const
  {
    return param_;
  }

  const Param & DefaultParamHandler::getDefaults() const
  {
    return defaults_;
  }

  const String & DefaultParamHandler::getName() const
  {
    return error_name_;
  }

  void DefaultParamHandler::setName(const String & name)
  {
    error_name_ = name;
  }

  const std::vector<String> & DefaultParamHandler::getSubsections() const
  {
    return subsections_;
  }

} // namespace OpenMS
