// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/MetaInfoDescription.h>

using namespace std;

namespace OpenMS
{
	
	MetaInfoDescription::MetaInfoDescription():
		MetaInfoInterface(),
		comment_(),
		name_(),
		data_processing_()
	{
	  
	}
	
	MetaInfoDescription::MetaInfoDescription(const MetaInfoDescription& source):
		MetaInfoInterface(source),
	  comment_(source.comment_),
	  name_(source.name_),
		data_processing_(source.data_processing_)
	{
	  
	}
	
	MetaInfoDescription::~MetaInfoDescription()
	{
	  
	}
	
	MetaInfoDescription& MetaInfoDescription::operator=(const MetaInfoDescription& source)
	{
	  if (&source == this) return *this;
	  
	  MetaInfoInterface::operator=(source);
	  comment_ = source.comment_;
	  name_ = source.name_;
		data_processing_ = source.data_processing_;
	  
	  return *this;
	}

  bool MetaInfoDescription::operator==(const MetaInfoDescription& rhs) const
  {
  	return 
		  comment_ == rhs.comment_ &&
		  name_ == rhs.name_ &&
		  data_processing_ == rhs.data_processing_ &&
  		MetaInfoInterface::operator==(rhs)
  		;
  }

	void MetaInfoDescription::setName(const String& name)
	{
	  name_ = name; 
	}

  const String& MetaInfoDescription::getName() const
  {
  	return name_;
  }
	const vector<DataProcessing>& MetaInfoDescription::getDataProcessing() const 
	{
	  return data_processing_; 
	}
	
	vector<DataProcessing>&  MetaInfoDescription::getDataProcessing()
	{
	  return data_processing_; 
	}
	
	void MetaInfoDescription::setDataProcessing(const vector<DataProcessing>& processing_method)
	{
	  data_processing_ = processing_method; 
	}

}


