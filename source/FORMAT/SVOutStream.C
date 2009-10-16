#include <OpenMS/FORMAT/SVOutStream.h>
#include <limits>

using namespace OpenMS;
using namespace std;

SVOutStream::SVOutStream(ostream& out, const String& sep,
												 const String& replacement,
												 String::QuotingMethod quoting):
	ostream(out.rdbuf()), sep_(sep), replacement_(replacement),	nan_("nan"),
	quoting_(quoting), modify_strings_(true), newline_(true)
{
	// use high decimal precision (appropriate for double):
	precision(std::numeric_limits<double>::digits10);
}


SVOutStream& SVOutStream::operator<<(String str)
{
	if (str.find('\n') != string::npos)
		throw Exception::IllegalArgument(
			__FILE__, __LINE__, __PRETTY_FUNCTION__,
			"argument must not contain newline characters");
	if (!newline_) (ostream&)*this << sep_;
	else newline_ = false;
	if (!modify_strings_)
		(ostream&)*this << str;
	else if (quoting_ != String::NONE)
	{
		(ostream&)*this << str.quote('"', quoting_);
	}
	else (ostream&)*this << str.substitute(sep_, replacement_);
	return *this;
}


SVOutStream& SVOutStream::operator<<(const char* c_str)
{
	return operator<<(String(c_str));
}


SVOutStream& SVOutStream::operator<<(const char c)
{
	return operator<<(String(c));
}


SVOutStream& SVOutStream::operator<<(const string& str)
{
	return operator<<((String&)str);
}


SVOutStream& SVOutStream::operator<<(ostream& (*fp)(ostream&))
{
	if (fp == static_cast<ostream&(*)(ostream&)>(endl))
		newline_ = true;
	(ostream&)*this << fp;
	return *this;
}


SVOutStream& SVOutStream::write(const String& str)
{
	ostream::write(str.c_str(), str.size());
	return *this;
}


bool SVOutStream::modifyStrings(bool modify)
{
	bool old = modify_strings_;
	modify_strings_ = modify;
	return old;
}
