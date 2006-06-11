// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: TimeStamp.sip,v 1.1.1.1 2005/04/04 13:03:04 marc_sturm Exp $
//


class PreciseTime
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/TimeStamp.h>
	using namespace OpenMS;
%End
	public:

	PreciseTime()	throw();
	PreciseTime(const PreciseTime& /time/) throw();
	PreciseTime(long /secs/, long /usecs/) throw();
	virtual ~PreciseTime() throw();

	static const PreciseTime ZERO;
	void set(long /secs/, long /usecs/) throw();
	void set(const PreciseTime& /time/) throw();
	//	const PreciseTime& operator = (const PreciseTime& /time/) throw();
	virtual void clear() throw();

	bool operator < (const PreciseTime& /time/) const throw();
	bool operator > (const PreciseTime& /time/) const throw();
	bool operator == (const PreciseTime& /time/) const throw();

	long getSeconds() const throw();
	long getMicroSeconds() const throw();
	static PreciseTime now() throw();
};

class TimeStamp
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/TimeStamp.h>
	using namespace OpenMS;
%End
	public:

	TimeStamp() throw();
	virtual ~TimeStamp() throw();

	bool isNewerThan(const PreciseTime& /time/) const throw();
	bool isOlderThan(const PreciseTime& /time/) const throw();
	bool isNewerThan(const TimeStamp& /stamp/) const throw();
	bool isOlderThan(const TimeStamp& /stamp/) const throw();
	bool operator == (const TimeStamp& /stamp/) const throw();

	bool operator < (const TimeStamp& /stamp/) const throw();
	bool operator > (const TimeStamp& /stamp/) const throw();

	virtual void stamp(const PreciseTime& /time/ = PreciseTime::ZERO) throw();
	const PreciseTime& getTime() const throw();

	//	const PreciseTime& operator = (const PreciseTime& /time/) throw();
	virtual void clear() throw();
};

// ??? __str__ __repr__ missing!
// std::ostream& operator << (std::ostream& os, const PreciseTime& /time/) throw();

// std::ostream& operator << (std::ostream& os, const TimeStamp& stamp) throw();

