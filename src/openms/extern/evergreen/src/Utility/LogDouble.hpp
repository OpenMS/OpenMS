#ifndef _LOGDOUBLE_H
#define _LOGDOUBLE_H

#include <math.h>
#include <limits>
#include <iostream>
#include <assert.h>

#include "ComparableMixin.hpp"
#include "sign.hpp"

// TODO: test for compatibility with Vector<LogDouble>, which will
// initialize sign to zero, log_absolute_value to zero.

class LogDouble : public ComparableMixin<LogDouble> {
private:
  signed char _sign;
  double _log_absolute_value;
  
  static double logaddexp(double log_a, double log_b) {
    // returns the log( exp(log_a) + exp(log_b) )

    // if both are infinite, taking the difference will result in a NaN;
    // simply return infinity
    if (isinf(log_a) && isinf(log_b))
      return log_a;

    if ( log_a > log_b )
      return logaddexp_first_larger(log_a, log_b);
    else
      return logaddexp_first_larger(log_b, log_a);
  }
  static double logaddexp_first_larger(double log_a, double log_b) {
    // Note: this check is for internal testing, and can be removed for greater speed
    assert( log_a >= log_b );
    if ( log_a == -std::numeric_limits<double>::infinity() )
      return log_b;
    return log1p( exp(log_b-log_a) ) + log_a;
  }
  static double logsubabsexp(double log_a, double log_b) {
    // returns the log( abs( exp(log_a) - exp(log_b) ) )
    if ( log_a > log_b )
      return logsubexp_first_larger(log_a, log_b);
    else
      return logsubexp_first_larger(log_b, log_a);
  }
  static double logsubexp_first_larger(double log_a, double log_b) {
    // Note: this check is for internal testing, and can be removed for greater speed
    assert(log_a >= log_b);
    if ( log_a == -std::numeric_limits<double>::infinity() )
      return log_b;
    return log1p( -exp(log_b-log_a) ) + log_a;
  }
public:
  LogDouble() {
    _log_absolute_value = std::numeric_limits<double>::quiet_NaN();
    _sign = 1;
  }
  explicit LogDouble(double x) {
    if ( x == 0.0 )
      _sign = 1;
    else
      _sign = (signed char) ::sign(x);
    
    _log_absolute_value = log(fabs(x));
  }

  static LogDouble create_from_log_absolute_value(double log_absolute_value_param) {
    LogDouble result;
    result._log_absolute_value = log_absolute_value_param;
    result._sign = 1;
    return result;
  }

  // +=, -=, *=, /=
  const LogDouble & operator +=(const LogDouble & rhs) {
    // if the signs are the same, simply use logaddexp
    if (_sign == rhs._sign)
      _log_absolute_value = logaddexp(_log_absolute_value, rhs._log_absolute_value);
    else
      {
	double new_log_absolute_value = logsubabsexp(_log_absolute_value, rhs._log_absolute_value);
	if ( _log_absolute_value < rhs._log_absolute_value )
	  // *this "loses"
	  _sign *= -1;
	_log_absolute_value = new_log_absolute_value;
      }
    return *this;
  }
  const LogDouble & operator -=(const LogDouble & rhs) {
    // if the signs are different, they will be the same after negation
    if (_sign != rhs._sign)
      _log_absolute_value = logaddexp(_log_absolute_value, rhs._log_absolute_value);
    else
      {
	double new_log_absolute_value = logsubabsexp(_log_absolute_value, rhs._log_absolute_value);
	if ( _log_absolute_value < rhs._log_absolute_value )
	  // *this "loses"
	  _sign *= -1;
	_log_absolute_value = new_log_absolute_value;
      }
    return *this;
  }
  const LogDouble & operator *=(const LogDouble & rhs) {
    _sign *= rhs._sign;
    _log_absolute_value += rhs._log_absolute_value;
    return *this;
  }
  const LogDouble & operator /=(const LogDouble & rhs) {
    _sign *= rhs._sign;
    _log_absolute_value -= rhs._log_absolute_value;
    return *this;
  }
  LogDouble operator -() const {
    LogDouble result(*this);
    result._sign = -result._sign;
    return result;
  }

  explicit operator double() const {
    return _sign * exp(_log_absolute_value);
  }

  double log_absolute_value() const {
    return _log_absolute_value;
  }

  double sign() const {
    return _sign;
  }

  bool operator <(LogDouble rhs) const {
    return (_sign < rhs._sign) || ( _sign == rhs._sign && ( (_sign == 1 && _log_absolute_value < rhs._log_absolute_value) || (_sign == -1 && _log_absolute_value > rhs._log_absolute_value) ) );
  }

  bool operator ==(LogDouble rhs) const {
    // magnitudes must be equal, and if nonzero, signs must be equal (if
    // it's zero, disregard sign)
    return _log_absolute_value == rhs._log_absolute_value && (_sign == rhs._sign || (_log_absolute_value == -std::numeric_limits<double>::infinity() ) );
  }

  static bool is_nan(LogDouble x) {
    return isnan(double(x));
  }

  static bool is_inf(LogDouble x) {
    return isinf(x.log_absolute_value());
  }

  friend LogDouble exp(LogDouble rhs) {
    rhs._log_absolute_value = double(rhs);
    rhs._sign = 1;
    return rhs;
  }
  
  friend std::ostream & operator <<(std::ostream & os, LogDouble rhs) {
    if (rhs._sign == -1)
      os << '-';
    os << "exp(" << rhs._log_absolute_value << ")~" << double(rhs);
    return os;
  }
};

LogDouble operator +(LogDouble lhs, LogDouble rhs) {
  lhs += rhs;
  return lhs;
}

LogDouble operator -(LogDouble lhs, LogDouble rhs) {
  lhs -= rhs;
  return lhs;
}

LogDouble operator *(LogDouble lhs, LogDouble rhs) {
  lhs *= rhs;
  return lhs;
}

LogDouble operator /(LogDouble lhs, LogDouble rhs) {
  lhs /= rhs;
  return lhs;
}

LogDouble pow(LogDouble lhs, LogDouble rhs) {
  assert(lhs.sign() >= 0);

  // for all x, x^0 --> 1
  if ( rhs.log_absolute_value() == -std::numeric_limits<double>::infinity() )
    return LogDouble(1.0);
  // for all y>0 (guaranteed by previous check), 0^y --> 0
  if ( lhs.log_absolute_value() == -std::numeric_limits<double>::infinity() )
    return LogDouble(0.0);
  return LogDouble::create_from_log_absolute_value(double(rhs) * lhs.log_absolute_value());
}

LogDouble fabs(LogDouble x) {
  return LogDouble::create_from_log_absolute_value(x.log_absolute_value());
}

#endif
