
/*
 Copyright (C) 2006, 2007 Eric Ehlers

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#ifndef qla_calc_conversions_hpp
#define qla_calc_conversions_hpp

#include <ql/time/date.hpp>
#include <ql/time/calendar.hpp>
#include <ql/math/matrix.hpp>
#include <Addins/Calc/qldefs.hpp>
#include <sal/types.h>
#include <vector>

// temp fix pending Calc implementation of CoerceQuoteHandle
#include <qlo/quotes.hpp>

void calcToScalar(QuantLib::Date &, const sal_Int32&);
void calcToScalar(QuantLib::Date &, const ANY&);
void calcToScalar(QuantLib::Calendar &, const STRING &id);
void calcToScalar(QuantLib::Period &, const STRING &id);
void calcToVector(std::vector<QuantLib::Date> &, const SEQSEQ(sal_Int32) &);
void calcToVector(QuantLib::Array &, const SEQSEQ(double) &);
void calcToVector(std::vector<std::string> &, const SEQSEQ(ANY) &);
void calcToVector(std::vector<long> &, const SEQSEQ(sal_Int32) &);
void calcToVector(std::vector<bool> &, const SEQSEQ(sal_Int32) &);
void calcToVector(std::vector<QuantLib::Period> &, const SEQSEQ(ANY) &);
void calcToVector(std::vector<boost::any> &, const SEQSEQ(ANY) &);
QuantLib::Matrix calcToQlMatrix(const SEQSEQ(double) &);

void scalarToCalc(sal_Int32 &, const QuantLib::Date &);
void scalarToCalc(double &, const QuantLib::Real &);
void vectorToCalc(SEQSEQ(sal_Int32) &, const std::vector<QuantLib::Date> &);
void vectorToCalc(SEQSEQ(double) &, const QuantLib::Array &);
void matrixToCalc(SEQSEQ(double) &, const QuantLib::Matrix &);

#endif
