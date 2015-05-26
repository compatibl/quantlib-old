/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2007 Marco Bianchetti
Copyright (C) 2007 Giorgio Facchinetti
Copyright (C) 2006 Chiara Fornarola
Copyright (C) 2005 StatPro Italia srl
Copyright (C) 2013 Peter Caspers
Copyright (C) 2015 CompatibL

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

#ifndef ad_hull_white_calibration_hpp
#define ad_hull_white_calibration_hpp

#include <boost/test/unit_test.hpp>

class AdjointShortRateModelsTest{
public:
    static bool testCachedHullWhite();
    static bool testCachedHullWhiteFixedReversion();
    static bool testFuturesConvexityBias();
    static boost::unit_test_framework::test_suite* suite();
};

#endif