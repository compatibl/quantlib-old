/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2008, 2009 StatPro Italia srl
Copyright (C) 2009 Ferdinando Ametrano
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

#ifndef cl_adjoint_default_curves_hpp
#define cl_adjoint_default_curves_hpp
#pragma once

#include <boost/test/unit_test.hpp>

class AdjointDefaultProbabilityCurveTest
{
public:
    static bool testFlatHazardRate();
    static bool testFlatHazardRateTime();
    static bool testFlatHazardConsistency();
    static bool testFlatDensityConsistency();
    static bool testLinearDensityConsistency();
    static bool testLogLinearSurvivalConsistency();
    static boost::unit_test_framework::test_suite* suite();
};


#endif