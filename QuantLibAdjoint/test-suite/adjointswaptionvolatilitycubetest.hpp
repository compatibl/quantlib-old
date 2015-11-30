/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2006 Ferdinando Ametrano
Copyright (C) 2006 Katiuscia Manzoni
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

#ifndef cl_adjoint_swaption_volatility_cube_hpp
#define cl_adjoint_swaption_volatility_cube_hpp
#pragma once

#include <boost/test/unit_test.hpp>

class AdjointSwaptionVolatilityCubeTest {
public:
    //test Implied volatility derivative on strike
    static bool testImpliedVolatilitiesDiffStrike();
    //test Implied volatility derivative on input to volatility cube volatility
    static bool testImpliedVolatilitiesDiffVols();
    //test Black Variance derivative on input to volatility cube volatility
    static bool testBlackVariancesDiffVols();
    //test Black Variance derivative on strike
    static bool testBlackVarianceDiffStrike();

    static boost::unit_test_framework::test_suite* suite();
};

#endif
