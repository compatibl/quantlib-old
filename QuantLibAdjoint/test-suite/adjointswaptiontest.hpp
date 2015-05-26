/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2003 RiskMap srl
Copyright (C) 2006, 2007 Ferdinando Ametrano
Copyright (C) 2006 Marco Bianchetti
Copyright (C) 2006 Cristina Duminuco
Copyright (C) 2007, 2008 StatPro Italia srl
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

#ifndef adjoint_swaption_hpp
#define adjoint_swaption_hpp

#include <boost/test/unit_test.hpp>

class AdjointSwaptionTest
{
public:
    static bool testSpreadDependency();
    static bool testStrikeDependency();
    static bool testSpreadTreatment();
    static bool testCachedValue();
    static bool testCashSettledSwaptions();
    static boost::unit_test_framework::test_suite* suite();
};

#endif 