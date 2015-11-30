/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
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

#ifndef cl_adjoint_complex_differentiation_hpp
#define cl_adjoint_complex_differentiation_hpp
#pragma once

#include <boost/test/unit_test.hpp>

class AdjointComplexTest
{
public:
    static bool testConstructor();
    static bool testAssign();
    static bool testImag();
    static bool testReal();
    static bool testConj();
    static bool testAdd();
    static bool testSub();
    static bool testMul();
    static bool testDiv();
    static bool testAddReal();
    static bool testSubReal();
    static bool testMulReal();
    static bool testDivReal();
    static bool testArg();
    static bool testPolar();
    static bool testPow2();
    static bool testPowReal();
    static bool testPowComplex();
    static bool testNorm();
    static bool testAbs();
    static bool testInverse();
    static bool testSqrt();
    static bool testExp();
    static bool testLog();
    static bool testLog10();
    static bool testSin();
    static bool testCos();
    static bool testTan();
    static bool testSinh();
    static bool testCosh();
    static bool testTanh();
    static boost::unit_test_framework::test_suite* suite();
};

#endif