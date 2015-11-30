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

#include "adjointrealmathtest.hpp"
#include "utilities.hpp"
#include <ql/types.hpp>

using namespace std;
using namespace QuantLib;


#define CL_TEST_UNARY_MATH_FUN(name, x)                 \
{                                                       \
    Real fd = Real(name(x));                            \
    Real fr = name(Real(x));                            \
    if (fr != fd)                                       \
    {                                                   \
        cout << "\n" << "arg = " << x << "\n"           \
        << #name << " double = " << fd << "\n"          \
        << #name << " Real   = " << fr << "\n";         \
        result = false;                                 \
    }                                                   \
}                                                       \


#define CL_TEST_BINARY_MATH_FUN(name, x, y)                 \
{                                                           \
    Real fd = Real(name(x, y));                             \
    Real fr0 = name(Real(x), y);                            \
    Real fr1 = name(x, Real(y));                            \
    Real fr2 = name(Real(x), Real(y));                      \
    if ((fd != fr0) || (fd != fr1) || (fd != fr2))          \
    {                                                       \
        cout << "\n" << "arg0 = " << x << "\n"              \
        << #name << " double = " << fd << "\n"              \
        << #name << " Real 0 = " << fr0 << "\n"             \
        << #name << " Real 1 = " << fr1 << "\n"             \
        << #name << " Real 2 = " << fr2 << "\n";            \
        result = false;                                     \
    }                                                       \
}                                                           \


namespace
{
    bool check(Real fr, Real fd, Real x, string fname)
    {
        const Real
            relErr = 3.0e-13,
            absErr = 0.0;

        bool result = true;
        if (fd != 0)
        {
            Real err = abs((fr - fd) / fd);
            if (err > relErr)
            {
                BOOST_ERROR("\narg = " << x << "\n" << fname << " relative error = " << err << "\n");
                result = false;
            }
        }
        else
        {
            Real err = abs(fr);
            if (err > absErr)
            {
                BOOST_ERROR("\narg = " << x << "\n" << fname << " absolute error = " << err << "\n");
                result = false;
            }
        }
        return result;
    }
}


bool AdjointRealMathTest::testNotCppadMathFuns(double x)
{
    bool result = true;

    Real fd = Real(asinh(x));
    Real fr = asinh(Real(x));
    result &= check(fr, fd, x, "asinh");

    if (x >= 1) {
        fd = Real(acosh(x));
        fr = acosh(Real(x));
        result &= check(fr, fd, x, "acosh");
    }
    if (x >= -1 && x <= 1) {
        fd = Real(atanh(x));
        fr = atanh(Real(x));
        result &= check(fr, fd, x, "atanh");
    }
    return result;
}


bool AdjointRealMathTest::testUnaryMathFuns(double x)
{
    bool result = true;

    CL_TEST_UNARY_MATH_FUN(fabs, x)
    CL_TEST_UNARY_MATH_FUN(abs, x)
    CL_TEST_UNARY_MATH_FUN(floor, x)
    CL_TEST_UNARY_MATH_FUN(ceil, x)
    CL_TEST_UNARY_MATH_FUN(exp, x)
    CL_TEST_UNARY_MATH_FUN(sin, x)
    CL_TEST_UNARY_MATH_FUN(cos, x)
    CL_TEST_UNARY_MATH_FUN(tan, x)
    CL_TEST_UNARY_MATH_FUN(atan, x)
    CL_TEST_UNARY_MATH_FUN(sinh, x)
    CL_TEST_UNARY_MATH_FUN(cosh, x)
    CL_TEST_UNARY_MATH_FUN(tanh, x)
    CL_TEST_UNARY_MATH_FUN(isnan, x)
    if (x >= 0) {
        CL_TEST_UNARY_MATH_FUN(sqrt, x)
        CL_TEST_UNARY_MATH_FUN(log, x)
    }
    if (x >= -1 && x <= 1) {
        CL_TEST_UNARY_MATH_FUN(asin, x)
        CL_TEST_UNARY_MATH_FUN(acos, x)
    }

    result &= testNotCppadMathFuns(x);

    return result;
}


bool AdjointRealMathTest::testBinaryMathFuns(double x, double y)
{
    bool result = true;
    CL_TEST_BINARY_MATH_FUN(min, x, y)
    CL_TEST_BINARY_MATH_FUN(max, x, y)
    //CL_TEST_BINARY_MATH_FUN(atan2, x, y)
    //CL_TEST_BINARY_MATH_FUN(fmod, x, y)
    if (x > 0) {
        CL_TEST_BINARY_MATH_FUN(pow, x, y)
    }
    return result;
}


bool AdjointRealMathTest::testMathFuns()
{
    bool result = true;

    const double
        mind = std::numeric_limits<double>::min(),
        maxd = std::numeric_limits<double>::max();

    for (double x = mind; x <= maxd; x *= 1.1)
    {
        result &= testUnaryMathFuns(x);
        result &= testUnaryMathFuns(-x);
    }
    result &= testUnaryMathFuns(0);

    for (double x = mind; x <= maxd; x *= 2718.28)
    {
        for (double y = mind; y <= maxd; y *= 3145.92)
        {
            result &= testBinaryMathFuns(x, y);
            result &= testBinaryMathFuns(-x, y);
            result &= testBinaryMathFuns(x, -y);
            result &= testBinaryMathFuns(-x, -y);
        }
        result &= testBinaryMathFuns(x, 0);
        result &= testBinaryMathFuns(-x, 0);
    }
    for (double y = mind; y <= maxd; y *= 2718.281)
    {
        result &= testBinaryMathFuns(0, y);
        result &= testBinaryMathFuns(0, -y);
    }
    result &= testBinaryMathFuns(0, 0);

    return result;
}


boost::unit_test_framework::test_suite* AdjointRealMathTest::suite()
{
    boost::unit_test_framework::test_suite* suite =
        BOOST_TEST_SUITE("AdjointRealMathTest");
    suite->add(QUANTLIB_TEST_CASE(&testMathFuns));
    return suite;
}


#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_real_math)
BOOST_AUTO_TEST_CASE(testRealMathFunctions)
{
    BOOST_CHECK(AdjointRealMathTest::testMathFuns());
}
BOOST_AUTO_TEST_SUITE_END()

#endif