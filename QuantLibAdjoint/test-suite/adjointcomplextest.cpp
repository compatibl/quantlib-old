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

#include "adjointcomplextest.hpp"
#include "adjointcomplexutils.hpp"
#include <ql/quantlib.hpp>

bool AdjointComplexTest::testConstructor()
{
    CL_TEST_MESSAGE("Testing complex construction...");

#if defined CL_TAPE_COMPLEX_ENABLED
    StdComplex sz = { 2, 7 };

    // Construct from doubles.
    std::vector<Complex> Z1 = odd_constructions(sz.real(), sz.imag());

    bool ok = true;
    for (Complex const& z : Z1)
    {
        ok = ok
            && check_consistency(z, sz)
            && check_mode(z, Complex::default_mode);
    }

    std::vector<Real> X = { sz.real(), sz.imag() };

    // Construct from parameter Re and Im.
    Z1 = odd_constructions(X[0], X[1]);

    for (Complex const& z : Z1)
    {
        ok = ok
            && check_consistency(z, sz)
            && check_mode(z, Complex::default_mode);
    }

    std::vector<Real> Re = { sz.real() };
    cl::Independent(Re);
    X[0] = Re[0];

    // Construct from variable Re and parameter Im
    std::vector<Complex> Z2 = odd_constructions(X[0], X[1]);

    for (Complex const& z : Z2)
    {
        ok = ok
            && check_consistency(z, sz)
            && check_mode(z, Complex::RealBase);
    }

    // Cleaning tape.
    cl::tape_function<double>(Re, X);
    std::vector<Real> Im = { sz.imag() };
    cl::Independent(Im);
    X[0] = sz.real();
    X[1] = Im[0];

    // Construct from parameter Re and variable Im
    Z2 = odd_constructions(X[0], X[1]);

    for (Complex const& z : Z2)
    {
        ok = ok
            && check_consistency(z, sz)
            && check_mode(z, Complex::RealBase);
    }

    // Cleaning tape.
    cl::tape_function<double>(Im, X);
    cl::Independent(X);

    // Construct from variable Re and variable Im
    Z2 = odd_constructions(X[0], X[1]);

    for (Complex const& z : Z2)
    {
        ok = ok
            && check_consistency(z, sz)
            && check_mode(z, Complex::RealBase);
    }

    // Cleaning tape.
    std::vector<Real> Y = { X[0], X[1] };
    cl::tape_function<double>(X, Y);

    // Constructions from std::complex<double>
    std::vector<Complex> Z3 = {
        Complex(sz)
        , Complex() = Complex(sz)
        , Complex(Complex(sz))
    };

    for (Complex const& z : Z3)
    {
        ok = ok
            && check_consistency(z, sz)
            && check_mode(z, Complex::default_mode);
    }

    cl::Independent(Z3);

    // Copy constructor from variable.
    std::vector<Complex> W = { Complex(Z3[0]) };

    ok &= check_mode(W[0], Complex::ComplBase);

    // Cleaning tape.
    cl::tape_function<Complex>(Z3, W);

    return ok;
#endif
    return false;
}

bool AdjointComplexTest::testAssign()
{
    CL_TEST_MESSAGE("Testing complex assignment...");

    auto std_f = [](StdComplex z1, StdComplex z2)
    {
        z1 = z2;
        return z1;
    };
    auto f = [](Complex z1, Complex z2)
    {
        z1 = z2;
        return z1;
    };
    auto df1 = [](Complex, Complex)
    {
        return Complex(0);
    };
    auto df2 = [](Complex, Complex)
    {
        return Complex(1);
    };

    auto std_g = [](StdComplex z, double x)
    {
        z = x;
        return z;
    };
    auto g = [](Complex z, Real x)
    {
        z = x;
        return z;
    };
    auto dg1 = [](Complex z, Real x)
    {
        return Complex(0);
    };

    return testBinFunction(f, std_f)
        && testBinFunctionDerivative<Real>(f, df1, df2)
        && testBinFunctionDerivative<Complex>(f, df1, df2)
        && testMixedFunction(g, std_g)
        && testMixedFunctionDerivative(g, dg1);
}

bool AdjointComplexTest::testImag()
{
    CL_TEST_MESSAGE("Testing complex imag() method...");

    auto std_f = [](StdComplex z)
    {
        return z.imag();
    };
    auto f = [](Complex z)
    {
        return z.imag();
    };
    auto df = [](Complex z)
    {
        return std::vector<Real>{ 0, 1 };
    };

    auto std_g = [](StdComplex z, double x)
    {
        z.imag(x);
        return z;
    };
    auto g = [](Complex z, Real x)
    {
        z.imag(x);
        return z;
    };

    return testFunction(f, std_f)
        && testFunctionDerivative<Real>(f, df)
        && testMixedFunctionRealBase(g, std_g);
}

bool AdjointComplexTest::testReal()
{
    CL_TEST_MESSAGE("Testing complex real() method...");

    auto std_f = [](StdComplex z)
    {
        return z.real();
    };
    auto f = [](Complex z)
    {
        return z.real();
    };
    auto df = [](Complex z)
    {
        return std::vector<Real>{ 1, 0 };
    };

    auto std_g = [](StdComplex z, double x)
    {
        z.imag(x);
        return z;
    };
    auto g = [](Complex z, Real x)
    {
        z.imag(x);
        return z;
    };

    return testFunction(f, std_f)
        && testFunctionDerivative<Real>(f, df)
        && testMixedFunctionRealBase(g, std_g);
}

bool AdjointComplexTest::testConj()
{
    CL_TEST_MESSAGE("Testing complex conjugation...");

    auto std_f = [](StdComplex z)
    {
        return conj(z);
    };
    auto f = [](Complex z)
    {
        return conj(z);
    };

    return testFunction<standart_sample, Real>(f, std_f);
}

bool AdjointComplexTest::testAdd()
{
    CL_TEST_MESSAGE("Testing complex addition...");

    auto std_f = [](StdComplex z1, StdComplex z2)
    {
        return z1 + z2;
    };
    auto f = [](Complex z1, Complex z2)
    {
        return z1 + z2;
    };
    auto df1 = [](Complex, Complex)
    {
        return Complex(1);
    };
    auto df2 = [](Complex, Complex)
    {
        return Complex(1);
    };

    return testBinFunction(f, std_f)
        && testBinFunctionDerivative<Real>(f, df1, df2)
        && testBinFunctionDerivative<Complex>(f, df1, df2);
}

bool AdjointComplexTest::testSub()
{
    CL_TEST_MESSAGE("Testing complex subtraction...");

    auto std_f = [](StdComplex z1, StdComplex z2)
    {
        return z1 - z2;
    };
    auto f = [](Complex z1, Complex z2)
    {
        return z1 - z2;
    };
    auto df1 = [](Complex, Complex)
    {
        return Complex(1);
    };
    auto df2 = [](Complex, Complex)
    {
        return Complex(-1);
    };

    return testBinFunction(f, std_f)
        && testBinFunctionDerivative<Real>(f, df1, df2)
        && testBinFunctionDerivative<Complex>(f, df1, df2);
}

bool AdjointComplexTest::testMul()
{
    CL_TEST_MESSAGE("Testing complex multiplication...");

    auto std_f = [](StdComplex z1, StdComplex z2)
    {
        return z1 * z2;
    };
    auto f = [](Complex z1, Complex z2)
    {
        return z1 * z2;
    };
    auto df1 = [](Complex, Complex z2)
    {
        return z2;
    };
    auto df2 = [](Complex z1, Complex)
    {
        return z1;
    };

    return testBinFunction(f, std_f)
        && testBinFunctionDerivative<Real>(f, df1, df2)
        && testBinFunctionDerivative<Complex>(f, df1, df2);
}

bool AdjointComplexTest::testDiv()
{
    CL_TEST_MESSAGE("Testing complex division...");

    auto std_f = [](StdComplex z1, StdComplex z2)
    {
        return z1 / z2;
    };
    auto f = [](Complex z1, Complex z2)
    {
        return z1 / z2;
    };
    auto df1 = [](Complex, Complex z2)
    {
        return 1. / z2;
    };
    auto df2 = [](Complex z1, Complex z2)
    {
        return -1. * z1 / (z2 * z2);
    };

    return testBinFunction(f, std_f, 1e-10)
        && testBinFunctionDerivative<Real>(f, df1, df2, 1e-10)
        && testBinFunctionDerivative<Complex>(f, df1, df2, 1e-10);
}

bool AdjointComplexTest::testAddReal()
{
    CL_TEST_MESSAGE("Testing addition of Complex and Real...");

    auto std_f = [](StdComplex z, double x)
    {
        return z + x;
    };
    auto f = [](Complex z, Real x)
    {
        return z + x;
    };
    auto df1 = [](Complex z, Real x)
    {
        return Complex(1.);
    };

    auto std_g = [](StdComplex z, double x)
    {
        return x + z;
    };
    auto g = [](Complex z, Real x)
    {
        return x + z;
    };
    auto dg1 = [](Complex z, Real x)
    {
        return Complex(1.);
    };

    return testMixedFunction(f, std_f)
        && testMixedFunctionDerivative(f, df1)
        && testMixedFunction(g, std_g)
        && testMixedFunctionDerivative(g, dg1);
}

bool AdjointComplexTest::testSubReal()
{
    CL_TEST_MESSAGE("Testing subtraction of Complex and Real...");

    auto std_f = [](StdComplex z, double x)
    {
        return z - x;
    };
    auto f = [](Complex z, Real x)
    {
        return z - x;
    };
    auto df1 = [](Complex z, Real x)
    {
        return Complex(1.);
    };

    auto std_g = [](StdComplex z, double x)
    {
        return x - z;
    };
    auto g = [](Complex z, Real x)
    {
        return x - z;
    };
    auto dg1 = [](Complex z, Real x)
    {
        return Complex(-1.);
    };

    return testMixedFunction(f, std_f)
        && testMixedFunctionDerivative(f, df1)
        && testMixedFunction(g, std_g)
        && testMixedFunctionDerivative(g, dg1);
}

bool AdjointComplexTest::testMulReal()
{
    CL_TEST_MESSAGE("Testing multiplication of Complex and Real...");

    auto std_f = [](StdComplex z, double x)
    {
        return z * x;
    };
    auto f = [](Complex z, Real x)
    {
        return z * x;
    };
    auto df1 = [](Complex z, Real x)
    {
        return Complex(x);
    };

    auto std_g = [](StdComplex z, double x)
    {
        return x * z;
    };
    auto g = [](Complex z, Real x)
    {
        return x * z;
    };
    auto dg1 = [](Complex z, Real x)
    {
        return Complex(x);
    };

    return testMixedFunction(f, std_f)
        && testMixedFunctionDerivative(f, df1)
        && testMixedFunction(g, std_g)
        && testMixedFunctionDerivative(g, dg1);
}

bool AdjointComplexTest::testDivReal()
{
    CL_TEST_MESSAGE("Testing division of Complex and Real...");

    auto std_f = [](StdComplex z, double x)
    {
        return z / x;
    };
    auto f = [](Complex z, Real x)
    {
        return z / x;
    };
    auto df1 = [](Complex z, Real x)
    {
        return Complex(1.) / x;
    };

    auto std_g = [](StdComplex z, double x)
    {
        return x / z;
    };
    auto g = [](Complex z, Real x)
    {
        return x / z;
    };
    auto dg1 = [](Complex z, Real x)
    {
        return -x / (z * z);
    };

    return testMixedFunction(f, std_f)
        && testMixedFunctionDerivative(f, df1)
        && testMixedFunction(g, std_g, 1e-10)
        && testMixedFunctionDerivative(g, dg1, 1e-10);
}

bool AdjointComplexTest::testArg()
{
    CL_TEST_MESSAGE("Testing complex argument function...");

    auto std_f = [](StdComplex z)
    {
        return arg(z);
    };
    auto f = [](Complex z)
    {
        return arg(z);
    };
    auto df = [](Complex z)
    {
        Real ratio = z.imag() / z.real();
        Real atanDeriv = 1 / (1 + ratio * ratio);
        return std::vector<Real>{
            atanDeriv * (-z.imag() / z.real() / z.real())
                , atanDeriv / z.real()
        };
    };

    return testFunction(f, std_f, 1e-10)
        && testFunctionDerivative<Real>(f, df, 1e-10);
}

bool AdjointComplexTest::testPolar()
{
    CL_TEST_MESSAGE("Testing complex polar function...");

    // Note, polar is (Real, Real) -> Complex function, but we use
    // (Complex, Real) -> Complex adaptation for testing.
    auto std_f = [](StdComplex z, double x)
    {
        return std::polar(z.real(), x);
    };
    auto f = [](Complex z, Real x)
    {
        return std::polar(z.real(), x);
    };

    return testMixedFunctionRealBase(f, std_f);
}

bool AdjointComplexTest::testPow2()
{
    CL_TEST_MESSAGE("Testing complex z^2 function...");

    bool ok;
    {
        std::vector<cl::tape_double> x = { 2, 7 };
        cl::Independent(x);
        std::complex<Real> z(x[0], x[1]);
        std::complex<Real> w = std::pow(z, 2);
        std::vector<cl::tape_double> y(2);
        y[0] = w.real();
        y[1] = w.imag();
        cl::tape_function<double> f(x, y);
        std::vector<std::vector<Real> > expected = {
            { 2 * x[0], -2 * x[1] }
            , { 2 * x[1], 2 * x[0] }
        };
        ok = calculate_and_check(f, expected);
    }

    auto std_f = [](StdComplex z)
    {
        return std::pow(z, 2);
    };
    auto f = [](Complex z)
    {
        return std::pow(z, 2);
    };
    auto df = [](Complex z)
    {
        return 2 * z;
    };

    return ok
        && testFunction(f, std_f)
        && testFunctionDerivative<Real>(f, df)
        && testFunctionDerivative<Complex>(f, df);
}

bool AdjointComplexTest::testPowReal()
{
    CL_TEST_MESSAGE("Testing complex pow(Complex, Real) and pow(Real, Complex) function...");

    auto std_f = [](StdComplex z, double x)
    {
        return std::pow(z, x);
    };
    auto f = [](Complex z, Real x)
    {
        return std::pow(z, x);
    };
    auto df1 = [](Complex z, Real x)
    {
        return x * std::pow(z, x - 1.);
    };

    auto std_g = [](StdComplex z, double x)
    {
        return std::pow(x, z);
    };
    auto g = [](Complex z, Real x)
    {
        return std::pow(x, z);
    };
    auto dg1 = [](Complex z, Real x)
    {
        return std::log(x) * std::pow(x, z);
    };

    return testMixedFunction(f, std_f, 1e-10)
        && testMixedFunctionDerivative<small_arguments>(f, df1, 1e-10)
        && testMixedFunction<small_arguments>(g, std_g, 1e-10)
        && testMixedFunctionDerivative<small_arguments>(g, dg1, 1e-10);
}

bool AdjointComplexTest::testPowComplex()
{
    CL_TEST_MESSAGE("Testing complex pow(Complex, Complex) function...");

    auto std_f = [](StdComplex z1, StdComplex z2)
    {
        return std::pow(z1, z2);
    };
    auto f = [](Complex z1, Complex z2)
    {
        return std::pow(z1, z2);
    };
    auto df1 = [](Complex z1, Complex z2)
    {
        return z2 * std::pow(z1, z2 - 1.);
    };
    auto df2 = [](Complex z1, Complex z2)
    {
        return std::log(z1) * std::pow(z1, z2);
    };

    return testBinFunction<small_arguments>(f, std_f, 1e-10)
        && testBinFunctionDerivative<Real, small_arguments>(f, df1, df2, 1e-10)
        && testBinFunctionDerivative<Complex, small_arguments>(f, df1, df2, 1e-10);
}

bool AdjointComplexTest::testNorm()
{
    CL_TEST_MESSAGE("Testing complex norm...");

    bool ok;
    {
        std::vector<cl::tape_double> x = { 2, 7 };
        cl::Independent(x);
        std::complex<Real> z(x[0], x[1]);
        std::vector<cl::tape_double> y(1, std::norm(z));
        cl::tape_function<double> f(x, y);
        std::vector<std::vector<Real> > expected = {
            { 2 * x[0], 2 * x[1] }
        };
        ok = calculate_and_check(f, expected);
    }

    auto std_f = [](StdComplex z)
    {
        return std::norm(z);
    };
    auto f = [](Complex z)
    {
        return std::norm(z);
    };
    auto df = [](Complex z)
    {
        return std::vector<Real>{ 2 * z.real(), 2 * z.imag() };
    };

    return ok
        && testFunction(f, std_f)
        && testFunctionDerivative<Real>(f, df);
}

bool AdjointComplexTest::testAbs()
{
    CL_TEST_MESSAGE("Testing complex absolute value...");

    bool ok;
    {
        std::vector<cl::tape_double> x = { 2, 7 };
        cl::Independent(x);
        std::complex<Real> z(x[0], x[1]);
        std::vector<cl::tape_double> y(1, std::abs(z));
        cl::tape_function<double> f(x, y);
        std::vector<std::vector<Real> > expected = {
            { x[0] / std::abs(z), x[1] / std::abs(z) }
        };
        ok = calculate_and_check(f, expected);
    }

    auto std_f = [](StdComplex z)
    {
        return std::abs(z);
    };
    auto f = [](Complex z)
    {
        return std::abs(z);
    };
    auto df = [](Complex z)
    {
        return std::vector<Real> {
            z.real() / std::abs(z)
                , z.imag() / std::abs(z)
        };
    };

    return ok
        && testFunction(f, std_f)
        && testFunctionDerivative<Real>(f, df, 1e-10);
}

bool AdjointComplexTest::testInverse()
{
    CL_TEST_MESSAGE("Testing complex 1/z function...");

    bool ok;
    {
        std::vector<cl::tape_double> x = { 2, 7 };
        cl::Independent(x);
        std::complex<Real> z(x[0], x[1]);
        std::complex<Real> w = 1 / z;
        std::vector<cl::tape_double> y(2);
        y[0] = w.real();
        y[1] = w.imag();
        cl::tape_function<double> f(x, y);
        std::vector<std::vector<Real> > expected = {
            { -(std::pow(x[0], 2) - std::pow(x[1], 2)) / std::pow((std::pow(x[0], 2) + std::pow(x[1], 2)), 2),
            -2 * x[0] * x[1] / std::pow((std::pow(x[0], 2) + std::pow(x[1], 2)), 2)
            }
            , { 2 * x[0] * x[1] / std::pow((std::pow(x[0], 2) + std::pow(x[1], 2)), 2),
                -(std::pow(x[0], 2) - std::pow(x[1], 2)) / std::pow((std::pow(x[0], 2) + std::pow(x[1], 2)), 2)
            }
        };
        ok = calculate_and_check(f, expected, 1e-10);
    }

    auto std_f = [](StdComplex z)
    {
        return 1. / z;
    };
    auto f = [](Complex z)
    {
        return 1 / z;
    };
    auto df = [](Complex z)
    {
        return -1 / (z * z);
    };

    return ok
        && testFunction(f, std_f, 1e-10)
        && testFunctionDerivative<Real>(f, df, 1e-10)
        && testFunctionDerivative<Complex>(f, df, 1e-10);
}

bool AdjointComplexTest::testSqrt()
{
    CL_TEST_MESSAGE("Testing complex square root...");

    auto std_f = [](StdComplex z)
    {
        return std::sqrt(z);
    };
    auto f = [](Complex z)
    {
        return std::sqrt(z);
    };
    auto df = [](Complex z)
    {
        return 0.5 / std::sqrt(z);
    };

    return testFunction(f, std_f, 1e-10)
        && testFunctionDerivative<Real>(f, df, 1e-10)
        && testFunctionDerivative<Complex>(f, df, 1e-10);
}

bool AdjointComplexTest::testExp()
{
    CL_TEST_MESSAGE("Testing complex exp function...");

    bool ok;
    {
        std::vector<cl::tape_double> x = { 2, 7 };
        cl::Independent(x);
        std::complex<Real> z(x[0], x[1]);
        std::complex<Real> w = exp(z);
        std::vector<cl::tape_double> y(2);
        y[0] = w.real();
        y[1] = w.imag();
        cl::tape_function<double> f(x, y);
        std::vector<std::vector<Real> > expected = {
            { std::exp(x[0])*std::cos(x[1]), -std::exp(x[0])*std::sin(x[1]) }
            , { std::exp(x[0])*std::sin(x[1]), std::exp(x[0])*std::cos(x[1]) }
        };
        ok = calculate_and_check(f, expected);
    }

    auto std_f = [](StdComplex z)
    {
        return std::exp(z);
    };
    auto f = [](Complex z)
    {
        return std::exp(z);
    };
    auto df = [](Complex z)
    {
        return std::exp(z);
    };

    return ok
        && testFunction(f, std_f, 1e-10)
        && testFunctionDerivative<Real>(f, df)
        && testFunctionDerivative<Complex>(f, df);
}

bool AdjointComplexTest::testLog()
{
    CL_TEST_MESSAGE("Testing complex log function...");

    bool ok;
    {
        std::vector<cl::tape_double> x = { 2, 7 };
        cl::Independent(x);
        std::complex<Real> z(x[0], x[1]);
        std::complex<Real> w = std::log(z);
        std::vector<cl::tape_double> y(2);
        y[0] = w.real();
        y[1] = w.imag();
        cl::tape_function<double> f(x, y);
        std::vector<std::vector<Real> > expected = {
            { x[0] / (std::pow(x[0], 2) + std::pow(x[1], 2)), x[1] / (std::pow(x[0], 2) + std::pow(x[1], 2)) }
            , { -x[1] / (std::pow(x[0], 2) + std::pow(x[1], 2)), x[0] / (std::pow(x[0], 2) + std::pow(x[1], 2)) }
        };
        ok = calculate_and_check(f, expected, 1e-10);
    }

    auto std_f = [](StdComplex z)
    {
        return std::log(z);
    };
    auto f = [](Complex z)
    {
        return std::log(z);
    };
    auto df = [](Complex z)
    {
        return 1 / z;
    };

    return ok
        && testFunction(f, std_f, 1e-10)
        && testFunctionDerivative<Real>(f, df, 1e-10)
        && testFunctionDerivative<Complex>(f, df, 1e-10);
}

bool AdjointComplexTest::testLog10()
{
    CL_TEST_MESSAGE("Testing complex log10 function...");

    auto std_f = [](StdComplex z)
    {
        return std::log10(z);
    };
    auto f = [](Complex z)
    {
        return std::log10(z);
    };
    auto df = [](Complex z)
    {
        return 0.43429448190325182765112891891660508L / z;
    };

    return testFunction(f, std_f, 1e-10)
        && testFunctionDerivative<Real>(f, df, 1e-10)
        && testFunctionDerivative<Complex>(f, df, 1e-10);
}

bool AdjointComplexTest::testSin()
{
    CL_TEST_MESSAGE("Testing complex sin function...");

    auto std_f = [](StdComplex z)
    {
        return std::sin(z);
    };
    auto f = [](Complex z)
    {
        return std::sin(z);
    };
    auto df = [](Complex z)
    {
        return std::cos(z);
    };

    return testFunction<small_arguments>(f, std_f, 1e-10)
        && testFunctionDerivative<Real, small_arguments>(f, df, 1e-10)
        && testFunctionDerivative<Complex, small_arguments>(f, df, 1e-10);
}

bool AdjointComplexTest::testCos()
{
    CL_TEST_MESSAGE("Testing complex cos function...");

    auto std_f = [](StdComplex z)
    {
        return std::cos(z);
    };
    auto f = [](Complex z)
    {
        return std::cos(z);
    };
    auto df = [](Complex z)
    {
        return -1. * std::sin(z);
    };

    return testFunction<small_arguments>(f, std_f, 1e-10)
        && testFunctionDerivative<Real, small_arguments>(f, df, 1e-10)
        && testFunctionDerivative<Complex, small_arguments>(f, df, 1e-10);
}

bool AdjointComplexTest::testTan()
{
    CL_TEST_MESSAGE("Testing complex tan function...");

    auto std_f = [](StdComplex z)
    {
        return std::tan(z);
    };
    auto f = [](Complex z)
    {
        return std::tan(z);
    };
    auto df = [](Complex z)
    {
        return 1. / pow(std::cos(z), 2);
    };

    return testFunction<small_arguments>(f, std_f, 1e-10)
        && testFunctionDerivative<Real, small_arguments>(f, df, 1e-10)
        && testFunctionDerivative<Complex, small_arguments>(f, df, 1e-10);
}

bool AdjointComplexTest::testSinh()
{
    CL_TEST_MESSAGE("Testing complex sinh function...");

    auto std_f = [](StdComplex z)
    {
        return std::sinh(z);
    };
    auto f = [](Complex z)
    {
        return std::sinh(z);
    };
    auto df = [](Complex z)
    {
        return std::cosh(z);
    };

    return testFunction<small_arguments>(f, std_f, 1e-10)
        && testFunctionDerivative<Real, small_arguments>(f, df, 1e-10)
        && testFunctionDerivative<Complex, small_arguments>(f, df, 1e-10);
}

bool AdjointComplexTest::testCosh()
{
    CL_TEST_MESSAGE("Testing complex cosh function...");

    auto std_f = [](StdComplex z)
    {
        return std::cosh(z);
    };
    auto f = [](Complex z)
    {
        return std::cosh(z);
    };
    auto df = [](Complex z)
    {
        return std::sinh(z);
    };

    return testFunction<small_arguments>(f, std_f, 1e-10)
        && testFunctionDerivative<Real, small_arguments>(f, df, 1e-10)
        && testFunctionDerivative<Complex, small_arguments>(f, df, 1e-10);
}

bool AdjointComplexTest::testTanh()
{
    CL_TEST_MESSAGE("Testing complex tanh function...");

    auto std_f = [](StdComplex z)
    {
        return std::tanh(z);
    };
    auto f = [](Complex z)
    {
        return std::tanh(z);
    };
    auto df = [](Complex z)
    {
        return 1. / pow(std::cosh(z), 2);
    };

    return testFunction<small_arguments>(f, std_f, 1e-10)
        && testFunctionDerivative<Real, small_arguments>(f, df, 1e-10)
        && testFunctionDerivative<Complex, small_arguments>(f, df, 1e-10);
}


test_suite* AdjointComplexTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Adjoint Complex Differentiation test");
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testConstructor));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testAssign));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testImag));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testReal));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testConj));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testAdd));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testSub));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testMul));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testDiv));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testAddReal));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testSubReal));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testMulReal));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testDivReal));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testArg));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testPolar));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testPow2));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testPowReal));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testPowComplex));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testNorm));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testAbs));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testInverse));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testSqrt));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testExp));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testLog));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testLog10));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testSin));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testCos));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testTan));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testSinh));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testCosh));
    suite->add(QUANTLIB_TEST_CASE(&AdjointComplexTest::testTanh));

    return suite;
}


#if defined CL_TAPE_COMPLEX_ENABLED

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_complex_differentation)

BOOST_AUTO_TEST_CASE(testComplexConstructor)
{
    BOOST_CHECK(AdjointComplexTest::testConstructor());
}
BOOST_AUTO_TEST_CASE(testComplexAssign)
{
    BOOST_CHECK(AdjointComplexTest::testAssign());
}
BOOST_AUTO_TEST_CASE(testComplexImag)
{
    BOOST_CHECK(AdjointComplexTest::testImag());
}
BOOST_AUTO_TEST_CASE(testComplexReal)
{
    BOOST_CHECK(AdjointComplexTest::testReal());
}
BOOST_AUTO_TEST_CASE(testComplexConj)
{
    BOOST_CHECK(AdjointComplexTest::testConj());
}
BOOST_AUTO_TEST_CASE(testComplexAdd)
{
    BOOST_CHECK(AdjointComplexTest::testAdd());
}
BOOST_AUTO_TEST_CASE(testComplexSub)
{
    BOOST_CHECK(AdjointComplexTest::testSub());
}
BOOST_AUTO_TEST_CASE(testComplexMul)
{
    BOOST_CHECK(AdjointComplexTest::testMul());
}
BOOST_AUTO_TEST_CASE(testComplexDiv)
{
    BOOST_CHECK(AdjointComplexTest::testDiv());
}
BOOST_AUTO_TEST_CASE(testComplexAddReal)
{
    BOOST_CHECK(AdjointComplexTest::testAddReal());
}
BOOST_AUTO_TEST_CASE(testComplexSubReal)
{
    BOOST_CHECK(AdjointComplexTest::testSubReal());
}
BOOST_AUTO_TEST_CASE(testComplexMulReal)
{
    BOOST_CHECK(AdjointComplexTest::testMulReal());
}
BOOST_AUTO_TEST_CASE(testComplexDivReal)
{
    BOOST_CHECK(AdjointComplexTest::testDivReal());
}
BOOST_AUTO_TEST_CASE(testComplexArg)
{
    BOOST_CHECK(AdjointComplexTest::testArg());
}
BOOST_AUTO_TEST_CASE(testComplexPolar)
{
    BOOST_CHECK(AdjointComplexTest::testPolar());
}
BOOST_AUTO_TEST_CASE(testComplexPow2)
{
    BOOST_CHECK(AdjointComplexTest::testPow2());
}
BOOST_AUTO_TEST_CASE(testComplexPowReal)
{
    BOOST_CHECK(AdjointComplexTest::testPowReal());
}
BOOST_AUTO_TEST_CASE(testComplexPowComplex)
{
    BOOST_CHECK(AdjointComplexTest::testPowComplex());
}
BOOST_AUTO_TEST_CASE(testComplexNorm)
{
    BOOST_CHECK(AdjointComplexTest::testNorm());
}
BOOST_AUTO_TEST_CASE(testComplexAbs)
{
    BOOST_CHECK(AdjointComplexTest::testAbs());
}
BOOST_AUTO_TEST_CASE(testComplexInverse)
{
    BOOST_CHECK(AdjointComplexTest::testInverse());
}
BOOST_AUTO_TEST_CASE(testComplexSqrt)
{
    BOOST_CHECK(AdjointComplexTest::testSqrt());
}
BOOST_AUTO_TEST_CASE(testComplexExp)
{
    BOOST_CHECK(AdjointComplexTest::testExp());
}
BOOST_AUTO_TEST_CASE(testComplexLog)
{
    BOOST_CHECK(AdjointComplexTest::testLog());
}
BOOST_AUTO_TEST_CASE(testComplexLog10)
{
    BOOST_CHECK(AdjointComplexTest::testLog10());
}
BOOST_AUTO_TEST_CASE(testComplexSin)
{
    BOOST_CHECK(AdjointComplexTest::testSin());
}
BOOST_AUTO_TEST_CASE(testComplexCos)
{
    BOOST_CHECK(AdjointComplexTest::testCos());
}
BOOST_AUTO_TEST_CASE(testComplexTan)
{
    BOOST_CHECK(AdjointComplexTest::testTan());
}
BOOST_AUTO_TEST_CASE(testComplexSinh)
{
    BOOST_CHECK(AdjointComplexTest::testSinh());
}
BOOST_AUTO_TEST_CASE(testComplexCosh)
{
    BOOST_CHECK(AdjointComplexTest::testCosh());
}
BOOST_AUTO_TEST_CASE(testComplexTanh)
{
    BOOST_CHECK(AdjointComplexTest::testTanh());
}

BOOST_AUTO_TEST_SUITE_END()

#endif

#endif
