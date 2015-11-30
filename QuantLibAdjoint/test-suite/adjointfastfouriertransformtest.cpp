/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2006 Joseph Wang
Copyright (C) 2009 Liquidnet Holdings, Inc.
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

#include "adjointfastfouriertransformtest.hpp"
#include "utilities.hpp"
#include <ql/math/fastfouriertransform.hpp>
#include <ql/math/array.hpp>
#include <complex>
#include <vector>
#include <functional>
#include "adjointtestutilities.hpp"
#include <ql/quantlib.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;
using namespace std;

typedef std::complex<cl::tape_double> Complex;
typedef std::complex<Real> ComplexReal;


std::vector<Complex> initVector(size_t size)
{
    std::vector<Complex> result;
    for (size_t i = 0; i < size; i++)
        result.push_back({i+1, 0});
    return result;
}

template <typename T>
void checkResults(std::vector<T> calc, std::vector<T> exp)
{
    if (calc.size() != exp.size())
        BOOST_ERROR("Size mismatch");
    for (size_t i = 0; i<calc.size(); i++)
    {
        if ((std::fabs(calc[i].real() - exp[i].real()) > 1.0e-2) ||
            (std::fabs(calc[i].imag() - exp[i].imag()) > 1.0e-2))
            BOOST_ERROR("Element(" << i << ")\n"
            << std::setprecision(4) << QL_SCIENTIFIC
            << "    calculated: " << calc[i] << "\n"
            << "    expected:   " << exp[i]);
    }
}

// A Fast Fourier Transform (FFT) is an algorithm to compute the discrete Fourier transform (DFT).
// The DFT is defined by the formula:
// Xk = sum{n=1}^(N-1)(xn*exp(-i*2*pi*k*n/N), k = 0..(N-1)
bool AdjointFastFourierTransformTest::testSimple()
{
    BOOST_TEST_MESSAGE("Testing complex direct FFT...");

    size_t m = 4;
    size_t N = pow(2, m);
    std::vector<Complex> Z = initVector(N);

    // Start tape recording.
    cl::Independent(Z);

    std::vector<Complex> W(N);
    FastFourierTransform fft(m);
    fft.transform(&Z[0], &Z[0]+N, &W[0]);

    // End of tape recording.
    cl::tape_function<Complex> f(Z, W);

    // Calculate derivatives using Forward mode.
    std::vector<std::complex<double> > dw, dX(N), result(N * N);
    for (size_t i = 0; i < N; i++)
    {
        dX[i] = 1;
        dw = f.Forward(1, dX);
        for (size_t j = 0; j < N; j++)
            result[i * N + j] = dw[j];
        dX[i] = 0;
    }

    // Calculate derivatives analytically.
    std::vector<std::complex<double>> analResult(N * N);
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
            analResult[i * N + j] = { cos(2*M_PI*i*j / N), -sin(2*M_PI*i*j / N) };
    }

    // Check analytical and adjoint results.
    checkResults<std::complex<double>>(result, analResult);

    return true;
}


// An Inverse Fast Fourier Transform (FFT) is an algorithm to compute the inverse discrete Fourier transform (DFT).
// The inverse DFT is defined by the formula:
// xk = sum{n=1}^(N-1)(Xn*exp(i*2*pi*k*n/N), k = 0..(N-1)
bool AdjointFastFourierTransformTest::testInverse()
{
    BOOST_TEST_MESSAGE("Testing complex inverse FFT...");

    size_t N = 8;
    std::vector<Complex> X = initVector(N);

    FastFourierTransform fft(N);
    size_t nFrq = fft.output_size();
    std::vector<Complex > Z(nFrq), ft(nFrq), res(N);

    fft.inverse_transform(X.begin(), X.end(), Z.begin());

    // Start tape recording.
    cl::Independent(Z);

    fft.inverse_transform(Z.begin(), Z.end(), ft.begin());

    // End of tape recording.
    cl::tape_function<Complex> f(Z, ft);

    // Calculate derivatives using Forward mode.
    std::vector<std::complex<double> > dw, dX(nFrq), result(nFrq*nFrq);
    for (size_t i = 0; i < nFrq; i++)
    {
        dX[i] = 1;
        dw = f.Forward(1, dX);
        for (size_t j = 0; j < nFrq; j++)
            result[i * nFrq + j] = dw[j];
        dX[i] = 0;
    }

    // Calculate derivatives analytically.
    std::vector<std::complex<double>> analResult(nFrq * nFrq);
    for (size_t i = 0; i < nFrq; i++)
    {
        for (size_t j = 0; j < nFrq; j++)
            analResult[i * nFrq + j] = { cos(2 * M_PI*i*j / nFrq), sin(2 * M_PI*i*j / nFrq) };
    }

    // Check analytical and adjoint results.
    checkResults<std::complex<double>>(result, analResult);

    return true;
}

test_suite* AdjointFastFourierTransformTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Adjoint Fast Fourier Transform test");
    suite->add(QUANTLIB_TEST_CASE(&AdjointFastFourierTransformTest::testSimple));
    suite->add(QUANTLIB_TEST_CASE(&AdjointFastFourierTransformTest::testInverse));
    return suite;
}

#if defined CL_TAPE_COMPLEX_ENABLED

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_fast_fourie_transform)

BOOST_AUTO_TEST_CASE(testSimple)
{
    BOOST_CHECK(AdjointFastFourierTransformTest::testSimple());
}

BOOST_AUTO_TEST_CASE(testInverse)
{
    BOOST_CHECK(AdjointFastFourierTransformTest::testInverse());
}

BOOST_AUTO_TEST_SUITE_END()

#endif

#endif