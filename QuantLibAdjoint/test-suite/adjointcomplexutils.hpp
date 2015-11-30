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

#ifndef ad_complex_differentiation_impl_hpp
#define ad_complex_differentiation_impl_hpp

#include "adjointtestutilities.hpp"
#include <ql/quantlib.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;

#define CL_ERROR(msg) BOOST_ERROR(msg)

#define CL_TEST_MESSAGE(msg) BOOST_TEST_MESSAGE(msg)

namespace
{
    typedef std::complex<cl::tape_double> Complex;
    typedef std::complex<double> StdComplex;

    // tol parameter is both relative and absolute tolerance.
    template <typename T>
    T applied_tol(T t, double tol)
    {
        return tol * std::max(1., std::abs(t));
    }

    // tol parameter is both relative and absolute tolerance.
    template <typename T1, typename T2>
    bool near_equal(T1 t1, T2 t2, double tol)
    {
        auto tolerance = applied_tol(t2, tol);
        return !(tolerance < std::abs(t1 - T1(t2)));
    }

    // tol parameter is both relative and absolute tolerance.
    template <typename Base>
    Base applied_tol(std::complex<Base> const& c, double tol)
    {
        return tol * std::max(1.
            , std::max(std::abs(c.imag())
            , std::abs(c.real())));
    }

    // tol parameter is both relative and absolute tolerance.
    template <typename Base1, typename Base2>
    bool near_equal(std::complex<Base1> const& c1, std::complex<Base2> const& c2, double tol)
    {
        auto tolerance = applied_tol(c2, tol);
        return !(tolerance < std::abs(c1.real() - c2.real()))
            && !(tolerance < std::abs(c1.imag() - c2.imag()));
    }

    bool check_mode(const Complex& z, Complex::Complex_Mode mode)
    {
        if (z.mode_ != mode)
        {
            CL_ERROR("Wrong complex mode:"
                << "\n\tCurrent:  " << z.mode_
                << "\n\tRequired: " << mode);
            return false;
        }
        return true;
    }

    template <typename T1, typename T2>
    bool check_consistency(T1 v1, T2 v2, double tol = 0)
    {
        if (!near_equal(v1, v2, tol))
        {
            CL_ERROR("Variable values mismatch:"
                << "\n\tFirst:  " << v1
                << "\n\tSecond: " << v2
                << "\n\tTolerance:  " << applied_tol(v2, tol));
            return false;
        }
        return true;
    }

    // tol parameter is both relative and absolute tolerance.
    template <typename T1, typename T2>
    bool check_derivative(T1 expected, T2 calculated, double tol)
    {
        if (!near_equal(expected, calculated, tol))
        {
            CL_ERROR("Derivative mismatch:"
                << "\n\tExpected:   " << expected
                << "\n\tCalculated: " << calculated
                << "\n\tTolerance:  " << applied_tol(calculated, tol));
            return false;
        }
        return true;
    }

    // tol parameter is both relative and absolute tolerance.
    template <typename T1, typename T2>
    bool check_value(T1 calculated, T2 expected, double tol)
    {
        if (!near_equal(calculated, expected, tol))
        {
            CL_ERROR("Function value mismatch:"
                << "\n\tExpected:   " << expected
                << "\n\tCalculated: " << calculated
                << "\n\tTolerance:  " << applied_tol(expected, tol));
            return false;
        }
        return true;
    }

    inline bool calculate_and_check(cl::tape_function<double>& f, const std::vector<std::vector<Real> >& expected, double tol = 0)
    {
        assert(f.Range() == expected.size());
        bool ok = true;
        std::vector<double> w(f.Range(), 0);
        for (int i = 0; i < expected.size(); i++)
        {
            w[i] = 1;
            std::vector<double> dw = f.Reverse(1, w);
            w[i] = 0;

            assert(dw.size() == expected[i].size());
            for (int j = 0; j < dw.size(); j++)
            {
                ok = ok && check_derivative(expected[i][j], dw[j], tol);
            }
        }
        return ok;
    }

    inline bool check_complex_derivative(cl::tape_function<double>& f, Complex dw_dz, double tol = 0)
    {
        std::vector<std::vector<Real> > expected = {
            { dw_dz.real(), -dw_dz.imag() }
            , { dw_dz.imag(), dw_dz.real() }
        };

        return calculate_and_check(f, expected, tol);
    }

    struct standart_sample
    {
        template <class ComplexType = StdComplex>
        static std::vector<ComplexType> test_sample()
        {
            std::vector<ComplexType> vz = {
                { 2, 7 }
                , { -5, 3 }
                , { 8, -0.7 }
                , { -0.04, -9 }
                , { 1e-15, 1e13 }
                , { -1e15, 3e-12 }
                , { 1e-8, -1e7 }
                , { -1e5, -3e-9 }
                , { 0.1, 50 }
                , { -30, 3.2 }
                , { 1.8, -40 }
                , { -60, -0.9 }
                , { 0.03, 0.09 }
                , { 0.06, -0.009 }
                , { -0.07, 0.12 }
                , { -0.02, -0.08 }
            };
            return vz;
        }
    };

    // Gives test sample with values < 100.
    struct small_arguments
    {
        template <class ComplexType = StdComplex>
        static std::vector<ComplexType> test_sample()
        {
            std::vector<ComplexType> vz = {
                { 2, 7 }
                , { -5, 3 }
                , { 8, -0.7 }
                , { -0.04, -9 }
                , { 0.1, 50 }
                , { -30, 3.2 }
                , { 1.8, -40 }
                , { -60, -0.9 }
                , { 0.03, 0.09 }
                , { 0.06, -0.009 }
                , { -0.07, 0.12 }
                , { -0.02, -0.08 }
            };
            return vz;
        }
    };

    struct positive_real
    {
        template <class ComplexType = StdComplex>
        static std::vector<ComplexType> test_sample()
        {
            std::vector<ComplexType> vz = {
                { 2, 7 }
                , { 5, 3 }
                , { 8, -0.7 }
                , { 0.04, -9 }
                , { 0.1, 50 }
                , { 30, 3.2 }
                , { 1.8, -40 }
                , { 60, -0.9 }
                , { 0.03, 0.09 }
                , { 0.06, -0.009 }
                , { 0.07, 0.12 }
                , { 0.02, -0.08 }
            };
            return vz;
        }
    };

    template <class Func, class StdFunc>
    struct TestImpl
    {
        typedef typename std::result_of<Func(Complex)>::type range_type;
        typedef typename std::result_of<StdFunc(StdComplex)>::type std_range_type;

        static_assert(std::is_same<range_type, Complex>::value
            && std::is_same<std_range_type, StdComplex>::value
            || std::is_same<range_type, Real>::value
            && std::is_same<std_range_type, double>::value
            , "Wrong functor return type");


        template <class Base, class RangeType = range_type>
        static bool testWhileTaping(StdComplex, Func, StdFunc, double);

        template <>
        static bool testWhileTaping<Real, Real>(StdComplex sz, Func f, StdFunc std_f, double tol)
        {
            std::vector<Real> X = { sz.real(), sz.imag() };
            cl::Independent(X);

            std::vector<Real> Y = { f(Complex(X[0], X[1])) };
            double sy = std_f(sz);

            bool ok = check_value(Y[0], sy, tol);

            cl::tape_function<double>(X, Y);

            return ok;
        }

        template <>
        static bool testWhileTaping<Real, Complex>(StdComplex sz, Func f, StdFunc std_f, double tol)
        {
            std::vector<Real> X = { sz.real(), sz.imag() };
            cl::Independent(X);

            std::vector<Complex> W = { f(Complex(X[0], X[1])) };
            StdComplex sw = std_f(sz);

            bool ok = check_value(W[0], sw, tol)
                && (W[0].mode_ == Complex::RealBase);

            std::vector<Real> Y = { W[0].real(), W[0].imag() };
            cl::tape_function<double>(X, Y);

            return ok;
        }

        template <>
        static bool testWhileTaping<Complex, Real>(StdComplex sz, Func f, StdFunc std_f, double tol) { return true; }

        template <>
        static bool testWhileTaping<Complex, Complex>(StdComplex sz, Func f, StdFunc std_f, double tol)
        {
            std::vector<Complex> Z = { Complex(sz) };
            cl::Independent(Z);

            bool ok = (Z[0].mode_ == Complex::ComplBase);

            std::vector<Complex> W = { f(Z[0]) };

            ok &= (W[0].mode_ == Complex::ComplBase);

            StdComplex sw = std_f(sz);
            cl::tape_function<Complex>(Z, W);

            return ok
                && check_value(W[0], sw, tol)
                && (W[0].mode_ == Complex::ComplBase);
        }

        static bool testWithoutTaping(StdComplex sz, Func f, StdFunc std_f, double tol = 0)
        {
            Complex z(sz);
            range_type y = f(z);
            std_range_type sy = std_f(sz);

            return check_value(y, sy, tol);
        }

        static bool test(StdComplex sz, Func f, StdFunc std_f, double tol = 0)
        {
            return testWhileTaping<Real>(sz, f, std_f, tol)
                && testWhileTaping<Complex>(sz, f, std_f, tol)
                && testWithoutTaping(sz, f, std_f, tol);
        }
    };

    // Test function consistency with standart realization.
    // f        - tested function,
    // std_f    - standart realization,
    // tol      - tolerance (both absolute and relative).
    // If f is Complex -> Complex function then f is tested with parameter,
    // real based variable and complex based variable argument.
    // If f is Complex -> Real function then f is tested with parameter and
    // real based variable argument.
    template <class Sample = standart_sample, class Func, class StdFunc>
    bool testFunction(Func f, StdFunc std_f, double tol = 0)
    {
#if defined CL_TAPE_COMPLEX_ENABLED
        auto vz = Sample::test_sample();

        bool ok = true;
        for (StdComplex z : vz)
        {
            ok = ok && TestImpl<Func, StdFunc>::test(z, f, std_f, tol);
        }
        return ok;
#endif
        return false;
    }

    // Test function consistency with standart realization.
    // f        - tested function,
    // std_f    - standart realization,
    // tol      - tolerance (both absolute and relative).
    // The function tested with Base based variable and parameter argument.
    template <class Sample, class Base, class Func, class StdFunc>
    bool testFunction(Func f, StdFunc std_f, double tol = 0)
    {
#if defined CL_TAPE_COMPLEX_ENABLED
        auto vz = Sample::test_sample();

        bool ok = true;
        for (StdComplex z : vz)
        {
            ok = ok
                && TestImpl<Func, StdFunc>::testWithoutTaping(z, f, std_f, tol)
                && TestImpl<Func, StdFunc>::testWhileTaping<Base>(z, f, std_f, tol);
        }
        return ok;
#endif
        return false;
    }

    template <class Base, class Func, class DFunc>
    struct DerivTestImpl;

    template <class Func, class DFunc>
    struct DerivTestImpl<Real, Func, DFunc>
    {
        typedef typename std::result_of<Func(Complex)>::type range_type;
        typedef typename std::result_of<DFunc(Complex)>::type d_range_type;

        template <class RangeType>
        static bool test(Complex, Func, DFunc, double);

        template <>
        static bool test<Real>(Complex zValue, Func f, DFunc df, double tol)
        {
            static_assert(std::is_same<d_range_type, std::vector<Real>>::value
                , "DFunc must return std::vector<Real> type with gradient of f(x+iy)");

            std::vector<Real> x = { zValue.real(), zValue.imag() };
            cl::Independent(x);
            Complex z(x[0], x[1]);
            std::vector<Real> y = { f(z) };
            cl::tape_function<double> tapeF(x, y);
            std::vector<std::vector<Real> > expected = { df(z) };
            return calculate_and_check(tapeF, expected, tol);
        }

        template <>
        static bool test<Complex>(Complex zValue, Func f, DFunc df, double tol)
        {
            static_assert(std::is_same<d_range_type, Complex>::value
                , "DFunc must return complex value of f'(z)");

            std::vector<Real> x = { zValue.real(), zValue.imag() };
            cl::Independent(x);
            Complex z(x[0], x[1]);
            Complex w = f(z);
            std::vector<Real> y = { w.real(), w.imag() };
            cl::tape_function<double> tapeF(x, y);
            return check_complex_derivative(tapeF, df(z), tol);
        }

        static bool test(Complex zValue, Func f, DFunc df, double tol = 0)
        {
            return test<range_type>(zValue, f, df, tol);
        }
    };

    template <class Func, class DFunc>
    struct DerivTestImpl<Complex, Func, DFunc>
    {
        typedef typename std::result_of<Func(Complex)>::type range_type;
        typedef typename std::result_of<DFunc(Complex)>::type d_range_type;

        static_assert(std::is_same<range_type, Complex>::value, "Func must be Complex -> Complex functor");
        static_assert(std::is_same<d_range_type, Complex>::value, "DFunc must be Complex -> Complex functor");

        static bool test(Complex zValue, Func f, DFunc df, double tol = 0)
        {
            std::vector<Complex> z = { Complex(zValue) };
            cl::Independent(z);
            std::vector<Complex> w = { f(z[0]) };
            cl::tape_function<Complex> tapeF(z, w);

            std::vector<std::complex<double> > dX = { 1 };
            StdComplex dw = tapeF.Forward(1, dX)[0];
            return check_derivative(df(z[0]), dw, tol);
        }
    };

    // Test function derivative consistency with expected results.
    // f        - tested function,
    // df       - expected derivative results,
    // tol      - tolerance (both absolute and relative).
    // f is tested for differentiation with respect to Base typed variable.
    // If f is Complex -> Complex function then df should return Complex value
    // f'(z) with f(z) derivative.
    // If f is Complex -> Real function then df should return std::vector<Real>
    // type with value { df/dx, df/dy } - gradient of f(x+i*y).
    template <class Base, class Sample = standart_sample, class Func, class DFunc>
    bool testFunctionDerivative(Func f, DFunc df, double tol = 0)
    {
#if defined CL_TAPE_COMPLEX_ENABLED

        auto vz = Sample::test_sample<Complex>();

        bool ok = true;
        for (Complex z : vz)
        {
            ok = ok && DerivTestImpl<Base, Func, DFunc>::test(z, f, df, tol);
        }
        return ok;
#endif
        return false;
    }

    // Test binary function consistency with standart realization.
    // f        - tested function,
    // std_f    - standart realization,
    // tol      - tolerance (both absolute and relative).
    template <class Sample = standart_sample, class Func, class StdFunc>
    bool testBinFunction(Func f, StdFunc std_f, double tol = 0)
    {
#if defined CL_TAPE_COMPLEX_ENABLED
        typedef typename std::result_of<StdFunc(StdComplex, StdComplex)>::type result_type;
        std::function<result_type(StdComplex, StdComplex)> sf(std_f);

        auto vz = Sample::test_sample();

        bool ok = true;
        for (StdComplex sz : vz)
        {
            using std::placeholders::_1;
            Complex z(sz);
            ok = ok
                && testFunction<Sample>(std::bind(f, _1, z), std::bind(sf, _1, sz), tol)
                && testFunction<Sample>(std::bind(f, z, _1), std::bind(sf, sz, _1), tol);
        }
        return ok;
#endif
        return false;
    }

    // Test binary function derivative consistency with expected results.
    // f        - tested function,
    // df1, df2 - expected derivatives results,
    // tol      - tolerance (both absolute and relative).
    // f is tested for differentiation with respect to Base typed variable.
    // If f is (Complex, Complex) -> Complex function f(z1, z2)
    // then df1(z1, z2) and df2(z1, z2) should return Complex value of df/dz1 and df/dz2 derivatives.
    template <class Base, class Sample = standart_sample, class Func, class D1Func, class D2Func>
    bool testBinFunctionDerivative(Func f, D1Func df1, D2Func df2, double tol = 0)
    {
#if defined CL_TAPE_COMPLEX_ENABLED
        auto vz = Sample::test_sample<Complex>();

        bool ok = true;
        for (Complex z : vz)
        {
            using std::placeholders::_1;
            ok = ok
                && testFunctionDerivative<Base, Sample>(std::bind(f, _1, z), std::bind(df1, _1, z), tol)
                && testFunctionDerivative<Base, Sample>(std::bind(f, z, _1), std::bind(df2, z, _1), tol);
        }
        return ok;
#endif
        return false;
    }

    // Test binary (Complex, Real) -> Complex function consistency with standart realization.
    // f        - tested function,
    // std_f    - standart realization,
    // tol      - tolerance (both absolute and relative).
    template <class Sample = standart_sample, class Func, class StdFunc>
    bool testMixedFunction(Func f, StdFunc std_f, double tol = 0)
    {
#if defined CL_TAPE_COMPLEX_ENABLED

        std::function<StdComplex(StdComplex, StdComplex)> sf = [&std_f](StdComplex z1, StdComplex z2)
        {
            return std_f(z1, z2.real());
        };

        std::function<Complex(Complex, Complex)> ff = [&f](Complex z1, Complex z2)
        {
            return f(z1, z2.real());
        };

        auto vz = Sample::test_sample();

        bool ok = true;
        for (StdComplex sz : vz)
        {
            using std::placeholders::_1;
            Complex z(sz);
            ok = ok
                && testFunction<Sample>(std::bind(ff, _1, z), std::bind(sf, _1, sz), tol)
                && testFunction<Sample, Real>(std::bind(ff, z, _1), std::bind(sf, sz, _1), tol);
        }

        return ok;
#endif
        return false;
    }

    // Test binary (Complex, Real) -> Complex function consistency with standart realization.
    // f        - tested function,
    // std_f    - standart realization,
    // tol      - tolerance (both absolute and relative).
    // The function tested with Real based variable only.
    template <class Sample = standart_sample, class Func, class StdFunc>
    bool testMixedFunctionRealBase(Func f, StdFunc std_f, double tol = 0)
    {
#if defined CL_TAPE_COMPLEX_ENABLED

        std::function<StdComplex(StdComplex, StdComplex)> sf = [&std_f](StdComplex z1, StdComplex z2)
        {
            return std_f(z1, z2.real());
        };

        std::function<Complex(Complex, Complex)> ff = [&f](Complex z1, Complex z2)
        {
            return f(z1, z2.real());
        };

        auto vz = Sample::test_sample();

        bool ok = true;
        for (StdComplex sz : vz)
        {
            using std::placeholders::_1;
            Complex z(sz);
            ok = ok
                && testFunction<Sample, Real>(std::bind(ff, _1, z), std::bind(sf, _1, sz), tol)
                && testFunction<Sample, Real>(std::bind(ff, z, _1), std::bind(sf, sz, _1), tol);
        }

        return ok;
#endif
        return false;
    }

    // Test binary (Complex, Real) -> Complex function derivative consistency with expected results.
    // f        - tested function,
    // df1      - expected derivative with respect to Complex argument,
    // tol      - tolerance (both absolute and relative).
    template <class Sample = standart_sample, class Func, class D1Func>
    bool testMixedFunctionDerivative(Func f, D1Func df1, double tol = 0)
    {
#if defined CL_TAPE_COMPLEX_ENABLED

        std::function<Complex(Complex, Complex)> d_f1 = [&df1](Complex z1, Complex z2)
        {
            return df1(z1, z2.real());
        };

        std::function<Complex(Complex, Complex)> ff = [&f](Complex z1, Complex z2)
        {
            return f(z1, z2.real());
        };

        auto vz = Sample::test_sample<Complex>();

        bool ok = true;
        for (Complex z : vz)
        {
            using std::placeholders::_1;
            ok = ok
                && testFunctionDerivative<Real, Sample>(std::bind(ff, _1, z), std::bind(d_f1, _1, z), tol)
                && testFunctionDerivative<Complex, Sample>(std::bind(ff, _1, z), std::bind(d_f1, _1, z), tol);
        }

        return ok;
#endif
        return false;
    }


    template <class R, class C = Complex>
    inline std::vector<C> odd_constructions(const R& real, const R& imag)
    {
        auto setImag = [](C z, R i)
        {
            z.imag(i);
            return z;
        };

        auto setReal = [](C z, R r)
        {
            z.real(r);
            return z;
        };

        const StdComplex i(0, 1);

        std::vector<R> X = { real, imag };

        std::vector<C> Z = {
            Complex(X[0], X[1])
            , setImag(Complex(X[0]), X[1])
            , Complex() = Complex(X[0], X[1])
            , Complex() = setImag(Complex(X[0]), X[1])
            , setImag(Complex() = Complex(X[0]), X[1])
            , Complex(X[1], X[0]) = Complex(X[0], X[1])
            , Complex(i) = Complex(X[0], X[1])
            , (Complex(42, 43)) = Complex(X[0], X[1])
            , (setImag(Complex(42), 43)) = Complex(X[0], X[1])
            , (Complex() = Complex(i * 43.)) = Complex(X[0], X[1])
            , setReal(Complex(42, X[1]), X[0])
            , setImag(setReal(Complex(), X[0]), X[1])
            , setReal(setImag(Complex(), X[1]), X[0])
        };

        return Z;
    }
}

#endif
