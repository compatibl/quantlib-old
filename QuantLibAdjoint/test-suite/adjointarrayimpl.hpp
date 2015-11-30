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


#ifndef cl_adjoint_array_optimization_impl_hpp
#define cl_adjoint_array_optimization_impl_hpp
#pragma once


#include <ql/quantlib.hpp>
#include <cl/tape/impl/ql/blackcalculator_t.hpp>

#include "adjointarraytest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"

using namespace QuantLib;
using namespace boost::unit_test_framework;

//#define BIG_FUNC
//#define IF_PROBLEM
//#define BLACK_SCHOLES_FUNC
//#define SIN_COS_EXP
#define SIMPLE_ARITHMETIC
//#define BLACK_GAMMA

//#define PRINT_CHECKPOINT

#define OUTPUT_FOLDER_NAME "AdjointArrayOptimization"


#if defined BLACK_SCHOLES_FUNC
#   define TARGET_FUNCTION(x) (black_func(x))
#   define TARGET_FUNCTION_NAME "BlackScholes price"
#elif defined BERMUDAN_SWAPTION_NVP
#   define TARGET_FUNCTION(x) (bermudan_nvp(x))
#   define TARGET_FUNCTION_NAME "Bermudan Swaption NVP"
#   define ENGINE tree
//#   define ENGINE fdHullWhite // If problem is not solved for fdHullWhite engine
#elif defined BIG_FUNC
#   define TARGET_FUNCTION(x) (big_func(x))
#   define TARGET_FUNCTION_NAME "Big function"
#elif defined IF_PROBLEM
#   define TARGET_FUNCTION(x) (if_problem(x))
#   define TARGET_FUNCTION_NAME "Function with if problem"
#elif defined SIN_COS_EXP
#   define TARGET_FUNCTION(x) (sin_cos_exp(x))
#   define TARGET_FUNCTION_NAME "sin & cos & exp"
#elif defined SIMPLE_ARITHMETIC
#   define TARGET_FUNCTION(x) (simple_arithmetic(x))
#   define TARGET_FUNCTION_NAME "Simple arithmetic"
#elif defined BLACK_GAMMA
#   define TARGET_FUNCTION(x) (black_gamma(x))
#   define TARGET_FUNCTION_NAME "Black-Scholes Gamma"
#else
#   define TARGET_FUNCTION(x) (some_func(x))
#   define TARGET_FUNCTION_NAME "Elementary function"
#endif

namespace std
{
    template <class T>
    inline CppAD::AD<T> max(const CppAD::AD<T>& x, const CppAD::AD<T>& y)
    {
        return CppAD::CondExpLt(x, y, y, x);
    }

    inline cl::tape_value max(const cl::tape_value& x, const cl::tape_value& y)
    {
        return CppAD::CondExpLt(x, y, y, x);
    }
}

namespace
{
    enum
    {
#if defined CL_GRAPH_GEN
#   if defined BERMUDAN_SWAPTION_NVP
        // Number of points for dependency graphics.
        pointNo = 20,
        // Number of points for performance graphic.
        iterNo = 20,
        // Step for portfolio size for performance testing .
        step = 1,
#   else
        // Number of points for dependency graphics.
        pointNo = 50,
        // Number of points for performance graphic.
        iterNo = 50,
        // Step for portfolio size for performance testing .
        step = 10,
#   endif
#else
        // Number of points for dependency graphics.
        pointNo = 0,
        // Number of points for performance graphic.
        iterNo = 0,
        // Step for portfolio size for performance testing .
        step = 1,
#endif
        // Defines performance accuracy. Its value is a minimum number
        // of calling of O(1) complexity methods per one performance test.
#if defined BLACK_SCHOLES_FUNC
        iterNumFactor = 40000,
#elif defined BERMUDAN_SWAPTION_NVP
        iterNumFactor = 0,
#elif defined BIG_FUNC
#   define iterNumFactor (50000 / big_func_count)
        //iterNumFactor = 3000,
#elif defined IF_PROBLEM || defined SIN_COS_EXP || defined SIMPLE_ARITHMETIC || defined BLACK_GAMMA
        iterNumFactor = 100000,
#else
        iterNumFactor = 500000,
#endif
    };


#if defined IF_PROBLEM
    const double finite_diff_step = 1e-9;
#else
    const double finite_diff_step = 1e-5;
#endif

    using std::max;
    using std::sin;
    using std::cos;
    using std::exp;
    using std::log;
    using std::pow;
    using std::sqrt;

    template <class T>
    inline T some_func(const T& x)
    {
        return x + x*x/50. + 10./(x + 1) + sin(x) * exp(x / 100.);
    }

    template <class T>
    inline T simple_arithmetic(const T& x)
    {
        return x
            + x / (1 + x * x)
            + 1 / (1 + x)
            + (x + 3) * (x + 5) * (x + 7) * (x + 11) / 1e6
            - x * (x * x + 8) / (x * x + 123)
                / (x * x + 2 / (x + 12))
                / (x + 4) / (x * x * x + 2 * x * x + 17)
                / (x * x + 23) / (x + 13 / (x * x + 9));
    }

    template <class T>
    inline T sin_cos_exp(const T& x)
    {
        return sin(cos(exp(cos(exp(sin(exp(cos(sin(exp(cos(sin(exp(cos(sin(exp(cos(x / 10.0)))))))))))))))));
    }

    template <class T>
    inline T black_gamma(const T& x)
    {
        T strike = 80 + 0.4 * x;
        T time = 1 + 0.01 * x;
        T spot = 110 - 0.2 * x;
        T rate = 0.02 + 0.0008 * x;
        T sigma = 0.2 + 0.003 * x;

        T stdDev = sigma / sqrt(time);
        T d1 = 1 / stdDev * (log(spot / strike) + (rate + pow(sigma, 2) / 2) * time);
        T fNorm = 1 / sqrt(2 * 3.14159265358979323846) * exp(-pow(d1, 2) / 2);
        T gamma = fNorm / spot / stdDev;
        return gamma;
    }

    static size_t big_func_count = 1000;

    template <class T>
    inline T big_func(T x)
    {
        T t = x + 10;
        for (int i = 0; i < big_func_count; i++)
        {
            t += some_func(x + 314 / (i + x + 1) + 420 / t + x*x / 157.8) / 1000.;
        }
        return t;
    }

    template <class T>
    inline T black_func(const T& x)
    {
        T strike = 100 + 0.2 * x;
        T discount = 0.5;
        T forward = 200;
        T stdDev = 0.5;

        return BlackCalculator_T<T>(Option::Call, strike, forward, stdDev, discount).value();
    }

    template <class T>
    inline T if_problem(const T& x)
    {
        return max(T(2.), T(3.))
            + max(x, -x + 0.0001)
            + max(T(0.1234), x * x / 100 - 60.1234)
            + max(T(50.1234), x)
            + max(T(90.1234), x)
            + max(T(-40.1234), -x)
            + max(T(-80.1234), -x)
            + max(T(-30.1234 - x / 2), -x)
            + max(T(-70.1234 - 0.1 * x), -x)
            + max(x, 100.1234 - x)
            + max(0.5 * x + 500.1234 / (x / 2 + 5), 80.1234 - x)
            + max(0.5 * x + x * x / 100.1234, 100.1234 - x)
            + max(0.5 * x + x * x / 100.1234, T(100.1234))
            + max(x, T(100.1234));
    }


    // Struct for graphics recording.
    struct OutStruct
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "X", "Adjoint/Tape results", "Fin. diff./Function results"
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type& operator<<(stream_type& stm, OutStruct& v)
        {
            stm << v.x_
                << ";" << v.adj_
                << ";" << v.thr_
                << std::endl;
            return stm;
        }

        Real x_;
        Real adj_;
        Real thr_;
    };


    struct OnePointTest
        : cl::AdjointTest<OnePointTest>
    {
        // Set it to non-zero to log the changes of compare operations result.
        static const size_t compare_change_count_ = 1;

        static const CalcMethod default_method = forward;

        explicit OnePointTest(Size size, cl::tape_empty_test_output* logger = nullptr)
            : AdjointTest()
            , size_(size)
            , input_(size)
            , output_(size)
            , X_()
            , Y_()
            , analyticalY_(size)
        {
            setLogger(logger);
            for (size_t i = 0; i < size_; i++)
            {
                input_[i] = 100 * (double(i) / size_);
            }
        }

        void recordTape()
        {
            size_t i = size_ / 2;
            X_ = { input_[i] };
            log_ << "Tape is recording for value   \t  [" << i << "] = " << input_[i]
                << "\nWhile the range is 0.." << size_ - 1 << ", with\t  [0] = " << input_[0]
                << "\n                              \t, [" << size_ - 1 << "] = " << input_[size_ - 1] << std::endl;
            cl::Independent(X_);
            calcY();
            f_ = std::make_unique<cl::tape_function<double>>(X_, Y_);
            f_->compare_change_count(compare_change_count_);
        }

        void calcY()
        {
            Y_ = { targetFunc(X_[0]) };
        }

        inline static Real targetFunc(const Real& x)
        {
            return TARGET_FUNCTION(x);
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            Real h = finite_diff_step;
            analyticalResults_.resize(size_);
            analyticalY_.resize(size_);
            for (Size i = 0; i < size_; i++)
            {
                Real left = targetFunc(input_[i] - h);
                Real right = targetFunc(input_[i] + h);
                analyticalResults_[i] = (right - left) / (2 * h);
                analyticalY_[i] = (right + left) / 2;
            }
        }

        void checkCompareChenges(size_t i)
        {
            if (f_->compare_change_number())
            {
                log_ << "Warning: Compare operation result chenged while calculating [" << i << "] derivative."
                    << "\n\t" << compare_change_count_ << "-th change index: " << f_->compare_change_op_index()
                    << "\n\tTotal changes number: " << f_->compare_change_number() << std::endl;
            }
        }

        void calcForward()
        {
            forwardResults_.resize(size_);
            std::vector<double> xq(1);
            std::vector<double> xq1(1, 1);
            for (Size i = 0; i < size_; i++)
            {
                xq[0] = input_[i];
                output_[i] = f_->Forward(0, xq)[0];
                checkCompareChenges(i);
                forwardResults_[i] = f_->Forward(1, xq1)[0];
            }
        }

        void calcReverse()
        {
            reverseResults_.resize(size_);
            std::vector<double> xq(1);
            std::vector<double> xy(1, 1);
            for (Size i = 0; i < size_; i++)
            {
                xq[0] = input_[i];
                output_[i] = f_->Forward(0, xq)[0];
                checkCompareChenges(i);
                reverseResults_[i] = f_->Reverse(1, xy)[0];
            }
        }

        Size indepVarNumber()
        {
            return size_;
        }

        size_t minPerfIteration() { return 70 * iterNumFactor; }

        size_t tapePerfIteration()
        {
            return iterNumFactor / 5 + 1;
        }

        double relativeTol() const { return 0.05; }

        double absTol() const { return 1e-5; }

        Size size_;
        std::vector<double> input_;
        std::vector<double> output_;
        std::vector<Real> X_;
        std::vector<Real> Y_;
        std::vector<Real> analyticalY_;
    };


    struct OnePointTestData
    {
        explicit OnePointTestData()
        : outPerform_(OUTPUT_FOLDER_NAME "//OnePoint " TARGET_FUNCTION_NAME, {
            { "title", TARGET_FUNCTION_NAME " differentiation performance with OnePoint optimization." }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "line_box_width", "-5" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        })
            , outAdjoint_(OUTPUT_FOLDER_NAME "//OnePoint " TARGET_FUNCTION_NAME, {
                { "title", TARGET_FUNCTION_NAME " Adjoint differentiation performance with OnePoint optimization." }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        })
            , outSize_(OUTPUT_FOLDER_NAME "//OnePoint " TARGET_FUNCTION_NAME, {
                { "title", "Tape size dependence on number of variables" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
            , { "cleanlog", "false" }
        })
            , outDeriv_(OUTPUT_FOLDER_NAME "//OnePoint " TARGET_FUNCTION_NAME "//output", {
                { "filename", "Derivative" }
            , { "not_clear", "Not" }
            , { "title", "Function derivative" }
            , { "ylabel", "Y" }
            , { "xlabel", "X" }
            , { "cleanlog", "false" }
        })
            , outValue_(OUTPUT_FOLDER_NAME "//OnePoint " TARGET_FUNCTION_NAME "//output", {
                { "filename", "FunctionValue" }
            , { "not_clear", "Not" }
            , { "title", "Function value" }
            , { "ylabel", "Y" }
            , { "xlabel", "X" }
            , { "cleanlog", "false" }
        })
        { }

        bool makeOutput()
        {
            bool ok = true;
            if (pointNo > 0)
            {
                recordDependencePlot(pointNo);
            }
            ok &= cl::recordPerformance(*this, iterNo, step);
            return ok;
        }

        std::shared_ptr<OnePointTest> getTest(size_t size)
        {
            return std::make_shared<OnePointTest>(size, &outPerform_);
        }

        void recordDependencePlot(size_t pointNo)
        {
            std::vector<OutStruct> outDeriv(pointNo);
            std::vector<OutStruct> outValue(pointNo);
            auto test = getTest(pointNo);
            test->recordTape();
            test->calcForward();
            test->calcAnalytical();
            for (Size i = 0; i < pointNo; i++)
            {
                outDeriv[i] = { test->input_[i], test->forwardResults_[i], test->analyticalResults_[i] };
                outValue[i] = { test->input_[i], test->output_[i], test->analyticalY_[i] };
            }
            outDeriv_ << outDeriv;
            outValue_ << outValue;
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output outDeriv_;
        cl::tape_empty_test_output outValue_;
    };


    struct NoOptTest
        : cl::AdjointTest<NoOptTest>
    {
        static const CalcMethod default_method = reverse;

        explicit NoOptTest(Size size, cl::tape_empty_test_output* logger = nullptr)
            : AdjointTest()
            , size_(size)
            , input_(size)
            , output_(size)
            , X_(size)
            , Y_()
            , analyticalY_(size)
        {
            setLogger(logger);
            for (size_t i = 0; i < size_; i++)
            {
                input_[i] = 100 * (double(i) / size_);
            }

            std::copy(input_.begin(), input_.end(), X_.begin());
        }

        void recordTape()
        {
            cl::Independent(X_);
            calcY();
            f_ = std::make_unique<cl::tape_function<double>>(X_, Y_);
        }

        void calcY()
        {
            std::transform(X_.begin(), X_.end(), output_.begin(), targetFunc);

            Real sum = std::accumulate(output_.begin(), output_.end(), Real());

            Y_ = { sum };
        }

        inline static Real targetFunc(const Real& x)
        {
            return TARGET_FUNCTION(x);
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            Real h = finite_diff_step;
            analyticalResults_.resize(size_);
            analyticalY_.resize(size_);
            for (Size i = 0; i < size_; i++)
            {
                Real left = targetFunc(input_[i] - h);
                Real right = targetFunc(input_[i] + h);
                analyticalResults_[i] = (right - left) / (2 * h);
                analyticalY_[i] = (right + left) / 2;
            }
        }

        Size indepVarNumber()
        {
            return size_;
        }

        size_t minPerfIteration() { return 2 * iterNumFactor; }

        double relativeTol() const { return 0.05; }

        double absTol() const { return 1e-5; }

        Size size_;
        std::vector<double> input_;
        std::vector<Real> output_;
        std::vector<Real> X_;
        std::vector<Real> Y_;
        std::vector<Real> analyticalY_;
    };


    struct NoOptTestData
    {
        explicit NoOptTestData()
        : outPerform_(OUTPUT_FOLDER_NAME "//NoOpt " TARGET_FUNCTION_NAME, {
            { "title", TARGET_FUNCTION_NAME " differentiation performance without optimization." }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "line_box_width", "-5" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        })
            , outAdjoint_(OUTPUT_FOLDER_NAME "//NoOpt " TARGET_FUNCTION_NAME, {
                { "title", TARGET_FUNCTION_NAME " Adjoint differentiation performance without optimization." }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        })
            , outSize_(OUTPUT_FOLDER_NAME "//NoOpt " TARGET_FUNCTION_NAME, {
                { "title", "Tape size dependence on number of variables" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
            , { "cleanlog", "false" }
        })
            , outDeriv_(OUTPUT_FOLDER_NAME "//NoOpt " TARGET_FUNCTION_NAME "//output", {
                { "filename", "Derivative" }
            , { "not_clear", "Not" }
            , { "title", "Function derivative" }
            , { "ylabel", "Y" }
            , { "xlabel", "X" }
            , { "cleanlog", "false" }
        })
            , outValue_(OUTPUT_FOLDER_NAME "//NoOpt " TARGET_FUNCTION_NAME "//output", {
                { "filename", "FunctionValue" }
            , { "not_clear", "Not" }
            , { "title", "Function value" }
            , { "ylabel", "Y" }
            , { "xlabel", "X" }
            , { "cleanlog", "false" }
        })
        { }

        bool makeOutput()
        {
            bool ok = true;
            if (pointNo > 0)
            {
                recordDependencePlot(pointNo);
            }
            ok &= cl::recordPerformance(*this, iterNo, step);
            return ok;
        }

        std::shared_ptr<NoOptTest> getTest(size_t size)
        {
            return std::make_shared<NoOptTest>(size, &outPerform_);
        }

        void recordDependencePlot(size_t pointNo)
        {
            std::vector<OutStruct> outDeriv(pointNo);
            std::vector<OutStruct> outValue(pointNo);
            auto test = getTest(pointNo);
            test->calcWithReverse();
            test->calcAnalytical();
            for (Size i = 0; i < pointNo; i++)
            {
                outDeriv[i] = { test->input_[i], test->reverseResults_[i], test->analyticalResults_[i] };
                outValue[i] = { test->input_[i], test->output_[i], test->analyticalY_[i] };
            }
            outDeriv_ << outDeriv;
            outValue_ << outValue;
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output outDeriv_;
        cl::tape_empty_test_output outValue_;
    };


    struct InnerArrayTest
        : cl::AdjointTest<InnerArrayTest>
    {
        typedef cl::tape_value inner_type;
        typedef cl::tape_wrapper<inner_type> tape_type;
        typedef cl::tape_function<inner_type> func_type;
        static const CalcMethod default_method = reverse;

        explicit InnerArrayTest(Size size, cl::tape_empty_test_output* logger = nullptr)
            : AdjointTest()
            , size_(size)
            , input_(size)
            , output_(size)
            , X_()
            , Y_()
            , analyticalY_(size)
            , adf_()
        {
            setLogger(logger);
            for (size_t i = 0; i < size_; i++)
            {
                input_[i] = 100 * (double(i) / size_);
            }

            tape_type x(inner_type(std::valarray<double>(input_.data(), size)));
            X_ = { std::move(x) };
        }

        void recordTape()
        {
            cl::Independent(X_);
            calcY();
            adf_ = std::make_unique<func_type>(X_, Y_);
            auto& res = CppAD::Value(Y_[0].value()).array_value_;
            output_.assign(std::begin(res), std::end(res));
        }

        void calcY()
        {
            Y_ = { targetFunc(X_[0]) };
        }

        template <class T>
        inline static T targetFunc(const T& x)
        {
            return TARGET_FUNCTION(x);
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            double h = finite_diff_step;
            analyticalResults_.resize(size_);
            analyticalY_.resize(size_);
            for (Size i = 0; i < size_; i++)
            {
                double left = targetFunc(input_[i] - h);
                double right = targetFunc(input_[i] + h);
                analyticalResults_[i] = (right - left) / (2 * h);
                analyticalY_[i] = (right + left) / 2;
            }
        }

        void calcForward()
        {
            typedef std::vector<inner_type> Vec;
            Vec dy = adf_->Forward(1, Vec{ 1 });
            setForwardResults(dy[0]);
        }

        void calcReverse()
        {
            typedef std::vector<inner_type> Vec;
            Vec dy = adf_->Reverse(1, Vec{ 1 });
            setReverseResults(dy[0]);
        }

        void setForwardResults(const inner_type& res)
        {
            auto& result = res.array_value_;
            forwardResults_.assign(std::begin(result), std::end(result));            
        }

        void setReverseResults(const inner_type& res)
        {
            auto& result = res.array_value_;
            reverseResults_.assign(std::begin(result), std::end(result));
        }

        size_t memory()
        {
            if (!adf_)
            {
                return 0;
            }
            return adf_->Memory();
        }

        Size indepVarNumber()
        {
            return size_;
        }

        size_t minPerfIteration() { return 5 * iterNumFactor; }

        size_t tapePerfIteration()
        {
            return iterNumFactor / size_ + 1;
        }

        double relativeTol() const { return 0.05; }

        double absTol() const { return 1e-5; }

        Size size_;
        std::vector<double> input_;
        std::vector<double> output_;
        std::vector<tape_type> X_;
        std::vector<tape_type> Y_;
        std::vector<double> analyticalY_;
        std::unique_ptr<func_type> adf_;
    };


    struct InnerArrayTestData
    {
        explicit InnerArrayTestData()
        : outPerform_(OUTPUT_FOLDER_NAME "//InnerArray " TARGET_FUNCTION_NAME, {
            { "title", TARGET_FUNCTION_NAME " differentiation performance using InnerArray." }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "line_box_width", "-5" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        })
            , outAdjoint_(OUTPUT_FOLDER_NAME "//InnerArray " TARGET_FUNCTION_NAME, {
                { "title", TARGET_FUNCTION_NAME " Adjoint differentiation performance using InnerArray." }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        })
            , outSize_(OUTPUT_FOLDER_NAME "//InnerArray " TARGET_FUNCTION_NAME, {
                { "title", "Tape size dependence on number of variables"
                "\\n (memory dynamically allocated for arrays is not included)" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
            , { "cleanlog", "false" }
        })
            , outDeriv_(OUTPUT_FOLDER_NAME "//InnerArray " TARGET_FUNCTION_NAME "//output", {
                { "filename", "Derivative" }
            , { "not_clear", "Not" }
            , { "title", "Function derivative" }
            , { "ylabel", "Y" }
            , { "xlabel", "X" }
            , { "cleanlog", "false" }
        })
            , outValue_(OUTPUT_FOLDER_NAME "//InnerArray " TARGET_FUNCTION_NAME "//output", {
                { "filename", "FunctionValue" }
            , { "not_clear", "Not" }
            , { "title", "Function value" }
            , { "ylabel", "Y" }
            , { "xlabel", "X" }
            , { "cleanlog", "false" }
        })
        { }

        bool makeOutput()
        {
            outPerform_.log() << "Start InnerArrayTestData output.\n";
            bool ok = true;
            if (pointNo > 0)
            {
                recordDependencePlot(pointNo);
            }
            ok &= cl::recordPerformance(*this, iterNo, step);
            return ok;
        }

        std::shared_ptr<InnerArrayTest> getTest(size_t size)
        {
            return std::make_shared<InnerArrayTest>(size, &outPerform_);
        }

        void recordDependencePlot(size_t pointNo)
        {
            std::vector<OutStruct> outDeriv(pointNo);
            std::vector<OutStruct> outValue(pointNo);
            auto test = getTest(pointNo);
            test->recordTape();
            test->calcReverse();
            test->calcAnalytical();
            for (Size i = 0; i < pointNo; i++)
            {
                outDeriv[i] = { test->input_[i], test->reverseResults_[i], test->analyticalResults_[i] };
                outValue[i] = { test->input_[i], test->output_[i], test->analyticalY_[i] };
            }
            outDeriv_ << outDeriv;
            outValue_ << outValue;
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output outDeriv_;
        cl::tape_empty_test_output outValue_;
    };


    struct MixedTest
        : cl::AdjointTest<MixedTest>
    {
        typedef cl::tape_value inner_type;
        typedef cl::tape_wrapper<inner_type> tape_type;
        typedef cl::tape_function<inner_type> func_type;
        static const CalcMethod default_method = reverse;

        explicit MixedTest(Size size, cl::tape_empty_test_output* logger = nullptr)
            : AdjointTest()
            , size_(size)
            , input_(size)
            , output_(size)
            , X_()
            , Y_()
            , analyticalY_(size)
            , adf_()
            , x_val_()
        {
            setLogger(logger);
            for (size_t i = 0; i < size_; i++)
            {
                input_[i] = 100 * (double(i) / size_);
            }

            x_val_ = std::valarray<double>(input_.data(), size);
            X_ = { tape_type() };
        }

        void recordTape()
        {
            cl::Independent(X_);
            calcY();
            adf_ = std::make_unique<func_type>();
            adf_->Dependent(X_, Y_);
            typedef std::vector<inner_type> Vec;
            auto res = adf_->Forward(0, Vec{ x_val_ })[0].array_value_;
            output_.assign(std::begin(res), std::end(res));
        }

        void calcY()
        {
            Y_ = { targetFunc(X_[0]) };
        }

        template <class T>
        inline static T targetFunc(const T& x)
        {
            return TARGET_FUNCTION(x);
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
#if defined BLACK_SCHOLES_FUNC
            // We can not use inner_type in QuantLib.
            // TODO: make real calculation here.
            auto const& adjRes = results<default_method>();
            analyticalResults_.assign(std::begin(adjRes), std::end(adjRes));
#else
            double h = finite_diff_step;
            inner_type left = targetFunc(x_val_ - h);
            inner_type right = targetFunc(x_val_ + h);
            auto findiff = ((right - left) / (2 * h)).array_value_;
            auto y = ((right + left) / 2).array_value_;
            analyticalResults_.assign(std::begin(findiff), std::end(findiff));
            analyticalY_.assign(std::begin(y), std::end(y));
#endif
        }

        void calcForward()
        {
            typedef std::vector<inner_type> Vec;
            Vec dy = adf_->Forward(1, Vec{ 1 });
            setForwardResults(dy[0]);
        }

        void calcReverse()
        {
            typedef std::vector<inner_type> Vec;
            Vec dy = adf_->Reverse(1, Vec{ 1 });
            setReverseResults(dy[0]);
        }

        void setForwardResults(const inner_type& res)
        {
            auto& result = res.array_value_;
            forwardResults_.assign(std::begin(result), std::end(result));
        }

        void setReverseResults(const inner_type& res)
        {
            auto& result = res.array_value_;
            reverseResults_.assign(std::begin(result), std::end(result));
        }

        size_t memory()
        {
            if (!adf_)
            {
                return 0;
            }
            return adf_->Memory();
        }

        Size indepVarNumber()
        {
            return size_;
        }

        size_t minPerfIteration() { return 10 * iterNumFactor; }

        size_t tapePerfIteration()
        {
            return iterNumFactor / size_ + 1;
        }

        double relativeTol() const { return 0.05; }

        double absTol() const { return 1e-5; }

        Size size_;
        std::vector<double> input_;
        std::vector<double> output_;
        std::vector<tape_type> X_;
        std::vector<tape_type> Y_;
        std::vector<double> analyticalY_;
        std::unique_ptr<func_type> adf_;
        inner_type x_val_;
    };


    struct MixedTestData
    {
        explicit MixedTestData()
        : outPerform_(OUTPUT_FOLDER_NAME "//Mixed " TARGET_FUNCTION_NAME, {
            { "title", TARGET_FUNCTION_NAME " differentiation performance with Mixed optimization." }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "line_box_width", "-5" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        })
            , outAdjoint_(OUTPUT_FOLDER_NAME "//Mixed " TARGET_FUNCTION_NAME, {
                { "title", TARGET_FUNCTION_NAME " Adjoint differentiation performance with Mixed optimization." }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        })
            , outSize_(OUTPUT_FOLDER_NAME "//Mixed " TARGET_FUNCTION_NAME, {
                { "title", "Tape size dependence on number of variables"
                "\\n (memory dynamically allocated for arrays is not included)" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
            , { "cleanlog", "false" }
        })
            , outDeriv_(OUTPUT_FOLDER_NAME "//Mixed " TARGET_FUNCTION_NAME "//output", {
                { "filename", "Derivative" }
            , { "not_clear", "Not" }
            , { "title", "Function derivative" }
            , { "ylabel", "Y" }
            , { "xlabel", "X" }
            , { "cleanlog", "false" }
        })
            , outValue_(OUTPUT_FOLDER_NAME "//Mixed " TARGET_FUNCTION_NAME "//output", {
                { "filename", "FunctionValue" }
            , { "not_clear", "Not" }
            , { "title", "Function value" }
            , { "ylabel", "Y" }
            , { "xlabel", "X" }
            , { "cleanlog", "false" }
        })
        { }

        bool makeOutput()
        {
            outPerform_.log() << "Start MixedTestData output.\n";
            bool ok = true;
            if (pointNo > 0)
            {
                recordDependencePlot(pointNo);
            }
            ok &= cl::recordPerformance(*this, iterNo, step);
            return ok;
        }

        std::shared_ptr<MixedTest> getTest(size_t size)
        {
            return std::make_shared<MixedTest>(size, &outPerform_);
        }

        void recordDependencePlot(size_t pointNo)
        {
            std::vector<OutStruct> outDeriv(pointNo);
            std::vector<OutStruct> outValue(pointNo);
            auto test = getTest(pointNo);
            test->recordTape();
            test->calcReverse();
            test->calcAnalytical();
            for (Size i = 0; i < pointNo; i++)
            {
                outDeriv[i] = { test->input_[i], test->reverseResults_[i], test->analyticalResults_[i] };
                outValue[i] = { test->input_[i], test->output_[i], test->analyticalY_[i] };
            }
            outDeriv_ << outDeriv;
            outValue_ << outValue;
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output outDeriv_;
        cl::tape_empty_test_output outValue_;
    };

#ifdef PRINT_CHECKPOINT
    struct CheckpointTest
        : cl::AdjointTest<CheckpointTest>
    {
        static const CalcMethod default_method = reverse;

        typedef CppAD::AD<double> AD;
        typedef CppAD::ADFun<double> ADFun;
        typedef std::vector<AD> ADVector;

        static inline void algo(const ADVector& x, ADVector& y)
        {
            AD val = TARGET_FUNCTION(x[0]);
            y[0] = val;
        }

        explicit CheckpointTest(Size size, cl::tape_empty_test_output* logger = nullptr)
            : AdjointTest()
            , size_(size)
            , input_(size)
            , output_(size)
            , X_(size)
            , Y_(1)
            , analyticalY_(size)
            , afun_("Checkpoint for target function", algo, ADVector(1), ADVector(1))
            , ff_()
        {
            setLogger(logger);
            for (size_t i = 0; i < size_; i++)
            {
                input_[i] = 100 * (double(i) / size_);
            }
            X_.assign(input_.begin(), input_.end());
        }

        void recordTape()
        {
            Independent(X_);
            calcY();
            ff_ = std::make_unique<ADFun>(X_, Y_);
        }

        void calcY()
        {
            Y_[0] = 0;
            for (size_t i = 0; i < size_; i++)
            {
                ADVector y(1);
                afun_(ADVector{ X_[i] }, y);
                Y_[0] += output_[i] = y[0];
            }
        }

        inline static Real targetFunc(const Real& x)
        {
            return TARGET_FUNCTION(x);
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            Real h = finite_diff_step;
            analyticalResults_.resize(size_);
            analyticalY_.resize(size_);
            for (Size i = 0; i < size_; i++)
            {
                Real left = targetFunc(input_[i] - h);
                Real right = targetFunc(input_[i] + h);
                analyticalResults_[i] = (right - left) / (2 * h);
                analyticalY_[i] = (right + left) / 2;
            }
        }

        void calcReverse()
        {
            reverseResults_ = ff_->Reverse(1, std::vector<double>(1, 1));
        }

        size_t memory()
        {
            if (!ff_)
            {
                return 0;
            }
            return ff_->Memory();
        }

        Size indepVarNumber()
        {
            return size_;
        }

        size_t minPerfIteration() { return iterNumFactor / 2; }

        size_t tapePerfIteration()
        {
            return minPerfIteration() / size_ / 5 + 1;
        }

        double relativeTol() const { return 0.05; }

        double absTol() const { return 1e-5; }

        Size size_;
        std::vector<double> input_;
        std::vector<AD> output_;
        std::vector<AD> X_;
        std::vector<AD> Y_;
        std::vector<Real> analyticalY_;
        CppAD::checkpoint<double> afun_;
        std::unique_ptr<ADFun> ff_;
    };


    struct CheckpointTestData
    {
        explicit CheckpointTestData()
        : outPerform_(OUTPUT_FOLDER_NAME "//Checkpoint " TARGET_FUNCTION_NAME, {
            { "title", TARGET_FUNCTION_NAME " differentiation performance with Checkpoint optimization." }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "line_box_width", "-5" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        })
            , outAdjoint_(OUTPUT_FOLDER_NAME "//Checkpoint " TARGET_FUNCTION_NAME, {
                { "title", TARGET_FUNCTION_NAME " Adjoint differentiation performance with Checkpoint optimization." }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        })
            , outSize_(OUTPUT_FOLDER_NAME "//Checkpoint " TARGET_FUNCTION_NAME, {
                { "title", "Tape size dependence on number of variables"
                "\\n (memory dynamically allocated for arrays is not included)" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
            , { "cleanlog", "false" }
        })
            , outDeriv_(OUTPUT_FOLDER_NAME "//Checkpoint " TARGET_FUNCTION_NAME "//output", {
                { "filename", "Derivative" }
            , { "not_clear", "Not" }
            , { "title", "Function derivative" }
            , { "ylabel", "Y" }
            , { "xlabel", "X" }
            , { "cleanlog", "false" }
        })
            , outValue_(OUTPUT_FOLDER_NAME "//Checkpoint " TARGET_FUNCTION_NAME "//output", {
                { "filename", "FunctionValue" }
            , { "not_clear", "Not" }
            , { "title", "Function value" }
            , { "ylabel", "Y" }
            , { "xlabel", "X" }
            , { "cleanlog", "false" }
        })
        { }

        bool makeOutput()
        {
            outPerform_.log() << "Start CheckpointTestData output.\n";
            bool ok = true;
            if (pointNo > 0)
            {
                recordDependencePlot(pointNo);
            }
            ok &= cl::recordPerformance(*this, iterNo, step);
            return ok;
        }

        std::shared_ptr<CheckpointTest> getTest(size_t size)
        {
            return std::make_shared<CheckpointTest>(size, &outPerform_);
        }

        void recordDependencePlot(size_t pointNo)
        {
            std::vector<OutStruct> outDeriv(pointNo);
            std::vector<OutStruct> outValue(pointNo);
            auto test = getTest(pointNo);
            test->recordTape();
            test->calcReverse();
            test->calcAnalytical();
            for (Size i = 0; i < pointNo; i++)
            {
                outDeriv[i] = { test->input_[i], test->reverseResults_[i], test->analyticalResults_[i] };
                outValue[i] = { test->input_[i], test->output_[i], test->analyticalY_[i] };
            }
            outDeriv_ << outDeriv;
            outValue_ << outValue;
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output outDeriv_;
        cl::tape_empty_test_output outValue_;
    };
#endif PRINT_CHECKPOINT


    struct MethodsCompare
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "N", "No Optimization", "OnePoint", "InnerArray", "Mixed"
#ifdef PRINT_CHECKPOINT
                , "Checkpoint"
#endif // PRINT_CHECKPOINT
                , "Finite diff. (double)", "Finite diff. (Real)", "Finite diff. (InnerArray)"
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type& operator<<(stream_type& stm, MethodsCompare& v)
        {
            stm << v.n_
                << ";" << v.NoOpt_
                << ";" << v.OnePoint_
                << ";" << v.InnerArrayOpt_
                << ";" << v.Mixed_
#ifdef PRINT_CHECKPOINT
                << ";" << v.Checkpoint_
#endif // PRINT_CHECKPOINT
                << ";" << v.double_
                << ";" << v.Real_
                << ";" << v.FinInner_
                << std::endl;
            return stm;
        }

        Real n_;
        Real NoOpt_;
        Real OnePoint_;
        Real InnerArrayOpt_;
        Real Mixed_;
        Real Checkpoint_;
        Real double_;
        Real Real_;
        Real FinInner_;
    };


    void compareSize()
    {
        cl::tape_empty_test_output outTotal(OUTPUT_FOLDER_NAME "//Compare " TARGET_FUNCTION_NAME, {
            { "title", TARGET_FUNCTION_NAME " Adjoint differentiation and tape recording"
            "\\n and function construction (including Forward(0) sweep) performance." }
            , { "filename", "Total" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "true" }
        });
        cl::tape_empty_test_output outTotalSmooth(OUTPUT_FOLDER_NAME "//Compare " TARGET_FUNCTION_NAME, {
            { "title", TARGET_FUNCTION_NAME " Adjoint differentiation and tape recording"
            "\\n and function construction (including Forward(0) sweep) performance." }
            , { "filename", "TotalSmooth" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        });
        cl::tape_empty_test_output outAdjoint(OUTPUT_FOLDER_NAME "//Compare " TARGET_FUNCTION_NAME, {
            { "title", TARGET_FUNCTION_NAME " Adjoint differentiation performance." }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "false" }
        });
        cl::tape_empty_test_output outAdjointSmooth(OUTPUT_FOLDER_NAME "//Compare " TARGET_FUNCTION_NAME, {
            { "title", TARGET_FUNCTION_NAME " Adjoint differentiation performance." }
            , { "filename", "AdjointSmooth" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        });
        cl::tape_empty_test_output outTape(OUTPUT_FOLDER_NAME "//Compare " TARGET_FUNCTION_NAME, {
            { "title", TARGET_FUNCTION_NAME " tape recording performance." }
            , { "filename", "Tape" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "false" }
            //, { "smooth", "default" }
        });
        cl::tape_empty_test_output outSize(OUTPUT_FOLDER_NAME "//Compare " TARGET_FUNCTION_NAME, {
            { "title", "Tape size dependence on number of variables"
            "\\n (memory dynamically allocated for arrays is not included)" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (KB)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "false" }
        });

        size_t graphSize = iterNo + step - 1;
        std::vector<MethodsCompare> total(graphSize)
            , adjoint(graphSize)
            , tape(graphSize)
            , mem(graphSize);
        for (size_t i = 0; i < graphSize; i++)
        {
            size_t n;
            if (i < step)
            {
                n = i + 1;
            }
            else
            {
                n = (i - step + 2) * step;
            }
            auto no = NoOptTest(n, &outTotal);
            auto one = OnePointTest(n, &outTotal);
            auto inn = InnerArrayTest(n, &outTotal);
            auto mix = MixedTest(n, &outTotal);
#ifdef PRINT_CHECKPOINT
            auto chp = CheckpointTest(n, &outTotal);
#else
            auto chp = MixedTest(n, &outTotal);
#endif

            tape[i] = {
                Real(n)
                , no.testTapePerformance()
                , one.testTapePerformance()
                , inn.testTapePerformance()
                , mix.testTapePerformance()
                , chp.testTapePerformance()
                , 0
                , 0
                , 0
            };
            adjoint[i] = {
                Real(n)
                , no.testPerformance<>()
                , one.testPerformance<>()
                , inn.testPerformance<>()
                , mix.testPerformance<>()
                , chp.testPerformance<>()
                , inn.testPerformance<InnerArrayTest::analytical>()
                , no.testPerformance<NoOptTest::analytical>()
                , mix.testPerformance<MixedTest::analytical>()
            };
            total[i] = {
                Real(n)
                , tape[i].NoOpt_            + adjoint[i].NoOpt_
                , tape[i].OnePoint_         + adjoint[i].OnePoint_
                , tape[i].InnerArrayOpt_    + adjoint[i].InnerArrayOpt_
                , tape[i].Mixed_            + adjoint[i].Mixed_
                , tape[i].Checkpoint_       + adjoint[i].Checkpoint_
                , tape[i].double_           + adjoint[i].double_
                , tape[i].Real_             + adjoint[i].Real_
                , tape[i].FinInner_         + adjoint[i].FinInner_
            };
            mem[i] = {
                Real(n)
                , no.memory() / 1024.0
                , one.memory() / 1024.0
                , inn.memory() / 1024.0
                , mix.memory() / 1024.0
                , chp.memory() / 1024.0
                , 0
                , 0
                , 0
            };
        }

        outTotal << total;
        outAdjoint << adjoint;
        outTotalSmooth << total;
        outAdjointSmooth << adjoint;
        outTape << tape;
        outSize << mem;
    }

    void compareFunc()
    {
#ifdef BIG_FUNC
        cl::tape_empty_test_output outTotal(OUTPUT_FOLDER_NAME "//Compare " "Functions", {
            { "title", TARGET_FUNCTION_NAME " Adjoint differentiation and tape recording"
            "\\n and function construction (including Forward(0) sweep) performance." }
            , { "filename", "Total" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Function compexity (number of operation)" }
            , { "cleanlog", "true" }
        });
        cl::tape_empty_test_output outTotalSmooth(OUTPUT_FOLDER_NAME "//Compare " "Functions", {
            { "title", TARGET_FUNCTION_NAME " Adjoint differentiation and tape recording"
            "\\n and function construction (including Forward(0) sweep) performance." }
            , { "filename", "TotalSmooth" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Function compexity (number of operation)" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        });
        cl::tape_empty_test_output outAdjoint(OUTPUT_FOLDER_NAME "//Compare " "Functions", {
            { "title", TARGET_FUNCTION_NAME " Adjoint differentiation performance." }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Function compexity (number of operation)" }
            , { "cleanlog", "false" }
        });
        cl::tape_empty_test_output outAdjointSmooth(OUTPUT_FOLDER_NAME "//Compare " "Functions", {
            { "title", TARGET_FUNCTION_NAME " Adjoint differentiation performance." }
            , { "filename", "AdjointSmooth" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Function compexity (number of operation)" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        });
        cl::tape_empty_test_output outTape(OUTPUT_FOLDER_NAME "//Compare " "Functions", {
            { "title", TARGET_FUNCTION_NAME " tape recording performance." }
            , { "filename", "Tape" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Function compexity (number of operation)" }
            , { "cleanlog", "false" }
            //, { "smooth", "default" }
        });
        cl::tape_empty_test_output outSize(OUTPUT_FOLDER_NAME "//Compare " "Functions", {
            { "title", "Tape size dependence on number of variables"
            "\\n (memory dynamically allocated for arrays is not included)" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
            , { "xlabel", "Function compexity (number of operation)" }
            , { "cleanlog", "false" }
        });

        const size_t max_count = 100;
        std::vector<MethodsCompare> total(max_count)
            , adjoint(max_count)
            , tape(max_count)
            , mem(max_count);
        size_t cont_saved = big_func_count;
        for (size_t i = 0; i < max_count; i++)
        {
            size_t n = 100;
            big_func_count = 1 + i + i*i/10;
            auto no = NoOptTest(n, &outTotal);
            auto one = OnePointTest(n, &outTotal);
            auto inn = InnerArrayTest(n, &outTotal);
            auto mix = MixedTest(n, &outTotal);
#ifdef PRINT_CHECKPOINT
            auto chp = CheckpointTest(n, &outTotal);
#else
            auto chp = MixedTest(n, &outTotal);
#endif

            tape[i] = {
                20 * big_func_count
                , no.testTapePerformance()
                , one.testTapePerformance()
                , inn.testTapePerformance()
                , mix.testTapePerformance()
                , chp.testTapePerformance()
                , 0
                , 0
                , 0
            };
            adjoint[i] = {
                20 * big_func_count
                , no.testPerformance<>()
                , one.testPerformance<>()
                , inn.testPerformance<>()
                , mix.testPerformance<>()
                , chp.testPerformance<>()
                , inn.testPerformance<InnerArrayTest::analytical>()
                , no.testPerformance<NoOptTest::analytical>()
                , mix.testPerformance<MixedTest::analytical>()
            };
            total[i] = {
                20 * big_func_count
                , tape[i].NoOpt_ + adjoint[i].NoOpt_
                , tape[i].OnePoint_ + adjoint[i].OnePoint_
                , tape[i].InnerArrayOpt_ + adjoint[i].InnerArrayOpt_
                , tape[i].Mixed_ + adjoint[i].Mixed_
                , tape[i].Checkpoint_ + adjoint[i].Checkpoint_
                , tape[i].double_ + adjoint[i].double_
                , tape[i].Real_ + adjoint[i].Real_
                , tape[i].FinInner_ + adjoint[i].FinInner_
            };
            mem[i] = {
                20 * big_func_count
                , no.memory() / 1048576.0
                , one.memory() / 1048576.0
                , inn.memory() / 1048576.0
                , mix.memory() / 1048576.0
                , chp.memory() / 1048576.0
                , 0
                , 0
                , 0
            };
        }
        big_func_count = cont_saved;

        outTotal << total;
        outAdjoint << adjoint;
        outTotalSmooth << total;
        outAdjointSmooth << adjoint;
        outTape << tape;
        outSize << mem;
#endif
    }
}

#endif // #ifndef cl_adjoint_array_optimization_impl_hpp
