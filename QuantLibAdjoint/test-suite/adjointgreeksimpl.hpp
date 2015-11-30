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


#ifndef cl_adjoint_greeks_test_impl_hpp
#define cl_adjoint_greeks_test_impl_hpp
#pragma once

#include <random>

#include <ql/quantlib.hpp>
#include <cl/tape/impl/ql/blackcalculator_t.hpp>

#include "adjointgreekstest.hpp"
#include "utilities.hpp"
#include "adjointtestbase.hpp"

using namespace QuantLib;
using namespace boost::unit_test_framework;

#define OUTPUT_FOLDER_NAME "AdjointGreeks"


namespace
{
    enum
    {
#if defined CL_GRAPH_GEN
        // Number of points for performance graphic.
        iterNo = 100,
        // Step for portfolio size for performance testing .
        step = 10,
#else
        // Number of points for performance graphic.
        iterNo = 0,
        // Step for portfolio size for performance testing .
        step = 1,
#endif
        // Defines performance accuracy. Its value is a minimum number
        // of calling of O(1) complexity methods per one performance test.
#if defined NDEBUG
        iterNumFactor = 500000,
#else
        iterNumFactor = 500,
#endif
    };

    // Parameters of single option.
    enum DataType
    {
        Strike
        , Stock
        , Sigma
        , Time
        , Rate
    };

    template <class T>
    inline T black_scholes_price(const std::vector<T>& data, size_t start = 0, size_t stride = 1)
    {
        T const& strike = data[start + stride * Strike];
        T const& stock  = data[start + stride * Stock];
        T const& sigma  = data[start + stride * Sigma];
        T const& time   = data[start + stride * Time];
        T const& rate   = data[start + stride * Rate];

        return black_scholes_price(strike, stock, sigma, time, rate);
    }

    template <class T>
    inline T black_scholes_price(T const& strike, T const& stock, T const& sigma, T const& time, T const& rate)
    {
        T inversed_discount_factor = std::exp(rate * time);
        T discount_factor = 1 / inversed_discount_factor;
        T forward = stock * inversed_discount_factor;
        T stdDev = sigma * std::sqrt(time);

        return BlackCalculator_T<T>(Option::Call, strike, forward, stdDev, discount_factor).value();
    }


    // Struct for graphics recording.
    struct OutStruct
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "X", "Adjoint (optimized)", "Brute force"
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

        double x_;
        double adj_;
        double thr_;
    };


    struct GreeksTest
        : cl::AdjointTest<GreeksTest>
    {
        typedef cl::tvalue inner_type;
        typedef cl::tape_wrapper<inner_type> tape_type;
        typedef cl::tfunc<inner_type> func_type;
        static const CalcMethod default_method = reverse;

        explicit GreeksTest(Size size, cl::tape_empty_test_output* logger = nullptr)
            : AdjointTest()
            , size_(size)
            , input_(5 * size)
            , X_(5)
            , Y_()
            , func_()
            , x_val_(5)
            , reverse_()
        {
            setLogger(logger);
            gen_data();
        }

        void gen_data()
        {            
            std::uniform_real_distribution<> distribution(0.98, 1.02);
            std::mt19937 gen;
            gen.seed(0);

            for (size_t i = 0; i < size_; i++)
            {
                double & strike = input_[i + size_ * Strike];
                double & stock  = input_[i + size_ * Stock];
                double & sigma  = input_[i + size_ * Sigma];
                double & time   = input_[i + size_ * Time];
                double & rate   = input_[i + size_ * Rate];

                strike = distribution(gen) * 120.0;
                stock  = distribution(gen) * 100.0;
                sigma  = distribution(gen) * 0.2;
                time   = distribution(gen) * 1.0;
                rate   = distribution(gen) * 0.05;
            }

            for (size_t i = 0; i < 5; i++)
            {
                X_[i] = input_[i * size_];
                x_val_[i] = std::valarray<double>(input_.data() + i * size_, size_);
            }
        }

        void recordTape()
        {
            cl::tape_start(X_);
            calcY();
            func_ = std::make_unique<func_type>();
            func_->dependent(X_, Y_);
            func_->forward(0, x_val_);
        }

        void calcY()
        {
            Y_ = { black_scholes_price(X_) };
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            std::vector<double> prices(size_);
            for (size_t i = 0; i < size_; i++)
            {
                prices[i] = black_scholes_price(input_, i, size_);
            }

            double h = 1e-7;
            analyticalResults_.resize(4 * size_);
            for (size_t i = size_; i < 5 * size_; i++)
            {
                double temp = input_[i];
                double bump = h * input_[i];
                input_[i] += bump;
                double bumped_price = black_scholes_price(input_, i % size_, size_);
                analyticalResults_[i - size_] = (bumped_price - prices[i % size_]) / bump;
                input_[i] = temp;
            }
        }

        void calcReverse()
        {
            typedef std::vector<inner_type> Vec;
            reverse_ = func_->reverse(1, Vec{ 1 });
            setReverseResults(reverse_);
        }

        void setReverseResults(std::vector<inner_type> const& reverse)
        {
            if (reverse.size() == 0)
                return;

            reverseResults_.resize(4 * size_);
            for (size_t i = 0; i < 4; i++)
            {
                auto& result = reverse[i + 1].array_value_;
                std::copy(std::begin(result), std::end(result)
                    , reverseResults_.begin() + i * size_);
            }
        }
        
        //template <CalcMethod method = default_method>
        //bool checkAdjoint()
        //{
        //    setReverseResults(reverse_);
        //    return cl::AdjointTest<GreeksTest>::checkAdjoint<reverse>();
        //}

        size_t memory()
        {
            if (!func_)
            {
                return 0;
            }
            return func_->Memory();
        }

        size_t indepVarNumber()
        {
            return size_;
        }

        size_t minPerfIteration() { return iterNumFactor; }

        size_t tapePerfIteration()
        {
            return iterNumFactor / size_ / 10 + 1;
        }

        double relativeTol() const { return 1e-5; }

        double absTol() const { return 1e-8; }

        Size size_;
        std::vector<double> input_;
        std::vector<tape_type> X_;
        std::vector<tape_type> Y_;
        std::unique_ptr<func_type> func_;
        std::vector<inner_type> x_val_;
        std::vector<inner_type> reverse_;
    };


    struct GreeksTestData
    {
        explicit GreeksTestData()
        : outPerform_(OUTPUT_FOLDER_NAME, {
            { "title", "Performance of calculation Delta, Vega, Theta and Rho greeks for portfolio." }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Size of portfolio" }
            , { "line_box_width", "-5" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        })
            , outAdjoint_(OUTPUT_FOLDER_NAME, {
                { "title", "Performance of calculation Delta, Vega, Theta and Rho greeks for portfolio." }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Size of portfolio" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        })
            , outTotal_(OUTPUT_FOLDER_NAME, {
                { "title", "Performance of calculation Delta, Vega, Theta and Rho greeks for portfolio." }
            , { "filename", "Total" }
            , { "not_clear", "Not" }
            , { "line_box_width", "-4" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Size of portfolio" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        })
            , outSize_(OUTPUT_FOLDER_NAME, {
                { "title", "Tape size dependence on size of portfolio"
                "\\n (memory dynamicly allocated for arrays is not included)" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
            , { "cleanlog", "false" }
        })
        { }

        bool makeOutput()
        {
            outPerform_.log() << "Start GreeksTestData output.\n";

            std::vector<PerformanceTime> perfTimes(iterNo);
            std::vector<AdjointTime> adjointTimes(iterNo);
            std::vector<TapeSize> tapeMemory(iterNo);
            std::vector<OutStruct> total(iterNo);
            auto p = perfTimes.begin();
            auto a = adjointTimes.begin();
            auto m = tapeMemory.begin();
            auto t = total.begin();
            bool ok = true;
            for (size_t i = 1; i <= iterNo; i++, ++p)
            {
                std::shared_ptr<GreeksTest> test(getTest(i * step));
                ok &= test->recordPerformance(*p, *a++);
                *m++ = { test->indepVarNumber(), test->memory() };
                *t++ = {
                    i * step
                    , p->timeAdjoint_ + p->timeTapeRecording_
                    , p->timeAnalytical_
                };
            }
            outPerform_ << perfTimes;
            outAdjoint_ << adjointTimes;
            outSize_ << tapeMemory;
            outTotal_ << total;
            return ok;
        }

        std::shared_ptr<GreeksTest> getTest(size_t size)
        {
            return std::make_shared<GreeksTest>(size, &outPerform_);
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outTotal_;
        cl::tape_empty_test_output outSize_;
    };

}

#endif
