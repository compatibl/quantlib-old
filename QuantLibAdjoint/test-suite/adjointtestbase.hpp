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

#ifndef cl_adjoint_test_base_hpp
#define cl_adjoint_test_base_hpp
#pragma once

#include "adjointtestutilities.hpp"

#define CL_ADJOINT_TEST_BASE_ERROR(msg)             \
    log_ << "Error in " __FUNCSIG__ << std::endl;   \
    log_ << msg << std::endl;                       \
    BOOST_ERROR(msg)

namespace cl
{
    using QuantLib::Real;
    using QuantLib::PerformanceTime;
    using QuantLib::AdjointTime;


    // Class for easier adjoint testing.
    template <typename Test>
    struct AdjointTestBase
    {
        typedef Test test_type;

        typedef std::vector<double> adjoint_result_type;

        enum CalcMethod
        {
            forward,
            reverse,
            analytical,
            other
        };

        template <class Ostr>
        friend Ostr& operator<<(Ostr& o, CalcMethod method)
        {
            switch (method)
            {
            case forward:
                o << "Forward mode";
                break;
            case reverse:
                o << "Reverse mode";
                break;
            case analytical:
                o << "Analytical";
                break;
            case other:
                o << "Adjoint";
                break;
            }
            return o;
        }

        static const CalcMethod default_method = reverse;

        // Encapsulates work with log.
        struct Log
        {
            Log(cl::tape_empty_test_output* p, const AdjointTestBase* base)
                : logger_(p)
                , that_(base->that())
                , timer_()
                , muted_(false)
            { }

            void setLogger(cl::tape_empty_test_output* p)
            {
                logger_ = p;
            }

            void mute()
            {
                muted_ = true;
            }
            void unmute()
            {
                muted_ = false;
            }

            void beforeTape()
            {
                if (muted_)
                    return;
                *this << "Start of tape recording: " << currentTime() << std::endl;
                timer_.restart();
            }
            void afterTape()
            {
                if (muted_)
                    return;
                *this << "End of tape recording." << std::endl;
                *this << "Time for tape recording: " << timer_.elapsed() << std::endl;
                timer_.restart();
            }

            template <CalcMethod method>
            void before()
            {
                if (muted_)
                    return;
                *this << "Start of differentiation (" << method << "): " << currentTime() << std::endl;
                timer_.restart();
            }

            template <CalcMethod method>
            void after()
            {
                if (muted_)
                    return;
                *this << "End of differentiation (" << method << ")." << std::endl;
                *this << "Time for " << method << ": " << timer_.elapsed() << std::endl;
                timer_.restart();
            }

            void beforeForward()
            {
                if (muted_)
                    return;
                *this << "Start of differentiation in Forward mode: " << currentTime() << std::endl;
                timer_.restart();
            }
            void afterForward()
            {
                if (muted_)
                    return;
                *this << "End of differentiation in Forward mode." << std::endl;
                *this << "Time for Forward mode: " << timer_.elapsed() << std::endl;
                timer_.restart();
            }

            void beforeReverse()
            {
                if (muted_)
                    return;
                *this << "Start of differentiation in Reverse mode: " << currentTime() << std::endl;
                timer_.restart();
            }
            void afterReverse()
            {
                if (muted_)
                    return;
                *this << "End of differentiation in Reverse mode." << std::endl;
                *this << "Time for Reverse mode: " << timer_.elapsed() << std::endl;
                timer_.restart();
            }

            void beforeAnalytical()
            {
                if (muted_)
                    return;
                *this << "Start of differentiation using " + that_->analyticalName() + " method: " << currentTime() << std::endl;
                timer_.restart();
            }
            void afterAnalytical()
            {
                if (muted_)
                    return;
                *this << "End of differentiation using " + that_->analyticalName() + " method." << std::endl;
                *this << "Time for " + that_->analyticalName() + " method: " << timer_.elapsed() << std::endl;
                timer_.restart();
            }

            template<class T>
            Log& operator<<(T&& t)
            {
                if (!muted_ && logger_)
                {
                    logger_->log() << std::forward<T>(t);
                }
                return *this;
            }

            // this is the type of std::cout
            typedef std::basic_ostream<char, std::char_traits<char> > cout_manip_type;

            // this is the function signature of std::endl
            typedef cout_manip_type& (*standard_manipulator)(cout_manip_type&);

            // define an operator<< to take in std::endl
            Log& operator<<(standard_manipulator mnp)
            {
                if (logger_ && !muted_)
                {
                    logger_->log() << mnp;
                }
                return *this;
            }

        private:
            std::string currentTime()
            {
                return QuantLib::currentTime();
            }

            cl::tape_empty_test_output* logger_;
            const test_type* that_;
            boost::timer timer_;
            bool muted_;
        };

        AdjointTestBase()
            : f_()
            , forwardResults_()
            , reverseResults_()
            , analyticalResults_()
            , log_(nullptr, this)
        { }

        test_type* that()
        {
            static_assert(std::is_base_of<AdjointTestBase<Test>, Test>::value
                , "template parameter Test must be"
                " AdjointTestBase<Test> derived class");

            return static_cast<Test*>(this);
        }

        const test_type* that() const
        {
            static_assert(std::is_base_of<AdjointTestBase<Test>, Test>::value
                , "template parameter Test must be"
                " AdjointTestBase<Test> derived class");

            return static_cast<const Test*>(this);
        }

        // Tests adjoint derivative calculation in forward and reverse mode.
        // The tape function and all derivatives will be recalculated.
        // Returns true if derivative consistency check passed, false otherwise.
        template<typename Class = void>
        bool test()
        {
            that()->doTape();
            that()->doCalculation<forward>();
            that()->doCalculation<reverse>();
            that()->doCalculation<analytical>();
            return that()->check();
        }

        // Tests adjoint derivative calculation with one method only.
        // The tape function, tested method derivatives and analytical derivatives will be recalculated.
        // Returns true if derivative consistency check passed, false otherwise.
        template<CalcMethod method = test_type::default_method>
        bool testAdjoint()
        {
            static_assert(method != analytical, "Use this method for testing Adjoint methods ONLY");
            that()->doTape();
            that()->doCalculation<method>();
            that()->doCalculation<analytical>();
            return that()->checkAdjoint<method>();
        }

        // Tests adjoint derivative calculation in reverse mode only.
        // The tape function, reverse mode and analytical derivatives will be recalculated.
        // Returns true if derivative consistency check passed, false otherwise.
        template<typename Class = void>
        bool testReverse()
        {
            return that()->testAdjoint<reverse>();
        }

        // Records performance for tape recording, reverse mode and analytical method (to calculate derivatives).
        // Return true if result consistency check passed.
        template<typename Class = void>
        bool recordPerformance(PerformanceTime& perfTime, AdjointTime& adjointTime)
        {
            log_ << "Performance recording." << std::endl;
            perfTime.indepVarNumber_
                = adjointTime.indepVarNumber_ = that()->indepVarNumber();
            perfTime.timeTapeRecording_ = that()->testTapePerformance();
            perfTime.timeAdjoint_
                = adjointTime.timeAdjoint_ = that()->testPerformance();
            perfTime.timeAnalytical_ = that()->testPerformance<analytical>();
            return that()->checkAdjoint<test_type::default_method>();
        }

        // Records performance for reverse mode (to calculate derivatives).
        template<typename Class = void>
        void recordAdjointTime(AdjointTime& adjointTime)
        {
            adjointTime.indepVarNumber_ = that()->indepVarNumber();
            that()->doTape();
            adjointTime.timeAdjoint_ = that()->testPerformance();
        }

        // Returns performance of tape recording.
        double testTapePerformance(size_t testNo = 0)
        {
            if (testNo == 0)
            {
                testNo = that()->tapePerfIteration();
            }
            log_ << "Tape recording performance testing." << std::endl;
            log_.beforeTape();
            log_.mute();
            boost::timer timer;
            for (size_t i = 0; i < testNo; i++)
            {
                that()->recordTape();
            }
            double elapsed = timer.elapsed();
            log_.unmute();
            log_.afterTape();
            return elapsed / testNo;
        }

        // Returns performance of derivatives calculations.
        template <CalcMethod method = test_type::default_method>
        double testPerformance(size_t testNo = 0)
        {
            if (!testNo)
            {
                testNo = that()->perfIteration<method>();
            }
            if (method != analytical && !f_)
            {
                that()->doTape();
            }
            log_ << method << " performance testing." << std::endl;
            log_.before<method>();
            log_.mute();
            boost::timer timer;
            for (size_t i = 0; i < testNo; i++)
            {
                that()->calculate<method>();
            }
            double elapsed = timer.elapsed();
            log_.unmute();
            log_.after<method>();
            return elapsed / testNo;
        }

        // Calculate Greeks with reverse mode only.
        template<typename Class = void>
        void calcWithReverse()
        {
            that()->doTape();
            that()->doCalculation<reverse>();
        }

        // Returns tape memory size (bytes).
        size_t memory()
        {
            if (!f_)
            {
                return 0;
            }
            return f_->Memory();
        }

        // Set object that provides log stream.
        void setLogger(cl::tape_empty_test_output* p)
        {
            log_.setLogger(p);
        }

        template <CalcMethod method>
        void doCalculation()
        {
            log_.before<method>();
            that()->calculate<method>();
            log_.after<method>();
        }

        // Records tape and sets new tape function.
        template<typename Class = void>
        void doTape()
        {
            log_.beforeTape();
            that()->recordTape();
            log_.afterTape();
        }

        template <CalcMethod method>
        void calculate()
        {
            switch (method)
            {
            case forward:
                that()->calcForward();
                break;
            case reverse:
                that()->calcReverse();
                break;
            case analytical:
                that()->calcAnalytical();
                break;
            default:
                BOOST_FAIL(__FUNCSIG__ " is not implemented.");
            }
        }

        template <>
        void calculate<other>()
        {
            that()->calcAdjoint();
        }

        // Checks forward mode, reverse mode and analytical result consistency.
        template<typename Class = void>
        bool check()
        {
            bool result = true;
            int n = analyticalResults_.size();
            if (n != forwardResults_.size())
            {
                CL_ADJOINT_TEST_BASE_ERROR("Forward mode and " << that()->analyticalName() << " vectors size mismatch.");
                return false;
            }
            if (n != reverseResults_.size())
            {
                CL_ADJOINT_TEST_BASE_ERROR("Reverse mode and " << that()->analyticalName() << " vectors size mismatch.");
                return false;
            }
            adjoint_result_type::const_iterator fwd = forwardResults_.begin();
            adjoint_result_type::const_iterator rev = reverseResults_.begin();
            std::vector<Real  >::const_iterator anl = analyticalResults_.begin();
            log_ << "Derivatives" << std::endl;
            log_ << "Forward mode || Reverse mode || " << that()->analyticalName() << std::endl;
            for (int i = 0; i < n; ++i, ++anl, ++fwd, ++rev)
            {
                log_ << std::setw(12) << std::right << *fwd << " ";
                log_ << std::setw(15) << std::right << *rev << " ";
                log_ << std::setw(15) << std::right << *anl << std::endl;

                // Maximum of absolut values of results.
                Real maxabs = std::max(std::abs(*anl), std::max(std::abs(*fwd), std::abs(*rev)));
                Real tol = std::max(that()->relativeTol() * maxabs, that()->absTol());
                if (std::abs(*fwd - *anl) > tol)
                {
                    result = false;
                    CL_ADJOINT_TEST_BASE_ERROR("Forward mode and " << that()->analyticalName() << " derivative[" << i << "] mismatch."
                        << std::setw(24) << std::left << "\n    Forward mode:" << *fwd
                        << std::setw(24) << std::left << "\n    " + that()->analyticalName() + ":" << *anl
                        << "\n    tolerance:  " << tol);
                }
                if (std::abs(*rev - *anl) > tol)
                {
                    result = false;
                    CL_ADJOINT_TEST_BASE_ERROR("Reverse mode and " << that()->analyticalName() << " derivative[" << i << "] mismatch."
                        << std::setw(24) << std::left << "\n    Reverse mode:" << *rev
                        << std::setw(24) << std::left << "\n    " + that()->analyticalName() + ": " << *anl
                        << "\n    tolerance:  " << tol);
                }
                if (std::abs(*rev - *fwd) > tol)
                {
                    result = false;
                    CL_ADJOINT_TEST_BASE_ERROR("Forward mode and Reverse mode derivative[" << i << "] mismatch."
                        << std::setw(24) << std::left << "\n    Forward mode:" << *fwd
                        << std::setw(24) << std::left << "\n    Reverse mode:" << *rev
                        << "\n    tolerance:  " << tol);
                }
            }
            return result;
        }

        template <CalcMethod method>
        adjoint_result_type& results();

        template <> adjoint_result_type& results<forward>() { return forwardResults_; }
        template <> adjoint_result_type& results<reverse>() { return reverseResults_; }
        template <> adjoint_result_type& results<other>() { return adjointResults_; }

        template <CalcMethod method = test_type::default_method>
        bool checkAdjoint()
        {
            bool ok = true;
            int n = analyticalResults_.size();
            if (n != results<method>().size())
            {
                CL_ADJOINT_TEST_BASE_ERROR(method << " and " << that()->analyticalName() << " vectors size mismatch.");
                return false;
            }
            adjoint_result_type::const_iterator adj = results<method>().begin();
            std::vector<Real  >::const_iterator anl = analyticalResults_.begin();
            log_ << "Derivatives" << std::endl;
            log_ << std::setw(12) << method << " || " << that()->analyticalName() << std::endl;
            for (int i = 0; i < n; ++i, ++anl, ++adj)
            {
                log_ << std::setw(12) << std::right << *adj << " ";
                log_ << std::setw(15) << std::right << *anl << std::endl;

                // Maximum of absolut values of results.
                Real maxabs = std::max(std::abs(*anl), std::abs(*adj));
                Real tol = std::max(that()->relativeTol() * maxabs, that()->absTol());
                if (std::abs(*adj - *anl) > tol)
                {
                    ok = false;
                    CL_ADJOINT_TEST_BASE_ERROR(method << " and " << that()->analyticalName() << " derivative[" << i << "] mismatch."
                        << std::setw(24) << std::left << "\n    Adjoint:" << *adj
                        << std::setw(24) << std::left << "\n    " + that()->analyticalName() + ": " << *anl
                        << "\n    tolerance:  " << tol);
                }
            }
            return ok;
        }

        // Checks reverse mode and analytical result consistency.
        template<typename Class = void>
        bool checkReverse()
        {
            bool result = true;
            int n = analyticalResults_.size();
            if (n != reverseResults_.size())
            {
                CL_ADJOINT_TEST_BASE_ERROR("Reverse mode and " << that()->analyticalName() << " vectors size mismatch.");
                return false;
            }
            adjoint_result_type::const_iterator rev = reverseResults_.begin();
            std::vector<Real  >::const_iterator anl = analyticalResults_.begin();
            log_ << "Derivatives" << std::endl;
            log_ << "Reverse mode || " << that()->analyticalName() << std::endl;
            for (int i = 0; i < n; ++i, ++anl, ++rev)
            {
                log_ << std::setw(12) << std::right << *rev << " ";
                log_ << std::setw(15) << std::right << *anl << std::endl;

                // Maximum of absolut values of results.
                Real maxabs = std::max(std::abs(*anl), std::abs(*rev));
                Real tol = std::max(that()->relativeTol() * maxabs, that()->absTol());
                if (std::abs(*rev - *anl) > tol)
                {
                    result = false;
                    CL_ADJOINT_TEST_BASE_ERROR("Reverse mode and " << that()->analyticalName() << " derivative[" << i << "] mismatch."
                        << std::setw(24) << std::left << "\n    Reverse mode:" << *rev
                        << std::setw(24) << std::left << "\n    " + that()->analyticalName() + ": " << *anl
                        << "\n    tolerance:  " << tol);
                }
            }
            return result;
        }

        std::unique_ptr<cl::tape_function<double> > f_;
        adjoint_result_type forwardResults_;
        adjoint_result_type reverseResults_;
        adjoint_result_type adjointResults_;
        std::vector<Real> analyticalResults_;
        Log log_;
    };

    // Extension of AdjointTestBase class with some default implementation of required methods.
    template <typename Test>
    struct AdjointTest
        : public AdjointTestBase<Test>
    {
        // Default name of analytical method.
        std::string analyticalName() const { return "Finite diff."; }

        // Default reverse calculation method implemetns case of one dimensional range space.
        template<typename Class = void>
        void calcReverse()
        {
            assert(f_->Range() == 1);
            reverseResults_ = f_->Reverse(1, std::vector<double>(1, 1));
        }

        // Default forward calculation method implemetns case of one dimensional range space.
        template<typename Class = void>
        void calcForward()
        {
            assert(f_->Range() == 1);
            forwardResults_.resize(that()->indepVarNumber());
            // Direction for directional derivative.
            std::vector<double> dX(f_->Domain(), 0);
            std::vector<double>::iterator px = dX.begin();
            std::generate(forwardResults_.begin(), forwardResults_.end(), [&]()
            {
                *px = 1;
                double temp = f_->Forward(1, dX)[0];
                *px++ = 0;
                return temp;
            });
        }

        // Relative tolerance for checking consistensy with analytical results.
        double relativeTol() const { return 1e-10; }

        // Absolute tolerance for checking consistensy with analytical results.
        double absTol() const { return 1e-10; }

        // Minimum iterations needed for one dimantional problem performance measuring.
        // It is used in other *PerfIteration() methods and scaled with indepVarNumber() and depVarNumber().
        // Overload this method to determinate the number of iterations for all performance measuring.
        size_t minPerfIteration() { return 0; }

        // Overload this method if your problem has more than one dependent variables.
        size_t depVarNumber()
        {
            assert(f_ ? (f_->Range() == 1) : true);
            return 1;
        }

        // Default number of tests for tape recording performance measuring.
        template<typename Class = void>
        size_t tapePerfIteration()
        {
            // recordTape() has at least O(depVarNumber + indepVarNumber) complexity.
            return  that()->minPerfIteration()
                / (that()->indepVarNumber() + that()->depVarNumber()) + 1;
        }

        // Default number of tests for performance measuring.
        template <CalcMethod method>
        size_t perfIteration();

        // Default number of tests for forward mode performance measuring.
        template<>
        size_t perfIteration<forward>()
        {
            // calcForward() has at least O(indepVarNumber *(indepVarNumber + depVarNumber)) complexity.
            return  that()->minPerfIteration() / that()->indepVarNumber()
                / (that()->indepVarNumber() + that()->depVarNumber()) + 1;
        }

        // Default number of tests for reverse mode performance measuring.
        template<>
        size_t perfIteration<reverse>()
        {
            // calcReverse() has at least O(depVarNumber *(indepVarNumber + depVarNumber)) complexity.
            return  that()->minPerfIteration() / that()->depVarNumber()
                / (that()->indepVarNumber() + that()->depVarNumber()) + 1;
        }

        // Default number of tests for analytical method performance measuring.
        template<>
        size_t perfIteration<analytical>()
        {
            // calcAnalytical() has at least O(depVarNumber + indepVarNumber) complexity.
            return  that()->minPerfIteration()
                / (that()->indepVarNumber() + that()->depVarNumber()) + 1;
        }

        // Default number of tests for other adjoint methods performance measuring.
        template<>
        size_t perfIteration<other>()
        {
            // calcReverse() has at least O(indepVarNumber + depVarNumber) complexity.
            return  that()->minPerfIteration()
                / (that()->indepVarNumber() + that()->depVarNumber()) + 1;
        }
    };


    struct SimpleTest
        : AdjointTest<SimpleTest>
    {};


    // Mesures performance for forward mode and analytical method.
    // Records plots for tape recording, adjoint and analytical method performance (to factory.outPerform_);
    // adjoint performance (to factory.outAdjoint_);
    // tape size (to factory.outSize_).
    // parameters:
    // factory - struct with output stream fields and test constructing method,
    // n - number of performanse mesuaring,
    // step - step for problem size (used for factory.getTest(i * step) call).
    // Returns true if all calculated derivatives checking pass.
    template <class TestFactory>
    bool recordPerformance(TestFactory& factory, size_t n, size_t step = 1)
    {
        typedef std::remove_reference_t<decltype(*factory.getTest(1))> test_type;
        static_assert(
            std::is_base_of<
                AdjointTestBase<test_type>
                , test_type
            >::value
            , "TestFactory class must have getTest(size_t) method that returns"
            " (smart) pointer to AdjointTestBase<Test> derived class");

        std::vector<PerformanceTime> perfTimes(n);
        std::vector<AdjointTime> adjointTimes(n);
        std::vector<TapeSize> tapeMemory(n);
        auto p = perfTimes.begin();
        auto a = adjointTimes.begin();
        auto m = tapeMemory.begin();
        bool ok = true;
        for (size_t i = 1; i <= n; i++)
        {
            std::shared_ptr<test_type> test(factory.getTest(i * step));
            ok &= test->recordPerformance(*p++, *a++);
            *m++ = { test->indepVarNumber(), test->memory() };
        }
        factory.outPerform_ << perfTimes;
        factory.outAdjoint_ << adjointTimes;
        factory.outSize_ << tapeMemory;
        return ok;
    }


    // Mesures performance for analytical method.
    // Records plot for adjoint performance (to factory.outAdjoint_).
    // parameters:
    // factory - struct with output stream field and test constructing method,
    // n - number of performanse mesuaring,
    // step - step for problem size (used for factory.getTest(i * step) call).
    template<typename TestFactory>
    void recordAdjointTime(TestFactory& factory, size_t n, size_t step = 1)
    {
        typedef std::remove_reference_t<decltype(*factory.getTest(1))> test_type;
        static_assert(
            std::is_base_of<
                AdjointTestBase<test_type>
                , test_type
            >::value
            , "TestFactory class must have getTest(size_t) method that returns"
            " (smart) pointer to AdjointTestBase<Test> derived class");

        std::vector<AdjointTime> adjointTimes(n);
        auto p = adjointTimes.begin();
        for (size_t i = 1; i <= n; i++, ++p)
        {
            std::shared_ptr<test_type> test(factory.getTest(i * step));
            test->recordAdjointTime(*p);
        }

        factory.outAdjoint_ << adjointTimes;
    }
}

#endif