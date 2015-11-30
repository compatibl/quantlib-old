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

#ifndef quantlib_adjoint_test_helper_hpp
#define quantlib_adjoint_test_helper_hpp


#include "adjointtestutilities.hpp"

namespace cl
{
    // Abstract class that define derivative caculation method.
    class DerivCalculator
    {
    public:
        ~DerivCalculator() {}
        virtual void calcForward(cl::TapeFunction<double>&, std::vector<double>&) = 0;
        virtual void calcReverse(cl::TapeFunction<double>&, std::vector<double>&) = 0;
    };

    // Calculates first derivatives.
    class FirstDerivCalculator
        : public DerivCalculator
    {
    public:
        FirstDerivCalculator(int signFactor = 1) : signFactor_(signFactor) {}

        void calcForward(cl::TapeFunction<double>& f, std::vector<double>& forwardResults)
        {
            forwardResults.resize(f.Domain());
            // Direction for directional derivative.
            std::vector<double> dX(f.Domain(), 0);
            std::vector<double>::iterator px = dX.begin();
            std::generate(forwardResults.begin(), forwardResults.end(), [&]()
            {
                *px = 1;
                double temp = signFactor_ * f.Forward(1, dX)[0];
                *px++ = 0;
                return temp;
            });
        }

        void calcReverse(cl::TapeFunction<double>& f, std::vector<double>& reverseResults)
        {
            reverseResults = f.Reverse(1, std::vector<double>(1, 1));
            for (double& res : reverseResults)
            {
                res *= signFactor_;
            }
        }

    private:
        int signFactor_;
    };

    // Calculates second derivatives, assuming that all mixed derivatives are zero.
    class SecondDerivCalculator
        : public DerivCalculator
    {
    public:
        SecondDerivCalculator(int signFactor = 1) : signFactor_(signFactor) {}

        void calcForward(cl::TapeFunction<double>& f, std::vector<double>& forwardResults)
        {
            forwardResults.resize(f.Domain());
            // Direction for directional derivative.
            std::vector<double> dX(f.Domain(), 0);
            std::vector<double>::iterator px = dX.begin();
            std::generate(forwardResults.begin(), forwardResults.end(), [&]()
            {
                *px = 1;
                f.Forward(1, dX);
                *px++ = 0;
                return signFactor_ * 2 * f.Forward(2, dX)[0];
            });
        }

        void calcReverse(cl::TapeFunction<double>& f, std::vector<double>& reverseResults)
        {
            const size_t n = f.Domain();
            f.Forward(1, std::vector<double>(n, 1));
            std::vector<double> dw = f.Reverse(2, std::vector<double>(1, 1));
            reverseResults.resize(n);
            for (size_t i = 0; i < n; i++)
            {
                reverseResults[i] = signFactor_ * dw[2 * i + 1];
            }
        }

    private:
        int signFactor_;
    };

    // Calculates third derivatives, assuming that all mixed derivatives are zero.
    class ThirdDerivCalculator
        : public DerivCalculator
    {
    public:
        ThirdDerivCalculator(int signFactor = 1) : signFactor_(signFactor) {}

        void calcForward(cl::TapeFunction<double>& f, std::vector<double>& forwardResults)
        {
            forwardResults.resize(f.Domain());
            // Direction for directional derivative.
            std::vector<double> dX(f.Domain(), 0);
            std::vector<double>::iterator px = dX.begin();
            std::generate(forwardResults.begin(), forwardResults.end(), [&]()
            {
                *px = 1;
                f.Forward(1, dX);
                *px++ = 0;
                f.Forward(2, dX);
                return signFactor_ * 6 * f.Forward(3, dX)[0];
            });
        }

        void calcReverse(cl::TapeFunction<double>& f, std::vector<double>& reverseResults)
        {
            const size_t n = f.Domain();
            f.Forward(1, std::vector<double>(n, 1));
            f.Forward(2, std::vector<double>(n, 0));
            std::vector<double> dw = f.Reverse(3, std::vector<double>(1, 1));
            reverseResults.resize(n);
            for (size_t i = 0; i < n; i++)
            {
                reverseResults[i] = signFactor_ * 2 * dw[3 * i + 2];
            }
        }

    private:
        int signFactor_;
    };

    // Calculates second order mixed derivatives with respect to i-th and (i+n)-th
    // variables from 2*n dimensional domain vector. For correct work all other
    // mixed derivatives should be zero.
    class SecondMixedCalculator
        : public DerivCalculator
    {
    public:
        SecondMixedCalculator(int signFactor = 1) : signFactor_(signFactor) {}

        void calcForward(cl::TapeFunction<double>& f, std::vector<double>& forwardResults)
        {
            size_t n = f.Domain() / 2;
            forwardResults.resize(n);
            // Direction for directional derivative.
            std::vector<double> dX(f.Domain(), 0);
            std::vector<double>::iterator px = dX.begin();
            std::vector<double>::iterator pxn = dX.begin() + n;
            std::generate(forwardResults.begin(), forwardResults.end(), [&]()
            {
                // Mixed derivative calculated from second order derivatives like this:
                // z = x + y
                // w = x - y
                // ddF/dzdz = ddF/dxdx + 2 * ddF/dxdy + ddF/dydy
                // ddF/dwdw = ddF/dxdx - 2 * ddF/dxdy + ddF/dydy
                // ddF/dxdy = (ddF/dzdz - ddF/dwdw) / 4
                *px = 1;
                *pxn = 1;
                f.Forward(1, dX);
                *px = 0;
                *pxn = 0;
                double dzdz = 2 * f.Forward(2, dX)[0];
                *px = 1;
                *pxn = -1;
                f.Forward(1, dX);
                *px++ = 0;
                *pxn++ = 0;
                double dwdw = 2 * f.Forward(2, dX)[0];
                return signFactor_ * (dzdz - dwdw) / 4;
            });
        }

        void calcReverse(cl::TapeFunction<double>& f, std::vector<double>& reverseResults)
        {
            const size_t n = f.Domain() / 2;
            std::vector<double> dX(f.Domain(), 0);
            std::fill(dX.begin(), dX.begin() + n, 1);
            f.Forward(1, dX);
            std::vector<double> dw = f.Reverse(2, std::vector<double>(1, 1));
            reverseResults.resize(n);
            for (size_t i = 0; i < n; i++)
            {
                reverseResults[i] = signFactor_ * dw[2 * (n + i) + 1];
            }
        }

    private:
        int signFactor_;
    };

    // Calculates third order mixed derivatives with respect to i-th variable
    // one time and (i+n)-th variable two times. Where domain vector is 2*n
    // dimensional. For correct work all other mixed derivatives should be zero.
    class ThirdMixedCalculator
        : public DerivCalculator
    {
    public:
        ThirdMixedCalculator(int signFactor = 1) : signFactor_(signFactor) {}

        void calcForward(cl::TapeFunction<double>& f, std::vector<double>& forwardResults)
        {
            size_t n = f.Domain() / 2;
            forwardResults.resize(n);
            // Direction for directional derivative.
            std::vector<double> dX(f.Domain(), 0);
            std::vector<double>::iterator px = dX.begin();
            std::vector<double>::iterator pxn = dX.begin() + n;
            std::generate(forwardResults.begin(), forwardResults.end(), [&]()
            {
                // Mixed derivative calculated from third order derivatives like this:
                // z = x + y
                // w = x - y
                // d3F/dz3 = d3F/dx3 + 3 * d3F/dx2dy + 3 * d3F/dxdy2 + d3F/dy3
                // d3F/dw3 = d3F/dx3 - 3 * d3F/dx2dy + 3 * d3F/dxdy2 - d3F/dy3
                // d3F/dx2dy = (d3F/dz3 - d3F/dw3 - 2 * d3F/dy3) / 6
                *px = 1;
                *pxn = 1;
                f.Forward(1, dX);
                *px = 0;
                *pxn = 0;
                f.Forward(2, dX);
                double dz3 = 6 * f.Forward(3, dX)[0];
                *px = -1;
                *pxn = 1;
                f.Forward(1, dX);
                *px = 0;
                *pxn = 0;
                f.Forward(2, dX);
                double dw3 = 6 * f.Forward(3, dX)[0];
                *px = 1;
                f.Forward(1, dX);
                *px = 0;
                f.Forward(2, dX);
                ++px;
                ++pxn;
                double dy3 = 6 * f.Forward(3, dX)[0];
                return signFactor_ * (dz3 - dw3 - 2 * dy3) / 6;
            });
        }

        void calcReverse(cl::TapeFunction<double>& f, std::vector<double>& reverseResults)
        {
            const size_t n = f.Domain() / 2;
            std::vector<double> dX(f.Domain(), 0);
            std::fill(dX.begin() + n, dX.end(), 1);
            f.Forward(1, dX);
            std::fill(dX.begin() + n, dX.end(), 0);
            f.Forward(2, dX);
            std::vector<double> dw = f.Reverse(3, std::vector<double>(1, 1));
            reverseResults.resize(n);
            for (size_t i = 0; i < n; i++)
            {
                reverseResults[i] = signFactor_ * 2 * dw[3 * i + 2];
            }
        }

    private:
        int signFactor_;
    };

    // Abstract class for easier adjoint testing.
    class AdjointTestHelper
    {
    protected:
        typedef QuantLib::Real Real;

        // Encapsulates work with log.
        class Log
        {
        public:
            Log(cl::AdjointTestOutput* p, const AdjointTestHelper* helper)
            : logger_(p)
            , helper_(helper)
            , timer_()
            , muted_(false)
            { }

            void setLogger(cl::AdjointTestOutput* p)
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
                *this << "Start of differentiation using " + helper_->analyticalName() + " method: " << currentTime() << std::endl;
                timer_.restart();
            }
            void afterAnalytical()
            {
                if (muted_)
                    return;
                *this << "End of differentiation using " + helper_->analyticalName() + " method." << std::endl;
                *this << "Time for " + helper_->analyticalName() + " method: " << timer_.elapsed() << std::endl;
                timer_.restart();
            }

            template<class T>
            Log& operator<<(T&& t)
            {
                if (logger_ && !muted_)
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

            cl::AdjointTestOutput* logger_;
            const AdjointTestHelper* helper_;
            boost::timer timer_;
            bool muted_;
        };

        AdjointTestHelper()
            : f_()
            , forwardResults_()
            , reverseResults_()
            , analyticalResults_()
            , log_(0, this)
            , calculator_(new FirstDerivCalculator(1))
        { }

        virtual ~AdjointTestHelper() {}

    public:
        // Returns number of independent variables.
        virtual size_t indepVarNumber() = 0;

        // Set derivative calculator.
        void setCalculator(std::shared_ptr<DerivCalculator> calculator)
        {
            calculator_ = calculator;
        }

        // Tests adjoint derivative calculation in forward and reverse mode.
        bool test()
        {
            recordTape();
            calcForward();
            calcReverse();
            calcAnalytical();
            return check();
        }

        // Tests adjoint derivative calculation in reverse mode only.
        bool testReverse()
        {
            recordTape();
            calcReverse();
            calcAnalytical();
            return checkReverse();
        }

        // Prints performance for all calculations.
        bool printPerformance(std::ostream& out = std::cout)
        {
            out << "\nNumber of independent variables: " << indepVarNumber() << std::endl;
            out << "Time for tape recording:\t" << testTapePerformance() << std::endl;
            out << "Time for Forward mode:  \t" << testForwardPerformance() << std::endl;
            out << "Time for Reverse mode:  \t" << testReversePerformance() << std::endl;
            out << "Time for " << analyticalName() << ":  \t" << testAnalyticalPerformance() << std::endl;
            out << std::endl;
            return check();
        }

        // Records performance for reverse mode and for analytical model (to calculate derivatives).
        // Return true if result check passed.
        bool recordPerformance(QuantLib::PerformanceTime& perfTime)
        {
            perfTime.indepVarNumber_ = indepVarNumber();
            perfTime.timeTapeRecording_ = testTapePerformance();
            perfTime.timeAdjoint_ = testReversePerformance();
            perfTime.timeAnalytical_ = testAnalyticalPerformance();
            return checkReverse();
        }

        // Returns performance of tape recording.
        double testTapePerformance(size_t testNo = 0)
        {
            if (testNo == 0)
            {
                // Setting default number of tests.
                // recordTape() has O(indepVarNumber) complexity, so divide by indepVarNumber.
                testNo = minPerfIteration() / indepVarNumber() + 1;
            }
            log_.beforeTape();
            log_.mute();
            boost::timer timer;
            for (size_t i = 0; i < testNo; i++)
            {
                recordTape();
            }
            double elapsed = timer.elapsed();
            log_.unmute();
            log_.afterTape();
            return elapsed / testNo;
        }

        // Returns performance of forward mode derivatives calculations.
        double testForwardPerformance(size_t testNo = 0)
        {
            if (!testNo)
            {
                // Setting default number of tests.
                // calcForward() has O(indepVarNumber^2) complexity, so divide by indepVarNumber^2.
                testNo = minPerfIteration() / indepVarNumber() / indepVarNumber() + 1;
            }
            if (!f_)
            {
                recordTape();
            }
            log_.beforeForward();
            log_.mute();
            boost::timer timer;
            for (size_t i = 0; i < testNo; i++)
            {
                calcForward();
            }
            double elapsed = timer.elapsed();
            log_.unmute();
            log_.afterForward();
            return elapsed / testNo;
        }

        // Returns performance of reverse mode derivatives calculations.
        double testReversePerformance(size_t testNo = 0)
        {
            if (!testNo)
            {
                // Setting default number of tests.
                // calcReverse() has O(indepVarNumber) complexity, so divide by indepVarNumber.
                testNo = minPerfIteration() / indepVarNumber() + 1;
            }
            if (!f_)
            {
                recordTape();
            }
            log_.beforeReverse();
            log_.mute();
            boost::timer timer;
            for (size_t i = 0; i < testNo; i++)
            {
                calcReverse();
            }
            double elapsed = timer.elapsed();
            log_.unmute();
            log_.afterReverse();
            return elapsed / testNo;
        }

        // Returns performance of Greek calculations with checking model.
        double testAnalyticalPerformance(size_t testNo = 0)
        {
            if (!testNo)
            {
                // Setting default number of tests.
                // calcAnalytical() has at least O(indepVarNumber) complexity, so divide by indepVarNumber.
                testNo = minPerfIteration() / indepVarNumber() + 1;
            }
            log_.beforeAnalytical();
            log_.mute();
            boost::timer timer;
            for (size_t i = 0; i < testNo; i++)
            {
                calcAnalytical();
            }
            double elapsed = timer.elapsed();
            log_.unmute();
            log_.afterAnalytical();
            return elapsed / testNo;
        }

        // Calculate Greeks with reverse mode only.
        void calcWithReverse()
        {
            recordTape();
            calcReverse();
        }

        // Returns Greeks calculated in reverse mode.
        std::vector<double> const& getReverse()
        {
            return reverseResults_;
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
        void setLogger(cl::AdjointTestOutput* p)
        {
            log_.setLogger(p);
        }

    protected:
        // Declare independent variables.
        virtual void independent() = 0;

        // Calculate dependent variables from independent.
        virtual void calcDependent() = 0;

        // Create new adjoint function.
        // AdjointTestHelper get ownership of it.
        virtual cl::TapeFunction<double>* newFunc() = 0;

        // Derived class could implement analytical calculation of derivatives.
        virtual std::vector<Real> analytical() { return std::vector<Real>(); }

        // Relative tolerance for checking consistensy with analytical results.
        virtual double relativeTol() const { return 1e-10; }

        // Absolute tolerance for checking consistensy with analytical results.
        virtual double absTol() const { return 1e-10; }

        // Name of analytical method.
        virtual std::string analyticalName() const { return "Finite diff."; }

        // Minimum iterations needed for performance mesuaring.
        // It is automaticly scaled using indepVarNumber() method.
        virtual size_t minPerfIteration() { return 1; }

        // Calculates derivatives in Forward mode.
        virtual void calcForward()
        {
            log_.beforeForward();
            calculator_->calcForward(*f_, forwardResults_);
            log_.afterForward();
        }

        // Calculates derivatives in Reverse mode.
        virtual void calcReverse()
        {
            log_.beforeReverse();
            calculator_->calcReverse(*f_, reverseResults_);
            log_.afterReverse();
        }

        // Records tape and sets new tape function.
        void recordTape()
        {
            log_.beforeTape();
            independent();
            calcDependent();
            f_.reset(newFunc());
            log_.afterTape();
        }

        // Calculates derivatives using analytical method.
        void calcAnalytical()
        {
            log_.beforeAnalytical();
            analyticalResults_ = analytical();
            log_.afterAnalytical();
        }

        // Checks forward mode, reverse mode and analytical result consistency.
        bool check()
        {
            bool result = true;
            int n = analyticalResults_.size();
            if (n != forwardResults_.size())
            {
                BOOST_ERROR("\nForward mode and " << analyticalName() << " vectors size mismatch.");
                return false;
            }
            if (n != reverseResults_.size())
            {
                BOOST_ERROR("\nReverse mode and " << analyticalName() << " vectors size mismatch.");
                return false;
            }
            std::vector<double>::const_iterator fwd = forwardResults_.begin();
            std::vector<double>::const_iterator rev = reverseResults_.begin();
            std::vector<Real  >::const_iterator anl = analyticalResults_.begin();
            log_ << "Derivatives" << std::endl;
            log_ << "Forward mode || Reverse mode || " << analyticalName() << std::endl;
            for (int i = 0; i < n; ++i, ++anl, ++fwd, ++rev)
            {
                log_ << std::setw(12) << *fwd << " ";
                log_ << std::setw(15) << *rev << " ";
                log_ << std::setw(15) << *anl << std::endl;

                // Maximum of absolut values of results.
                Real maxabs = std::max(std::abs(*anl), std::max(std::abs(*fwd), std::abs(*rev)));
                Real tol = std::max(relativeTol() * maxabs, absTol());
                if (std::abs(*fwd - *anl) > tol)
                {
                    result = false;
                    BOOST_ERROR("\nForward mode and " << analyticalName() << " derivative[" << i << "] mismatch."
                        << std::setw(24) << std::left << "\n    Forward mode:" << *fwd
                        << std::setw(24) << std::left << "\n    " + analyticalName() + ":" << *anl
                        << "\n    tolerance:  " << tol);
                }
                if (std::abs(*rev - *anl) > tol)
                {
                    result = false;
                    BOOST_ERROR("\nReverse mode and " << analyticalName() << " derivative[" << i << "] mismatch."
                        << std::setw(24) << std::left << "\n    Reverse mode:" << *rev
                        << std::setw(24) << std::left << "\n    " + analyticalName() + ": " << *anl
                        << "\n    tolerance:  " << tol);
                }
                if (std::abs(*rev - *fwd) > tol)
                {
                    result = false;
                    BOOST_ERROR("\nForward mode and Reverse mode derivative[" << i << "] mismatch."
                        << std::setw(24) << std::left << "\n    Forward mode:" << *fwd
                        << std::setw(24) << std::left << "\n    Reverse mode:" << *rev
                        << "\n    tolerance:  " << tol);
                }
            }
            return result;
        }

        // Checks reverse mode and analytical result consistency.
        bool checkReverse()
        {
            bool result = true;
            int n = analyticalResults_.size();
            if (n != reverseResults_.size())
            {
                BOOST_ERROR("\nReverse mode and " << analyticalName() << " vectors size mismatch.");
                return false;
            }
            std::vector<double>::const_iterator rev = reverseResults_.begin();
            std::vector<Real  >::const_iterator anl = analyticalResults_.begin();
            log_ << "Derivatives" << std::endl;
            log_ << "Reverse mode || " << analyticalName() << std::endl;
            for (int i = 0; i < n; ++i, ++anl, ++rev)
            {
                log_ << std::setw(12) << *rev << " ";
                log_ << std::setw(15) << *anl << std::endl;

                // Maximum of absolut values of results.
                Real maxabs = std::max(std::abs(*anl), std::abs(*rev));
                Real tol = std::max(relativeTol() * maxabs, absTol());
                if (std::abs(*rev - *anl) > tol)
                {
                    result = false;
                    BOOST_ERROR("\nReverse mode and " << analyticalName() << " derivative[" << i << "] mismatch."
                        << std::setw(24) << std::left << "\n    Reverse mode:" << *rev
                        << std::setw(24) << std::left << "\n    " + analyticalName() + ": " << *anl
                        << "\n    tolerance:  " << tol);
                }
            }
            return result;
        }

        std::unique_ptr<cl::TapeFunction<double> > f_;
        std::vector<double> forwardResults_;
        std::vector<double> reverseResults_;
        std::vector<Real> analyticalResults_;
        Log log_;
        std::shared_ptr<DerivCalculator> calculator_;
    };


    // Mesures performance for forward mode and analytical method.
    // outPerform, outSize - output for performanse time and memory size,
    // factory - functor that provides AdjointTestHelper classes for test,
    // n - number of performanse mesuaring,
    // step - step for problem size (used for factory(size_t) call).
    template <class HelperFactory>
    bool recordPerformance(cl::AdjointTestOutput& outPerform, cl::AdjointTestOutput& outSize, HelperFactory factory, size_t n, size_t step = 1)
    {
        std::vector<PerformanceTime> perfTimes(n);
        std::vector<TapeSize> tapeMemory(n);
        auto p = perfTimes.begin();
        auto m = tapeMemory.begin();
        bool result = true;
        for (size_t i = 1; i <= n; i++, ++p)
        {
#ifndef CL_GRAPH_GEN
            i = n;
#endif
            std::shared_ptr<AdjointTestHelper> helper = factory(step * i);
            helper->setLogger(&outPerform);
            result &= helper->recordPerformance(*p);
            *m++ = { helper->indepVarNumber(), helper->memory() };
        }
        outPerform << perfTimes;
        outSize << tapeMemory;
        return result;
    }
}

#endif