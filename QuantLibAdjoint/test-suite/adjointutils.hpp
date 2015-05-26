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

#ifndef quantlib_utilities_adjoint_hpp
#define quantlib_utilities_adjoint_hpp

#include <ql/math/matrix.hpp>
#include <boost/timer.hpp>
#include <boost/date_time.hpp>

namespace QuantLib
{
    struct PerformanceTime
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "Number of independent variables", "Tape recording time", "Adjoint calculations time", "No Adjoint calculations time"
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, PerformanceTime& v)
        {
                stm << v.indepVarNumber_
                    << ";" << v.timeTapeRecording_
                    << ";" << v.timeAdjoint_
                    << ";" << v.timeAnalytical_ << std::endl;

                return stm;
            }

        double timeTapeRecording_;
        double timeAdjoint_;
        double timeAnalytical_;
        Size indepVarNumber_;
    };

    struct AdjointTime
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "Number of independent variables", ""
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, AdjointTime& v)
        {
                stm << v.indepVarNumber_
                    << ";" << v.timeAdjoint_
                    << std::endl;
                return stm;
            }

        double timeAdjoint_;
        Size indepVarNumber_;
    };

    struct TapeSize
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "Number of independent variables", ""
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, TapeSize& v)
        {
                stm << v.indepVarNumber_
                    << ";" << v.memory_ / 1048576.0
                    << std::endl;
                return stm;
            }

        Size indepVarNumber_;
        Size memory_;
    };

    inline std::string currentTime()
    {
        boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();
        std::stringstream ss;
        ss << now.time_of_day().hours() << "h  " << now.time_of_day().minutes() << "m  "
            << now.time_of_day().seconds() << "s  " << now.time_of_day().fractional_seconds() << "mcs";
        return ss.str();
    }

    // checking consistency of elements of three vectors
    inline bool checkConsistencyWithBC(
        const std::vector<double>& forwardMode,
        const std::vector<double>& reverseMode,
        const std::vector<Real>& theoretical,
        double relativeTol, double absTol = 1e-10)
    {
        bool result = true;
        int n = theoretical.size();
        if (theoretical.size() != forwardMode.size())
        {
            std::cout << "Forward mode and Theoretical results vectors size mismatch.\n";
            return false;
        }
        if (theoretical.size() != reverseMode.size())
        {
            std::cout << "Reverse mode and Theoretical results vectors size mismatch.\n";
            return false;
        }
        std::vector<double>::const_iterator fwd = forwardMode.begin();
        std::vector<double>::const_iterator rev = reverseMode.begin();
        std::vector<Real  >::const_iterator thr = theoretical.begin();
        for (int i = 0; i < n; ++i, ++thr, ++fwd, ++rev)
        {
            Real maxabs = std::max(std::abs(*thr), std::max(std::abs(*fwd), std::abs(*rev)));
            Real tol = relativeTol * maxabs + absTol;
            if (std::abs(*fwd - *thr) > tol)
            {
                result = false;
                std::cout << "Forward mode derivative[" << i << "] and Theoretical results mismatch."
                    << "\n    forward mode:      " << *fwd
                    << "\n    theoretical:       " << *thr
                    << "\n    tolerance:  " << tol << std::endl;
            }
            if (std::abs(*rev - *thr) > tol)
            {
                result = false;
                std::cout << "Reverse mode derivative[" << i << "] and Theoretical results mismatch."
                    << "\n    reverse mode:      " << *rev
                    << "\n    theoretical:       " << *thr
                    << "\n    tolerance:  " << tol << std::endl;
            }
        }
        return result;
    }


    // checking consistency of elements of three vectors
    inline bool checkConsistencyWithFD(
        const std::vector<double>& forwardMode,
        const std::vector<double>& reverseMode,
        const std::vector<Real>& finiteDiff,
        double relativeTol, double absTol = 1e-10)
    {
        bool result = true;
        int n = finiteDiff.size();
        if (finiteDiff.size() != forwardMode.size())
        {
            std::cout << "Forward mode and Finite difference vectors size mismatch.\n";
            return false;
        }
        if (finiteDiff.size() != reverseMode.size())
        {
            std::cout << "Reverse mode and Finite difference vectors size mismatch.\n";
            return false;
        }
        std::vector<double>::const_iterator fwd = forwardMode.begin();
        std::vector<double>::const_iterator rev = reverseMode.begin();
        std::vector<Real  >::const_iterator fin = finiteDiff.begin();
        for (int i = 0; i < n; ++i, ++fin, ++fwd, ++rev)
        {
            Real maxabs = std::max(std::abs(*fin), std::max(std::abs(*fwd), std::abs(*rev)));
            Real tol = relativeTol * maxabs + absTol;
            if (std::abs(*fwd - *fin) > tol)
            {
                result = false;
                std::cout << "Forward mode and Finite difference derivative[" << i << "] mismatch."
                    << "\n    forward mode:      " << *fwd
                    << "\n    finite difference: " << *fin
                    << "\n    tolerance:  " << tol << std::endl;
            }
            if (std::abs(*rev - *fin) > tol)
            {
                result = false;
                std::cout << "Reverse mode and Finite difference derivative[" << i << "] mismatch."
                    << "\n    reverse mode:      " << *rev
                    << "\n    finite difference: " << *fin
                    << "\n    tolerance:  " << tol << std::endl;
            }
            if (std::abs(*rev - *fwd) > tol)
            {
                result = false;
                std::cout << "Forward mode and Reverse mode derivative[" << i << "] mismatch."
                    << "\n    forward mode:      " << *fwd
                    << "\n    reverse mode:      " << *rev
                    << "\n    tolerance:  " << tol << std::endl;
            }
        }
        return result;
    }


    // computes gradient in Forward AD mode and outputs results
    inline double gradForward(cl::TapeFunction<double>& f,
        std::vector<double>& grad,
        bool timeOutput,
        bool gradOutput)
    {
        using namespace std;
        size_t n = f.Domain();
        size_t m = f.Range();
        grad.resize(m*n);
        vector<double> dx(n), dy(m);


        const size_t iterNum = 10;
        boost::timer timer;
        timer.restart();
        for (size_t k = 0; k < iterNum; k++)
        {
            for (size_t i = 0; i < n; i++)
            {
                dx[i] = 1;
                dy = f.Forward(1, dx);
                dx[i] = 0;
                for (size_t j = 0; j < m; j++)
                    grad[j*n + i] = dy[j];
            }
        }
        double timeForward = timer.elapsed() / iterNum;
        if (timeOutput)
            cout << "Time for Forward AD mode : " << timeForward << endl;


        if (gradOutput)
        {
            cout << "Forward AD mode : " << endl;
            for (size_t i = 0; i < grad.size(); i++)
                cout << "dy" << i / n << "/dx" << i % n << "= " << grad[i] << endl;
        }
        return timeForward;
    }


    // computes gradient in Reverse AD mode and outputs results
    inline double gradReverse(cl::TapeFunction<double>& f,
        std::vector<double>& grad,
        bool timeOutput,
        bool gradOutput)
    {
        using namespace std;
        size_t n = f.Domain();
        size_t m = f.Range();
        grad.resize(m*n);
        vector<double> dw(m), dy(n);


        const size_t iterNum = 10;
        boost::timer timer;
        timer.restart();
        for (size_t k = 0; k < iterNum; k++)
        for (size_t i = 0; i < m; i++)
        {
            dw[i] = 1;
            dy = f.Reverse(1, dw);
            dw[i] = 0;
            for (size_t j = 0; j < n; j++)
                grad[i*n + j] = dy[j];
        }
        double timeReverse = timer.elapsed() / iterNum;
        if (timeOutput)
            cout << "Time for Reverse AD mode : " << timer.elapsed() / iterNum << endl;

        if (gradOutput)
        {
            cout << "Reverse AD mode : " << endl;
            for (size_t i = 0; i < grad.size(); i++)
                cout << "dy" << i / n << "/dx" << i % n << "= " << grad[i] << endl;
        }
        return timeReverse;
    }

    // computes gradient in Forward AD mode and outputs results
    inline double gradForward(cl::TapeFunction<double>& f,
        std::vector<double>& grad,
        cl::AdjointTestOutput& out,
        bool timeOutput,
        bool gradOutput,
        size_t iterNum = 10)
    {
#ifndef CL_GRAPH_GEN
        iterNum = 1;
#endif
        using namespace std;
        out.log() << "Start of differentiation in Forward mode:  " << currentTime() << endl;
        size_t n = f.Domain();
        size_t m = f.Range();
        grad.resize(m*n);
        vector<double> dx(n), dy(m);

        try{
            boost::timer timer;
            timer.restart();
            for (size_t k = 0; k < iterNum; k++)
            {
                for (size_t i = 0; i < n; i++)
                {
                    dx[i] = 1;
                    dy = f.Forward(1, dx);
                    dx[i] = 0;
                    for (size_t j = 0; j < m; j++)
                        grad[j*n + i] = dy[j];
                }
            }
            double timeForward = timer.elapsed() / iterNum;
            if (timeOutput)
                cout << "Time for Forward AD mode : " << timeForward << endl;


            if (gradOutput)
            {
                cout << "Forward AD mode : " << endl;
                for (size_t i = 0; i < grad.size(); i++)
                    cout << "dy" << i / n << "/dx" << i % n << "= " << grad[i] << endl;
            }
            out.log() << "Forward Mode Derivatives calculated successfully" << std::endl;
            out.log() << "Time in Forward mode:   " << timeForward << "s"<< endl;

            return timeForward;
        }
        catch (std::exception& ex)
        {
            out.log() << ex.what() << std::endl;
            return 0;
        }
    }


    // computes gradient in Reverse AD mode and outputs results
    inline double gradReverse(cl::TapeFunction<double>& f,
        std::vector<double>& grad,
        cl::AdjointTestOutput& out,
        bool timeOutput,
        bool gradOutput,
        size_t iterNum = 10)
    {
#ifndef CL_GRAPH_GEN
        iterNum = 1;
#endif
        using namespace std;
        out.log() << "Start of differentiation in Reverse mode:   " << currentTime() << endl;
        size_t n = f.Domain();
        size_t m = f.Range();
        grad.resize(m*n);
        vector<double> dw(m), dy(n);

        try{
            boost::timer timer;
            timer.restart();
            for (size_t k = 0; k < iterNum; k++)
            for (size_t i = 0; i < m; i++)
            {
                dw[i] = 1;
                dy = f.Reverse(1, dw);
                dw[i] = 0;
                for (size_t j = 0; j < n; j++)
                    grad[i*n + j] = dy[j];
            }
            double timeReverse = timer.elapsed() / iterNum;
            if (timeOutput)
                cout << "Time for Reverse AD mode : " << timer.elapsed() / iterNum << endl;

            if (gradOutput)
            {
                cout << "Reverse AD mode : " << endl;
                for (size_t i = 0; i < grad.size(); i++)
                    cout << "dy" << i / n << "/dx" << i % n << "= " << grad[i] << endl;
            }

            out.log() << "Reverse Mode Derivatives calculated successfully" << std::endl;
            out.log() << "Time in Reverse mode:   " << timeReverse <<" s"<< endl;

            return timeReverse;
        }
        catch (std::exception& ex)
        {
            out.log() << ex.what() << std::endl;
            return 0;
        }
    }

    // vector-matrix converters

    template<class C>
    inline Matrix vector2matrix(const C& vec, int rows, int columns)
    {
        Matrix A(rows, columns);
        copy(vec.begin(), vec.end(), A.begin());
        return A;
    }

    inline std::vector<cl::TapeDouble>  matrix2vector(const Matrix& A)
    {
        return std::vector<cl::TapeDouble>(A.begin(), A.end());
    }


    inline PerformanceTime& operator+=(PerformanceTime& left, PerformanceTime const& right)
    {
        left.timeAdjoint_ += right.timeAdjoint_;
        left.timeAnalytical_ += right.timeAnalytical_;
        left.timeTapeRecording_ += right.timeTapeRecording_;
        return left;
    }
    inline PerformanceTime& operator-=(PerformanceTime& left, PerformanceTime const& right)
    {
        left.timeAdjoint_ -= right.timeAdjoint_;
        left.timeAnalytical_ -= right.timeAnalytical_;
        left.timeTapeRecording_ -= right.timeTapeRecording_;
        return left;
    }
    inline PerformanceTime& operator*=(PerformanceTime& left, double right)
    {
        left.timeAdjoint_ *= right;
        left.timeAnalytical_ *= right;
        left.timeTapeRecording_ *= right;
        return left;
    }
    inline PerformanceTime operator+(PerformanceTime const& left, PerformanceTime const& right)
    {
        PerformanceTime temp = left;
        return temp += right;
    }
    inline PerformanceTime operator-(PerformanceTime const& left, PerformanceTime const& right)
    {
        PerformanceTime temp = left;
        return temp -= right;
    }
    inline PerformanceTime operator*(double left, PerformanceTime const& right)
    {
        PerformanceTime temp = right;
        return temp *= left;
    }
    // what = max(what, where);
    inline void upTo(PerformanceTime & what, PerformanceTime const& where)
    {
        what.timeAdjoint_ = std::max(what.timeAdjoint_, where.timeAdjoint_);
        what.timeAnalytical_ = std::max(what.timeAnalytical_, where.timeAnalytical_);
        what.timeTapeRecording_ = std::max(what.timeTapeRecording_, where.timeTapeRecording_);
    }
    // what = max(what, where);
    inline void upTo(PerformanceTime & what, double where)
    {
        upTo(what, PerformanceTime{ where, where, where, 0 });
    }

    // Dummy realization of heat equation.
    template <class RanIt>
    void heatbackward(RanIt first, RanIt last, double a, size_t n = 1)
    {
        if (first - last < 3)
            return;

        for (size_t i = 0; i < n; i++)
        {
            std::reverse_iterator<RanIt> prev(last);
            auto curr = prev + 1;
            auto next = curr + 1;
            for (; next != std::reverse_iterator<RanIt>(first); ++curr, ++next, ++prev)
            {
                typename RanIt::value_type deriv = *next + *prev - 2 * (*curr);
                *curr += 0.33 * a * deriv;
                *prev -= 0.165 * a * deriv;
                *next -= 0.165 * a * deriv;
                //upTo(*next, 0);
            }
            //upTo(*curr, 0.5 * (*prev));
        }
    }

    template <class RanIt>
    inline void makeSmooth(RanIt first, RanIt last, double sigma = 10.0)
    {
        size_t n = size_t(sigma*sigma * 10) + 1;
        double a = sigma*sigma / n;
        heatbackward(first, last, a, n);
    }

    // Makes data in c smooth.
    template <class Cont>
    inline void makeSmooth(Cont& c, double sigma = 10.0)
    {
        makeSmooth(c.begin(), c.end(), sigma);
    }
}

#endif