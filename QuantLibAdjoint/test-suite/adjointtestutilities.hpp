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

#ifndef cl_adjoint_utilities_hpp
#define cl_adjoint_utilities_hpp
#pragma once

#include "utilities.hpp"
#include "adjointutils.hpp"

namespace QuantLib {

    // checking consistency of elements of three vectors
    inline bool checkWithFiniteDiff(
        const std::vector<double>& forwardMode,
        const std::vector<double>& reverseMode,
        const std::vector<Real  >& finiteDiff,
        double relativeTol,
        double absTol = 1e-10)
    {
        bool result = true;
        int n = finiteDiff.size();
        if (finiteDiff.size() != forwardMode.size()) {
            BOOST_ERROR("\nForward mode and Finite difference vectors size mismatch.");
            return false;
        }
        if (finiteDiff.size() != reverseMode.size()) {
            BOOST_ERROR("\nReverse mode and Finite difference vectors size mismatch.");
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
                BOOST_ERROR("\nForward mode and Finite difference derivative[" << i << "] mismatch."
                    << "\n    forward mode:      " << *fwd
                    << "\n    finite difference: " << *fin
                    << "\n    tolerance:  " << tol);
            }
            if (std::abs(*rev - *fin) > tol)
            {
                result = false;
                BOOST_ERROR("\nReverse mode and Finite difference derivative[" << i << "] mismatch."
                    << "\n    reverse mode:      " << *rev
                    << "\n    finite difference: " << *fin
                    << "\n    tolerance:  " << tol);
            }
            if (std::abs(*rev - *fwd) > tol)
            {
                result = false;
                BOOST_ERROR("\nForward mode and Reverse mode derivative[" << i << "] mismatch."
                    << "\n    forward mode:      " << *fwd
                    << "\n    reverse mode:      " << *rev
                    << "\n    tolerance:  " << tol);
            }
        }
        return result;
    }

    // checking consistency of elements of three vectors
    inline bool checkWithFiniteDiff(
        const std::vector<double>& forwardMode,
        const std::vector<double>& reverseMode,
        const std::vector<Real  >& finiteDiff,
        cl::tape_empty_test_output& out,
        double relativeTol,
        double absTol = 1e-10)
    {
        bool result = true;
        int n = finiteDiff.size();
        if (finiteDiff.size() != forwardMode.size()) {
            BOOST_ERROR("\nForward mode and Finite difference vectors size mismatch.");
            return false;
        }
        if (finiteDiff.size() != reverseMode.size()) {
            BOOST_ERROR("\nReverse mode and Finite difference vectors size mismatch.");
            return false;
        }
        std::vector<double>::const_iterator fwd = forwardMode.begin();
        std::vector<double>::const_iterator rev = reverseMode.begin();
        std::vector<Real  >::const_iterator fin = finiteDiff.begin();

        out.log() << "Forward, Reverse and Finite diference derivatives:" << std::endl;

        for (int i = 0; i < n; ++i, ++fin, ++fwd, ++rev)
        {
            out.log() << "[ "<< i << "th ] derivatives:"<< * fwd << "  \t" << *rev << "  \t" << *fin << std::endl;
            Real maxabs = std::max(std::abs(*fin), std::max(std::abs(*fwd), std::abs(*rev)));
            Real tol = relativeTol * maxabs + absTol;
            if (std::abs(*fwd - *fin) > tol)
            {
                result = false;
                BOOST_ERROR("\nForward mode and Finite difference derivative[" << i << "] mismatch."
                    << "\n    forward mode:      " << *fwd
                    << "\n    finite difference: " << *fin
                    << "\n    tolerance:  " << tol);
                out.log() << "\nForward mode and Finite difference derivative[" << i << "] mismatch:" << "forward mode" << *fwd << "finite difference" << *fin << std::endl;

            }
            if (std::abs(*rev - *fin) > tol)
            {
                result = false;
                BOOST_ERROR("\nReverse mode and Finite difference derivative[" << i << "] mismatch."
                    << "\n    reverse mode:      " << *rev
                    << "\n    finite difference: " << *fin
                    << "\n    tolerance:  " << tol);
                out.log() << "\n Reverse mode and Finite difference derivative[" << i << "] mismatch:" << "reverse mode" << *rev << "finite difference" << *fin << std::endl;
            }
            if (std::abs(*rev - *fwd) > tol)
            {
                result = false;
                BOOST_ERROR("\nForward mode and Reverse mode derivative[" << i << "] mismatch."
                    << "\n    forward mode:      " << *fwd
                    << "\n    reverse mode:      " << *rev
                    << "\n    tolerance:  " << tol);
                out.log() << "\nForward mode and Reverse mode derivative[" << i << "] mismatch:" << "forward mode" << *fwd << "reverse mode" << *rev << std::endl;
            }
        }
        return result;
    }


    // checking consistency of vectors elements
    template <class T>
    inline bool checkWithFiniteDiff(
        const std::vector<T>& adjointDeriv,
        const std::vector<Real>& finiteDiffDeriv,
        Real relativeTol,
        Real absoluteTol)
    {
        bool result = true;
        if (adjointDeriv.size() != finiteDiffDeriv.size())
        {
            result = false;
            BOOST_ERROR("\nAn adjoint derivatives vector and a finite "
                "difference derivatives vector have different sizes.");
        }
        for (size_t i = 0, n = adjointDeriv.size(); i < n; ++i)
        {
            Real err = std::abs(adjointDeriv[i] - finiteDiffDeriv[i]);
            Real toler = std::max(absoluteTol,
                relativeTol * std::max(std::abs(adjointDeriv[i]), std::abs(finiteDiffDeriv[i])));
            if (err > toler)
            {
                result = false;
                BOOST_ERROR("\nAdjoint derivative and finite difference derivative at position " << i << " mismatch."
                    << "\n  adjoint: " << adjointDeriv[i]
                    << "\n  finite difference: " << finiteDiffDeriv[i]
                    << "\n  tolerance: " << toler);
            }
        }
        return result;
    }


    // checking consistency of vectors elements
    template <class T>
    inline bool checkWithFiniteDiff(
        const std::vector<T>& adjointDeriv,
        const std::vector<Real>& finiteDiffDeriv,
        cl::tape_empty_test_output& out,
        Real relativeTol,
        Real absoluteTol)
    {
        bool result = true;
        if (adjointDeriv.size() != finiteDiffDeriv.size())
        {
            result = false;
            BOOST_ERROR("\nAn adjoint derivatives vector and a finite "
                "difference derivatives vector have different sizes.");
            out.log() << "\nAn adjoint derivatives vector and a finite "
                << "difference derivatives vector have different sizes." << std::endl;
        }
        for (size_t i = 0, n = adjointDeriv.size(); i < n; ++i)
        {
            Real err = std::abs(adjointDeriv[i] - finiteDiffDeriv[i]);
            Real toler = std::max(absoluteTol,
                relativeTol * std::max(std::abs(adjointDeriv[i]), std::abs(finiteDiffDeriv[i])));
            if (err > toler)
            {
                result = false;
                BOOST_ERROR("\nAdjoint derivative and finite difference derivative at position " << i << " mismatch."
                    << "\n  adjoint: " << adjointDeriv[i]
                    << "\n  finite difference: " << finiteDiffDeriv[i]
                    << "\n  tolerance: " << toler);
                out.log() << "\nAdjoint derivative and finite difference derivative at position " << i << " mismatch."
                    << "\n  adjoint: " << adjointDeriv[i]
                    << "\n  finite difference: " << finiteDiffDeriv[i]
                    << "\n  tolerance: " << toler << std::endl;
            }
        }
        return result;
    }


    // adapter for Jacobian
    inline std::vector<Real> jacobian(cl::tape_function<double>& f, std::vector<Real> const& x)
    {
        std::vector<double> xd, yd;
        std::for_each(x.cbegin(), x.cend(),
            [&](Real const& a){ xd.push_back(double(a)); }
        );
        yd = f.Jacobian(xd);
        std::vector<Real> y;
        std::for_each(yd.cbegin(), yd.cend(),
            [&](double const& a){ y.push_back(Real(a)); }
        );
        return y;
    }


    // adapter with a timer for Jacobian
    inline std::vector<Real> jacobian(cl::tape_function<double>& f, std::vector<Real> const& x, size_t iterNum, double& time)
    {
        std::vector<double> xd, yd;
        std::for_each(x.cbegin(), x.cend(),
            [&](Real const& a){ xd.push_back(double(a)); }
        );

        boost::timer timer;
        for (size_t i = 0; i < iterNum; ++i)
            yd = f.Jacobian(xd);
        time = timer.elapsed() / iterNum;

        std::vector<Real> y;
        std::for_each(yd.cbegin(), yd.cend(),
            [&](double const& a){ y.push_back(Real(a)); }
        );
        return y;
    }

}

#endif