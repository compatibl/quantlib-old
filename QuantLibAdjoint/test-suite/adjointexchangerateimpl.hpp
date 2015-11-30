/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2004 StatPro Italia srl
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

//based on exchangerate.cpp from test-suite

#ifndef cl_adjoint_exchange_rate_impl_hpp
#define cl_adjoint_exchange_rate_impl_hpp
#pragma once

#include "adjointexchangeratetest.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include "utilities.hpp"
#include <ql/exchangerate.hpp>
#include <ql/currencies/europe.hpp>
#include <ql/currencies/america.hpp>
#include <ql/currencies/asia.hpp>
#include <ql/currencies/exchangeratemanager.hpp>
#include <boost/make_shared.hpp>
#include <ql/quantlib.hpp>

using namespace QuantLib;
using namespace std;
using namespace boost::unit_test_framework;

namespace
{
   struct ExchangeTest
        : cl::AdjointTest<ExchangeTest>
    {
        explicit ExchangeTest(Size size)
            : AdjointTest()
            , size_(size)
            , rates_(size)
            , derivedAmount_(1)
            , EUR_(EURCurrency())
            , USD_(USDCurrency())
            , GBP_(GBPCurrency())
            , m1_(50000.0 * GBP_)
        {
            Money::conversionType = Money::NoConversion;
            rates_ = { 1.2042, 0.6612 };
        }

        Size indepVarNumber() { return size_; }

        Size depVarNumber() { return 1; }

        Size minPerfIteration() { return 0; }

        void recordTape()
        {
            cl::Independent(rates_);
            calculateAmount();
            f_ = std::make_unique<cl::tape_function<double>>(rates_, derivedAmount_);
        }

        void calculateAmount()
        {
            ExchangeRate derived = ExchangeRate::chain(ExchangeRate(EUR_, USD_, rates_[0]), ExchangeRate(EUR_, GBP_, rates_[1]));
            derivedAmount_[0] = derived.exchange(m1_).value();
        }

        // Calculates derivatives using analytical formula.
        void calcAnalytical()
        {
            analyticalResults_.resize(size_, 0.0);
            analyticalResults_[0] = m1_.value() / 0.6612;
            analyticalResults_[1] = -m1_.value()*rates_[0] / std::pow(rates_[1], 2);
        }

        double relativeTol() const { return 1e-3; }

        double absTol() const { return 1e-5; }

        Size size_;
        std::vector<cl::tape_double> rates_;
        std::vector<cl::tape_double> derivedAmount_;
        Currency EUR_, USD_, GBP_;
        Money m1_;
    };
}

#endif
