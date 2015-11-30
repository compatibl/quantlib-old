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

#ifndef cl_adjoint_bond_portfolio_impl_hpp
#define cl_adjoint_bond_portfolio_impl_hpp
#pragma once

#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include <boost/timer.hpp>
#include <iostream>
#include <ql/quantlib.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;

#define OUTPUT_FOLDER_NAME "AdjointBondPortfolio"

namespace
{
    enum
    {
#if defined CL_GRAPH_GEN
        // Number of points for dependency plots.
        pointNo = 100,
        // Number of points for performance plot.
        iterNo = 1000,
        // Step for portfolio size for performance testing .
        step = 10,
#else
        // Number of points for dependency plots.
        pointNo = 10,
        // Number of points for performance plot.
        iterNo = 30,
        // Step for portfolio size for performance testing .
        step = 10,
#endif
        // Defines performance accuracy. Its value is a minimum number
        // of calling of O(1) complexity methods per one performance test.
        iterNumFactor = 0,
    };

    struct RateVariation
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "Rate", ""
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, RateVariation& v)
        {
                stm << v.inputRate_
                    << ";" << v.bondPortfolioPrice_ << std::endl;

                return stm;
            }

        Real inputRate_;
        Real bondPortfolioPrice_;
    };

    struct BondPortfolioData
    {
        BondPortfolioData() :
          calendar_(TARGET())
        , today_(calendar_.adjust(Date::todaysDate()))
        , frequency_(Semiannual)
        , bondDayCount_(Thirty360())
        , compounding_(Compounded)
        // Value of bond at maturity
        , faceAmount_(1000000.0)
        , issueMonth_(6)
        , length_(3)
        , settlementDays_(3)
        , coupon_(0.02)
        , accrualConvention_(Unadjusted)
        , paymentConvention_(ModifiedFollowing)
        , redemption_(100.0)
        {
        }

        // Method to create a bond object of a Quantlib FixedRateBond class.
        // It returns an instance of FixedRateBond class.
        FixedRateBond bondCreate()
        {
            Date issue = calendar_.advance(today_,
                                           issueMonth_, Months);

            //Maturity time for a bond.
            Date maturity = calendar_.advance(issue, length_, Years);
            Schedule sch(issue, maturity,
                         Period(frequency_), calendar_,
                         accrualConvention_, accrualConvention_,
                         DateGeneration::Backward, false);
            return FixedRateBond(settlementDays_, faceAmount_, sch,
                                  std::vector<Rate>(1, coupon_),
                                  bondDayCount_, paymentConvention_,
                                  redemption_, issue);
        }

        Real calculatePrice(Real rate)
        {
            return BondFunctions::cleanPrice(bondCreate(), rate,
                                                       bondDayCount_,
                                                       compounding_,
                                                       frequency_);
        }

        // Finite difference calculation of derivatives.
        //  - indep - vector of independent variables.
        //  - sfFinite - vector to store derivative values,
        // calculated by finite difference scheme.
        //  - h - step of approximation.
        void calculateCentralFinDiff(std::vector<cl::tape_double>& indep, std::vector<Real>& sfFinite, double h)
        {
            Size size = indep.size();
            sfFinite.resize(size);

            // Evaluate derivatives using step-forward and step-backward values using the formula:
            // Derivative = (stepforward - stepbackward) / (2 * h), h - step size.
            for (Size i = 0; i < size; i++)
            {
                // Calculate approximate derivative:
                sfFinite[i] = (calculatePrice(indep[i] + h) - calculatePrice(indep[i] - h)) / (2 * h);
            }
        }

        Calendar calendar_;
        Date today_;
        Real faceAmount_;
        Frequency frequency_;
        DayCounter bondDayCount_;
        Compounding compounding_;
        Integer issueMonth_;
        Integer length_;
        Natural settlementDays_;
        Real coupon_;
        BusinessDayConvention accrualConvention_;
        BusinessDayConvention paymentConvention_;
        Real redemption_;
    };


    struct TestData
        : public BondPortfolioData
    {
        struct Test
        : public cl::AdjointTest<Test>
        {
            Test(Size size, TestData* data)
            : size_(size)
            , data_(data)
            , rate_(size)
            , portfolioPrice_()
            {
                setLogger(&data_->outPerform_);

                if (size_ != 1)
                {
                    for (Size i = 0; i < size_; i++)
                    {
                        rate_[i] = 0.03 + 0.001*i;
                    }
                }
                else
                {
                    rate_[0] = 0.030;
                }
            }

            Size indepVarNumber() { return size_; }

            Size depVarNumber() { return 1; }

            Size minPerfIteration() { return iterNumFactor; }

            void recordTape()
            {
                cl::Independent(rate_);
                calculatePortfolioPrice();
                f_ = std::make_unique<cl::tape_function<double>>(rate_, portfolioPrice_);
            }

            // Calculates price of portfolio and each option.
            void calculatePortfolioPrice()
            {
                portfolioPrice_.resize(1, 0);
                for (Size i = 0; i < size_; i++)
                    portfolioPrice_[0]+=data_->calculatePrice(rate_[i]);
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                double h = 1.0e-3;  // shift for finite diff. method
                analyticalResults_.resize(size_);
                data_->calculateCentralFinDiff(rate_, analyticalResults_, h);
            }

            double relativeTol() const { return 1e-2; }

            double absTol() const { return 1e-5; }

            Size size_;
            TestData* data_;
            std::vector<cl::tape_double> rate_;
            std::vector<cl::tape_double> portfolioPrice_;
        };

        TestData()
            : BondPortfolioData()
            , outPerform_(OUTPUT_FOLDER_NAME "//",
            { { "filename", "AdjointPerformance" }
            , { "not_clear", "Not" }
            , { "title", "Bond portfolio value differentiation performance with respect to yield rates " }
            , { "xlabel", "Number of bonds in a portfolio" }
            , { "ylabel", "Time (s)" }
            , { "line_box_width", "-5" }
            , { "smooth", "50" }
            , { "cleanlog", "true" }
            })
            , outAdjoint_(OUTPUT_FOLDER_NAME "//",
            { { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "title", "Bond portfolio value adjoint differentiation performance with respect to yield rates" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of yield rates" }
            , { "smooth", "50" }
            , { "cleanlog", "false" }
            })
            , outSize_(OUTPUT_FOLDER_NAME "//",
            { { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "title", "Tape size dependence on number of yield rates" }
            , { "ylabel", "Memory (MB)" }
            , { "smooth", "50" }
            , { "cleanlog", "false" }
            })
            , out_(OUTPUT_FOLDER_NAME "//output",
            { { "filename", "BondPriceOnRate" }
            , { "not_clear", "Not" }
            , { "title", "Bond price dependence on yield rate" }
            , { "ylabel", "Bond value" }
            , { "smooth", "50" }
            , { "cleanlog", "false" }
            })
            {
            }

        bool makeOutput()
        {
            bool ok = true;
            if (pointNo > 0)
            {
                ok &= recordDependencePlot();
            }
            ok &= cl::recordPerformance(*this, iterNo, step);
            return ok;
        }

        std::shared_ptr<Test> getTest(size_t size)
        {
            return std::make_shared<Test>(size, this);
        }

        // Makes plots for strike sensitivity dependence.
        bool recordDependencePlot()
        {
            std::vector<RateVariation> outData(pointNo);
            auto test = getTest(pointNo);

            for (Size i = 0; i < pointNo; i++)
            {
                outData[i] = { test->rate_[i], test->data_->calculatePrice(test->rate_[i]) };
            }
            out_ << outData;
            return true;
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output out_;
    };


    typedef TestData::Test BondPortfolioTest;
}
#endif