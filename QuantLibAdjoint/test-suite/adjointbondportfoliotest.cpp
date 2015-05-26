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

#define CL_OPEN_GENERAL_OUPUT

#include "adjointbondportfoliotest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include <boost/timer.hpp>
#include <iostream>
#include <ql/quantlib.hpp>

// Based on bonds.cpp file in testsuite
// Bond portfolio consists of a sum of bonds with ranging interest rates: r_1,...,r_n and fixed face amount.
// The portfolio value is BondPortfolioValue= \sum_{i=1}^n e^{-r_i*t_i}.
// In this test we calculate derivatives of BondPortfolioValue with respect to r_i, i=1,...,n; n-- number of  bonds in the portfolio.

using namespace QuantLib;
using namespace std;
using namespace boost::unit_test_framework;

struct CommonVars
{
    // Setup initial parameters for a bond and finite-differences calculation of derivatives.
    CommonVars(Size size): calendar_(TARGET())
        , today_(calendar_.adjust(Date::todaysDate()))
        , frequency_(Semiannual)
        , bondDayCount_(Thirty360())
        , compounding_(Compounded)
        // Value of bond at maturity
        , faceAmount_(1000000.0)
        , independent_(size)
        , portfolioPrice_(1)
        , h_(1.0e-8)
        , relTolerance_(1e-2)
        , absTolerance_(1e-4)
    {
        Settings::instance().evaluationDate() = today_;
        prepare();
    }

    // Initialize  values for interest rates in a bond portfolio
    void prepare()
    {
        std::vector<cl::TapeDouble>::iterator it = independent_.begin();
        //Init yield curve values

        for (Size i = 0; it != independent_.end(); it++, i++)
            *it = 0.03 + 0.001*i;
    }

    // Method to create a bond object of a Quantlib FixedRateBond class.
    // It returns an instance of  FixedRateBond class.

    FixedRateBond  bondCreate()
    {
        Integer issueMonth = 6;
        Integer length = 3;
        Natural settlementDays = 3;
        Real coupon = .02;
        BusinessDayConvention accrualConvention = Unadjusted;
        BusinessDayConvention paymentConvention = ModifiedFollowing;
        Real redemption = 100.0;
        Date dated = calendar_.advance(today_,
            issueMonth, Months);
        Date issue = dated;

        //Maturity time for a bond.
        Date maturity = calendar_.advance(issue,
            length, Years);
        Schedule sch(dated, maturity,
            Period(frequency_), calendar_,
            accrualConvention, accrualConvention,
            DateGeneration::Backward, false);
        return  FixedRateBond(settlementDays, faceAmount_, sch,
            std::vector<Rate>(1, coupon),
            bondDayCount_, paymentConvention,
            redemption, issue);
    }

    // Common data.
    Calendar calendar_;
    Date today_;
    Real faceAmount_;
    Frequency frequency_;
    DayCounter bondDayCount_;
    Compounding compounding_;

    // Adjoint routine vectors
    std::vector<cl::TapeDouble> independent_;
    std::vector<cl::TapeDouble>  portfolioPrice_;

    // Finite-differences scheme step parameter.
    double h_;

    // Parameters to check derivative values by comparison with finite-differences approach.
    double relTolerance_;
    double absTolerance_;

    // Performance results to be stored in vectors
    std::vector<PerformanceTime> performanceTime_;

    // Performance results for adjoint mode only
    std::vector<AdjointTime> adjointTime_;

    // Finite-differences calculation of derivatives.
    // std::vector<cl::TapeDouble>& indep --- vector of independent variables.
    // vector<Real>& sfFinite --- vector to store derivative values,
    // calculated by finite-differences scheme.
    // h --- step of approximation.
    double  finiteDiffCalculate(std::vector<cl::TapeDouble>& indep, vector<Real>& sfFinite, double h, cl::AdjointTestOutput& out)
    {
        Size size = indep.size();
        sfFinite.resize(size);
        out.log() << "Start differentiation using finite-differences:\t " << currentTime() << std::endl;
        //Finite differences (central formula)
        boost::timer timerDiff;

        //Creates bond object.
        FixedRateBond bond = bondCreate();

        // Evaluate derivatives using step-forward and step-backward values by the formula:
        // Derivative=(stepforward - stepbackward) / (2 * h), h-step.
        std::vector<cl::TapeDouble>::iterator it, it_fin;
        for (it = indep.begin(), it_fin = sfFinite.begin(); it != indep.end(); it++, it_fin++)
        {
            Real price = 0.0;

            // Evaluate bond value for a rate decreased by step h.
            price -= BondFunctions::cleanPrice(bond, *it,
                bondDayCount_,
                compounding_,
                frequency_);
            *it -= h;
            price += BondFunctions::cleanPrice(bond, *it,
                bondDayCount_,
                compounding_,
                frequency_);
            Real stepbackward = price;

            // Evaluate bond value for a rate increased by step h.
            price -= BondFunctions::cleanPrice(bond, *it,
                bondDayCount_,
                compounding_,
                frequency_);
            *it += 2 * h;
            price += BondFunctions::cleanPrice(bond, *it,
                bondDayCount_,
                compounding_,
                frequency_);
            Real stepforward = price;

            // Calculate approximated derivative:
            *it_fin = (stepforward - stepbackward) / (2 * h);
            *it -= h;
        }
        double timeCalculated_= timerDiff.elapsed();
        out.log() << "Time for differentiation using finite-differences" << timeCalculated_ << std::endl;
        // Return time  to calculate all derivatives in sfFinite vector.
        return timeCalculated_;
    }
};

// Struct to plot dependence of bond value on rate.
struct RateDependence
{
    static std::deque<std::string > get_columns()
    {
        static std::deque<std::string > columns =
        {
            "Rate", "bondPortfolioPrice"
        };

        return columns;
    }

    template <typename stream_type>
    friend inline stream_type&
        operator << (stream_type& stm, RateDependence& v)
    {
            stm << v.inputRate_
                << ";" << v.bondPortfolioPrice_ << std::endl;

            return stm;
        }

    Real inputRate_;
    Real bondPortfolioPrice_;
};

// Plot dependence of bond value on rate.
// CommonVars vars --- initial settings struct.
// pointsNum --- number of points on a plot.
// cl::AdjointTestOutput --- output stream to store plot generation info to a log file.
void rateDependencePlot(CommonVars vars, Size pointsNum, cl::AdjointTestOutput& output)
{

    std::vector<RateDependence> rateDependence_;
    //Init range for the  rate and output portfolio price
    for (Size i = 0; i<pointsNum; i++)
    {
        Real rate=0.025 + 0.002*i;

        // Calculate bond value.
        Real price = BondFunctions::cleanPrice(vars.bondCreate(), rate,
            vars.bondDayCount_,
            vars.compounding_,
            vars.frequency_);
        rateDependence_.push_back(RateDependence{rate, price});
    }
    output << rateDependence_;
    output.log() << "Plot Bond on rate dependence successfully generated"<<std::endl;
}

// Bond portfolio price is calculated using the formula:
// BondPortfolioPrice=\sum_{i=1}^N FaceVal*exp^{-r_i*t_i},
// where $r_i$ is a required return of the i-th  bond,
// FaceVal is a face value of i-th bond, $t_i$ is a maturity time; N is a number of bonds.

bool AdjointBondPortfolioTest::testBondPortfolio()
{
    bool result = false;
#ifdef CL_TAPE_CPPAD

    // Number of bonds in a portfolio.
    Size size = 1100;

    // Step to increase number of bonds on each consequtive iteration.
    Size step = 10;

    // Starting number of bonds
    Size bondsStart= size;

#ifdef CL_GRAPH_GEN
    // Starting number of bonds
     bondsStart =100;
#endif
    // Initialize settings for bonds
    CommonVars vars(size);

    // Vector to store derivatives calculated in Forward, Reverse modes and using finite-differences method.
    std::vector<double>sfForward;
    std::vector<double> sfReverse;
    std::vector<Real> sfFinite;

    // Plots streams.
    cl::AdjointTestOutput out("AdjointBondPortfolio\\AdjointBondPortfolioTest\\output",{
         { "filename", "BondonRate"}
        ,{ "not_clear", "Not" }
        ,{ "title", "Bond  price dependence on interest rate"}
        ,{ "ylabel", "Bond value" }
    });
    cl::AdjointTestOutput outSize("AdjointBondPortfolio\\AdjointBondPortfolioTest",{
         { "filename", "TapeSize"}
        ,{ "not_clear", "Not"}
        ,{ "title", "Tape size dependence on  number of interest rates"}
        ,{ "ylabel", "Memory (MB)"}
    });
    cl::AdjointTestOutput outPerform("AdjointBondPortfolio\\AdjointBondPortfolioTest",{
         { "filename", "BondPortfolio"}
        ,{ "not_clear", "Not" }
        ,{ "title", "Bond portfolio value differentiation performance with respect to interest rates "}
        ,{ "xlabel", "Number of bonds in a portfolio"}
        ,{ "ylabel", "Time (s)" }
    });
    cl::AdjointTestOutput outAdjoint("AdjointBondPortfolio\\AdjointBondPortfolioTest",{
         { "filename", "Adjoint"}
        ,{ "not_clear", "Not"}
        ,{ "title", "Bond portfolio value adjoint differentiation with respect to interest rates"}
        ,{ "ylabel", "Time (s)"}
    });

    // TapeSize on Number of bonds in a portfolio vector.
    std::vector<TapeSize> tapeSize;


    // Running Adjoint Differentiation for various numbers of bonds in a portfolio
    for (Size pos = bondsStart; pos <= size;)
    {
        // Beginning of tape recording
        boost::timer timer;
        vars.portfolioPrice_[0] = 0;

        // Create independent variables vector.
        std::vector<cl::TapeDouble> indepVar(vars.independent_.begin(), vars.independent_.begin()+pos);
        outPerform.log() << "Start of tape recording:\t " << currentTime() << std::endl;
        Independent(indepVar);
        for (Size i = 0; i<pos; i++)
        {
            //Calculate bond price
            Real price = BondFunctions::cleanPrice(vars.bondCreate(), indepVar[i],
                vars.bondDayCount_,
                vars.compounding_,
                vars.frequency_);

            // Adds bond to a portfolio.
            vars.portfolioPrice_[0] += price;
        }

        cl::TapeFunction<double> f(indepVar, vars.portfolioPrice_);

        outPerform.log() << "End of tape recording" << std::endl;

        // Memory for tape calculation.
        tapeSize.push_back(TapeSize{pos, f.Memory()});

        outPerform.log() << "Tape Memory in bytes:" << f.Memory() << std::endl;

        // End of tape recording. Differentiaion wiil be held with respect to  the independent variables vector.
        double timeTapeRecording = timer.elapsed();
        outPerform.log() << "Time for tape recording:\t " << timeTapeRecording << std::endl;

        // Start differentiation in Forward mode.
        gradForward(f, sfForward, outPerform, false, false);

        // Start differentiation in Reverse mode.
        double timeAdjoint = gradReverse(f, sfReverse, outPerform, false, false);

        // Central-difference scheme.
        double timeAnalytical = vars.finiteDiffCalculate(indepVar, sfFinite, vars.h_,outPerform);

        // Adding new data to the performace result vector.
        vars.performanceTime_.push_back(PerformanceTime{ timeTapeRecording, timeAdjoint, timeAnalytical, pos });
        vars.adjointTime_.push_back({ timeAdjoint, pos });

        outPerform.log() << std::endl;
        pos += step;
    }

    // Central finite-difference scheme check.
    // cl::AdjointTestOutput outPerform is an output stream
    // to write sfForward, sfReverse, sfFinite vectors of derivatives to a  log file.
    result = checkWithFiniteDiff(sfForward, sfReverse, sfFinite, outPerform, vars.relTolerance_, vars.absTolerance_);

    // output into CSV files.
    outSize << tapeSize;
    outPerform << vars.performanceTime_;
    outAdjoint << vars.adjointTime_;

    // Number of points on an output dependence plot
    Size pointsNumber = 100;
    // Plot the Bond  price dependence on interest rate plot.
    // cl::AdjointTestOutput output --- output stream for Bond on Rate dependence plot
    rateDependencePlot(vars, pointsNumber, out);
#endif
    return result;
}

test_suite*  AdjointBondPortfolioTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("CppADBond Portfolio  test");
    suite->add(QUANTLIB_TEST_CASE(&AdjointBondPortfolioTest::testBondPortfolio));
    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_BondPortfolioTest)

BOOST_AUTO_TEST_CASE(testBondPortfolio)
{
    BOOST_CHECK(AdjointBondPortfolioTest::testBondPortfolio());
}

BOOST_AUTO_TEST_SUITE_END()

#endif