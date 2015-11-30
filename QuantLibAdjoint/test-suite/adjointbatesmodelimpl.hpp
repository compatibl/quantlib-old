/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2005, 2008 Klaus Spanderen
Copyright (C) 2007 StatPro Italia srl
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

// Based on batesmodel.cpp file from test-suite.


#ifndef cl_adjoint_bates_model_impl_hpp
#define cl_adjoint_bates_model_impl_hpp
#pragma once


#include "adjointbatesmodeltest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include <ql/time/calendars/target.hpp>
#include <ql/processes/batesprocess.hpp>
#include <ql/processes/merton76process.hpp>
#include <ql/instruments/europeanoption.hpp>
#include <ql/time/daycounters/actualactual.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/yield/zerocurve.hpp>
#include <ql/pricingengines/blackformula.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
#include <ql/pricingengines/vanilla/batesengine.hpp>
#include <ql/pricingengines/vanilla/jumpdiffusionengine.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
#include <ql/pricingengines/vanilla/mceuropeanhestonengine.hpp>
#include <ql/pricingengines/vanilla/fdbatesvanillaengine.hpp>
#include <ql/models/equity/batesmodel.hpp>
#include <ql/models/equity/hestonmodelhelper.hpp>
#include <ql/time/period.hpp>
#include <ql/quotes/simplequote.hpp>
#include <boost/make_shared.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;

// Based on QuantLib\test-suite\batesmodel.cpp file

#define OUTPUT_FOLDER_NAME "AdjointBates"

/*!
AdjointBatesModelsTest is tested using adjoint differentiation for the class BatesProcess.

The Bates model (1996) combines the Merton and Heston approaches and proposes a model
with stochastic volatility and jumps governed by:
\f[
\begin{array}{rcl}
dS(t)  &=& (r-d-\lambda m) S dt +\sqrt{v} S dW_1 + (\exp^J - 1) S dN \\
dv(t)  &=& \kappa (\theta - v) dt + \sigma \sqrt{v} dW_2 \\
dW_1 dW_2 &=& \rho dt \\
\omega(J) &=& \frac{1}{\sqrt{2\pi \delta^2}}
\exp\left[-\frac{(J-\nu)^2}{2\delta^2}\right]
\end{array}
\f]
Here:
\f$ S(t) \f$                    is an asset (underlying) price,
\f$ v(t) \f$                    is volatility,
\f$ J \f$ and \f$ \omega(J) \f$ are a jump size and its PDF
\f$ r \f$                       is a risk-free (domestic) interest rate,
\f$ d \f$                       is a dividend (foreign interest) rate,
\f$ m \f$                       is an average jump amplitude
\f$ \sigma \f$                  is a volatility of volatility,
\f$ \kappa \f$                  is a mean-reverting rate,
\f$ \theta \f$                  is a long-term variance,
\f$ \rho \f$                    is correlation between price and variance (Wiener processes \f$ dW_1 \f$ and \f$ dW_2 \f$),
\f$ \lambda \f$, \f$ \nu \f$, \f$ \delta \f$ are jump intensity, mean and volatility for the Merton jump-diffusion, respectively.

References:
A. Sepp, Pricing European-Style Options under Jump Diffusion
Processes with Stochastic Volatility: Applications of Fourier
Transform (<http://math.ut.ee/~spartak/papers/stochjumpvols.pdf>)

In this test, the dependence of the NPV of a portfolio consisted of several options
on option strike prices is calculated using adjoint and finite difference differentiation.
In addition, the NPV is checked with the Black formula
(the parameters of the Bates model are chosen to produce results identical to the Black formula).
*/

namespace {
    enum
    {
#if defined CL_GRAPH_GEN
        // Number of points for dependency plots.
        pointNo = 100,
        // Number of points for performance plot.
        iterNo = 30,
        // Step for portfolio size for performance testing .
        step = 1,
#else
        // Number of points for dependency plots.
        pointNo = 0,
        // Number of points for performance plot.
        iterNo = 0,
        // Step for portfolio size for performance testing .
        step = 1,
#endif
        // Defines performance accuracy. Its value is a minimum number
        // of calling of O(1) complexity methods per one performance test.
        iterNumFactor = 0,
    };

    // This structure contains all constant variables
    // Parameters of the Bates model are chosen to produce results identical to the Black formula
    // (based on BatesModelTest::testAnalyticVsBlack() from QuantLib\test-suite\batesmodel.cpp)
    struct CommonVars
    {
        CommonVars()
        : settlementDate_(Date(1, September, 2005))
        , exerciseDate_(settlementDate_ + 6 * Months)
        , riskfreeRate_(0.10)
        , dividentRate_(0.04)
        , spotPrice_(30.0)
        , volatility_(0.05)
        , kappa_(5.0)
        , theta_(0.05)
        , sigma_(1.0e-4)
        , rho_(0.0)
        , lambda_(0.0001)
        , nu_(0.0)
        , delta_(0.0001)
        , exercise_(boost::make_shared<EuropeanExercise>(exerciseDate_))
        , strikePrice_(18.0)
        {
        }

        // global data
        const Date settlementDate_;
        const Date exerciseDate_;
        const Real riskfreeRate_;
        const Real dividentRate_;
        const Real spotPrice_;
        const Real volatility_;
        // Bates model specific data
        const Real kappa_;
        const Real theta_;
        const Real sigma_;
        const Real rho_;
        const Real lambda_;
        const Real nu_;
        const Real delta_;
        // variables and objects needed for calculations
        boost::shared_ptr<Exercise> exercise_;
        const Real strikePrice_;
        Handle<YieldTermStructure> riskFreeTS_;
    };

    // This structure contains all constant variables from CommonVars + strike price
    // In addition it provides methods for calcualtions
    struct BatesData
        : public CommonVars
    {
        explicit BatesData()
        : CommonVars()
        {
        }

        // Provides payoffs that are different only by their strike price
        boost::shared_ptr<StrikedTypePayoff> makePayoff(Rate strikePrice)
        {
            boost::shared_ptr<StrikedTypePayoff> payoff(
                new PlainVanillaPayoff(Option::Put, strikePrice));
            return payoff;
        }

        // Calculate forward price for given divident rate
        Real CalculateForwardPrice(Real dividentRate)
        {
            Settings::instance().evaluationDate() = settlementDate_;
            DayCounter dayCounter = ActualActual();
            Handle<Quote> s0(boost::make_shared<SimpleQuote>(spotPrice_));
            Real yearFraction = dayCounter.yearFraction(settlementDate_, exerciseDate_);
            Real forwardPrice = s0->value()*std::exp((riskfreeRate_ - dividentRate)*yearFraction);
            return forwardPrice;
        }

        // Provides options that are different only by their strike price
        boost::shared_ptr<VanillaOption> makeOptionForStrikePrice(Real strikePrice)
        {
            Settings::instance().evaluationDate() = settlementDate_;
            DayCounter dayCounter = ActualActual();

            boost::shared_ptr<StrikedTypePayoff> payoff = makePayoff(strikePrice);
            boost::shared_ptr<VanillaOption> option(new VanillaOption(payoff, exercise_));
            Handle<Quote> s0(boost::make_shared<SimpleQuote>(spotPrice_));

            Handle<YieldTermStructure> riskFreeTS(flatRate(riskfreeRate_, dayCounter));
            Handle<YieldTermStructure> dividendTS(flatRate(dividentRate_, dayCounter));
            boost::shared_ptr<BatesProcess> process(
                new BatesProcess(riskFreeTS, dividendTS, s0, volatility_,
                kappa_, theta_, sigma_, rho_, lambda_, nu_, delta_));

            boost::shared_ptr<PricingEngine> engine(new BatesEngine(
                boost::shared_ptr<BatesModel>(new BatesModel(process)), 64));

            option->setPricingEngine(engine);

            return option;
        }

        // Returns net present value (NPV) of option with given fixed rate strike and other
        // parameters determined from class instance.
        Real optionNPVForStrikePrice(Real strikePrice)
        {
            auto option = makeOptionForStrikePrice(strikePrice);
            return option->NPV();
        }

        // Black-Scholes result (for cross check)
        Real optionNPVForStrikePrice_BS(Real strikePrice)
        {
            Settings::instance().evaluationDate() = settlementDate_;
            DayCounter dayCounter = ActualActual();
            Real yearFraction = dayCounter.yearFraction(settlementDate_, exerciseDate_);
            boost::shared_ptr<StrikedTypePayoff> payoff = makePayoff(strikePrice);
            Real forwardPrice = CalculateForwardPrice(dividentRate_);
            Real NPV = blackFormula(payoff->optionType(), payoff->strike(),
                forwardPrice, std::sqrt(volatility_*yearFraction)) *
                std::exp(-1 * riskfreeRate_*yearFraction);
            return NPV;
        }

        // Provides options that are different only by the divident rate
        boost::shared_ptr<VanillaOption> makeOptionForDividentRate(Real dividentRate)
        {
            Settings::instance().evaluationDate() = settlementDate_;
            DayCounter dayCounter = ActualActual();

            boost::shared_ptr<StrikedTypePayoff> payoff = makePayoff(strikePrice_);
            boost::shared_ptr<VanillaOption> option(new VanillaOption(payoff, exercise_));
            Handle<Quote> s0(boost::make_shared<SimpleQuote>(spotPrice_));

            Handle<YieldTermStructure> riskFreeTS(flatRate(riskfreeRate_, dayCounter));
            Handle<YieldTermStructure> dividendTS(flatRate(dividentRate, dayCounter));
            boost::shared_ptr<BatesProcess> process(
                new BatesProcess(riskFreeTS, dividendTS, s0, volatility_,
                kappa_, theta_, sigma_, rho_, lambda_, nu_, delta_));

            boost::shared_ptr<PricingEngine> engine(new BatesEngine(
                boost::shared_ptr<BatesModel>(new BatesModel(process)), 64));

            option->setPricingEngine(engine);

            return option;
        }

        // Returns net present value (NPV) of option with given fixed rate strike and other
        // parameters determined from class instance.
        Real optionNPVForDividentRate(Real dividentRate)
        {
            auto option = makeOptionForDividentRate(dividentRate);
            return option->NPV();
        }

        // Black-Scholes result (for cross check)
        Real optionNPVForDividentRate_BS(Real dividentRate)
        {
            Settings::instance().evaluationDate() = settlementDate_;
            DayCounter dayCounter = ActualActual();
            Real yearFraction = dayCounter.yearFraction(settlementDate_, exerciseDate_);
            boost::shared_ptr<StrikedTypePayoff> payoff = makePayoff(strikePrice_);
            Real forwardPrice = CalculateForwardPrice(dividentRate);
            Real NPV = blackFormula(payoff->optionType(), payoff->strike(),
                forwardPrice, std::sqrt(volatility_*yearFraction)) *
                std::exp(-1*riskfreeRate_*yearFraction);
            return NPV;
        }
    };

    // Struct for graphics recording of NPV dependence on anything
    // (taken from adjointbermudanswaptionimpl.hpp)
    struct NPVonSmth
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "Independent variable", ""
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, NPVonSmth& v)
        {
                stm << v.var_
                    << ";" << v.sensitivity_
                    << std::endl;
                return stm;
            }

        Real var_;
        Real sensitivity_;
    };

    // Bates test structure (dependence on strike rate)
    struct BatesTestStrikePrice
        : public cl::AdjointTest<BatesTestStrikePrice>
    {
        // Create test with a portfolio given by the sum of size options with strike prices in the range:
        //                                                    data_->strikePrice_ + 0.5 * i, 0 < i < size
        BatesTestStrikePrice(Size size, BatesData* data, cl::tape_empty_test_output* logger = nullptr)
        : AdjointTest()
        , size_(size)
        , data_(data)
        , strikes_(size)
        , NPV_(size)
        , totalNPV_()
        {
            setLogger(logger);

            for (Size i = 0; i < size_; i++)
            {
                strikes_[i] = data_->strikePrice_ + 0.5 * i;
            }
        }

        // Do tape recording for adjoint calculation
        void recordTape()
        {
            cl::Independent(strikes_);
            calculateTotalNPV();
            f_ = std::make_unique<cl::tape_function<double>>(strikes_, totalNPV_);
        }

        // Calculates price of  each option and portfolio (the sum)
        Real calculateTotalNPV()
        {
            totalNPV_.resize(1, 0);
            NPV_.resize(size_);
            for (Size i = 0; i < size_; i++)
            {
                NPV_[i] = data_->optionNPVForStrikePrice(strikes_[i]);
                totalNPV_.front() += NPV_[i];
            }
            return totalNPV_.front();
        }

        // Calculates price of  each option and portfolio (the sum) using Black formula (for cross check)
        Real calculateTotalNPV_BS()
        {
            Real sum = 0.0;
            for (Size i = 0; i < size_; i++)
            {
                sum += data_->optionNPVForStrikePrice_BS(strikes_[i]);
            }
            return sum;
        }

        // Check with BS formula
        bool CheckBatesVsBlack()
        {
            Real calculated = this->calculateTotalNPV();
            Real expected = this->calculateTotalNPV_BS();
            Real maxabs = std::max(std::abs(calculated), std::abs(expected));
            Real tol = std::max(this->relativeTol() * maxabs, this->absTol());
            return (std::abs(expected - calculated) < tol);
        }

        // Calculates derivatives using finite difference method
        void calcAnalytical()
        {
            const Real diff = 1e-6;
            Real h = diff * data_->strikePrice_;  // shift for finite diff. method
            analyticalResults_.resize(size_);
            for (Size i = 0; i < size_; i++)
            {
                analyticalResults_[i] = (data_->optionNPVForStrikePrice(strikes_[i] + h)
                    - data_->optionNPVForStrikePrice(strikes_[i] - h)) / (2 * h);
            }
        }

        // get relative tolerance
        double relativeTol() const { return 1e-5; }

        // get absolute tolerance
        double absTol() const { return 1e-7; }

        // get number of independent variables
        Size indepVarNumber() { return size_; }

        Size size_;
        BatesData* data_;
        std::vector<cl::tape_double> strikes_;
        std::vector<Real> NPV_;
        std::vector<cl::tape_double> totalNPV_;
    };


    // Bates test structure (dependence on divident rate)
    struct BatesTestDividentRate
        : public cl::AdjointTest<BatesTestDividentRate>
    {
        // Create test with a portfolio given by the sum of options on assets with divident rates in the range:
        //                                                    data_->dividentRate_ + 0.01 * i, 0 < i < size
        BatesTestDividentRate(Size size, BatesData* data, cl::tape_empty_test_output* logger = nullptr)
        : AdjointTest()
        , size_(size)
        , data_(data)
        , dividents_(size)
        , NPV_(size)
        , totalNPV_()
        {
            setLogger(logger);

            for (Size i = 0; i < size_; i++)
            {
                dividents_[i] = data_->dividentRate_ + 0.005 * i;
            }
        }

        // Do tape recording for adjoint calculation
        void recordTape()
        {
            cl::Independent(dividents_);
            calculateTotalNPV();
            f_ = std::make_unique<cl::tape_function<double>>(dividents_, totalNPV_);
        }

        // Calculates price of  each option and portfolio (the sum)
        Real calculateTotalNPV()
        {
            totalNPV_.resize(1, 0);
            NPV_.resize(size_);
            for (Size i = 0; i < size_; i++)
            {
                NPV_[i] = data_->optionNPVForDividentRate(dividents_[i]);
                totalNPV_.front() += NPV_[i];
            }
            return totalNPV_.front();
        }

        // Calculates price of each option and portfolio (the sum) using Black formula (for cross check)
        Real calculateTotalNPV_BS()
        {
            Real sum = 0.0;
            for (Size i = 0; i < size_; i++)
            {
                sum += data_->optionNPVForDividentRate_BS(dividents_[i]);
            }
            return sum;
        }

        // Check with BS formula
        bool CheckBatesVsBlack()
        {
            Real calculated = this->calculateTotalNPV();
            Real expected = this->calculateTotalNPV_BS();
            Real maxabs = std::max(std::abs(calculated), std::abs(expected));
            Real tol = std::max(this->relativeTol() * maxabs, this->absTol());
            return (std::abs(expected - calculated) < tol);
        }

        // Calculates derivatives using finite difference method
        void calcAnalytical()
        {
            const Real diff = 1e-4;
            Real h = diff * data_->dividentRate_;  // shift for finite diff. method
            analyticalResults_.resize(size_);
            for (Size i = 0; i < size_; i++)
            {
                analyticalResults_[i] = (data_->optionNPVForDividentRate(dividents_[i] + h)
                    - data_->optionNPVForDividentRate(dividents_[i] - h)) / (2 * h);
            }
        }

        // get relative tolerance
        double relativeTol() const { return 1e-5; }

        // get absolute tolerance
        double absTol() const { return 1e-7; }

        // get number of independent variables
        Size indepVarNumber() { return size_; }

        Size size_;
        BatesData* data_;
        std::vector<cl::tape_double> dividents_;
        std::vector<Real> NPV_;
        std::vector<cl::tape_double> totalNPV_;
    };

    // BatesData structure extended for making plots with dependence on the divident rate
    struct TestDataDividentRate
        : public BatesData
    {
        explicit TestDataDividentRate()
        : outPerform_(OUTPUT_FOLDER_NAME "//" + std::string("DividentRate"), {
            { "title", "Bates NPV differentiation performance with respect to divident rate" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of divident rates" }
            , { "line_box_width", "-5" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        })
            , outAdjoint_(OUTPUT_FOLDER_NAME "//" + std::string("DividentRate"), {
                { "title", "Bates NPV adjoint differentiation with respect to fixed divident rate" }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of divident rates" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        })
            , outSize_(OUTPUT_FOLDER_NAME "//" + std::string("DividentRate"), {
                { "title", "Tape size dependence on number of divident rates" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
            , { "cleanlog", "false" }
        })
            , out_(OUTPUT_FOLDER_NAME "//" + std::string("DividentRate") + "//output", {
                { "filename", "NPVonDividentRate" }
            , { "not_clear", "Not" }
            , { "title", "Option NPV dependence on fixed divident rate" }
            , { "ylabel", "Option NPV" }
            , { "xlabel", "Divident rate" }
            , { "cleanlog", "false" }
        })
        { }

        // main method to produce output
        bool makeOutput()
        {
            bool ok = true;
            if (pointNo > 0)
            {
                ok &= recordDependencePlot(pointNo);
                ok &= cl::recordPerformance(*this, iterNo, step);
            }
            return ok;
        }

        std::shared_ptr<BatesTestDividentRate> getTest(size_t size)
        {
            return std::make_shared<BatesTestDividentRate>(size, this, &outPerform_);
        }

        // Makes graphics for divident sensitivity dependence
        bool recordDependencePlot(size_t pointNo)
        {
            bool ok = true;
            std::vector<NPVonSmth> outData(pointNo);
            auto test = getTest(pointNo);
            test->recordTape();
            for (Size i = 0; i < pointNo; i++)
            {
                outData[i] = { test->dividents_[i], test->NPV_[i] };
            }
            out_ << outData;
            return ok;
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output out_;
    }; 
}


#endif
