/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2005, 2006, 2007, 2008, 2009 StatPro Italia srl
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

//based on piecewiseyieldcurve.cpp from test-suite

#include "adjointpiecewiseyieldcurvetest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include <ql/quantlib.hpp>


using namespace QuantLib;
using namespace boost::unit_test_framework;

namespace
{
    struct Datum
    {
        Integer n_;
        TimeUnit units_;
        Rate rate_;
    };

    struct BondDatum
    {
        Integer n_;
        TimeUnit units_;
        Integer length_;
        Frequency frequency_;
        Rate coupon_;
        Real price_;
    };

    struct CommonVars
    {
        // Constructor.
        CommonVars()
        {
            calendar_ = TARGET();
            settlementDays_ = 2;
            today_ = calendar_.adjust(Date::todaysDate());
            Settings::instance().evaluationDate() = today_;
            settlement_ = calendar_.advance(today_, settlementDays_, Days);
            fixedLegConvention_ = Unadjusted;
            fixedLegFrequency_ = Annual;
            fixedLegDayCounter_ = Thirty360();
            bondSettlementDays_ = 3;
            bondDayCounter_ = ActualActual();
            bondConvention_ = Following;
            bondRedemption_ = 100.0;
            bmaFrequency_ = Quarterly;
            bmaConvention_ = Following;
            bmaDayCounter_ = ActualActual();

            depositData_ = {
                { 1, Weeks, 4.559 },
                { 1, Months, 4.581 },
                { 2, Months, 4.573 },
                { 3, Months, 4.557 },
                { 6, Months, 4.496 },
                { 9, Months, 4.490 }
            };

            Size maxSize = 50;
            Rate startRate = 4.5;
            Rate maxRate = 6;
            Real step = (maxRate - startRate) / maxSize;

            for (Size i = 0; i < maxSize; i++)
            {
                Datum swapData;
                swapData.n_ = i + 1;
                swapData.units_ = Years;
                swapData.rate_ = startRate + i*step;
                swapData_.push_back(swapData);
            }

            deposits_ = depositData_.size();
            swaps_ = swapData_.size();

            // Market elements.
            rates_ =
                std::vector<boost::shared_ptr<SimpleQuote> >(deposits_ + swaps_);

            for (Size i = 0; i < deposits_; i++)
            {
                rates_[i] = boost::shared_ptr<SimpleQuote>(
                    new SimpleQuote(depositData_[i].rate_ / 100));
            }
            for (Size i = 0; i < swaps_; i++)
            {
                rates_[i + deposits_] = boost::shared_ptr<SimpleQuote>(
                    new SimpleQuote(swapData_[i].rate_ / 100));
            }


            // Rate helpers.
            instruments_ =
                std::vector<boost::shared_ptr<RateHelper> >(deposits_ + swaps_);

            boost::shared_ptr<IborIndex> euribor6m(new Euribor6M);
            for (Size i = 0; i < deposits_; i++)
            {
                Handle<Quote> r(rates_[i]);
                instruments_[i] = boost::shared_ptr<RateHelper>(new
                                                                DepositRateHelper(r, depositData_[i].n_*depositData_[i].units_,
                                                                euribor6m->fixingDays(), calendar_,
                                                                euribor6m->businessDayConvention(),
                                                                euribor6m->endOfMonth(),
                                                                euribor6m->dayCounter()));
            }
            for (Size i = 0; i < swaps_; i++)
            {
                Handle<Quote> r(rates_[i + deposits_]);
                instruments_[i + deposits_] = boost::shared_ptr<RateHelper>(new
                                                                            SwapRateHelper(r, swapData_[i].n_*swapData_[i].units_,
                                                                            calendar_,
                                                                            fixedLegFrequency_, fixedLegConvention_,
                                                                            fixedLegDayCounter_, euribor6m));
            }



        }

        // Calculate total NPV of swaps based on input rates using Japanese Yen LIBOR rate.
        // Japanese Yen LIBOR rate can be considered as the interbank cost of borrowing funds in Japanese yens.
        Real calculateSwapNPV(std::vector<Real> swapRate)
        {
            Real NPV = 0;

            // Create vector of simple quotes
            rates_ = std::vector < boost::shared_ptr<SimpleQuote> >(swapRate.size());
            for (Size i = 0; i < swapRate.size(); i++)
            {
                rates_[i] = boost::shared_ptr<SimpleQuote>(new SimpleQuote(swapRate[i]));
            }

            // Create vector of rate helpers.
            instruments_ = std::vector<boost::shared_ptr<RateHelper> >(swapRate.size());

            // Create ibor index using 6 Month JPY Libor. 
            boost::shared_ptr<IborIndex> index(new JPYLibor(6 * Months));
            for (Size i = 0; i < swapRate.size(); i++)
            {
                Handle<Quote> r(rates_[i]);

                // Add new swap rate helper to instruments.
                instruments_[i] = boost::shared_ptr<RateHelper>(
                    new SwapRateHelper(r, swapData_[i].n_*swapData_[i].units_,
                    calendar_,
                    fixedLegFrequency_, fixedLegConvention_,
                    fixedLegDayCounter_, index));
            }

            // Initialize yield term structure.
            termStructure_ = boost::shared_ptr<YieldTermStructure>(
                new PiecewiseYieldCurve<Discount, LogLinear>(
                settlement_, instruments_,
                Actual360(),
                1.0e-12));


            RelinkableHandle<YieldTermStructure> curveHandle;
            curveHandle.linkTo(termStructure_);

            // Create ibor index using 6 Month JPY Libor. 
            boost::shared_ptr<IborIndex> jpylibor6m(new JPYLibor(6 * Months, curveHandle));
            for (Size i = 0; i < swapRate.size(); i++)
            {
                Period tenor = swapData_[i].n_*swapData_[i].units_;

                // Create vanilla swap
                VanillaSwap swap = MakeVanillaSwap(tenor, jpylibor6m, 0.0)
                    .withEffectiveDate(settlement_)
                    .withFixedLegDayCount(fixedLegDayCounter_)
                    .withFixedLegTenor(Period(fixedLegFrequency_))
                    .withFixedLegConvention(fixedLegConvention_)
                    .withFixedLegTerminationDateConvention(fixedLegConvention_)
                    .withFixedLegCalendar(calendar_)
                    .withFloatingLegCalendar(calendar_);

                NPV += swap.NPV();
            }
            return NPV;
        }

        // Global variables.
        Calendar calendar_;
        Natural settlementDays_;
        Date today_, settlement_;
        BusinessDayConvention fixedLegConvention_;
        Frequency fixedLegFrequency_;
        DayCounter fixedLegDayCounter_;
        Natural bondSettlementDays_;
        DayCounter bondDayCounter_;
        BusinessDayConvention bondConvention_;
        Real bondRedemption_;
        Frequency bmaFrequency_;
        BusinessDayConvention bmaConvention_;
        DayCounter bmaDayCounter_;

        Size deposits_, fras_, swaps_, bonds_, bmas_;
        std::vector<boost::shared_ptr<SimpleQuote> > rates_, fraRates_,
            prices_, fractions_;
        std::vector<boost::shared_ptr<RateHelper> > instruments_, fraHelpers_,
            bondHelpers_, bmaHelpers_;
        std::vector<Schedule> schedules_;
        boost::shared_ptr<YieldTermStructure> termStructure_;

        std::vector<Datum> depositData_;
        std::vector<Datum> swapData_;


        // Cleanup.
        SavedSettings backup_;
        IndexHistoryCleaner cleaner_;
    };

    struct Variation
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
            operator << (stream_type& stm, Variation& v)
        {
                stm << v.rate_
                    << ";" << v.swapnNPV_
                    << std::endl;
                return stm;
            }

        Real rate_;
        Real swapnNPV_;
    };
}

// Test bootstrap over Japanese Yen LIBOR swaps.
bool AdjointPiecewiseYieldCurveTest::testJpyLibor()
{
    BOOST_TEST_MESSAGE("Testing bootstrap over JPY LIBOR swaps...");

    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVars vars;

    // Create output streams for plots.
    cl::AdjointTestOutput outAdjointPerformance("AdjointPiecewiseYieldCurve//JpyLibor"
                                                , { { "filename", "AdjointPerformance" }
    , { "not_clear", "Not" }
    , { "line_box_width", "-5" }
    , { "title", "Swap NPV differentiation performance with respect to rate" }
    , { "ylabel", "Time (s)" }
    , { "xlabel", "Number of rates" } });

    cl::AdjointTestOutput outAdjoint("AdjointPiecewiseYieldCurve//JpyLibor"
                                     , { { "filename", "Adjoint" }
    , { "not_clear", "Not" }
    , { "title", "Swap NPV adjoint differentiation with respect to rate" }
    , { "ylabel", "Time (s)" }
    , { "xlabel", "Number of rates" } });

    cl::AdjointTestOutput outTapeSize("AdjointPiecewiseYieldCurve//JpyLibor"
                                      , { { "filename", "TapeSize" }
    , { "not_clear", "Not" }
    , { "title", "Tape size dependence on number of rates" }
    , { "ylabel", "Size (MB)" }
    , { "xlabel", "Number of rates" } });

    cl::AdjointTestOutput outVariation("AdjointPiecewiseYieldCurve//JpyLibor//output"
                                       , { { "filename", "SwapNPVonRate" }
    , { "ylabel", "Swap NPV" }
    , { "not_clear", "Not" }
    , { "title", "Swap NPV dependence on rate" }
    , { "xlabel", "Rate" } });

    // Set today date.
    vars.today_ = Date(4, October, 2007);
    Settings::instance().evaluationDate() = vars.today_;

    // Set calendar type.
    vars.calendar_ = Japan();

    // Set settlement date.
    vars.settlement_ =
        vars.calendar_.advance(vars.today_, vars.settlementDays_, Days);

    // Create vectors for store performance results.
    std::vector<PerformanceTime> performResults;
    std::vector<AdjointTime> adjointResults;
    std::vector<TapeSize> tapeSizeResults;
    boost::timer timer;

    Size sizeof_indep;


#if defined CL_GRAPH_GEN
    Size startSize = 0;
    Size maxSize = vars.swapData_.size();
#else
    Size startSize = 5;
    Size maxSize = startSize + 1;
#endif

    // Variate number of independent variables.
    for (Size s = startSize; s < maxSize; s += 2)
    {
        // Create temporary structures for store performance results.
        PerformanceTime perfomTime;
        AdjointTime adTime;
        TapeSize tSize;

        // Start timing of tape recording.
        outAdjointPerformance.log() << "Start taping : " << currentTime() << std::endl;
        timer.restart();
        sizeof_indep = s + 1;

        // Store number of independent variables
        perfomTime.indepVarNumber_ = sizeof_indep;
        adTime.indepVarNumber_ = sizeof_indep;
        tSize.indepVarNumber_ = sizeof_indep;

        // Create vector of independent variables.
        std::vector<cl::TapeDouble> swapRate(sizeof_indep);

        // Initialize vector of independent variables
        for (Size i = 0; i < sizeof_indep; i++)
            swapRate[i] = vars.swapData_[i].rate_ / 100;

        // Start taping. Declare rates as independent variables.
        Independent(swapRate);

        // Create vector of dependent variables.
        std::vector<cl::TapeDouble> Y(1);

        // Calculate sum of swap Net Present Values.
        Y[0] = vars.calculateSwapNPV(swapRate);

        // End of tape recording.Declare sum of swap NPVs as dependent variable.
        // Differentiaion will be held with respect to the independent variables vector.
        cl::TapeFunction<double> f(swapRate, Y);

        // Store time of tape recording.
        perfomTime.timeTapeRecording_ = timer.elapsed();
        outAdjointPerformance.log() << "End of tape recording. Time for tape recording : " << perfomTime.timeTapeRecording_ << std::endl;
        // Store size of tape.
        tSize.memory_ = f.Memory();

        std::vector<double> sf_Forward, sf_Reverse;

        //Start differentiation in Forward mode.
        double tf = gradForward(f, sf_Forward, outAdjointPerformance, false, false);

        //Start differentiation in Reverse mode.
        double tr = gradReverse(f, sf_Reverse, outAdjointPerformance, false, false);

        // Store time of adjoint differntiation.
        perfomTime.timeAdjoint_ = std::min(tf, tr);
        adTime.timeAdjoint_ = perfomTime.timeAdjoint_;

        // Start timing for calculating derivatives by finite difference.
        outAdjointPerformance.log() << "Start differentiation using central finite differences : " << currentTime() << std::endl;
        timer.restart();

        // Start differentiation using finite differences method.
        double h = 1.0e-4;
        double tol = 1.0e-2;
        //  Create vector for derivatives
        std::vector<Real> sf_Finite(sizeof_indep);

        // Calculate sum of swap NPVs.
        Real YF = vars.calculateSwapNPV(swapRate);

        for (Size i = 0; i < sizeof_indep; i++)
        {
            swapRate[i] += h;
            //Evaluate derivative using finite difference.
            sf_Finite[i] = (vars.calculateSwapNPV(swapRate) - YF) / h;
            swapRate[i] -= h;
        }

        // Store time of finite difference differntiation.
        perfomTime.timeAnalytical_ = timer.elapsed();
        outAdjointPerformance.log() << "Time for differentiation using central finite differences : " << perfomTime.timeAnalytical_ << " s" << std::endl;
        //Adding new data to the result vectors.
        performResults.push_back(perfomTime);
        adjointResults.push_back(adTime);
        tapeSizeResults.push_back(tSize);

        // Check derivatives calculated by adjoint and finite differences methods.
        result = checkWithFiniteDiff(sf_Forward, sf_Reverse, sf_Finite, tol, tol);
    }


    std::vector<Variation> variationResults;
#ifdef CL_GRAPH_GEN
    // Calculate swaption NPV dependence on variation of volatility.
    for (Size i = 0; i < sizeof_indep; i++)
    {
        std::vector<Real> swapRate(1);
        swapRate[0] = vars.swapData_[i].rate_ / 100;
        Variation var;
        var.rate_ = swapRate[0];
        var.swapnNPV_ = vars.calculateSwapNPV(swapRate);
        variationResults.push_back(var);

    }
#endif

    //Output results to .csv and .plt files.
    outAdjointPerformance << performResults;
    outAdjoint << adjointResults;
    outTapeSize << tapeSizeResults;
    outVariation << variationResults;
#endif
    return result;
}


test_suite* AdjointPiecewiseYieldCurveTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Piecewise yield curve tests");
    suite->add(QUANTLIB_TEST_CASE(&AdjointPiecewiseYieldCurveTest::testJpyLibor));
    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_piecewise)

BOOST_AUTO_TEST_CASE(testJpyLiborSwapNpvs)
{
    BOOST_CHECK(AdjointPiecewiseYieldCurveTest::testJpyLibor());
}

BOOST_AUTO_TEST_SUITE_END()

#endif