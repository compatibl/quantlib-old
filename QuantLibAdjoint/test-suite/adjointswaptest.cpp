/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2003, 2004, 2007 StatPro Italia srl
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

// based on swap.cpp file from test-suite

#include <ql/quantlib.hpp>
#include "utilities.hpp"
#include "adjointswaptest.hpp"
#include "adjointtestutilities.hpp"

using namespace QuantLib;
using namespace boost::unit_test_framework;

namespace
{
    struct CommonVars
    {
        // global data
        Date today_, settlement_;
        VanillaSwap::Type type_;
        Real nominal_;
        Calendar calendar_;
        BusinessDayConvention fixedConvention_, floatingConvention_;
        Frequency fixedFrequency_, floatingFrequency_;
        DayCounter fixedDayCount_;
        boost::shared_ptr<IborIndex> index_;
        Natural settlementDays_;
        RelinkableHandle<YieldTermStructure> termStructure_;
        std::vector<Integer> lengths_;
        std::vector<Rate> rates_;
        std::vector<Spread> spreads_;

        // cleanup
        SavedSettings backup_;

        //constructor
        CommonVars()
        {
            type_ = VanillaSwap::Payer;
            settlementDays_ = 2;
            nominal_ = 100.0;
            fixedConvention_ = Unadjusted;
            floatingConvention_ = ModifiedFollowing;
            fixedFrequency_ = Annual;
            floatingFrequency_ = Semiannual;
            fixedDayCount_ = Thirty360();
            index_ = boost::shared_ptr<IborIndex>(new
                                                 Euribor(Period(floatingFrequency_), termStructure_));
            calendar_ = index_->fixingCalendar();
            today_ = calendar_.adjust(Date::todaysDate());
            Settings::instance().evaluationDate() = today_;
            settlement_ = calendar_.advance(today_, settlementDays_, Days);
            termStructure_.linkTo(flatRate(settlement_, 0.05, Actual365Fixed()));
            lengths_ = { 1, 2, 5, 10, 20 };
            rates_ = { 0.04, 0.05, 0.06, 0.07};
            spreads_ = { -0.01, -0.001, 0.001, 0.01 };
        }

        // utilities
        boost::shared_ptr<VanillaSwap> makeSwap(Integer length
                                                , Rate fixedRate
                                                , Spread floatingSpread)
        {
            Date maturity = calendar_.advance(settlement_, length, Years,
                                             floatingConvention_);
            Schedule fixedSchedule(settlement_, maturity, Period(fixedFrequency_),
                                   calendar_, fixedConvention_, fixedConvention_,
                                   DateGeneration::Forward, false);
            Schedule floatSchedule(settlement_, maturity,
                                   Period(floatingFrequency_),
                                   calendar_, floatingConvention_,
                                   floatingConvention_,
                                   DateGeneration::Forward, false);
            boost::shared_ptr<VanillaSwap> swap(
                new VanillaSwap(type_, nominal_,
                fixedSchedule, fixedRate, fixedDayCount_,
                floatSchedule, index_, floatingSpread,
                index_->dayCounter()));
            swap->setPricingEngine(boost::shared_ptr<PricingEngine>(
                new DiscountingSwapEngine(termStructure_)));
            return swap;
        }


        Real FairRateNPV(std::vector<Real> spread, std::vector<Spread> spreads)
        {
            Real Y = 0;
            for (Size i = 0; i < lengths_.size(); i++)
            {
                for (Size j = 0; j < spread.size(); j++)
                {
                    boost::shared_ptr<VanillaSwap> swap = makeSwap(lengths_[i], 0.0, spreads[j]);

                    swap = makeSwap(lengths_[i], swap->fairRate(), spread[j]);

                    if (std::fabs(swap->NPV()) > 1.0e-10)
                    {
                        BOOST_ERROR("recalculating with implied rate:\n"
                                    << std::setprecision(2)
                                    << "    length: " << lengths_[i] << " years\n"
                                    << "    floating spread: "
                                    << io::rate(spreads[j]) << "\n"
                                    << "    swap value: " << swap->NPV());
                    }

                    Y += swap->NPV();
                }
            }

            return Y;
        }

        Real FairSpreadNPV(std::vector<Real> rate, std::vector<Rate> rates)
        {
            Real Y = 0;

            for (Size i = 0; i < lengths_.size(); i++)
            {
                for (Size j = 0; j < rate.size(); j++)
                {
                    boost::shared_ptr<VanillaSwap> swap = makeSwap(lengths_[i], rates[j], 0.0);

                    swap = makeSwap(lengths_[i], rate[j], swap->fairSpread());

                    if (std::fabs(swap->NPV()) > 1.0e-10)
                    {
                        BOOST_ERROR("recalculating with implied spread:\n"
                                    << std::setprecision(2)
                                    << "    length: " << lengths_[i] << " years\n"
                                    << "    fixed rate: " << io::rate(rate[j]) << "\n"
                                    << "    swap value: " << swap->NPV());
                    }

                    Y += swap->NPV();
                }
            }
            return Y;
        }

        Real RateDependencyNPV(std::vector<Real> rate)
        {
            Real Y = 0;

            for (Size i = 0; i < lengths_.size(); i++)
            {
                for (Size j = 0; j < spreads_.size(); j++)
                {
                    // store the results for different rates...
                    std::vector<Real> swap_values;
                    for (Size k = 0; k < rate.size(); k++)
                    {
                        boost::shared_ptr<VanillaSwap> swap =
                            makeSwap(lengths_[i], rate[k], spreads_[j]);
                        swap_values.push_back(swap->NPV());

                        Y += swap->NPV();
                    }
                    // and check that they go the right way
                    std::vector<Real>::iterator it =
                        std::adjacent_find(swap_values.begin(), swap_values.end(),
                        std::less<Real>());
                    if (it != swap_values.end())
                    {
                        Size n = it - swap_values.begin();
                        BOOST_ERROR("NPV is increasing with the fixed rate in a swap: \n"
                                    << "    length: " << lengths_[i] << " years\n"
                                    << "    value:  " << swap_values[n]
                                    << " paying fixed rate: " << io::rate(rate[n]) << "\n"
                                    << "    value:  " << swap_values[n + 1]
                                    << " paying fixed rate: " << io::rate(rate[n + 1]));
                    }
                }
            }

            return Y;
        }

        Real RateVariationNPV(Rate rate)
        {
            boost::shared_ptr<VanillaSwap> swap = makeSwap(10, rate, 0.01);
            return swap->NPV();
        }

        Real SpreadDependencyNPV(std::vector<Real> spread)
        {
            Real Y = 0;

            for (Size i = 0; i < lengths_.size(); i++)
            {
                for (Size j = 0; j < rates_.size(); j++)
                {
                    // store the results for different spreads...
                    std::vector<Real> swap_values;
                    for (Size k = 0; k < spread.size(); k++)
                    {
                        boost::shared_ptr<VanillaSwap> swap =
                            makeSwap(lengths_[i], rates_[j], spread[k]);
                        swap_values.push_back(swap->NPV());
                        Y += swap->NPV();
                    }
                    // and check that they go the right way
                    std::vector<Real>::iterator it =
                        std::adjacent_find(swap_values.begin(), swap_values.end(),
                        std::greater<Real>());
                    if (it != swap_values.end())
                    {
                        Size n = it - swap_values.begin();
                        BOOST_ERROR("NPV is decreasing with the floating spread in a swap: \n"
                                    << "    length: " << lengths_[i] << " years\n"
                                    << "    value:  " << swap_values[n]
                                    << " receiving spread: " << io::rate(spread[n]) << "\n"
                                    << "    value:  " << swap_values[n + 1]
                                    << " receiving spread: " << io::rate(spread[n + 1]));
                    }
                }
            }

            return Y;
        }

        Real InArrearsNPV(std::vector<Real> capletVolatility)
        {
            Date maturity = today_ + 5 * Years;
            Calendar calendar = NullCalendar();
            Schedule schedule(today_, maturity, Period(Annual), calendar,
                                Following, Following,
                                DateGeneration::Forward, false);
            DayCounter dayCounter = SimpleDayCounter();
            std::vector<Real> nominals(1, 100000000.0);
            boost::shared_ptr<IborIndex> index_(new IborIndex("dummy", 1 * Years, 0,
                EURCurrency(), calendar,
                Following, false, dayCounter,
                termStructure_));
            Rate oneYear = 0.05;
            Rate r = std::log(1.0 + oneYear);
            termStructure_.linkTo(flatRate(today_, r, dayCounter));


            std::vector<Rate> coupons(1, oneYear);
            Leg fixedLeg = FixedRateLeg(schedule)
                .withNotionals(nominals)
                .withCouponRates(coupons, dayCounter);

            std::vector<Real> gearings;
            std::vector<Rate> spreads;
            Natural fixingDays = 0;

            Real Y = 0;
            for (Size i = 0; i < capletVolatility.size(); i++)
            {
                Handle<OptionletVolatilityStructure> vol(
                    boost::shared_ptr<OptionletVolatilityStructure>(new
                    ConstantOptionletVolatility(today_, NullCalendar(), Following,
                    capletVolatility[i], dayCounter)));

                boost::shared_ptr<IborCouponPricer> pricer(new
                                                           BlackIborCouponPricer(vol));

                Leg floatingLeg = IborLeg(schedule, index_)
                    .withNotionals(nominals)
                    .withPaymentDayCounter(dayCounter)
                    .withFixingDays(fixingDays)
                    .withGearings(gearings)
                    .withSpreads(spreads)
                    .inArrears();
                setCouponPricer(floatingLeg, pricer);

                Swap swap(floatingLeg, fixedLeg);
                swap.setPricingEngine(boost::shared_ptr<PricingEngine>(
                    new DiscountingSwapEngine(termStructure_)));

                Y += swap.NPV();
            }
            return Y;
        }
    };

    struct Variation
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "Independent Variable", "SwapNPV"
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, Variation& v)
        {
                stm << v.param_
                    << ";" << v.swapNpv_
                    << std::endl;
                return stm;
            }

        Real param_;
        Real swapNpv_;
    };
}


bool AdjointSwapTest::testFairRate()
{
    BOOST_MESSAGE("Testing vanilla-swap calculation of fair rate...");

    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVars vars;

    std::vector<Spread> spreads;

    Size maxSize = 100;
    Rate startSpread = -0.01;
    Rate maxSpread = 0.01;
    Real step = (maxSpread  - startSpread) / maxSize;

    for (Size i = 0; i < maxSize; i++)
        spreads.push_back(startSpread + i * step);

    std::vector<PerformanceTime> performResults;
    std::vector<AdjointTime> adjointResults;
    std::vector<TapeSize> tapeSizeResults;

    boost::timer timer;
    Size sizeof_indep;

    for (Size s = 0; s < spreads.size(); s += 4)
    {
        PerformanceTime perfomTime;
        AdjointTime adTime;
        TapeSize tSize;

        timer.restart();
        sizeof_indep = s + 1;

        perfomTime.indepVarNumber_ = sizeof_indep;
        adTime.indepVarNumber_ = sizeof_indep;
        tSize.indepVarNumber_ = sizeof_indep;

        std::vector<cl::TapeDouble> spread(sizeof_indep);
        for (Size i = 0; i < sizeof_indep; i++)
        {
            spread[i] = spreads[i];
        }

        Independent(spread);

        std::vector<cl::TapeDouble> Y(1);

        Y[0] = vars.FairRateNPV(spread, spreads);

        cl::TapeFunction<double> f(spread, Y);

        perfomTime.timeTapeRecording_ = timer.elapsed();
        tSize.memory_ = f.Memory();

        std::vector<double> sf_Forward(sizeof_indep), sf_Reverse(sizeof_indep);

        //Start differentiation in Forward mode
        double tf = gradForward(f, sf_Forward, false, false);

        //Start differentiation in Reverse mode
        double tr =  gradReverse(f, sf_Reverse, false, false);

        perfomTime.timeAdjoint_ = std::min(tf, tr);
        adTime.timeAdjoint_ = perfomTime.timeAdjoint_;

        //Finite differences
        timer.restart();
        double h = 1.0e-14;
        std::vector<Real> sf_Finite(sizeof_indep);

        Real YF = vars.FairRateNPV(spread, spreads);

        for (Size i = 0; i < sizeof_indep; i++)
        {
            spread[i] += h;
            sf_Finite[i] = (vars.FairRateNPV(spread, spreads) - YF) / h;
            spread[i] -= h;
        }

        perfomTime.timeAnalytical_ = timer.elapsed();

        performResults.push_back(perfomTime);
        adjointResults.push_back(adTime);
        tapeSizeResults.push_back(tSize);

        result = checkWithFiniteDiff(sf_Forward, sf_Reverse, sf_Finite, 1e-2, 1e-10);
    }

    cl::AdjointTestOutput outAdjointPerformance("Swap//FairRate", { { "filename", "AdjointPerformance" }, { "not_clear", "Not" }, { "title", "Fair Rate" }, { "ylabel", "Time" } });
    outAdjointPerformance << performResults;

    cl::AdjointTestOutput outAdjoint("Swap//FairRate", { { "filename", "Adjoint" }, { "not_clear", "Not" }, { "title", "Fair Rate" }, { "ylabel", "Time" } });
    outAdjoint << adjointResults;

    cl::AdjointTestOutput outTapeSize("Swap//FairRate", { { "filename", "TapeSize" }, { "not_clear", "Not" }, { "title", "Fair Rate" }, { "ylabel", "Size" } });
    outTapeSize << tapeSizeResults;

    std::vector<Real> spread;
    for (Size i = 0; i < spreads.size(); i++)
    {
        spread.push_back(spreads[i]);
    }

    Real stepv = step / maxSize;
    std::vector<Variation> variationResults;
    for (Size i = 0; i < maxSize; i++)
    {
        spread[0] += i * stepv;

        Variation var;
        var.param_ = spread[0];
        var.swapNpv_ = vars.FairRateNPV(spread, spread);
        variationResults.push_back(var);

        spread[0] -= i * stepv;
    }

    cl::AdjointTestOutput outVariation("Swap//FairRate//output", { { "filename", "SwapNPVonSpreads" }, { "ylabel", "SwapNPV" }, { "not_clear", "Not" }, { "title", "SwapNPV on Spreads" } });
    outVariation << variationResults;

#endif
    return result;
}

bool AdjointSwapTest::testFairSpread()
{
    BOOST_MESSAGE("Testing vanilla-swap calculation of fair spread...");

    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVars vars;

    std::vector<Rate> rates;

    Size maxSize = 100;
    Rate startRate = 0.01;
    Rate maxRate = 0.06;
    Real step = (maxRate - startRate) / maxSize;

    for (Size i = 0; i < maxSize; i++)
        rates.push_back(startRate + i*step);

    std::vector<PerformanceTime> performResults;
    std::vector<AdjointTime> adjointResults;
    std::vector<TapeSize> tapeSizeResults;

    boost::timer timer;
    Size sizeof_indep;

    for (Size s = 0; s < rates.size() ; s += 4)
    {
        PerformanceTime perfomTime;
        AdjointTime adTime;
        TapeSize tSize;

        timer.restart();
        sizeof_indep = s + 1;

        perfomTime.indepVarNumber_ = sizeof_indep;
        adTime.indepVarNumber_ = sizeof_indep;
        tSize.indepVarNumber_ = sizeof_indep;

        std::vector<cl::TapeDouble> rate(sizeof_indep);
        for (Size i = 0; i < sizeof_indep; i++)
            rate[i] = rates[i];

        Independent(rate);

        std::vector<cl::TapeDouble> Y(1);
        Y[0] = vars.FairSpreadNPV(rate, rates);

        cl::TapeFunction<double> f(rate, Y);

        perfomTime.timeTapeRecording_ = timer.elapsed();
        tSize.memory_ = f.Memory();

        std::vector<double> sf_Forward(sizeof_indep), sf_Reverse(sizeof_indep);

        //Start differentiation in Forward mode
        double tf = gradForward(f, sf_Forward, false, false);

        //Start differentiation in Reverse mode
        double tr = gradReverse(f, sf_Reverse, false, false);

        perfomTime.timeAdjoint_ = std::min(tf, tr);
        adTime.timeAdjoint_ = perfomTime.timeAdjoint_;

        //Finite differences
        timer.restart();
        double h = 1.0e-14;
        std::vector<Real> sf_Finite(sizeof_indep);

        Real YF = vars.FairSpreadNPV(rate, rates);

        for (Size i = 0; i < sizeof_indep; i++)
        {
            rate[i] += h;
            sf_Finite[i] = (vars.FairSpreadNPV(rate, rates) - YF) / h;
            rate[i] -= h;
        }

        perfomTime.timeAnalytical_ = timer.elapsed();

        performResults.push_back(perfomTime);
        adjointResults.push_back(adTime);
        tapeSizeResults.push_back(tSize);

        result = checkWithFiniteDiff(sf_Forward, sf_Reverse, sf_Finite, 1e-2, 1e-10);
    }

    cl::AdjointTestOutput outAdjointPerformance("Swap//FairSpread", { { "filename", "AdjointPerformance" }, { "not_clear", "Not" }, { "title", "Fair Spread" }, { "ylabel", "Time" } });
    outAdjointPerformance << performResults;

    cl::AdjointTestOutput outAdjoint("Swap//FairSpread", { { "filename", "Adjoint" }, { "not_clear", "Not" }, { "title", "Fair Spread" }, { "ylabel", "Time" } });
    outAdjoint << adjointResults;

    cl::AdjointTestOutput outTapeSize("Swap//FairSpread", { { "filename", "TapeSize" }, { "not_clear", "Not" }, { "title", "Fair Spread" }, { "ylabel", "Size" } });
    outTapeSize << tapeSizeResults;

    std::vector<Real> rate;
    for (Size i = 0; i < rates.size(); i++)
    {
        rate.push_back(rates[i]);
    }

    Real stepv = step / maxSize;
    std::vector<Variation> variationResults;
    for (Size i = 0; i < maxSize; i++)
    {
        rate[0] += i * stepv;

        Variation var;
        var.param_ = rate[0];
        var.swapNpv_ = vars.FairSpreadNPV(rate, rate);
        variationResults.push_back(var);

        rate[0] -= i * stepv;
    }

    cl::AdjointTestOutput outVariation("Swap//FairSpread//output", { { "filename", "SwapNPVonRates" }, { "ylabel", "SwapNPV" }, { "not_clear", "Not" }, { "title", "SwapNPV on Rates" } });
    outVariation << variationResults;
#endif
    return result;
}

bool AdjointSwapTest::testRateDependency()
{
    BOOST_MESSAGE("Testing vanilla-swap dependency on rate...");

    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVars vars;

    std::vector<Rate> rates;

    Size maxSize = 10;
    Rate startRate = 0.01;
    Rate maxRate = 0.06;
    Real step = (maxRate - startRate) / maxSize;

    for (Size i = 0; i < maxSize; i++)
        rates.push_back(startRate + i*step);

    std::vector<PerformanceTime> performResults;
    std::vector<AdjointTime> adjointResults;
    std::vector<TapeSize> tapeSizeResults;

    boost::timer timer;
    Size sizeof_indep;

    for (Size s = 0; s < rates.size(); s++)
    {
        PerformanceTime perfomTime;
        AdjointTime adTime;
        TapeSize tSize;

        timer.restart();
        sizeof_indep = s + 1;

        perfomTime.indepVarNumber_ = sizeof_indep;
        adTime.indepVarNumber_ = sizeof_indep;
        tSize.indepVarNumber_ = sizeof_indep;

        std::vector<cl::TapeDouble> rate(sizeof_indep);
        for (Size i = 0; i < sizeof_indep; i++)
            rate[i] = rates[i];

        Independent(rate);

        std::vector<cl::TapeDouble> Y(1);
        Y[0] = vars.RateDependencyNPV(rate);

        cl::TapeFunction<double> f(rate, Y);

        perfomTime.timeTapeRecording_ = timer.elapsed();
        tSize.memory_ = f.Memory();

        std::vector<double> sf_Forward(sizeof_indep), sf_Reverse(sizeof_indep);

        //Start differentiation in Forward mode
        double tf = gradForward(f, sf_Forward, false, false);

        //Start differentiation in Reverse mode
        double tr = gradReverse(f, sf_Reverse, false, false);

        perfomTime.timeAdjoint_ = std::min(tf, tr);
        adTime.timeAdjoint_ = perfomTime.timeAdjoint_;

        //Finite differences
        timer.restart();
        double h = 1.0e-14;
        std::vector<Real> sf_Finite(sizeof_indep);

        Real YF = vars.RateDependencyNPV(rate);

        for (Size i = 0; i < sizeof_indep; i++)
        {
            rate[i] += h;
            sf_Finite[i] = (vars.RateDependencyNPV(rate) - YF) / h;
            rate[i] -= h;
        }

        perfomTime.timeAnalytical_ = timer.elapsed();

        performResults.push_back(perfomTime);
        adjointResults.push_back(adTime);
        tapeSizeResults.push_back(tSize);

        result = checkWithFiniteDiff(sf_Forward, sf_Reverse, sf_Finite, 1e-2, 1e-10);
    }

    cl::AdjointTestOutput outAdjointPerformance("Swap//RateDependency", { { "filename", "AdjointPerformance" }, { "not_clear", "Not" }, { "title", "Rate Dependency" }, { "ylabel", "Time" } });
    outAdjointPerformance << performResults;

    cl::AdjointTestOutput outAdjoint("Swap//RateDependency", { { "filename", "Adjoint" }, { "not_clear", "Not" }, { "title", "Rate Dependency" }, { "ylabel", "Time" } });
    outAdjoint << adjointResults;

    cl::AdjointTestOutput outTapeSize("Swap//RateDependency", { { "filename", "TapeSize" }, { "not_clear", "Not" }, { "title", "Rate Dependency" }, { "ylabel", "Size" } });
    outTapeSize << tapeSizeResults;


    std::vector<Real> rate;
    for (Size i = 0; i < rates.size(); i++)
    {
        rate.push_back(rates[i]);
    }

    Real stepv = step / maxSize;
    std::vector<Variation> variationResults;
    for (Size i = 0; i < maxSize; i++)
    {
        rate[0] += i * stepv;

        Variation var;
        var.param_ = rate[0];
        var.swapNpv_ = vars.RateDependencyNPV(rate);
        variationResults.push_back(var);

        rate[0] -= i * stepv;
    }

    cl::AdjointTestOutput outVariation("Swap//RateDependency//output", { { "filename", "SwapNPVonRates" }, { "ylabel", "SwapNPV" }, { "not_clear", "Not" }, { "title", "SwapNPV on Rates" } });
    outVariation << variationResults;
#endif
    return result;
}

bool AdjointSwapTest::testSpreadDependency()
{
    BOOST_MESSAGE("Testing vanilla-swap dependency on floating spread...");

    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVars vars;

    std::vector<Spread> spreads;

    Size maxSize = 20;
    Rate startSpread = -0.01;
    Rate maxSpread = 0.01;
    Real step = (maxSpread - startSpread) / maxSize;

    for (Size i = 0; i < maxSize; i++)
        spreads.push_back(startSpread + i * step);

    std::vector<PerformanceTime> performResults;
    std::vector<AdjointTime> adjointResults;
    std::vector<TapeSize> tapeSizeResults;

    boost::timer timer;
    Size sizeof_indep;

    for (Size s = 0; s < spreads.size(); s += 2)
    {
        PerformanceTime perfomTime;
        AdjointTime adTime;
        TapeSize tSize;

        timer.restart();
        sizeof_indep = s + 1;

        perfomTime.indepVarNumber_ = sizeof_indep;
        adTime.indepVarNumber_ = sizeof_indep;
        tSize.indepVarNumber_ = sizeof_indep;

        std::vector<cl::TapeDouble> spread(sizeof_indep);
        for (Size i = 0; i < sizeof_indep; i++)
        {
            spread[i] = spreads[i];
        }

        Independent(spread);

        std::vector<cl::TapeDouble> Y(1);

        Y[0] = vars.SpreadDependencyNPV(spread);

        cl::TapeFunction<double> f(spread, Y);

        perfomTime.timeTapeRecording_ = timer.elapsed();
        tSize.memory_ = f.Memory();

        std::vector<double> sf_Forward(sizeof_indep), sf_Reverse(sizeof_indep);

        //Start differentiation in Forward mode
        double tf = gradForward(f, sf_Forward, false, false);

        //Start differentiation in Reverse mode
        double tr = gradReverse(f, sf_Reverse, false, false);

        perfomTime.timeAdjoint_ = std::min(tf, tr);
        adTime.timeAdjoint_ = perfomTime.timeAdjoint_;

        //Finite differences
        timer.restart();
        double h = 1.0e-14;
        std::vector<Real> sf_Finite(sizeof_indep);

        Real YF = vars.SpreadDependencyNPV(spread);

        for (Size i = 0; i < sizeof_indep; i++)
        {
            spread[i] += h;
            sf_Finite[i] = (vars.SpreadDependencyNPV(spread) - YF) / h;
            spread[i] -= h;
        }

        perfomTime.timeAnalytical_ = timer.elapsed();

        performResults.push_back(perfomTime);
        adjointResults.push_back(adTime);
        tapeSizeResults.push_back(tSize);

        result = checkWithFiniteDiff(sf_Forward, sf_Reverse, sf_Finite, 1e-2, 1e-10);
    }

    cl::AdjointTestOutput outAdjointPerformance("Swap//SpreadDependency", { { "filename", "AdjointPerformance" }, { "not_clear", "Not" }, { "title", "Spread Dependency" }, { "ylabel", "Time" } });
    outAdjointPerformance << performResults;

    cl::AdjointTestOutput outAdjoint("Swap//SpreadDependency", { { "filename", "Adjoint" }, { "not_clear", "Not" }, { "title", "Spread Dependency" }, { "ylabel", "Time" } });
    outAdjoint << adjointResults;

    cl::AdjointTestOutput outTapeSize("Swap//SpreadDependency", { { "filename", "TapeSize" }, { "not_clear", "Not" }, { "title", "Spread Dependency" }, { "ylabel", "Size" } });
    outTapeSize << tapeSizeResults;

    std::vector<Real> spread;
    for (Size i = 0; i < spreads.size(); i++)
    {
        spread.push_back(spreads[i]);
    }

    Real stepv = step / maxSize;
    std::vector<Variation> variationResults;
    for (Size i = 0; i < maxSize; i++)
    {
        spread[0] += i * stepv;

        Variation var;
        var.param_ = spread[0];
        var.swapNpv_ = vars.SpreadDependencyNPV(spread);
        variationResults.push_back(var);

        spread[0] -= i * stepv;
    }

    cl::AdjointTestOutput outVariation("Swap//SpreadDependency//output", { { "filename", "SwapNPVonSpreads" }, { "ylabel", "SwapNPV" }, { "not_clear", "Not" }, { "title", "SwapNPV on Spreads" } });
    outVariation << variationResults;
#endif
    return result;
}

bool AdjointSwapTest::testInArrears()
{
    BOOST_MESSAGE("Testing in-arrears swap calculation...");

    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVars vars;

    std::vector<Volatility> capletVolatilities;

    Size maxSize = 201;
    Rate startVol = 0.18;
    Rate maxVol = 0.22;
    Real step = (maxVol - startVol) / maxSize;

    for (Size i = 0; i < maxSize; i++)
        capletVolatilities.push_back(startVol + i * step);

    std::vector<PerformanceTime> performResults;
    std::vector<AdjointTime> adjointResults;
    std::vector<TapeSize> tapeSizeResults;

    boost::timer timer;
    Size sizeof_indep;
    for (Size s = 0; s < capletVolatilities.size(); s += 20)
    {
        PerformanceTime perfomTime;
        AdjointTime adTime;
        TapeSize tSize;

        timer.restart();
        sizeof_indep = s + 1;

        perfomTime.indepVarNumber_ = sizeof_indep;
        adTime.indepVarNumber_ = sizeof_indep;
        tSize.indepVarNumber_ = sizeof_indep;

        std::vector<cl::TapeDouble> capletVol(sizeof_indep);
        for (Size i = 0; i < sizeof_indep; i++)
        {
            capletVol[i] = capletVolatilities[i];
        }

        Independent(capletVol);

        std::vector<cl::TapeDouble> Y(1);

        Y[0] = vars.InArrearsNPV(capletVol);

        cl::TapeFunction<double> f(capletVol, Y);

        perfomTime.timeTapeRecording_ = timer.elapsed();
        tSize.memory_ = f.Memory();

        std::vector<double> sf_Forward(sizeof_indep), sf_Reverse(sizeof_indep);

        //Start differentiation in Forward mode
        double tf = gradForward(f, sf_Forward, false, false);

        //Start differentiation in Reverse mode
        double tr = gradReverse(f, sf_Reverse, false, false);

        perfomTime.timeAdjoint_ = std::min(tf, tr);
        adTime.timeAdjoint_ = perfomTime.timeAdjoint_;

        //Finite differences
        timer.restart();
        double h = 1.0e-12;
        std::vector<Real> sf_Finite(sizeof_indep);

        Real YF = vars.InArrearsNPV(capletVol);

        for (Size i = 0; i < sizeof_indep; i++)
        {
            capletVol[i] += h;
            sf_Finite[i] = (vars.InArrearsNPV(capletVol) - YF) / h;
            capletVol[i] -= h;
        }

        perfomTime.timeAnalytical_ = timer.elapsed();

        performResults.push_back(perfomTime);
        adjointResults.push_back(adTime);
        tapeSizeResults.push_back(tSize);

        result = checkWithFiniteDiff(sf_Forward, sf_Reverse, sf_Finite, 1e-2, 1e-10);
    }

    cl::AdjointTestOutput outAdjointPerformance("Swap//InArrears", { { "filename", "AdjointPerformance" }, { "not_clear", "Not" }, { "title", "Cached Value" }, { "ylabel", "Time" } });
    outAdjointPerformance << performResults;

    cl::AdjointTestOutput outAdjoint("Swap//InArrears", { { "filename", "Adjoint" }, { "not_clear", "Not" }, { "title", "Cached Value" }, { "ylabel", "Time" } });
    outAdjoint << adjointResults;

    cl::AdjointTestOutput outTapeSize("Swap//InArrears", { { "filename", "TapeSize" }, { "not_clear", "Not" }, { "title", "Cached Value" }, { "ylabel", "Size" } });
    outTapeSize << tapeSizeResults;


    std::vector<Real> volatility;
    for (Size i = 0; i < capletVolatilities.size(); i++)
    {
        volatility.push_back(capletVolatilities[i]);
    }

    Real stepv = step / 10;
    std::vector<Variation> variationResults;
    for (Size i = 0; i < 10; i++)
    {
        volatility[0] += i * stepv;

        Variation var;
        var.param_ = volatility[0];
        var.swapNpv_ = vars.InArrearsNPV(volatility);
        variationResults.push_back(var);

        volatility[0] -= i * stepv;
    }

    cl::AdjointTestOutput outVariation("Swap//InArrears//output", { { "filename", "SwapNPVonVol" }, { "ylabel", "SwapNPV" }, { "not_clear", "Not" }, { "title", "SwapNPV on Volatilities" } });
    outVariation << variationResults;

#endif
    return result;
}

test_suite* AdjointSwapTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Adjoint swap tests");
    suite->add(QUANTLIB_TEST_CASE(&AdjointSwapTest::testFairRate));
    suite->add(QUANTLIB_TEST_CASE(&AdjointSwapTest::testFairSpread));
    suite->add(QUANTLIB_TEST_CASE(&AdjointSwapTest::testRateDependency));
    suite->add(QUANTLIB_TEST_CASE(&AdjointSwapTest::testSpreadDependency));
    suite->add(QUANTLIB_TEST_CASE(&AdjointSwapTest::testInArrears));
    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_swap)

BOOST_AUTO_TEST_CASE(testSwapCalcOfFairRate)
{
    BOOST_CHECK(AdjointSwapTest::testFairRate());
}

BOOST_AUTO_TEST_CASE(testSwapCalcOfFairSpread)
{
    BOOST_CHECK(AdjointSwapTest::testFairSpread());
}

BOOST_AUTO_TEST_CASE(testSwapDependencyOnRate)
{
    BOOST_CHECK(AdjointSwapTest::testRateDependency());
}

BOOST_AUTO_TEST_CASE(testSwapDependencyOnSpread)
{
    BOOST_CHECK(AdjointSwapTest::testSpreadDependency());
}

BOOST_AUTO_TEST_CASE(testInArrearsSwapCalc)
{
    BOOST_CHECK(AdjointSwapTest::testInArrears());
}

BOOST_AUTO_TEST_SUITE_END()

#endif


