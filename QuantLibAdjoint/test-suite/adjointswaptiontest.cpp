/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2003 RiskMap srl
Copyright (C) 2006, 2007 Ferdinando Ametrano
Copyright (C) 2006 Marco Bianchetti
Copyright (C) 2006 Cristina Duminuco
Copyright (C) 2007, 2008 StatPro Italia srl
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

#include <ql/quantlib.hpp>
#include "utilities.hpp"
#include "adjointswaptiontest.hpp"
#include <test-suite/adjointtestutilities.hpp>

// based on swaption.cpp file from test-suite
using namespace QuantLib;
using namespace boost::unit_test_framework;

namespace
{
    struct CommonVars
    {
        // Global data.
        Date today_, settlement_;
        Real nominal_;
        Calendar calendar_;

        BusinessDayConvention fixedConvention_;
        Frequency fixedFrequency_;
        DayCounter fixedDayCount_;

        BusinessDayConvention floatingConvention_;
        Period floatingTenor_;
        boost::shared_ptr<IborIndex> index_;

        Natural settlementDays_;
        RelinkableHandle<YieldTermStructure> termStructure_;

        std::vector<Period> exercises_;
        std::vector<Period> lengths_;
        VanillaSwap::Type type_;

        // Cleanup
        SavedSettings backup_;

        // Make swaption with given parameters.
        boost::shared_ptr<Swaption> makeSwaption(
            const boost::shared_ptr<VanillaSwap>& swap,
            const Date& exercise,
            Volatility volatility,
            Settlement::Type settlementType = Settlement::Physical)
        {
            // Set volatility from simple quote.
            Handle<Quote> vol(boost::shared_ptr<Quote>(
                new SimpleQuote(volatility)));

            // Create pricing engine.
            boost::shared_ptr<PricingEngine> engine(
                new BlackSwaptionEngine(termStructure_, vol));

            // Create swaption using vanilla swap.
            boost::shared_ptr<Swaption> result(new
                                               Swaption(swap,
                                               boost::shared_ptr<Exercise>(
                                               new EuropeanExercise(exercise)),
                                               settlementType));

            // Set pricing engine in swaption.
            result->setPricingEngine(engine);
            return result;
        }

        // Constructor
        CommonVars()
        {
            settlementDays_ = 2;
            nominal_ = 1000000.0;
            fixedConvention_ = Unadjusted;
            fixedFrequency_ = Annual;
            fixedDayCount_ = Thirty360();

            index_ = boost::shared_ptr<IborIndex>(new Euribor6M(termStructure_));
            floatingConvention_ = index_->businessDayConvention();
            floatingTenor_ = index_->tenor();
            calendar_ = index_->fixingCalendar();
            today_ = calendar_.adjust(Date::todaysDate());
            Settings::instance().evaluationDate() = today_;
            settlement_ = calendar_.advance(today_, settlementDays_, Days);
            termStructure_.linkTo(flatRate(settlement_, 0.05, Actual365Fixed()));

            exercises_ = { 1 * Years, 3 * Years, 7 * Years, 10 * Years };
            lengths_ = {1 * Years, 2*Years, 5 * Years, 10 * Years, 20 * Years };
            type_ = VanillaSwap::Payer;
        };

        // Calculate sum of swaption NPV for given strike rates.
        std::vector<Real> StrikeDependency(std::vector<Real> strikes)
        {
            std::vector<Real> Y(2, 0);

            for (Size i = 0; i < exercises_.size(); i++)
            {
                for (Size j = 0; j < lengths_.size(); j++)
                {
                    // Set date of exercize.
                    Date exerciseDate = calendar_.advance(today_, exercises_[i]);

                    // Set start date.
                    Date startDate = calendar_.advance(exerciseDate, settlementDays_, Days);

                    // Store the results for different rates.
                    std::vector<Real> values;
                    std::vector<Real> values_cash;
                    Volatility vol = 0.20;
                    for (Size l = 0; l < strikes.size(); l++)
                    {
                        // Create vanilla swap.
                        boost::shared_ptr<VanillaSwap> swap =
                            MakeVanillaSwap(lengths_[j], index_, strikes[l])
                            .withEffectiveDate(startDate)
                            .withFixedLegTenor(1 * Years)
                            .withFixedLegDayCount(fixedDayCount_)
                            .withFloatingLegSpread(0.0)
                            .withType(type_);

                        // Create swaption using vanilla swap.
                        boost::shared_ptr<Swaption> swaption =
                            makeSwaption(swap, exerciseDate, vol);

                        values.push_back(swaption->NPV());
                        Y[0] += swaption->NPV();

                        // Create swaption with cash settlement.
                        boost::shared_ptr<Swaption> swaption_cash =
                            makeSwaption(swap, exerciseDate, vol,
                            Settlement::Cash);

                        values_cash.push_back(swaption_cash->NPV());
                        Y[1] += swaption_cash->NPV();

                    }

                    // Check that they go the right way.
                    std::vector<Real>::iterator it =
                        std::adjacent_find(values.begin(), values.end(),
                        std::less<Real>());
                    if (it != values.end())
                    {
                        Size n = it - values.begin();
                        BOOST_ERROR("NPV of Payer swaption with delivery settlement_"
                                    "is increasing with the strike:" <<
                                    "\noption tenor: " << exercises_[i] <<
                                    "\noption date:  " << exerciseDate <<
                                    "\nvolatility:   " << io::rate(vol) <<
                                    "\nswap tenor:   " << lengths_[j] <<
                                    "\nvalue:        " << values[n] << " at strike: " << io::rate(strikes[n]) <<
                                    "\nvalue:        " << values[n + 1] << " at strike: " << io::rate(strikes[n + 1]));
                    }
                    std::vector<Real>::iterator it_cash =
                        std::adjacent_find(values_cash.begin(), values_cash.end(),
                        std::less<Real>());
                    if (it_cash != values_cash.end())
                    {
                        Size n = it_cash - values_cash.begin();
                        BOOST_ERROR("NPV of Payer swaption with cash settlement_"
                                    "is increasing with the strike:" <<
                                    "\noption tenor: " << exercises_[i] <<
                                    "\noption date:  " << exerciseDate <<
                                    "\nvolatility:   " << io::rate(vol) <<
                                    "\nswap tenor:   " << lengths_[j] <<
                                    "\nvalue:        " << values_cash[n] << " at strike: " << io::rate(strikes[n]) <<
                                    "\nvalue:        " << values_cash[n + 1] << " at strike: " << io::rate(strikes[n + 1]));
                    }

                }

            }
            return Y;
        }

        // Calculate sum of swaption NPV for given spreads.
        std::vector<Real> SpreadDependency(std::vector<Real> spreads)
        {
            std::vector<Real> Y(2, 0);

            for (Size i = 0; i < exercises_.size(); i++)
            {
                for (Size j = 0; j < lengths_.size(); j++)
                {
                    // Set date of exercize.
                    Date exerciseDate = calendar_.advance(today_, exercises_[i]);

                    // Set start date.
                    Date startDate = calendar_.advance(exerciseDate, settlementDays_, Days);

                    // Store the results for different rates.
                    std::vector<Real> values;
                    std::vector<Real> values_cash;
                    for (Size l = 0; l< spreads.size(); l++)
                    {
                        // Create vanilla swap.
                        boost::shared_ptr<VanillaSwap> swap =
                            MakeVanillaSwap(lengths_[j], index_, 0.06)
                            .withFixedLegTenor(1 * Years)
                            .withFixedLegDayCount(fixedDayCount_)
                            .withEffectiveDate(startDate)
                            .withFloatingLegSpread(spreads[l])
                            .withType(type_);

                        // Create swaption using vanilla swap.
                        boost::shared_ptr<Swaption> swaption =
                            makeSwaption(swap, exerciseDate, 0.20);

                        values.push_back(swaption->NPV());
                        Y[0] += swaption->NPV();

                        // Create swaption with cash settlement.
                        boost::shared_ptr<Swaption> swaption_cash =
                            makeSwaption(swap, exerciseDate, 0.20,
                            Settlement::Cash);

                        values_cash.push_back(swaption_cash->NPV());
                        Y[1] += swaption_cash->NPV();
                    }
                    // Check that they go the right way.

                    std::vector<Real>::iterator it =
                        std::adjacent_find(values.begin(), values.end(),
                        std::greater<Real>());
                    if (it != values.end())
                    {
                        Size n = it - values.begin();
                        BOOST_ERROR("NPV is decreasing with the spread " <<
                                    "in a payer swaption (physical delivered):" <<
                                    "\nexercise date: " << exerciseDate <<
                                    "\nlength:        " << lengths_[j] <<
                                    "\nvalue:         " << values[n] << " for spread: " << io::rate(spreads[n]) <<
                                    "\nvalue:         " << values[n + 1] << " for spread: " << io::rate(spreads[n + 1]));
                    }
                    std::vector<Real>::iterator it_cash =
                        std::adjacent_find(values_cash.begin(), values_cash.end(),
                        std::greater<Real>());
                    if (it_cash != values_cash.end())
                    {
                        Size n = it_cash - values_cash.begin();
                        BOOST_ERROR("NPV is decreasing with the spread " <<
                                    "in a payer swaption (cash delivered):" <<
                                    "\nexercise date: " << exerciseDate <<
                                    "\nlength: " << lengths_[j] <<
                                    "\nvalue:  " << values_cash[n] << " for spread: " << io::rate(spreads[n]) <<
                                    "\nvalue:  " << values_cash[n + 1] << " for spread: " << io::rate(spreads[n + 1]));
                    }
                }

            }
            return Y;
        }

        // Calculate sum of swaption NPV for given spreads.
        std::vector<Real> SpreadTreatment(std::vector<Real> spreads)
        {
            std::vector<Real> Y(4, 0);

            for (Size i = 0; i < exercises_.size(); i++)
            {
                for (Size j = 0; j < lengths_.size(); j++)
                {
                    // Set date of exercize.
                    Date exerciseDate = calendar_.advance(today_, exercises_[i]);

                    // Set start date.
                    Date startDate = calendar_.advance(exerciseDate, settlementDays_, Days);
                    for (Size l = 0; l < spreads.size(); l++)
                    {
                        // Create vanilla swap.
                        boost::shared_ptr<VanillaSwap> swap =
                            MakeVanillaSwap(lengths_[j], index_, 0.06)
                            .withFixedLegTenor(1 * Years)
                            .withFixedLegDayCount(fixedDayCount_)
                            .withEffectiveDate(startDate)
                            .withFloatingLegSpread(spreads[l])
                            .withType(type_);

                        // Calculate spread correction.
                        Spread correction = spreads[l] *
                            swap->floatingLegBPS() /
                            swap->fixedLegBPS();

                        // Create equivalent vanilla swap.
                        boost::shared_ptr<VanillaSwap> equivalentSwap =
                            MakeVanillaSwap(lengths_[j], index_, 0.06 + correction)
                            .withFixedLegTenor(1 * Years)
                            .withFixedLegDayCount(fixedDayCount_)
                            .withEffectiveDate(startDate)
                            .withFloatingLegSpread(0.0)
                            .withType(type_);

                        // Create swaption using vanilla swap.
                        boost::shared_ptr<Swaption> swaption1 =
                            makeSwaption(swap, exerciseDate, 0.20);

                        Y[0] += swaption1->NPV();

                        // Create swaption using equivalent vanilla swap.
                        boost::shared_ptr<Swaption> swaption2 =
                            makeSwaption(equivalentSwap, exerciseDate, 0.20);

                        Y[1] += swaption2->NPV();

                        // Create swaption with cash settlement using vanilla swap.
                        boost::shared_ptr<Swaption> swaption1_cash =
                            makeSwaption(swap, exerciseDate, 0.20,
                            Settlement::Cash);

                        Y[2] += swaption1_cash->NPV();

                        // Create swaption with cash settlement using equivalent vanilla swap.
                        boost::shared_ptr<Swaption> swaption2_cash =
                            makeSwaption(equivalentSwap, exerciseDate, 0.20,
                            Settlement::Cash);

                        Y[3] += swaption2_cash->NPV();

                        // Check swaption NPV and equivalent swaption NPV
                        if (std::fabs(swaption1->NPV() - swaption2->NPV()) > 1.0e-6)
                            BOOST_ERROR("wrong spread treatment:" <<
                            "\nexercise: " << exerciseDate <<
                            "\nlength:   " << lengths_[j] <<
                            "\ntype      " << type_ <<
                            "\nspread:   " << io::rate(spreads[l]) <<
                            "\noriginal swaption value:   " << swaption1->NPV() <<
                            "\nequivalent swaption value: " << swaption2->NPV());

                        if (std::fabs(swaption1_cash->NPV() - swaption2_cash->NPV()) > 1.0e-6)
                            BOOST_ERROR("wrong spread treatment:" <<
                            "\nexercise date: " << exerciseDate <<
                            "\nlength: " << lengths_[j] <<
                            "\npay " << (type_ ? "fixed" : "floating") <<
                            "\nspread: " << io::rate(spreads[l]) <<
                            "\nvalue of original swaption:   " << swaption1_cash->NPV() <<
                            "\nvalue of equivalent swaption: " << swaption2_cash->NPV());
                    }
                }
            }
            return Y;
        }

        // Calculate sum of swaption NPV for given volatilities.
        Real CachedValue(std::vector<Real> volatility)
        {
            today_ = Date(13, March, 2002);
            settlement_ = Date(15, March, 2002);
            Settings::instance().evaluationDate() = today_;
            termStructure_.linkTo(flatRate(settlement_, 0.05, Actual365Fixed()));

            // Set date of exercize.
            Date exerciseDate = calendar_.advance(settlement_, 5 * Years);

            // Set start date.
            Date startDate = calendar_.advance(exerciseDate,
                                                    settlementDays_, Days);

            // Create vanilla swap.
            boost::shared_ptr<VanillaSwap> swap =
                MakeVanillaSwap(10 * Years, index_, 0.06)
                .withEffectiveDate(startDate)
                .withFixedLegTenor(1 * Years)
                .withFixedLegDayCount(fixedDayCount_);

            Real Y = 0;
            for (Size i = 0; i < volatility.size(); i++)
            {
                // Create swaption using vanilla swap.
                boost::shared_ptr<Swaption> swaption =
                    makeSwaption(swap, exerciseDate, volatility[i]);
                Y += swaption->NPV();
            }
            return Y;
        }

        // Calculate sum of swaption NPV for given strike rates.
        std::vector<Real> CashSettledSwaptions(std::vector<Real> strike)
        {
            std::vector<Real> Y(8);
            for (Size s = 0; s < strike.size(); s++)
            {
                for (Size i = 0; i < exercises_.size(); i++)
                {
                    for (Size j = 0; j < lengths_.size(); j++)
                    {
                        // Set date of exercise.
                        Date exerciseDate = calendar_.advance(today_, exercises_[i]);

                        // Set start date.
                        Date startDate = calendar_.advance(exerciseDate,
                                                           settlementDays_, Days);
                        // set date of maturity.
                        Date maturity =
                            calendar_.advance(startDate, lengths_[j],
                            floatingConvention_);

                        // Create floating schedule.
                        Schedule floatSchedule(startDate, maturity, floatingTenor_,
                                               calendar_, floatingConvention_,
                                               floatingConvention_,
                                               DateGeneration::Forward, false);

                        // Swap with fixed leg conventions: Business Days = Unadjusted, DayCount = 30/360.
                        Schedule fixedSchedule_u(startDate, maturity,
                                                 Period(fixedFrequency_),
                                                 calendar_, Unadjusted, Unadjusted,
                                                 DateGeneration::Forward, true);
                        boost::shared_ptr<VanillaSwap> swap_u360(
                            new VanillaSwap(type_, nominal_,
                            fixedSchedule_u, strike[s], Thirty360(),
                            floatSchedule, index_, 0.0,
                            index_->dayCounter()));

                        // Swap with fixed leg conventions: Business Days = Unadjusted, DayCount = Act/365.
                        boost::shared_ptr<VanillaSwap> swap_u365(
                            new VanillaSwap(type_, nominal_,
                            fixedSchedule_u, strike[s], Actual365Fixed(),
                            floatSchedule, index_, 0.0,
                            index_->dayCounter()));

                        // Swap with fixed leg conventions: Business Days = Modified Following, DayCount = 30/360.
                        Schedule fixedSchedule_a(startDate, maturity,
                                                 Period(fixedFrequency_),
                                                 calendar_, ModifiedFollowing,
                                                 ModifiedFollowing,
                                                 DateGeneration::Forward, true);
                        boost::shared_ptr<VanillaSwap> swap_a360(
                            new VanillaSwap(type_, nominal_,
                            fixedSchedule_a, strike[s], Thirty360(),
                            floatSchedule, index_, 0.0,
                            index_->dayCounter()));

                        // Swap with fixed leg conventions: Business Days = Modified Following, DayCount = Act/365.
                        boost::shared_ptr<VanillaSwap> swap_a365(
                            new VanillaSwap(type_, nominal_,
                            fixedSchedule_a, strike[s], Actual365Fixed(),
                            floatSchedule, index_, 0.0,
                            index_->dayCounter()));

                        // Create pricing engine.
                        boost::shared_ptr<PricingEngine> swapEngine(
                            new DiscountingSwapEngine(termStructure_));

                        // Set pricing engine in swaps.
                        swap_u360->setPricingEngine(swapEngine);
                        swap_a360->setPricingEngine(swapEngine);
                        swap_u365->setPricingEngine(swapEngine);
                        swap_a365->setPricingEngine(swapEngine);

                        const Leg& swapFixedLeg_u360 = swap_u360->fixedLeg();
                        const Leg& swapFixedLeg_a360 = swap_a360->fixedLeg();
                        const Leg& swapFixedLeg_u365 = swap_u365->fixedLeg();
                        const Leg& swapFixedLeg_a365 = swap_a365->fixedLeg();

                        // FlatForward curves.
                        Handle<YieldTermStructure> termStructure_u360(
                            boost::shared_ptr<YieldTermStructure>(
                            new FlatForward(settlement_, swap_u360->fairRate(),
                            Thirty360(), Compounded,
                            fixedFrequency_)));
                        Handle<YieldTermStructure> termStructure_a360(
                            boost::shared_ptr<YieldTermStructure>(
                            new FlatForward(settlement_, swap_a360->fairRate(),
                            Thirty360(), Compounded,
                            fixedFrequency_)));
                        Handle<YieldTermStructure> termStructure_u365(
                            boost::shared_ptr<YieldTermStructure>(
                            new FlatForward(settlement_, swap_u365->fairRate(),
                            Actual365Fixed(), Compounded,
                            fixedFrequency_)));
                        Handle<YieldTermStructure> termStructure_a365(
                            boost::shared_ptr<YieldTermStructure>(
                            new FlatForward(settlement_, swap_a365->fairRate(),
                            Actual365Fixed(), Compounded,
                            fixedFrequency_)));

                        // Annuity calculated by swap method fixedLegBPS().
                        // Fixed leg conventions: Unadjusted, 30/360.
                        Real annuity_u360 = swap_u360->fixedLegBPS() / 0.0001;
                        annuity_u360 = swap_u360->type() == VanillaSwap::Payer ?
                            -annuity_u360 : annuity_u360;
                        // Fixed leg conventions: ModifiedFollowing, act/365.
                        Real annuity_a365 = swap_a365->fixedLegBPS() / 0.0001;
                        annuity_a365 = swap_a365->type() == VanillaSwap::Payer ?
                            -annuity_a365 : annuity_a365;
                        // Fixed leg conventions: ModifiedFollowing, 30/360.
                        Real annuity_a360 = swap_a360->fixedLegBPS() / 0.0001;
                        annuity_a360 = swap_a360->type() == VanillaSwap::Payer ?
                            -annuity_a360 : annuity_a360;
                        // Fixed leg conventions: Unadjusted, act/365.
                        Real annuity_u365 = swap_u365->fixedLegBPS() / 0.0001;
                        annuity_u365 = swap_u365->type() == VanillaSwap::Payer ?
                            -annuity_u365 : annuity_u365;

                        // Calculation of Modified Annuity (cash settlement_).
                        // Fixed leg conventions of swap: unadjusted, 30/360.
                        Real cashannuity_u360 = 0.;
                        Size i;
                        for (i = 0; i < swapFixedLeg_u360.size(); i++)
                        {
                            cashannuity_u360 += swapFixedLeg_u360[i]->amount() / strike[s]
                                * termStructure_u360->discount(
                                swapFixedLeg_u360[i]->date());
                        }
                        // Fixed leg conventions of swap: unadjusted, act/365.
                        Real cashannuity_u365 = 0.;
                        for (i = 0; i < swapFixedLeg_u365.size(); i++)
                        {
                            cashannuity_u365 += swapFixedLeg_u365[i]->amount() / strike[s]
                                * termStructure_u365->discount(
                                swapFixedLeg_u365[i]->date());
                        }
                        // Fixed leg conventions of swap: modified following, 30/360.
                        Real cashannuity_a360 = 0.;
                        for (i = 0; i < swapFixedLeg_a360.size(); i++)
                        {
                            cashannuity_a360 += swapFixedLeg_a360[i]->amount() / strike[s]
                                * termStructure_a360->discount(
                                swapFixedLeg_a360[i]->date());
                        }
                        // Fixed leg conventions of swap: modified following, act/365.
                        Real cashannuity_a365 = 0.;
                        for (i = 0; i < swapFixedLeg_a365.size(); i++)
                        {
                            cashannuity_a365 += swapFixedLeg_a365[i]->amount() / strike[s]
                                * termStructure_a365->discount(
                                swapFixedLeg_a365[i]->date());
                        }

                        // Swaptions: underlying swap fixed leg conventions:
                        // unadjusted, 30/360.

                        // Physical settled swaption.
                        boost::shared_ptr<Swaption> swaption_p_u360 =
                            makeSwaption(swap_u360, exerciseDate, 0.20);
                        Real value_p_u360 = swaption_p_u360->NPV();

                        Y[0] += swaption_p_u360->NPV();

                        // Cash settled swaption.
                        boost::shared_ptr<Swaption> swaption_c_u360 =
                            makeSwaption(swap_u360, exerciseDate, 0.20,
                            Settlement::Cash);
                        Real value_c_u360 = swaption_c_u360->NPV();

                        Y[1] += swaption_c_u360->NPV();

                        // the NPV's ratio must be equal to annuities ratio.
                        Real npv_ratio_u360 = value_c_u360 / value_p_u360;
                        Real annuity_ratio_u360 = cashannuity_u360 / annuity_u360;

                        // Swaptions: underlying swap fixed leg conventions:
                        // modified following, act/365.

                        // Physical settled swaption.
                        boost::shared_ptr<Swaption> swaption_p_a365 =
                            makeSwaption(swap_a365, exerciseDate, 0.20);
                        Real value_p_a365 = swaption_p_a365->NPV();

                        Y[2] += swaption_p_a365->NPV();

                        // Cash settled swaption.
                        boost::shared_ptr<Swaption> swaption_c_a365 =
                            makeSwaption(swap_a365, exerciseDate, 0.20,
                            Settlement::Cash);
                        Real value_c_a365 = swaption_c_a365->NPV();

                        Y[3] = swaption_c_a365->NPV();

                        // the NPV's ratio must be equal to annuities ratio.
                        Real npv_ratio_a365 = value_c_a365 / value_p_a365;
                        Real annuity_ratio_a365 = cashannuity_a365 / annuity_a365;

                        // Swaptions: underlying swap fixed leg conventions:
                        // modified following, 30/360.

                        // Physical settled swaption.
                        boost::shared_ptr<Swaption> swaption_p_a360 =
                            makeSwaption(swap_a360, exerciseDate, 0.20);
                        Real value_p_a360 = swaption_p_a360->NPV();

                        Y[4] += swaption_p_a360->NPV();

                        // Cash settled swaption.
                        boost::shared_ptr<Swaption> swaption_c_a360 =
                            makeSwaption(swap_a360, exerciseDate, 0.20,
                            Settlement::Cash);
                        Real value_c_a360 = swaption_c_a360->NPV();

                        Y[5] += swaption_c_a360->NPV();

                        // the NPV's ratio must be equal to annuities ratio.
                        Real npv_ratio_a360 = value_c_a360 / value_p_a360;
                        Real annuity_ratio_a360 = cashannuity_a360 / annuity_a360;

                        // Swaptions: underlying swap fixed leg conventions:
                        // unadjusted, act/365.

                        // Physical settled swaption.
                        boost::shared_ptr<Swaption> swaption_p_u365 =
                            makeSwaption(swap_u365, exerciseDate, 0.20);
                        Real value_p_u365 = swaption_p_u365->NPV();

                        Y[6] += swaption_p_u365->NPV();

                        // Cash settled swaption.
                        boost::shared_ptr<Swaption> swaption_c_u365 =
                            makeSwaption(swap_u365, exerciseDate, 0.20,
                            Settlement::Cash);
                        Real value_c_u365 = swaption_c_u365->NPV();

                        Y[7] += swaption_c_u365->NPV();

                        // the NPV's ratio must be equal to annuities ratio.
                        Real npv_ratio_u365 = value_c_u365 / value_p_u365;
                        Real annuity_ratio_u365 = cashannuity_u365 / annuity_u365;

                        if (std::fabs(annuity_ratio_u360 - npv_ratio_u360) > 1e-10)
                        {
                            BOOST_ERROR("\n" <<
                                        "    The npv's ratio must be equal to " <<
                                        " annuities ratio" << "\n"
                                        "    Swaption " <<
                                        exercises_[i].units() << "y x " << lengths_[j].units() << "y" <<
                                        " (underlying swap fixed leg Unadjusted, 30/360)" << "\n" <<
                                        "    today_           : " <<
                                        today_ << "\n" <<
                                        "    Settlement date : " <<
                                        settlement_ << "\n" <<
                                        "    Exercise date   : " <<
                                        exerciseDate << "\n" <<
                                        "    Swap start date : " <<
                                        startDate << "\n" <<
                                        "    Swap end date   : " <<
                                        maturity << "\n" <<
                                        "    physical delivered swaption npv : " <<
                                        value_p_u360 << "\t\t\t" <<
                                        "    annuity : " <<
                                        annuity_u360 << "\n" <<
                                        "    cash delivered swaption npv :     " <<
                                        value_c_u360 << "\t\t\t" <<
                                        "    annuity : " <<
                                        cashannuity_u360 << "\n" <<
                                        "    npv ratio : " <<
                                        npv_ratio_u360 << "\n" <<
                                        "    annuity ratio : " <<
                                        annuity_ratio_u360 << "\n" <<
                                        "    difference : " <<
                                        (annuity_ratio_u360 - npv_ratio_u360));
                        }
                        if (std::fabs(annuity_ratio_a365 - npv_ratio_a365) > 1e-10)
                        {
                            BOOST_ERROR("\n" <<
                                        "    The npv's ratio must be equal to " <<
                                        " annuities ratio" << "\n"
                                        "    Swaption " <<
                                        exercises_[i].units() << "y x " << lengths_[j].units() << "y" <<
                                        " (underlying swap fixed leg Modified Following, act/365" << "\n" <<
                                        "    today_           : " <<
                                        today_ << "\n" <<
                                        "    Settlement date : " <<
                                        settlement_ << "\n" <<
                                        "    Exercise date   : " <<
                                        exerciseDate << "\n" <<
                                        "    Swap start date : " <<
                                        startDate << "\n" <<
                                        "    Swap end date   : " <<
                                        maturity << "\n" <<
                                        "    physical delivered swaption npv : " <<
                                        value_p_a365 << "\t\t\t" <<
                                        "    annuity : " <<
                                        annuity_a365 << "\n" <<
                                        "    cash delivered swaption npv :     " <<
                                        value_c_a365 << "\t\t\t" <<
                                        "    annuity : " <<
                                        cashannuity_a365 << "\n" <<
                                        "    npv ratio : " <<
                                        npv_ratio_a365 << "\n" <<
                                        "    annuity ratio : " <<
                                        annuity_ratio_a365 << "\n" <<
                                        "    difference : " <<
                                        (annuity_ratio_a365 - npv_ratio_a365));
                        }
                        if (std::fabs(annuity_ratio_a360 - npv_ratio_a360) > 1e-10)
                        {
                            BOOST_ERROR("\n" <<
                                        "    The npv's ratio must be equal to " <<
                                        " annuities ratio" << "\n"
                                        "    Swaption " <<
                                        exercises_[i].units() << "y x " << lengths_[j].units() << "y" <<
                                        " (underlying swap fixed leg Unadjusted, 30/360)" << "\n" <<
                                        "    today_           : " <<
                                        today_ << "\n" <<
                                        "    Settlement date : " <<
                                        settlement_ << "\n" <<
                                        "    Exercise date   : " <<
                                        exerciseDate << "\n" <<
                                        "    Swap start date : " <<
                                        startDate << "\n" <<
                                        "    Swap end date   : " <<
                                        maturity << "\n" <<
                                        "    physical delivered swaption npv : " <<
                                        value_p_a360 << "\t\t\t" <<
                                        "    annuity : " <<
                                        annuity_a360 << "\n" <<
                                        "    cash delivered swaption npv :     " <<
                                        value_c_a360 << "\t\t\t" <<
                                        "    annuity : " <<
                                        cashannuity_a360 << "\n" <<
                                        "    npv ratio : " <<
                                        npv_ratio_a360 << "\n" <<
                                        "    annuity ratio : " <<
                                        annuity_ratio_a360 << "\n" <<
                                        "    difference : " <<
                                        (annuity_ratio_a360 - npv_ratio_a360));
                        }
                        if (std::fabs(annuity_ratio_u365 - npv_ratio_u365) > 1e-10)
                        {
                            BOOST_ERROR("\n" <<
                                        "    The npv's ratio must be equal to " <<
                                        " annuities ratio" << "\n"
                                        "    Swaption " <<
                                        exercises_[i].units() << "y x " << lengths_[j].units() << "y" <<
                                        " (underlying swap fixed leg Unadjusted, act/365)" << "\n" <<
                                        "    today_           : " <<
                                        today_ << "\n" <<
                                        "    Settlement date : " <<
                                        settlement_ << "\n" <<
                                        "    Exercise date   : " <<
                                        exerciseDate << "\n" <<
                                        "    Swap start date : " <<
                                        startDate << "\n" <<
                                        "    Swap end date   : " <<
                                        maturity << "\n" <<
                                        "    physical delivered swaption npv : " <<
                                        value_p_u365 << "\t\t\t" <<
                                        "    annuity : " <<
                                        annuity_u365 << "\n" <<
                                        "    cash delivered swaption npv :     " <<
                                        value_c_u365 << "\t\t\t" <<
                                        "    annuity : " <<
                                        cashannuity_u365 << "\n" <<
                                        "    npv ratio : " <<
                                        npv_ratio_u365 << "\n" <<
                                        "    annuity ratio : " <<
                                        annuity_ratio_u365 << "\n" <<
                                        "    difference : " <<
                                        (annuity_ratio_u365 - npv_ratio_u365));
                        }
                    }
                }
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
                "param", ""
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, Variation& v)
        {
                stm << v.param_
                    << ";" << v.swaptionNpv_
                    << std::endl;
                return stm;
            }

        Real param_;
        Real swaptionNpv_;
    };
}


//Test swaption dependency on strike.
bool AdjointSwaptionTest::testStrikeDependency()
{
    BOOST_MESSAGE("Testing swaption dependency on strike...");

    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVars vars;

    // Create output streams for plots.
    cl::AdjointTestOutput outPerform("AdjointSwaption//StrikeDependency"
                                                , { { "filename", "AdjointPerformance" }
                                                  , { "not_clear", "Not" }
                                                  , { "line_box_width", "-5" }
                                                  , { "title", "Swaption NPV differentiation performance with respect to strike rate" }
                                                  , { "ylabel", "Time (s)" }
                                                  , { "xlabel", "Number of strike rates" } });

    cl::AdjointTestOutput outAdjoint("AdjointSwaption//StrikeDependency"
                                     , { { "filename", "Adjoint" }
                                       , { "not_clear", "Not" }
                                       , { "title", "Swaption NPV adjoint differentiation with respect to strike rate" }
                                       , { "cleanlog", "false" }
                                       , { "ylabel", "Time (s)" }
                                       , { "xlabel", "Number of strike rates" } });

    cl::AdjointTestOutput outSize("AdjointSwaption//StrikeDependency"
                                      , { { "filename", "TapeSize" }
                                        , { "not_clear", "Not" }
                                        , { "title", "Tape size dependence on number of strike rates" }
                                        , { "cleanlog", "false" }
                                        , { "ylabel", "Size (MB)" }
                                        , { "xlabel", "Number of strike rates" } });

    cl::AdjointTestOutput out("AdjointSwaption//StrikeDependency//output"
                                       , { { "filename", "SwaptNPVonStrikes" }
                                         , { "ylabel", "Swaption NPV" }
                                         , { "not_clear", "Not" }
                                         , { "title", "Swaption NPV dependence on strike rate" }
                                         , { "cleanlog", "false" }
                                         , { "xlabel", "Strike rate" } });

    // Set limit values.
    Size maxSize = 50;
    double startStrike = 0.01;
    double maxStrike = 0.06;
    double step = (maxStrike - startStrike) / maxSize;

#if defined CL_GRAPH_GEN
    Size startSize = 5;
#else
    Size startSize = maxSize - 1;
#endif


    // Create vector of strike rates.
    std::vector<double> strikes;
    for (Size i = 0; i < maxSize; i++)
    {
        strikes.push_back(startStrike + i*step);
    }

    // Create vectors for store performance results.
    std::vector<PerformanceTime> performResults;
    std::vector<AdjointTime> adjointResults;
    std::vector<TapeSize> tapeSizeResults;

    // Create timer.
    boost::timer timer;
    Size sizeof_indep;

    outPerform.log() << "Testing swaption dependency on strike." << std::endl << std::endl;

    // Variate number of independent variables.
    for (Size s = startSize; s <= maxSize; s += 5)
    {
        // Create temporary structures for store performance results.
        PerformanceTime performTime;
        AdjointTime adTime;
        TapeSize tSize;

       
        sizeof_indep = s;

        // Store number of independent variables
        performTime.indepVarNumber_ = sizeof_indep;
        adTime.indepVarNumber_ = sizeof_indep;
        tSize.indepVarNumber_ = sizeof_indep;

        outPerform.log() << "Number of independent variables: " << sizeof_indep << std::endl;

        //Set number of dependent vaiables
        Size sizeof_dep = 2;

        outPerform.log() << "Number of dependent variables: " << sizeof_dep << std::endl;

        outPerform.log() << "Start of tape recording: " << currentTime() << std::endl;

        // Start timing of tape recording.
        timer.restart();

        // Create vector of independent variables.
        std::vector<cl::TapeDouble> strike(sizeof_indep);

        // Create vector of doubles for Jacobian calculation.
        std::vector<double> strikeD(sizeof_indep);

        // Initialize vector of independent variables
        // and vector of doubles for Jacobian calculation.
        for (Size i = 0; i < sizeof_indep; i++)
        {
            strike[i] = strikes[i];
            strikeD[i] = strikes[i];
        }

        // Create vector of dependent variables.
        std::vector<cl::TapeDouble> Y(sizeof_dep);

        // Start taping. Declare strikes as independent variables.
        Independent(strike);

        // Calculate sum of swaption Net Present Values.
        Y = vars.StrikeDependency(strike);

        // End of tape recording. Declare sum of swaption NPVs as dependent variable.
        // Differentiaion will be held with respect to the independent variables vector.
        cl::TapeFunction<double> f(strike, Y);

        // Store time of tape recording.
        performTime.timeTapeRecording_ = timer.elapsed();

        outPerform.log() << "End of tape recording. " << std::endl;
        outPerform.log() << "Time for tape recording: " << performTime.timeTapeRecording_ << std::endl;

        // Store size of tape.
        tSize.memory_ = f.Memory();

        outPerform.log() << "Start of differentiation using Jacobian: " << currentTime() << std::endl;

        // Start timing for calculating derivatives by adjoint Jacobian.
        timer.restart();
        std::vector<double> Jacobian;

        //Compute derivatives using Jacobian.
        Jacobian = f.Jacobian(strikeD);

        // Store time of adjoint differntiation.
        performTime.timeAdjoint_ = timer.elapsed();
        adTime.timeAdjoint_ = performTime.timeAdjoint_;

        outPerform.log() << "End of differentiation." << std::endl;
        outPerform.log() << "Time for Jacobian: " << performTime.timeAdjoint_ << std::endl;

        outPerform.log() << "Start of  differentiation using finite differences method: " << currentTime() << std::endl;
        // Start timing for calculating derivatives by finite difference.
        timer.restart();

        //Start differentiation using finite differences method.
        double h = 1.0e-8;
        //  Create vector for Jacobian
        std::vector<Real> Jacobian_Finite(sizeof_dep*sizeof_indep);

        // Calculate sum of swaption NPVs.
        std::vector<Real> YF = vars.StrikeDependency(strike);
        for (Size i = 0; i < sizeof_indep; i++)
        {
            strike[i] += h;
            std::vector<Real> FiniteY(sizeof_dep);

            // Calculate sum of swaption NPVs with shifting
            // current strike rate with a step size of + h.
            FiniteY = vars.StrikeDependency(strike);
            strike[i] -= h;

            for (Size j = 0; j < sizeof_dep; j++)
            {
                //Evaluate derivative using finite difference.
                Jacobian_Finite[sizeof_indep*j + i] = (FiniteY[j] - YF[j]) / h;
            }
        }

        // Store time of finite difference differntiation.
        performTime.timeAnalytical_ = timer.elapsed();

        outPerform.log() << "End of differentiation using finite  differences method." << std::endl;
        outPerform.log() << "Time for finite differences method: " << performTime.timeAnalytical_ << std::endl;

        //Adding new data to the result vectors.
        performResults.push_back(performTime);
        adjointResults.push_back(adTime);
        tapeSizeResults.push_back(tSize);

        // Check derivatives calculated by adjoint and finite differences methods.
        result = checkWithFiniteDiff(Jacobian, Jacobian_Finite, 1e-2, 1e-10);

        outPerform.log() << "Derivatives Jacobian adjoint and Jacobian finite differences: " << std::endl;
        if (result)
        {
            for (int i = 0; i < sizeof_indep; i++)
            {
                for (int j = 0; j < sizeof_dep; j++)
                {
                    outPerform.log() 
                        << "dy[" << j << "] / dx[" << i << "] : "
                        << Jacobian[sizeof_indep*j + i] << " \t"
                        << Jacobian_Finite[sizeof_indep*j + i] << std::endl;
                }
            }
        }

        outPerform.log() << std::endl << std::endl;
    }

std::vector<Variation> variationResults;
#if defined CL_GRAPH_GEN
    // Set new limit values.
    startStrike = 0.07;
    maxStrike = 0.08;
    step = (maxStrike - startStrike) / maxSize;

    std::vector<Real> strike;
    strike.push_back(0.01);

    // Recalculate input strikes vector.
    for (Size i = 0; i < maxSize; i++)
    {
        strike.push_back(startStrike + i * step);
    }

    Real stepv = (strike[1] - strike[0]) / maxSize;

    // Calculate swaption NPV on variation of strike.
    for (Size i = 0; i < maxSize; i++)
    {
        strike[0] += i * stepv;

        Variation var;
        var.param_ = strike[0];
        var.swaptionNpv_ = vars.StrikeDependency(strike)[0];
        variationResults.push_back(var);

        strike[0] -= i * stepv;
    }
#endif

    
    //Output results to .csv and .plt files.
    outPerform << performResults;
    outAdjoint << adjointResults;
    outSize << tapeSizeResults;
    out << variationResults;

#endif
    return result;
}

// Test swaption dependency on spread.
bool AdjointSwaptionTest::testSpreadDependency()
{
    BOOST_MESSAGE("Testing swaption dependency on spread...");

    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVars vars;

    // Create output streams for plots.
    cl::AdjointTestOutput outPerform("AdjointSwaption//SpreadDependency"
                                                , { { "filename", "AdjointPerformance" }
                                                  , { "not_clear", "Not" }
                                                  , { "line_box_width", "-5" }
                                                  , { "title", "Swaption NPV differentiation performance with respect to spread" }
                                                  , { "ylabel", "Time (s)" }
                                                  , { "xlabel", "Number of spreads" } });

    cl::AdjointTestOutput outAdjoint("AdjointSwaption//SpreadDependency"
                                     , { { "filename", "Adjoint" }
                                       , { "not_clear", "Not" }
                                       , { "title", "Swaption NPV adjoint differentiation with respect to spread" }
                                       , { "cleanlog", "false" }
                                       , { "ylabel", "Time (s)" }
                                       , { "xlabel", "Number of spreads" } });

    cl::AdjointTestOutput outSize("AdjointSwaption//SpreadDependency"
                                      , { { "filename", "TapeSize" }
                                        , { "not_clear", "Not" }
                                        , { "title", "Tape size dependence on number of spreads" }
                                        , { "cleanlog", "false" }
                                        , { "ylabel", "Size (MB)" }
                                        , { "xlabel", "Number of spreads" } });

    cl::AdjointTestOutput out("AdjointSwaption//SpreadDependency//output"
                                       , { { "filename", "SwaptNPVonSpreads" }
                                         , { "ylabel", "Swaption NPV" }
                                         , { "not_clear", "Not" }
                                         , { "title", "Swaption NPV dependence on spread" }
                                         , { "cleanlog", "false" }
                                         , { "xlabel", "Spread" } });


    // Set limit values.
    Size maxSize = 51;
    double startSpread = -0.01;
    double maxSpread = 0.01;
    double step = (maxSpread - startSpread) / maxSize;

#if defined CL_GRAPH_GEN
    Size startSize = 0;
#else
    Size startSize = maxSize - 1;
#endif

    // Create vector of spreads.
    std::vector<double> spreads;
    double value;
    for (Size i = 0; i < maxSize; i++)
    {
        value = startSpread + i * step;
        value == 0 ? value += 0.001 : value;
        spreads.push_back(value);
    }

    // Create vectors for store performance results.
    std::vector<PerformanceTime> performResults;
    std::vector<AdjointTime> adjointResults;
    std::vector<TapeSize> tapeSizeResults;

    // Create timer.
    boost::timer timer;
    Size sizeof_indep;

    outPerform.log() << "Testing swaption dependency on spread." << std::endl << std::endl;

    // Variate number of independent variables.
    for (Size s = startSize; s < maxSize; s += 5)
    {
        // Create temporary structures for store performance results.
        PerformanceTime performTime;
        AdjointTime adTime;
        TapeSize tSize;

        s == 50 ? s = 49 : s;
        sizeof_indep = s + 1;

        // Store number of independent variables.
        performTime.indepVarNumber_ = sizeof_indep;
        adTime.indepVarNumber_ = sizeof_indep;
        tSize.indepVarNumber_ = sizeof_indep;

        outPerform.log() << "Number of independent variables: " << sizeof_indep << std::endl;

        //Set number of dependent vaiables.
        Size sizeof_dep = 2;

        outPerform.log() << "Number of dependent variables: " << sizeof_dep << std::endl;

        outPerform.log() << "Start of tape recording: " << currentTime() << std::endl;

        // Start timing of tape recording.
        timer.restart();

        // Create vector of independent variables.
        std::vector<cl::TapeDouble> spread(sizeof_indep);

        // Create vector of doubles for Jacobian calculation.
        std::vector<double> spreadD(sizeof_indep);

        // Initialize vector of independent variables
        // and vector of doubles for Jacobian calculation.
        for (Size i = 0; i < sizeof_indep; i++)
        {
            spread[i] = spreads[i];
            spreadD[i] = spreads[i];
        }

        // Create vector of dependent variables.
        std::vector<cl::TapeDouble> Y(sizeof_dep);

        // Start taping. Declare spreads as independent variables.
        Independent(spread);

        // Calculate sum of swaption Net Present Values.
        Y = vars.SpreadDependency(spread);

        // End of tape recording. Declare sum of swaption NPVs as dependent variable.
        // Differentiaion will be held with respect to the independent variables vector.
        cl::TapeFunction<double> f(spread, Y);

        // Store time of tape recording.
        performTime.timeTapeRecording_ = timer.elapsed();

        outPerform.log() << "End of tape recording. " << std::endl;
        outPerform.log() << "Time for tape recording: " << performTime.timeTapeRecording_ << std::endl;

        // Store size of tape.
        tSize.memory_ = f.Memory();

        outPerform.log() << "Start of differentiation using Jacobian: " << currentTime() << std::endl;

        // Start timing for calculating derivatives by adjoint Jacobian.
        timer.restart();
        std::vector<double> Jacobian;

        // Compute derivatives using Jacobian.
        Jacobian = f.Jacobian(spreadD);

        // Store time of adjoint differntiation.
        performTime.timeAdjoint_ = timer.elapsed();
        adTime.timeAdjoint_ = performTime.timeAdjoint_;

        outPerform.log() << "End of differentiation." << std::endl;
        outPerform.log() << "Time for Jacobian: " << performTime.timeAdjoint_ << std::endl;

        outPerform.log() << "Start of  differentiation using finite differences method: " << currentTime() << std::endl;

        // Start timing for calculating derivatives by finite difference.
        timer.restart();

        // Start differentiation using finite differences method.
        double h = 1.0e-8;
        // Create vector for Jacobian
        std::vector<Real> Jacobian_Finite(sizeof_dep*sizeof_indep);

        // Calculate sum of swaption NPVs.
        std::vector<Real> YF = vars.SpreadDependency(spread);
        for (Size i = 0; i < sizeof_indep; i++)
        {
            spread[i] += h;
            std::vector<Real> FiniteY(sizeof_dep);

            // Calculate sum of swaption NPVs with shifting
            // current spread with a step size of + h.
            FiniteY = vars.SpreadDependency(spread);
            spread[i] -= h;

            for (Size j = 0; j < sizeof_dep; j++)
            {
                //Evaluate derivative using finite difference.
                Jacobian_Finite[sizeof_indep*j + i] = (FiniteY[j] - YF[j]) / h;
            }
        }

        // Store time of finite difference differntiation.
        performTime.timeAnalytical_ = timer.elapsed();

        outPerform.log() << "End of differentiation using finite  differences method." << std::endl;
        outPerform.log() << "Time for finite differences method: " << performTime.timeAnalytical_ << std::endl;

        //Adding new data to the result vectors.
        performResults.push_back(performTime);
        adjointResults.push_back(adTime);
        tapeSizeResults.push_back(tSize);

        // Check derivatives calculated by adjoint and finite differences methods.
        result = checkWithFiniteDiff(Jacobian, Jacobian_Finite, 1e-2, 1e-10);

        outPerform.log() << "Derivatives Jacobian adjoint and Jacobian finite differences: " << std::endl;
        if (result)
        {
            for (int i = 0; i < sizeof_indep; i++)
            {
                for (int j = 0; j < sizeof_dep; j++)
                {
                    outPerform.log()
                        << "dy[" << j << "] / dx[" << i << "] : "
                        << Jacobian[sizeof_indep*j + i] << " \t"
                        << Jacobian_Finite[sizeof_indep*j + i] << std::endl;
                }
            }
        }

        outPerform.log() << std::endl << std::endl;
    }
    std::vector<Variation> variationResults;
#if defined CL_GRAPH_GEN
    // Calculate dependence on spread variation
    for (Size i = 0; i < maxSize; i++)
    {
        std::vector<Real> spread;
        spread.push_back(spreads[i]);
        Variation var;
        var.param_ = spread[0];
        var.swaptionNpv_ = vars.SpreadDependency(spread)[0];
        variationResults.push_back(var);
    }
#endif

    //Output results to .csv and .plt files.
    out << variationResults;
    outPerform << performResults;
    outAdjoint << adjointResults;
    outSize << tapeSizeResults;

#endif
    return result;
}

// Test swaption treatment of spread.
bool AdjointSwaptionTest::testSpreadTreatment()
{
    BOOST_MESSAGE("Testing swaption treatment of spread...");

    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVars vars;

    // Create output streams for plots.
    cl::AdjointTestOutput outPerform("AdjointSwaption//SpreadTreatment"
                                                , { { "filename", "AdjointPerformance" }
                                                  , { "not_clear", "Not" }
                                                  , { "line_box_width", "-5" }
                                                  , { "title", "Swaption NPV differentiation performance with respect to spread" }
                                                  , { "ylabel", "Time (s)" }
                                                  , { "xlabel", "Number of spreads" } });

    cl::AdjointTestOutput outAdjoint("AdjointSwaption//SpreadTreatment"
                                     , { { "filename", "Adjoint" }
                                       , { "not_clear", "Not" }
                                       , { "title", "Swaption NPV adjoint differentiation with respect to spread" }
                                       , { "cleanlog", "false" }
                                       , { "ylabel", "Time (s)" }
                                       , { "xlabel", "Number of spreads" } });

    cl::AdjointTestOutput outSize("AdjointSwaption//SpreadTreatment"
                                      , { { "filename", "TapeSize" }
                                        , { "not_clear", "Not" }
                                        , { "title", "Tape size dependence on number of spreads" }
                                        , { "cleanlog", "false" }
                                        , { "ylabel", "Size (MB)" }
                                        , { "xlabel", "Number of spreads" } });

    cl::AdjointTestOutput out("AdjointSwaption//SpreadTreatment//output"
                                       , { { "filename", "SwaptNPVonSpreads" }
                                         , { "ylabel", "Swaption NPV" }
                                         , { "not_clear", "Not" }
                                         , { "title", "Swaption NPV dependence on spread" }
                                         , { "cleanlog", "false" }
                                         , { "xlabel", "Spread" } });

    // Set limit values.
    Size maxSize = 20;
    double startSpread = -0.01;
    double maxSpread = 0.01;
    double step = (maxSpread - startSpread) / maxSize;

#if defined CL_GRAPH_GEN
    Size startSize = 0;
#else
    Size startSize = maxSize - 1;
#endif

    // Create vector of spreads.
    std::vector<double> spreads;
    double value;
    for (Size i = 0; i < maxSize; i++)
    {
        value = startSpread + i * step;
        value == 0 ? value += 0.00001 : value;
        spreads.push_back(value);
    }

    // Create vectors for store performance results.
    std::vector<PerformanceTime> performResults;
    std::vector<AdjointTime> adjointResults;
    std::vector<TapeSize> tapeSizeResults;

    // Create timer.
    boost::timer timer;
    Size sizeof_indep;

    outPerform.log() << "Testing swaption treatment of spread." << std::endl << std::endl;

    // Variate number of independent variables.
    for (Size s = startSize; s < maxSize; s++)
    {
        // Create temporary structures for store performance results.
        PerformanceTime performTime;
        AdjointTime adTime;
        TapeSize tSize;

        
        sizeof_indep = s + 1;

        // Store number of independent variables.
        performTime.indepVarNumber_ = sizeof_indep;
        adTime.indepVarNumber_ = sizeof_indep;
        tSize.indepVarNumber_ = sizeof_indep;

        outPerform.log() << "Number of independent variables: " << sizeof_indep << std::endl;

        //Set number of dependent vaiables.
        Size sizeof_dep = 4;

        outPerform.log() << "Number of dependent variables: " << sizeof_dep << std::endl;

        outPerform.log() << "Start of tape recording: " << currentTime() << std::endl;

        // Start timing of tape recording.
        timer.restart();

        // Create vector of independent variables.
        std::vector<cl::TapeDouble> spread(sizeof_indep);

        // Create vector of doubles for Jacobian calculation.
        std::vector<double> spreadD(sizeof_indep);

        // Initialize vector of independent variables
        // and vector of doubles for Jacobian calculation.
        for (Size i = 0; i < sizeof_indep; i++)
        {
            spread[i] = spreads[i];
            spreadD[i] = spreads[i];
        }

        // Create vector of dependent variables.
        std::vector<cl::TapeDouble> Y(sizeof_dep);

        // Start taping. Declare spreads as independent variables.
        Independent(spread);

        // Calculate sum of swaption Net Present Values.
        Y = vars.SpreadTreatment(spread);

        // End of tape recording. Declare sum of swaption NPVs as dependent variable.
        // Differentiaion will be held with respect to the independent variables vector.
        cl::TapeFunction<double> f(spread, Y);

        // Store time of tape recording.
        performTime.timeTapeRecording_ = timer.elapsed();

        outPerform.log() << "End of tape recording. " << std::endl;
        outPerform.log() << "Time for tape recording: " << performTime.timeTapeRecording_ << std::endl;

        // Store size of tape.
        tSize.memory_ = f.Memory();

        outPerform.log() << "Start of differentiation using Jacobian: " << currentTime() << std::endl;

        // Start timing for calculating derivatives by adjoint Jacobian.
        timer.restart();
        std::vector<double> Jacobian;

        // Compute derivatives using Jacobian.
        Jacobian = f.Jacobian(spreadD);

        // Store time of finite difference differntiation
        performTime.timeAdjoint_ = timer.elapsed();
        adTime.timeAdjoint_ = performTime.timeAdjoint_;

        outPerform.log() << "End of differentiation." << std::endl;
        outPerform.log() << "Time for Jacobian: " << performTime.timeAdjoint_ << std::endl;

        outPerform.log() << "Start of  differentiation using finite differences method: " << currentTime() << std::endl;
        // Start timing for calculating derivatives by finite difference.
        timer.restart();

        // Start differentiation using finite differences method.
        double h = 1.0e-8;
        // Create vector for Jacobian.
        std::vector<Real> Jacobian_Finite(sizeof_dep*sizeof_indep);

        // Calculate sum of swaption NPVs.
        std::vector<Real> YF = vars.SpreadTreatment(spread);
        for (Size i = 0; i < sizeof_indep; i++)
        {
            spread[i] += h;
            std::vector<Real> FiniteY(sizeof_dep);

            // Calculate sum of swaption NPVs with shifting
            // current spread with a step size of + h.
            FiniteY = vars.SpreadTreatment(spread);
            spread[i] -= h;

            for (Size j = 0; j < sizeof_dep; j++)
            {
                //Evaluate derivative using finite difference.
                Jacobian_Finite[sizeof_indep*j + i] = (FiniteY[j] - YF[j]) / h;
            }
        }

        // Store time of finite difference differntiation.
        performTime.timeAnalytical_ = timer.elapsed();

        outPerform.log() << "End of differentiation using finite  differences method." << std::endl;
        outPerform.log() << "Time for finite differences method: " << performTime.timeAnalytical_ << std::endl;

        //Adding new data to the result vectors.
        performResults.push_back(performTime);
        adjointResults.push_back(adTime);
        tapeSizeResults.push_back(tSize);

        // Check derivatives calculated by adjoint and finite differences methods.
        result = checkWithFiniteDiff(Jacobian, Jacobian_Finite, 1e-2, 1e-10);

        outPerform.log() << "Derivatives Jacobian adjoint and Jacobian finite differences: " << std::endl;
        if (result)
        {
            for (int i = 0; i < sizeof_indep; i++)
            {
                for (int j = 0; j < sizeof_dep; j++)
                {
                    outPerform.log()
                        << "dy[" << j << "] / dx[" << i << "] : "
                        << Jacobian[sizeof_indep*j + i] << " \t"
                        << Jacobian_Finite[sizeof_indep*j + i] << std::endl;
                }
            }
        }

        outPerform.log() << std::endl << std::endl;
    }
std::vector<Variation> variationResults;
#if defined CL_GRAPH_GEN
    // Calculate dependence on spread variation
    for (Size i = 0; i < maxSize; i++)
    {
        std::vector<Real> spread;
        spread.push_back(spreads[i]);
        Variation var;
        var.param_ = spread[0];
        var.swaptionNpv_ = vars.SpreadTreatment(spread)[0];
        variationResults.push_back(var);
    }
#endif

    //Output results to .csv and .plt files.
    out << variationResults;
    outPerform << performResults;
    outAdjoint << adjointResults;
    outSize << tapeSizeResults;

#endif
    return result;
}

// Test swaption value against cached value.
bool AdjointSwaptionTest::testCachedValue()
{
    BOOST_MESSAGE("Testing swaption value against cached value...");

    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVars vars;

    // Create output streams for plots.
    cl::AdjointTestOutput outPerform("AdjointSwaption//CachedValue"
                                                , { { "filename", "AdjointPerformance" }
                                                  , { "not_clear", "Not" }
                                                  , { "line_box_width", "-5" }
                                                  , { "title", "Swaption NPV differentiation performance with respect to volatility" }
                                                  , { "ylabel", "Time (s)" }
                                                  , { "xlabel", "Number of volatilities" } });

    cl::AdjointTestOutput outAdjoint("AdjointSwaption//CachedValue"
                                     , { { "filename", "Adjoint" }
                                       , { "not_clear", "Not" }
                                       , { "title", "Swaption NPV adjoint differentiation with respect to volatility" }
                                       , { "cleanlog", "false" }
                                       , { "ylabel", "Time (s)" }
                                       , { "xlabel", "Number of volatilities" } });

    cl::AdjointTestOutput outSize("AdjointSwaption//CachedValue"
                                      , { { "filename", "TapeSize" }
                                        , { "not_clear", "Not" }
                                        , { "title", "Tape size dependence on number of volatilities" }
                                        , { "cleanlog", "false" }
                                        , { "ylabel", "Size (MB)" }
                                        , { "xlabel", "Number of volatilities" } });

    cl::AdjointTestOutput out("AdjointSwaption//CachedValue//output"
                                       , { { "filename", "SwaptNPVonVol" }
                                         , { "ylabel", "SwaptionNPV" }
                                         , { "not_clear", "Not" }
                                         , { "title", "Swaption NPV dependence on volatilities" }
                                         , { "cleanlog", "false" }
                                         , { "xlabel", "Volatility" } });


    // Set limit values.
    Size maxSize = 401;
    Rate startVol = 0.15;
    Rate maxVol = 0.25;
    Real step = (maxVol - startVol) / maxSize;

#if defined CL_GRAPH_GEN
    Size startSize = 0;
#else
    Size startSize = maxSize - 1;
#endif

    // Create vector of volatilities.
    std::vector<Volatility> volatilities;
    for (Size i = 0; i < maxSize; i++)
        volatilities.push_back(startVol + i * step);

    // Create vectors for store performance results.
    std::vector<PerformanceTime> performResults;
    std::vector<AdjointTime> adjointResults;
    std::vector<TapeSize> tapeSizeResults;

    // Create timer.
    boost::timer timer;
    Size sizeof_indep;

    outPerform.log() << "Testing swaption value against cached value." << std::endl << std::endl;

    // Variate number of independent variables.
    for (Size s = startSize; s < maxSize; s += 25)
    {
        // Create temporary structures for store performance results.
        PerformanceTime performTime;
        AdjointTime adTime;
        TapeSize tSize;

        

        s == 400 ? s = 399 : s;
        sizeof_indep = s + 1;

        // Store number of independent variables.
        performTime.indepVarNumber_ = sizeof_indep;
        adTime.indepVarNumber_ = sizeof_indep;
        tSize.indepVarNumber_ = sizeof_indep;

        outPerform.log() << "Number of independent variables: " << sizeof_indep << std::endl;
        outPerform.log() << "Number of dependent variables: 1" << std::endl;

        outPerform.log() << "Start of tape recording: " << currentTime() << std::endl;

        // Start timing of tape recording.
        timer.restart();

        // Create vector of independent variables.
        std::vector<cl::TapeDouble> volatility(sizeof_indep);
        for (Size i = 0; i < sizeof_indep; i++)
        {
            volatility[i] = volatilities[i];
        }

        // Start taping. Declare volatilities as independent variables.
        Independent(volatility);

        // Create vector of dependent variables.
        std::vector<cl::TapeDouble> Y(1);

        // Calculate sum of swaption Net Present Values.
        Y[0] = vars.CachedValue(volatility);

        // End of tape recording. Declare sum of swaption NPVs as dependent variable.
        // Differentiaion will be held with respect to the independent variables vector.
        cl::TapeFunction<double> f(volatility, Y);

        // Store time of tape recording.
        performTime.timeTapeRecording_ = timer.elapsed();

        outPerform.log() << "End of tape recording. " << std::endl;
        outPerform.log() << "Time for tape recording: " << performTime.timeTapeRecording_ << std::endl;

        // Store size of tape.
        tSize.memory_ = f.Memory();

        std::vector<double> sf_Forward(sizeof_indep), sf_Reverse(sizeof_indep);
        
        //Start differentiation in Forward mode.
        double tf = gradForward(f, sf_Forward, outPerform, false, false);

        //Start differentiation in Reverse mode.
        double tr = gradReverse(f, sf_Reverse, outPerform, false, false);

        // Store time of finite difference differntiation
        performTime.timeAdjoint_ = std::min(tf, tr);
        adTime.timeAdjoint_ = performTime.timeAdjoint_;

        outPerform.log() << "Start of  differentiation using finite differences method: " << currentTime() << std::endl;
        // Start timing for calculating derivatives by finite difference.
        timer.restart();

        // Start differentiation using finite differences method.
        double h = 1.0e-8;
        // Create vector for derivatives.
        std::vector<Real> sf_Finite(sizeof_indep);

        // Calculate sum of swaption NPVs.
        Real YF = vars.CachedValue(volatility);
        for (Size i = 0; i < sizeof_indep; i++)
        {
            volatility[i] += h;
            //Evaluate derivative using finite difference.
            sf_Finite[i] = (vars.CachedValue(volatility) - YF) / h;
            volatility[i] -= h;
        }

        // Store time of finite difference differntiation.
        performTime.timeAnalytical_ = timer.elapsed();

        outPerform.log() << "End of differentiation using finite  differences method." << std::endl;
        outPerform.log() << "Time for finite differences method: " << performTime.timeAnalytical_ << std::endl;

        //Adding new data to the result vectors.
        performResults.push_back(performTime);
        adjointResults.push_back(adTime);
        tapeSizeResults.push_back(tSize);

        // Check derivatives calculated by adjoint and finite differences methods.
        result = checkWithFiniteDiff(sf_Forward, sf_Reverse, sf_Finite, outPerform, 1e-2, 1e-10);

        outPerform.log() << std::endl << std::endl;
    }

    std::vector<Variation> variationResults;
#if defined CL_GRAPH_GEN
    // Set new limit values.
    startVol = 0.23;
    maxVol = 0.25;
    step = (maxVol - startVol) / 100;

    std::vector<Real> volatility;
    volatility.push_back(0.15);

    // Recalculate input  volatilities vector.
    for (Size i = 0; i < 100; i++)
    {
         volatility.push_back(startVol + i * step);
    }

    Real stepv = (volatility[1] - volatility[0]) /10;

    // Calculate swaption NPV on variation of volatility.
    for (Size i = 0; i < 10; i++)
    {
        volatility[0] += i * stepv;

        Variation var;
        var.param_ = volatility[0];
        var.swaptionNpv_ = vars.StrikeDependency(volatility)[0];
        variationResults.push_back(var);

        volatility[0] -= i * stepv;
    }
#endif

    //Output results to .csv and .plt files.
    outPerform << performResults;
    outAdjoint << adjointResults;
    outSize << tapeSizeResults;
    out << variationResults;

#endif
    return result;
}

//Test cash settled swaptions modified annuity.
bool AdjointSwaptionTest::testCashSettledSwaptions()
{
    BOOST_MESSAGE("Testing cash settled swaptions modified annuity...");

    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVars vars;

    // Create output streams for plots.
    cl::AdjointTestOutput outPerform("AdjointSwaption//CashSettledSwaptions"
                                                , { { "filename", "AdjointPerformance" }
                                                  , { "not_clear", "Not" }
                                                  , { "line_box_width", "-5" }
                                                  , { "title", "Swaption NPV differentiation performance with respect to strike" }
                                                  , { "ylabel", "Time (s)" }
                                                  , { "xlabel", "Number of strikes" } });

    cl::AdjointTestOutput outAdjoint("AdjointSwaption//CashSettledSwaptions"
                                     , { { "filename", "Adjoint" }
                                       , { "not_clear", "Not" }
                                       , { "title", "Swaption NPV adjoint differentiation with respect to strike" }
                                       , { "cleanlog", "false" }
                                       , { "ylabel", "Time (s)" }
                                       , { "xlabel", "Number of strikes" } });

    cl::AdjointTestOutput outSize("AdjointSwaption//CashSettledSwaptions"
                                      , { { "filename", "TapeSize" }
                                        , { "not_clear", "Not" }
                                        , { "title", "Tape size dependence on number of strikes" }
                                        , { "cleanlog", "false" }
                                        , { "ylabel", "Size (MB)" }
                                        , { "xlabel", "Number of strikes" } });

    cl::AdjointTestOutput out("AdjointSwaption//CashSettledSwaptions//output"
                                       , { { "filename", "SwaptNPVonStrikes" }
                                       , { "ylabel", "SwaptionNPV" }
                                       , { "not_clear", "Not" }
                                       , { "title", "Swaption NPV dependence on strike" }
                                       , { "cleanlog", "false" }
                                       , { "xlabel", "Strike" } });

    // Set limit values.
    Size maxSize = 20;
    double startStrike = 0.01;
    double maxStrike = 0.06;
    double step = (maxStrike - startStrike) / maxSize;

#if defined CL_GRAPH_GEN
    Size startSize = 0;
#else
    Size startSize = maxSize - 1;
#endif

    // Create vector of strikes.
    std::vector<double> strikes;
    for (Size i = 0; i < maxSize; i++)
    {
        strikes.push_back(startStrike + i*step);
    }

    // Create vectors for store performance results.
    std::vector<PerformanceTime> performResults;
    std::vector<AdjointTime> adjointResults;
    std::vector<TapeSize> tapeSizeResults;

    // Create timer.
    boost::timer timer;
    Size sizeof_indep;

    outPerform.log() << "Testing cash settled swaptions modified annuity." << std::endl << std::endl;

    // Variate number of independent variables.
    for (Size s = startSize; s < maxSize; s++)
    {
        // Create temporary structures for store performance results.
        PerformanceTime performTime;
        AdjointTime adTime;
        TapeSize tSize;

        sizeof_indep = s + 1;

        // Store number of independent variables.
        performTime.indepVarNumber_ = sizeof_indep;
        adTime.indepVarNumber_ = sizeof_indep;
        tSize.indepVarNumber_ = sizeof_indep;
        
        outPerform.log() << "Number of independent variables: " << sizeof_indep << std::endl;
        outPerform.log() << "Start of tape recording: " << currentTime() << std::endl;

        // Start timing of tape recording.
        timer.restart();
                
        // Create vector of independent variables.
        std::vector<cl::TapeDouble> strike(sizeof_indep);

        // Create vector of doubles for Jacobian calculation.
        std::vector<double> strikeD(sizeof_indep);

        // Initialize vector of independent variables
        // and vector of doubles for Jacobian calculation.
        for (Size i = 0; i < sizeof_indep; i++)
        {
            strike[i] = strikes[i];
            strikeD[i] = strikes[i];
        }

        // Create vector of dependent variables.
        std::vector<cl::TapeDouble> Y;

        // Start taping. Declare strikes as independent variables.
        Independent(strike);

        // Calculate sum of swaption Net Present Values.
        Y = vars.CashSettledSwaptions(strike);

        // End of tape recording. Declare sum of swaption NPVs as dependent variable.
        // Differentiaion will be held with respect to the independent variables vector.
        cl::TapeFunction<double> f(strike, Y);

        // Store time of tape recording.
        performTime.timeTapeRecording_ = timer.elapsed();

        outPerform.log() << "End of tape recording. " << std::endl;
        outPerform.log() << "Time for tape recording: " << performTime.timeTapeRecording_ << std::endl;

        // Store size of tape.
        tSize.memory_ = f.Memory();

        //Set number of dependent vaiables.
        Size sizeof_dep = Y.size();

        outPerform.log() << "Number of dependent variables: " << sizeof_dep << std::endl;
        outPerform.log() << "Start of differentiation using Jacobian: " << currentTime() << std::endl;

        // Start timing for calculating derivatives by adjoint Jacobian.
        timer.restart();
        std::vector<double> Jacobian;

        // Compute derivatives using Jacobian.
        Jacobian = f.Jacobian(strikeD);

        // Store time of finite difference differntiation
        performTime.timeAdjoint_ = timer.elapsed();
        adTime.timeAdjoint_ = performTime.timeAdjoint_;

        outPerform.log() << "End of differentiation." << std::endl;
        outPerform.log() << "Time for Jacobian: " << performTime.timeAdjoint_ << std::endl;

        outPerform.log() << "Start of  differentiation using finite differences method: " << currentTime() << std::endl;
        // Start timing for calculating derivatives by finite difference.
        timer.restart();

        // Start differentiation using finite differences method.
        double h = 1.0e-8;
        // Create vector for Jacobian.
        std::vector<Real> Jacobian_Finite(sizeof_dep*sizeof_indep);

        // Calculate sum of swaption NPVs.
        std::vector<Real> YF = vars.CashSettledSwaptions(strike);
        for (Size i = 0; i < sizeof_indep; i++)
        {
            strike[i] += h;
            std::vector<Real> FiniteY(sizeof_dep);

            // Calculate sum of swaption NPVs with shifting
            // current spread with a step size of + h.
            FiniteY = vars.CashSettledSwaptions(strike);
            strike[i] -= h;

            for (Size j = 0; j < sizeof_dep; j++)
            {
                //Evaluate derivative using finite difference.
                Jacobian_Finite[sizeof_indep*j + i] = (FiniteY[j] - YF[j]) / h;
            }
        }

        // Store time of finite difference differntiation.
        performTime.timeAnalytical_ = timer.elapsed();

        outPerform.log() << "End of differentiation using finite  differences method." << std::endl;
        outPerform.log() << "Time for finite differences method: " << performTime.timeAnalytical_ << std::endl;

        //Adding new data to the result vectors.
        performResults.push_back(performTime);
        adjointResults.push_back(adTime);
        tapeSizeResults.push_back(tSize);

        // Check derivatives calculated by adjoint and finite differences methods.
        result = checkWithFiniteDiff(Jacobian, Jacobian_Finite, 1e-2, 1e-10);

        outPerform.log() << "Derivatives Jacobian adjoint and Jacobian finite differences: " << std::endl;
        if (result)
        {
            for (int i = 0; i < sizeof_indep; i++)
            {
                for (int j = 0; j < sizeof_dep; j++)
                {
                    outPerform.log()
                        << "dy[" << j << "] / dx[" << i << "] : "
                        << Jacobian[sizeof_indep*j + i] << " \t"
                        << Jacobian_Finite[sizeof_indep*j + i] << std::endl;
                }
            }
        }

        outPerform.log() << std::endl << std::endl;
    }

    std::vector<Variation> variationResults;
#if defined CL_GRAPH_GEN
    // Set new limit values.
    startStrike = 0.07;
    maxStrike = 0.08;
    step = (maxStrike - startStrike) / maxSize;

    std::vector<Real> strike;
    strike.push_back(0.01);

    // Recalculate input strikes vector.
    for (Size i = 0; i < maxSize; i++)
    {
        strike.push_back(startStrike + i * step);
    }

    Real stepv = (strike[1] - strike[0]) / maxSize;

    // Calculate swaption NPV on variation of strike.
    for (Size i = 0; i < maxSize; i++)
    {
        strike[0] += i * stepv;

        Variation var;
        var.param_ = strike[0];
        var.swaptionNpv_ = vars.StrikeDependency(strike)[0];
        variationResults.push_back(var);

        strike[0] -= i * stepv;
    }
#endif

    //Output results to .csv and .plt files.
    outPerform << performResults;
    outAdjoint << adjointResults;
    outSize << tapeSizeResults;
    out << variationResults;
#endif
    return result;
}



test_suite* AdjointSwaptionTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Adjoint swaption test");
    suite->add(QUANTLIB_TEST_CASE(&AdjointSwaptionTest::testStrikeDependency));
    suite->add(QUANTLIB_TEST_CASE(&AdjointSwaptionTest::testSpreadDependency));
    suite->add(QUANTLIB_TEST_CASE(&AdjointSwaptionTest::testSpreadTreatment));
    suite->add(QUANTLIB_TEST_CASE(&AdjointSwaptionTest::testCachedValue));
    suite->add(QUANTLIB_TEST_CASE(&AdjointSwaptionTest::testCashSettledSwaptions));
    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(AdjointSwaption)

BOOST_AUTO_TEST_CASE(testSwaptDependencyOnStrike)
{
    BOOST_CHECK(AdjointSwaptionTest::testStrikeDependency());
}

BOOST_AUTO_TEST_CASE(testSwaptDependencyOnSpread)
{
    BOOST_CHECK(AdjointSwaptionTest::testSpreadDependency());
}

BOOST_AUTO_TEST_CASE(testSwaptCorrectionOfSpread)
{
    BOOST_CHECK(AdjointSwaptionTest::testSpreadTreatment());
}

BOOST_AUTO_TEST_CASE(testSwaptCachedValue)
{
    BOOST_CHECK(AdjointSwaptionTest::testCachedValue());
}

BOOST_AUTO_TEST_CASE(testCashSettledSwaptions)
{
    BOOST_CHECK(AdjointSwaptionTest::testCashSettledSwaptions());
}

BOOST_AUTO_TEST_SUITE_END()

#endif
