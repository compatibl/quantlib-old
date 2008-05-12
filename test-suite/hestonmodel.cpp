/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2005, 2007 Klaus Spanderen

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

#include "hestonmodel.hpp"
#include "utilities.hpp"
#include <ql/processes/hestonprocess.hpp>
#include <ql/models/equity/hestonmodel.hpp>
#include <ql/models/equity/hestonmodelhelper.hpp>
#include <ql/pricingengines/vanilla/analytichestonengine.hpp>
#include <ql/pricingengines/vanilla/mceuropeanhestonengine.hpp>
#include <ql/pricingengines/blackformula.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/calendars/nullcalendar.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/daycounters/actualactual.hpp>
#include <ql/termstructures/yield/zerocurve.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
#include <ql/time/period.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;

void HestonModelTest::testBlackCalibration() {
    BOOST_MESSAGE(
       "Testing Heston model calibration using a flat volatility surface...");

    SavedSettings backup;

    /* calibrate a Heston model to a constant volatility surface without
       smile. expected result is a vanishing volatility of the volatility.
       In addition theta and v0 should be equal to the constant variance */

    Date today = Date::todaysDate();
    Settings::instance().evaluationDate() = today;

    DayCounter dayCounter = Actual360();
    Calendar calendar = NullCalendar();

    Handle<YieldTermStructure> riskFreeTS(flatRate(0.04, dayCounter));
    Handle<YieldTermStructure> dividendTS(flatRate(0.50, dayCounter));

    std::vector<Period> optionMaturities;
    optionMaturities.push_back(Period(1, Months));
    optionMaturities.push_back(Period(2, Months));
    optionMaturities.push_back(Period(3, Months));
    optionMaturities.push_back(Period(6, Months));
    optionMaturities.push_back(Period(9, Months));
    optionMaturities.push_back(Period(1, Years));
    optionMaturities.push_back(Period(2, Years));

    std::vector<boost::shared_ptr<CalibrationHelper> > options;
    Handle<Quote> s0(boost::shared_ptr<Quote>(new SimpleQuote(1.0)));
    Handle<Quote> vol(boost::shared_ptr<Quote>(new SimpleQuote(0.1)));
    Volatility volatility = vol->value();

    for (Size i = 0; i < optionMaturities.size(); ++i) {
        for (Real moneyness = -1.0; moneyness < 2.0; moneyness += 1.0) {
            // FLOATING_POINT_EXCEPTION
            const Time tau = dayCounter.yearFraction(
                                 riskFreeTS->referenceDate(),
                                 calendar.advance(riskFreeTS->referenceDate(),
                                                  optionMaturities[i]));
        const Real fwdPrice = s0->value()*dividendTS->discount(tau)
                            / riskFreeTS->discount(tau);
        const Real strikePrice = fwdPrice * std::exp(-moneyness * volatility
                                                     * std::sqrt(tau));

        options.push_back(boost::shared_ptr<CalibrationHelper>(
                          new HestonModelHelper(optionMaturities[i], calendar,
                                                s0->value(), strikePrice, vol,
                                                riskFreeTS, dividendTS)));
        }
    }

    for (Real sigma = 0.1; sigma < 0.7; sigma += 0.2) {
        const Real v0=0.01;
        const Real kappa=0.2;
        const Real theta=0.02;
        const Real rho=-0.75;

        boost::shared_ptr<HestonProcess> process(
            new HestonProcess(riskFreeTS, dividendTS,
                              s0, v0, kappa, theta, sigma, rho));

        boost::shared_ptr<HestonModel> model(new HestonModel(process));
        boost::shared_ptr<PricingEngine> engine(
                                         new AnalyticHestonEngine(model, 96));

        for (Size i = 0; i < options.size(); ++i)
            options[i]->setPricingEngine(engine);

        LevenbergMarquardt om(1e-8, 1e-8, 1e-8);
        model->calibrate(options, om, EndCriteria(400, 40, 1.0e-8, 
                                                  1.0e-8, 1.0e-8));

        Real tolerance = 3.0e-3;

        if (model->sigma() > tolerance) {
            BOOST_ERROR("Failed to reproduce expected sigma"
                        << "\n    calculated: " << model->sigma()
                        << "\n    expected:   " << 0.0
                        << "\n    tolerance:  " << tolerance);
        }

        if (std::fabs(model->kappa()
                  *(model->theta()-volatility*volatility)) > tolerance) {
            BOOST_ERROR("Failed to reproduce expected theta"
                        << "\n    calculated: " << model->theta()
                        << "\n    expected:   " << volatility*volatility);
        }

        if (std::fabs(model->v0()-volatility*volatility) > tolerance) {
            BOOST_ERROR("Failed to reproduce expected v0"
                        << "\n    calculated: " << model->v0()
                        << "\n    expected:   " << volatility*volatility);
        }
    }
}


void HestonModelTest::testDAXCalibration() {
    /* this example is taken from A. Sepp
       Pricing European-Style Options under Jump Diffusion Processes
       with Stochstic Volatility: Applications of Fourier Transform
       http://math.ut.ee/~spartak/papers/stochjumpvols.pdf
    */

    BOOST_MESSAGE(
             "Testing Heston model calibration using DAX volatility data...");

    SavedSettings backup;

    Date settlementDate(5, July, 2002);
    Settings::instance().evaluationDate() = settlementDate;

    DayCounter dayCounter = Actual365Fixed();
    Calendar calendar = TARGET();

    Integer t[] = { 13, 41, 75, 165, 256, 345, 524, 703 };
    Rate r[] = { 0.0357,0.0349,0.0341,0.0355,0.0359,0.0368,0.0386,0.0401 };

    std::vector<Date> dates;
    std::vector<Rate> rates;
    dates.push_back(settlementDate);
    rates.push_back(0.0357);
    Size i;
    for (i = 0; i < 8; ++i) {
        dates.push_back(settlementDate + t[i]);
        rates.push_back(r[i]);
    }
    // FLOATING_POINT_EXCEPTION
    Handle<YieldTermStructure> riskFreeTS(
                       boost::shared_ptr<YieldTermStructure>(
                                    new ZeroCurve(dates, rates, dayCounter)));

    Handle<YieldTermStructure> dividendTS(
                                   flatRate(settlementDate, 0.0, dayCounter));

    Volatility v[] =
      { 0.6625,0.4875,0.4204,0.3667,0.3431,0.3267,0.3121,0.3121,
        0.6007,0.4543,0.3967,0.3511,0.3279,0.3154,0.2984,0.2921,
        0.5084,0.4221,0.3718,0.3327,0.3155,0.3027,0.2919,0.2889,
        0.4541,0.3869,0.3492,0.3149,0.2963,0.2926,0.2819,0.2800,
        0.4060,0.3607,0.3330,0.2999,0.2887,0.2811,0.2751,0.2775,
        0.3726,0.3396,0.3108,0.2781,0.2788,0.2722,0.2661,0.2686,
        0.3550,0.3277,0.3012,0.2781,0.2781,0.2661,0.2661,0.2681,
        0.3428,0.3209,0.2958,0.2740,0.2688,0.2627,0.2580,0.2620,
        0.3302,0.3062,0.2799,0.2631,0.2573,0.2533,0.2504,0.2544,
        0.3343,0.2959,0.2705,0.2540,0.2504,0.2464,0.2448,0.2462,
        0.3460,0.2845,0.2624,0.2463,0.2425,0.2385,0.2373,0.2422,
        0.3857,0.2860,0.2578,0.2399,0.2357,0.2327,0.2312,0.2351,
        0.3976,0.2860,0.2607,0.2356,0.2297,0.2268,0.2241,0.2320 };

    Handle<Quote> s0(boost::shared_ptr<Quote>(new SimpleQuote(4468.17)));
    Real strike[] = { 3400,3600,3800,4000,4200,4400,
                      4500,4600,4800,5000,5200,5400,5600 };

    std::vector<boost::shared_ptr<CalibrationHelper> > options;

    for (Size s = 0; s < 13; ++s) {
        for (Size m = 0; m < 8; ++m) {
            Handle<Quote> vol(boost::shared_ptr<Quote>(
                                                  new SimpleQuote(v[s*8+m])));

            Period maturity((int)((t[m]+3)/7.), Weeks); // round to weeks
            options.push_back(boost::shared_ptr<CalibrationHelper>(
                        new HestonModelHelper(maturity, calendar,
                                              s0->value(), strike[s], vol,
                                              riskFreeTS, dividendTS, true)));
        }
    }

    const Real v0=0.1;
    const Real kappa=1.0;
    const Real theta=0.1;
    const Real sigma=0.5;
    const Real rho=-0.5;

    boost::shared_ptr<HestonProcess> process(new HestonProcess(
              riskFreeTS, dividendTS, s0, v0, kappa, theta, sigma, rho));

    boost::shared_ptr<HestonModel> model(new HestonModel(process));

    boost::shared_ptr<PricingEngine> engine(
                                         new AnalyticHestonEngine(model, 64));

    for (i = 0; i < options.size(); ++i)
        options[i]->setPricingEngine(engine);

    LevenbergMarquardt om(1e-8, 1e-8, 1e-8);
    model->calibrate(options, om, EndCriteria(400, 40, 1.0e-8, 1.0e-8, 1.0e-8));

    Real sse = 0;
    for (i = 0; i < 13*8; ++i) {
        const Real diff = options[i]->calibrationError()*100.0;
        sse += diff*diff;
    }
    Real expected = 177.2; //see article by A. Sepp.
    if (std::fabs(sse - expected) > 1.0) {
        BOOST_FAIL("Failed to reproduce calibration error"
                   << "\n    calculated: " << sse
                   << "\n    expected:   " << expected);
    }
}

void HestonModelTest::testAnalyticVsBlack() {
    BOOST_MESSAGE("Testing analytic Heston engine against Black formula...");

    SavedSettings backup;

    Date settlementDate = Date::todaysDate();
    Settings::instance().evaluationDate() = settlementDate;
    DayCounter dayCounter = ActualActual();
    Date exerciseDate = settlementDate + 6*Months;

    boost::shared_ptr<StrikedTypePayoff> payoff(
                                     new PlainVanillaPayoff(Option::Put, 30));
    boost::shared_ptr<Exercise> exercise(new EuropeanExercise(exerciseDate));

    Handle<YieldTermStructure> riskFreeTS(flatRate(0.1, dayCounter));
    Handle<YieldTermStructure> dividendTS(flatRate(0.04, dayCounter));

    Handle<Quote> s0(boost::shared_ptr<Quote>(new SimpleQuote(32.0)));

    const Real v0=0.05;
    const Real kappa=5.0;
    const Real theta=0.05;
    const Real sigma=1.0e-4;
    const Real rho=0.0;

    boost::shared_ptr<HestonProcess> process(new HestonProcess(
                   riskFreeTS, dividendTS, s0, v0, kappa, theta, sigma, rho));

    VanillaOption option(payoff, exercise);
    // FLOATING_POINT_EXCEPTION
    boost::shared_ptr<PricingEngine> engine(new AnalyticHestonEngine(
              boost::shared_ptr<HestonModel>(new HestonModel(process)), 144));

    option.setPricingEngine(engine);
    Real calculated = option.NPV();

    Real yearFraction = dayCounter.yearFraction(settlementDate, exerciseDate);
    Real forwardPrice = 32*std::exp((0.1-0.04)*yearFraction);
    Real expected = blackFormula(payoff->optionType(), payoff->strike(),
        forwardPrice, std::sqrt(0.05*yearFraction)) *
                                            std::exp(-0.1*yearFraction);
    Real error = std::fabs(calculated - expected);
    Real tolerance = 2.0e-7;
    if (error > tolerance) {
        BOOST_FAIL("failed to reproduce Black price"
                   << "\n    calculated: " << calculated
                   << "\n    expected:   " << expected
                   << "\n    error:      " << QL_SCIENTIFIC << error);
    }
}


void HestonModelTest::testAnalyticVsCached() {
    BOOST_MESSAGE("Testing analytic Heston engine against cached values...");

    SavedSettings backup;

    Date settlementDate(27, December, 2004);
    Settings::instance().evaluationDate() = settlementDate;
    DayCounter dayCounter = ActualActual();
    Date exerciseDate(28, March, 2005);

    boost::shared_ptr<StrikedTypePayoff> payoff(
                                  new PlainVanillaPayoff(Option::Call, 1.05));
    boost::shared_ptr<Exercise> exercise(new EuropeanExercise(exerciseDate));

    Handle<YieldTermStructure> riskFreeTS(flatRate(0.0225, dayCounter));
    Handle<YieldTermStructure> dividendTS(flatRate(0.02, dayCounter));

    Handle<Quote> s0(boost::shared_ptr<Quote>(new SimpleQuote(1.0)));
    const Real v0 = 0.1;
    const Real kappa = 3.16;
    const Real theta = 0.09;
    const Real sigma = 0.4;
    const Real rho = -0.2;

    boost::shared_ptr<HestonProcess> process(new HestonProcess(
                   riskFreeTS, dividendTS, s0, v0, kappa, theta, sigma, rho));

    VanillaOption option(payoff, exercise);

    boost::shared_ptr<AnalyticHestonEngine> engine(new AnalyticHestonEngine(
               boost::shared_ptr<HestonModel>(new HestonModel(process)), 64));

    option.setPricingEngine(engine);

    Real expected1 = 0.0404774515;
    Real calculated1 = option.NPV();
    Real tolerance = 1.0e-8;

    if (std::fabs(calculated1 - expected1) > tolerance) {
        BOOST_ERROR("Failed to reproduce cached analytic price"
                    << "\n    calculated: " << calculated1
                    << "\n    expected:   " << expected1);
    }


    // reference values from www.wilmott.com, technical forum
    // search for "Heston or VG price check"

    Real K[] = {0.9,1.0,1.1};
    Real expected2[] = { 0.1330371,0.0641016, 0.0270645 };
    Real calculated2[6];

    Size i;
    for (i = 0; i < 6; ++i) {
        Date exerciseDate(8+i/3, September, 2005);

        boost::shared_ptr<StrikedTypePayoff> payoff(
                                new PlainVanillaPayoff(Option::Call, K[i%3]));
        boost::shared_ptr<Exercise> exercise(
                                          new EuropeanExercise(exerciseDate));

        Handle<YieldTermStructure> riskFreeTS(flatRate(0.05, dayCounter));
        Handle<YieldTermStructure> dividendTS(flatRate(0.02, dayCounter));

        Real s = riskFreeTS->discount(0.7)/dividendTS->discount(0.7);
        Handle<Quote> s0(boost::shared_ptr<Quote>(new SimpleQuote(s)));

        boost::shared_ptr<HestonProcess> process(new HestonProcess(
                   riskFreeTS, dividendTS, s0, 0.09, 1.2, 0.08, 1.8, -0.45));

        VanillaOption option(payoff, exercise);

        boost::shared_ptr<PricingEngine> engine(new AnalyticHestonEngine(
                   boost::shared_ptr<HestonModel>(new HestonModel(process))));

        option.setPricingEngine(engine);
        calculated2[i] = option.NPV();
    }

    // we are after the value for T=0.7
    Time t1 = dayCounter.yearFraction(settlementDate, Date(8, September,2005));
    Time t2 = dayCounter.yearFraction(settlementDate, Date(9, September,2005));

    for (i = 0; i < 3; ++i) {
        const Real interpolated =
            calculated2[i]+(calculated2[i+3]-calculated2[i])/(t2-t1)*(0.7-t1);

        if (std::fabs(interpolated - expected2[i]) > 100*tolerance) {
            BOOST_ERROR("Failed to reproduce cached analytic prices:"
                        << "\n    calculated: " << interpolated
                        << "\n    expected:   " << expected2[i] );
        }
    }
}


void HestonModelTest::testMcVsCached() {
    BOOST_MESSAGE(
                "Testing Monte Carlo Heston engine against cached values...");

    SavedSettings backup;

    Date settlementDate(27, December, 2004);
    Settings::instance().evaluationDate() = settlementDate;

    DayCounter dayCounter = ActualActual();
    Date exerciseDate(28, March, 2005);

    boost::shared_ptr<StrikedTypePayoff> payoff(
                                   new PlainVanillaPayoff(Option::Put, 1.05));
    boost::shared_ptr<Exercise> exercise(new EuropeanExercise(exerciseDate));

    Handle<YieldTermStructure> riskFreeTS(flatRate(0.7, dayCounter));
    Handle<YieldTermStructure> dividendTS(flatRate(0.4, dayCounter));

    Handle<Quote> s0(boost::shared_ptr<Quote>(new SimpleQuote(1.05)));

    boost::shared_ptr<HestonProcess> process(new HestonProcess(
                   riskFreeTS, dividendTS, s0, 0.3, 1.16, 0.2, 0.8, 0.8));

    VanillaOption option(payoff, exercise);

    boost::shared_ptr<PricingEngine> engine;
    engine = MakeMCEuropeanHestonEngine<PseudoRandom>(process)
        .withStepsPerYear(91)
        .withAntitheticVariate()
        .withSamples(50000)
        .withSeed(1234);

    option.setPricingEngine(engine);

    Real expected = 0.0632851308977151;
    Real calculated = option.NPV();
    Real errorEstimate = option.errorEstimate();
    Real tolerance = 7.5e-4;

    if (std::fabs(calculated - expected) > 2.34*errorEstimate) {
        BOOST_ERROR("Failed to reproduce cached price"
                    << "\n    calculated: " << calculated
                    << "\n    expected:   " << expected
                    << " +/- " << errorEstimate);
    }

    if (errorEstimate > tolerance) {
        BOOST_ERROR("failed to reproduce error estimate"
                    << "\n    calculated: " << errorEstimate
                    << "\n    expected:   " << tolerance);
    }
}


void HestonModelTest::testKahlJaeckelCase() {
    BOOST_MESSAGE(
          "Testing Monte Carlo Heston engine for the Kahl-Jaeckel example...");

    /* Example taken from Wilmott mag (Sept. 2005).
       "Not-so-complex logarithms in the Heston model",
       Example was also discussed within the Wilmott thread
       "QuantLib code is very high quatlity"
    */

    SavedSettings backup;

    Date settlementDate(30, March, 2007);
    Settings::instance().evaluationDate() = settlementDate;

    DayCounter dayCounter = ActualActual();
    Date exerciseDate(30, March, 2017);

    boost::shared_ptr<StrikedTypePayoff> payoff(
                                   new PlainVanillaPayoff(Option::Call, 200));
    boost::shared_ptr<Exercise> exercise(new EuropeanExercise(exerciseDate));

    Handle<YieldTermStructure> riskFreeTS(flatRate(0.0, dayCounter));
    Handle<YieldTermStructure> dividendTS(flatRate(0.0, dayCounter));

    Handle<Quote> s0(boost::shared_ptr<Quote>(new SimpleQuote(100)));

    boost::shared_ptr<HestonProcess> process(new HestonProcess(
                   riskFreeTS, dividendTS, s0, 0.16, 1.0, 0.16, 2.0, -0.8,
                   HestonProcess::ExactVariance));

    VanillaOption option(payoff, exercise);

    const Real tolerance = 0.1;

    boost::shared_ptr<PricingEngine> engine =
        MakeMCEuropeanHestonEngine<PseudoRandom>(process)
        .withSteps(10)
        .withAntitheticVariate()
        .withTolerance(tolerance)
        .withSeed(1234);
    option.setPricingEngine(engine);

    const Real expected = 4.95212;
    const Real calculated = option.NPV();
    const Real errorEstimate = option.errorEstimate();

    if (std::fabs(calculated - expected) > 2.34*errorEstimate) {
        BOOST_ERROR("Failed to reproduce cached price"
                    << "\n    calculated: " << calculated
                    << "\n    expected:   " << expected
                    << " +/- " << errorEstimate);
    }

    if (errorEstimate > tolerance) {
        BOOST_ERROR("failed to reproduce error estimate"
                    << "\n    calculated: " << errorEstimate
                    << "\n    expected:   " << tolerance);
    }
}

namespace {
    struct HestonParameter {
        Real v0, kappa, theta, sigma, rho; };
}

void HestonModelTest::testDifferentIntegrals() {
    BOOST_MESSAGE(
       "Testing different numerical Heston integration algorithms...");

    SavedSettings backup;

    const Date settlementDate(27, December, 2004);
    Settings::instance().evaluationDate() = settlementDate;

    const DayCounter dayCounter = ActualActual();

    Handle<YieldTermStructure> riskFreeTS(flatRate(0.05, dayCounter));
    Handle<YieldTermStructure> dividendTS(flatRate(0.03, dayCounter));

    const Real strikes[] = { 0.5, 0.7, 1.0, 1.25, 1.5, 2.0 };
    const Integer maturities[] = { 1, 2, 3, 12, 60, 120, 360};
    const Option::Type types[] ={ Option::Put, Option::Call };

    const HestonParameter equityfx      = { 0.07, 2.0, 0.04, 0.55, -0.8 };
    const HestonParameter highCorr      = { 0.07, 1.0, 0.04, 0.55,  0.995 };
    const HestonParameter lowVolOfVol   = { 0.07, 1.0, 0.04, 0.025, -0.75 };
    const HestonParameter highVolOfVol  = { 0.07, 1.0, 0.04, 5.0, -0.75 };
    const HestonParameter kappaEqSigRho = { 0.07, 0.4, 0.04, 0.5, 0.8 };

    std::vector<HestonParameter> params;
    params.push_back(equityfx);
    params.push_back(highCorr);
    params.push_back(lowVolOfVol);
    params.push_back(highVolOfVol);
    params.push_back(kappaEqSigRho);

    const Real tol[] = { 1e-3, 1e-3, 0.2, 0.01, 1e-3 };

    for (std::vector<HestonParameter>::const_iterator iter = params.begin();
         iter != params.end(); ++iter) {

        Handle<Quote> s0(boost::shared_ptr<Quote>(new SimpleQuote(1.0)));
        boost::shared_ptr<HestonProcess> process(new HestonProcess(
            riskFreeTS, dividendTS, 
            s0, iter->v0, iter->kappa, 
            iter->theta, iter->sigma, iter->rho));

        boost::shared_ptr<HestonModel> model(new HestonModel(process));

        boost::shared_ptr<AnalyticHestonEngine> lobattoEngine(
                              new AnalyticHestonEngine(model, 1e-10,
                                                       1000000));
        boost::shared_ptr<AnalyticHestonEngine> laguerreEngine(
                                        new AnalyticHestonEngine(model, 192));
        boost::shared_ptr<AnalyticHestonEngine> legendreEngine(
            new AnalyticHestonEngine(
                model, AnalyticHestonEngine::Gatheral,
                AnalyticHestonEngine::Integration::gaussLegendre(512)));
        boost::shared_ptr<AnalyticHestonEngine> chebyshevEngine(
            new AnalyticHestonEngine(
                model, AnalyticHestonEngine::Gatheral,
                AnalyticHestonEngine::Integration::gaussChebyshev(512)));
        boost::shared_ptr<AnalyticHestonEngine> chebyshev2thEngine(
            new AnalyticHestonEngine(
                model, AnalyticHestonEngine::Gatheral,
                AnalyticHestonEngine::Integration::gaussChebyshev2th(512)));

        Real maxLegendreDiff    = 0.0;
        Real maxChebyshevDiff   = 0.0;
        Real maxChebyshev2thDiff= 0.0;
        Real maxLaguerreDiff    = 0.0;

        for (Size i=0; i < LENGTH(maturities); ++i) {
            boost::shared_ptr<Exercise> exercise(
                new EuropeanExercise(settlementDate
                                     + Period(maturities[i], Months)));

            for (Size j=0; j < LENGTH(strikes); ++j) {
                for (Size k=0; k < LENGTH(types); ++k) {

                    boost::shared_ptr<StrikedTypePayoff> payoff(
                        new PlainVanillaPayoff(types[k], strikes[j]));
                    
                    VanillaOption option(payoff, exercise);

                    option.setPricingEngine(lobattoEngine);
                    const Real lobattoNPV = option.NPV();

                    option.setPricingEngine(laguerreEngine);
                    const Real laguerre = option.NPV();

                    option.setPricingEngine(legendreEngine);
                    const Real legendre = option.NPV();

                    option.setPricingEngine(chebyshevEngine);
                    const Real chebyshev = option.NPV();

                    option.setPricingEngine(chebyshev2thEngine);
                    const Real chebyshev2th = option.NPV();

                    maxLaguerreDiff 
                        = std::max(maxLaguerreDiff,
                                   std::fabs(lobattoNPV-laguerre));
                    maxLegendreDiff
                        = std::max(maxLegendreDiff,
                                   std::fabs(lobattoNPV-legendre));
                    maxChebyshevDiff
                        = std::max(maxChebyshevDiff,
                                   std::fabs(lobattoNPV-chebyshev));
                    maxChebyshev2thDiff
                        = std::max(maxChebyshev2thDiff,
                                   std::fabs(lobattoNPV-chebyshev2th));

                }
            }
        }
        const Real maxDiff = std::max(std::max(
            std::max(maxLaguerreDiff,maxLegendreDiff), 
                                     maxChebyshevDiff), maxChebyshev2thDiff);

        const Real tr = tol[iter - params.begin()];
        if (maxDiff > tr) {
            BOOST_ERROR("Failed to reproduce Heston pricing values "
                        "within given tolerance"
                        << "\n    maxDifference: " << maxDiff
                        << "\n    tolerance:     " << tr);
        }
    }
}

test_suite* HestonModelTest::suite() {
    test_suite* suite = BOOST_TEST_SUITE("Heston model tests");

    // FLOATING_POINT_EXCEPTION
    suite->add(QUANTLIB_TEST_CASE(&HestonModelTest::testBlackCalibration));
    // FLOATING_POINT_EXCEPTION
    suite->add(QUANTLIB_TEST_CASE(&HestonModelTest::testDAXCalibration));
    // FLOATING_POINT_EXCEPTION
    suite->add(QUANTLIB_TEST_CASE(&HestonModelTest::testAnalyticVsBlack));
    suite->add(QUANTLIB_TEST_CASE(&HestonModelTest::testAnalyticVsCached));
    suite->add(QUANTLIB_TEST_CASE(&HestonModelTest::testKahlJaeckelCase));
    suite->add(QUANTLIB_TEST_CASE(&HestonModelTest::testDifferentIntegrals));

    // this passes but takes way too long
    // suite->add(QUANTLIB_TEST_CASE(&HestonModelTest::testEngines));
    //suite->add(QUANTLIB_TEST_CASE(&HestonModelTest::testMcVsCached));
    
    return suite;
}
