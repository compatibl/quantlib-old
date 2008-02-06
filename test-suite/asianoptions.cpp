/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003, 2004 Ferdinando Ametrano
 Copyright (C) 2005, 2007 StatPro Italia srl

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

#include "asianoptions.hpp"
#include "utilities.hpp"
#include <ql/time/daycounters/actual360.hpp>
#include <ql/instruments/asianoption.hpp>
#include <ql/pricingengines/asian/analytic_discr_geom_av_price.hpp>
#include <ql/pricingengines/asian/analytic_cont_geom_av_price.hpp>
#include <ql/pricingengines/asian/mc_discr_geom_av_price.hpp>
#include <ql/pricingengines/asian/mc_discr_arith_av_price.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>
#include <ql/utilities/dataformatters.hpp>
#include <boost/progress.hpp>
#include <map>

using namespace QuantLib;
using namespace boost::unit_test_framework;

#define REPORT_FAILURE(greekName, averageType, \
                       runningAccumulator, pastFixings, \
                       fixingDates, payoff, exercise, s, q, r, today, v, \
                       expected, calculated, tolerance) \
    BOOST_ERROR( \
        exerciseTypeToString(exercise) \
        << " Asian option with " \
        << averageTypeToString(averageType) << " and " \
        << payoffTypeToString(payoff) << " payoff:\n" \
        << "    running variable: " \
        << io::checknull(runningAccumulator) << "\n" \
        << "    past fixings:     " \
        << io::checknull(pastFixings) << "\n" \
        << "    future fixings:   " << fixingDates.size() << "\n" \
        << "    underlying value: " << s << "\n" \
        << "    strike:           " << payoff->strike() << "\n" \
        << "    dividend yield:   " << io::rate(q) << "\n" \
        << "    risk-free rate:   " << io::rate(r) << "\n" \
        << "    reference date:   " << today << "\n" \
        << "    maturity:         " << exercise->lastDate() << "\n" \
        << "    volatility:       " << io::volatility(v) << "\n\n" \
        << "    expected   " << greekName << ": " << expected << "\n" \
        << "    calculated " << greekName << ": " << calculated << "\n"\
        << "    error:            " << std::fabs(expected-calculated) \
        << "\n" \
        << "    tolerance:        " << tolerance);

namespace {

    std::string averageTypeToString(Average::Type averageType) {

        if (averageType == Average::Geometric)
            return "Geometric Averaging";
        else if (averageType == Average::Arithmetic)
            return "Arithmetic Averaging";
        else
            QL_FAIL("unknown averaging");
    }

}


void AsianOptionTest::testAnalyticContinuousGeometricAveragePrice() {

    BOOST_MESSAGE("Testing analytic continuous geometric average-price "
                  "Asians...");

    // data from "Option Pricing Formulas", Haug, pag.96-97

    DayCounter dc = Actual360();
    Date today = Date::todaysDate();

    boost::shared_ptr<SimpleQuote> spot(new SimpleQuote(80.0));
    boost::shared_ptr<SimpleQuote> qRate(new SimpleQuote(-0.03));
    boost::shared_ptr<YieldTermStructure> qTS = flatRate(today, qRate, dc);
    boost::shared_ptr<SimpleQuote> rRate(new SimpleQuote(0.05));
    boost::shared_ptr<YieldTermStructure> rTS = flatRate(today, rRate, dc);
    boost::shared_ptr<SimpleQuote> vol(new SimpleQuote(0.20));
    boost::shared_ptr<BlackVolTermStructure> volTS = flatVol(today, vol, dc);

    boost::shared_ptr<BlackScholesMertonProcess> stochProcess(new
        BlackScholesMertonProcess(Handle<Quote>(spot),
                                  Handle<YieldTermStructure>(qTS),
                                  Handle<YieldTermStructure>(rTS),
                                  Handle<BlackVolTermStructure>(volTS)));

    boost::shared_ptr<PricingEngine> engine(new
            AnalyticContinuousGeometricAveragePriceAsianEngine(stochProcess));

    Average::Type averageType = Average::Geometric;
    Option::Type type = Option::Put;
    Real strike = 85.0;
    Date exerciseDate = today + 90;

    Size pastFixings = Null<Size>();
    Real runningAccumulator = Null<Real>();

    boost::shared_ptr<StrikedTypePayoff> payoff(
                                        new PlainVanillaPayoff(type, strike));

    boost::shared_ptr<Exercise> exercise(new EuropeanExercise(exerciseDate));

    ContinuousAveragingAsianOption option(averageType, payoff, exercise);
    option.setPricingEngine(engine);

    Real calculated = option.NPV();
    Real expected = 4.6922;
    Real tolerance = 1.0e-4;
    if (std::fabs(calculated-expected) > tolerance) {
        REPORT_FAILURE("value", averageType, runningAccumulator, pastFixings,
                       std::vector<Date>(), payoff, exercise, spot->value(),
                       qRate->value(), rRate->value(), today,
                       vol->value(), expected, calculated, tolerance);
    }

    // trying to approximate the continuous version with the discrete version
    runningAccumulator = 1.0;
    pastFixings = 0;
    std::vector<Date> fixingDates(exerciseDate-today+1);
    for (Size i=0; i<fixingDates.size(); i++) {
        fixingDates[i] = today + i;
    }
    boost::shared_ptr<PricingEngine> engine2(new
              AnalyticDiscreteGeometricAveragePriceAsianEngine(stochProcess));
    DiscreteAveragingAsianOption option2(averageType,
                                         runningAccumulator, pastFixings,
                                         fixingDates,
                                         payoff,
                                         exercise);
    option2.setPricingEngine(engine2);

    calculated = option2.NPV();
    tolerance = 3.0e-3;
    if (std::fabs(calculated-expected) > tolerance) {
        REPORT_FAILURE("value", averageType, runningAccumulator, pastFixings,
                       fixingDates, payoff, exercise, spot->value(),
                       qRate->value(), rRate->value(), today,
                       vol->value(), expected, calculated, tolerance);
    }

}


void AsianOptionTest::testAnalyticContinuousGeometricAveragePriceGreeks() {

    BOOST_MESSAGE("Testing analytic continuous geometric average-price Asian "
                  "greeks...");

    SavedSettings backup;

    std::map<std::string,Real> calculated, expected, tolerance;
    tolerance["delta"]  = 1.0e-5;
    tolerance["gamma"]  = 1.0e-5;
    tolerance["theta"]  = 1.0e-5;
    tolerance["rho"]    = 1.0e-5;
    tolerance["divRho"] = 1.0e-5;
    tolerance["vega"]   = 1.0e-5;

    Option::Type types[] = { Option::Call, Option::Put };
    Real underlyings[] = { 100.0 };
    Real strikes[] = { 90.0, 100.0, 110.0 };
    Rate qRates[] = { 0.04, 0.05, 0.06 };
    Rate rRates[] = { 0.01, 0.05, 0.15 };
    Integer lengths[] = { 1, 2 };
    Volatility vols[] = { 0.11, 0.50, 1.20 };

    DayCounter dc = Actual360();
    Date today = Date::todaysDate();
    Settings::instance().evaluationDate() = today;

    boost::shared_ptr<SimpleQuote> spot(new SimpleQuote(0.0));
    boost::shared_ptr<SimpleQuote> qRate(new SimpleQuote(0.0));
    Handle<YieldTermStructure> qTS(flatRate(qRate, dc));
    boost::shared_ptr<SimpleQuote> rRate(new SimpleQuote(0.0));
    Handle<YieldTermStructure> rTS(flatRate(rRate, dc));
    boost::shared_ptr<SimpleQuote> vol(new SimpleQuote(0.0));
    Handle<BlackVolTermStructure> volTS(flatVol(vol, dc));

    boost::shared_ptr<BlackScholesMertonProcess> process(
         new BlackScholesMertonProcess(Handle<Quote>(spot), qTS, rTS, volTS));

    for (Size i=0; i<LENGTH(types); i++) {
      for (Size j=0; j<LENGTH(strikes); j++) {
        for (Size k=0; k<LENGTH(lengths); k++) {

            boost::shared_ptr<EuropeanExercise> maturity(
                              new EuropeanExercise(today + lengths[k]*Years));

            boost::shared_ptr<PlainVanillaPayoff> payoff(
                                new PlainVanillaPayoff(types[i], strikes[j]));

            boost::shared_ptr<PricingEngine> engine(new
                 AnalyticContinuousGeometricAveragePriceAsianEngine(process));

            ContinuousAveragingAsianOption option(Average::Geometric,
                                                  payoff, maturity);
            option.setPricingEngine(engine);

            Size pastFixings = Null<Size>();
            Real runningAverage = Null<Real>();

            for (Size l=0; l<LENGTH(underlyings); l++) {
              for (Size m=0; m<LENGTH(qRates); m++) {
                for (Size n=0; n<LENGTH(rRates); n++) {
                  for (Size p=0; p<LENGTH(vols); p++) {

                      Real u = underlyings[l];
                      Rate q = qRates[m],
                           r = rRates[n];
                      Volatility v = vols[p];
                      spot->setValue(u);
                      qRate->setValue(q);
                      rRate->setValue(r);
                      vol->setValue(v);

                      Real value = option.NPV();
                      calculated["delta"]  = option.delta();
                      calculated["gamma"]  = option.gamma();
                      calculated["theta"]  = option.theta();
                      calculated["rho"]    = option.rho();
                      calculated["divRho"] = option.dividendRho();
                      calculated["vega"]   = option.vega();

                      if (value > spot->value()*1.0e-5) {
                          // perturb spot and get delta and gamma
                          Real du = u*1.0e-4;
                          spot->setValue(u+du);
                          Real value_p = option.NPV(),
                               delta_p = option.delta();
                          spot->setValue(u-du);
                          Real value_m = option.NPV(),
                               delta_m = option.delta();
                          spot->setValue(u);
                          expected["delta"] = (value_p - value_m)/(2*du);
                          expected["gamma"] = (delta_p - delta_m)/(2*du);

                          // perturb rates and get rho and dividend rho
                          Spread dr = r*1.0e-4;
                          rRate->setValue(r+dr);
                          value_p = option.NPV();
                          rRate->setValue(r-dr);
                          value_m = option.NPV();
                          rRate->setValue(r);
                          expected["rho"] = (value_p - value_m)/(2*dr);

                          Spread dq = q*1.0e-4;
                          qRate->setValue(q+dq);
                          value_p = option.NPV();
                          qRate->setValue(q-dq);
                          value_m = option.NPV();
                          qRate->setValue(q);
                          expected["divRho"] = (value_p - value_m)/(2*dq);

                          // perturb volatility and get vega
                          Volatility dv = v*1.0e-4;
                          vol->setValue(v+dv);
                          value_p = option.NPV();
                          vol->setValue(v-dv);
                          value_m = option.NPV();
                          vol->setValue(v);
                          expected["vega"] = (value_p - value_m)/(2*dv);

                          // perturb date and get theta
                          Time dT = dc.yearFraction(today-1, today+1);
                          Settings::instance().evaluationDate() = today-1;
                          value_m = option.NPV();
                          Settings::instance().evaluationDate() = today+1;
                          value_p = option.NPV();
                          Settings::instance().evaluationDate() = today;
                          expected["theta"] = (value_p - value_m)/dT;

                          // compare
                          std::map<std::string,Real>::iterator it;
                          for (it = calculated.begin();
                               it != calculated.end(); ++it) {
                              std::string greek = it->first;
                              Real expct = expected  [greek],
                                   calcl = calculated[greek],
                                   tol   = tolerance [greek];
                              Real error = relativeError(expct,calcl,u);
                              if (error>tol) {
                                  REPORT_FAILURE(greek, Average::Geometric,
                                                 runningAverage, pastFixings,
                                                 std::vector<Date>(),
                                                 payoff, maturity,
                                                 u, q, r, today, v,
                                                 expct, calcl, tol);
                              }
                          }
                      }
                  }
                }
              }
            }
        }
      }
    }
}


void AsianOptionTest::testAnalyticDiscreteGeometricAveragePrice() {

    BOOST_MESSAGE("Testing analytic discrete geometric average-price "
                  "Asians...");

    // data from "Implementing Derivatives Model",
    // Clewlow, Strickland, p.118-123

    DayCounter dc = Actual360();
    Date today = Date::todaysDate();

    boost::shared_ptr<SimpleQuote> spot(new SimpleQuote(100.0));
    boost::shared_ptr<SimpleQuote> qRate(new SimpleQuote(0.03));
    boost::shared_ptr<YieldTermStructure> qTS = flatRate(today, qRate, dc);
    boost::shared_ptr<SimpleQuote> rRate(new SimpleQuote(0.06));
    boost::shared_ptr<YieldTermStructure> rTS = flatRate(today, rRate, dc);
    boost::shared_ptr<SimpleQuote> vol(new SimpleQuote(0.20));
    boost::shared_ptr<BlackVolTermStructure> volTS = flatVol(today, vol, dc);

    boost::shared_ptr<BlackScholesMertonProcess> stochProcess(new
        BlackScholesMertonProcess(Handle<Quote>(spot),
                                  Handle<YieldTermStructure>(qTS),
                                  Handle<YieldTermStructure>(rTS),
                                  Handle<BlackVolTermStructure>(volTS)));

    boost::shared_ptr<PricingEngine> engine(
          new AnalyticDiscreteGeometricAveragePriceAsianEngine(stochProcess));

    Average::Type averageType = Average::Geometric;
    Real runningAccumulator = 1.0;
    Size pastFixings = 0;
    Size futureFixings = 10;
    Option::Type type = Option::Call;
    Real strike = 100.0;
    boost::shared_ptr<StrikedTypePayoff> payoff(
                                        new PlainVanillaPayoff(type, strike));

    Date exerciseDate = today + 360;
    boost::shared_ptr<Exercise> exercise(new EuropeanExercise(exerciseDate));

    std::vector<Date> fixingDates(futureFixings);
    Integer dt = Integer(360/futureFixings+0.5);
    fixingDates[0] = today + dt;
    for (Size j=1; j<futureFixings; j++)
        fixingDates[j] = fixingDates[j-1] + dt;

    DiscreteAveragingAsianOption option(averageType, runningAccumulator,
                                        pastFixings, fixingDates,
                                        payoff, exercise);
    option.setPricingEngine(engine);

    Real calculated = option.NPV();
    Real expected = 5.3425606635;
    Real tolerance = 1e-10;
    if (std::fabs(calculated-expected) > tolerance) {
        REPORT_FAILURE("value", averageType, runningAccumulator, pastFixings,
                       fixingDates, payoff, exercise, spot->value(),
                       qRate->value(), rRate->value(), today,
                       vol->value(), expected, calculated, tolerance);
    }
}


void AsianOptionTest::testMCDiscreteGeometricAveragePrice() {

    BOOST_MESSAGE("Testing Monte Carlo discrete geometric average-price "
                  "Asians...");

    // data from "Implementing Derivatives Model",
    // Clewlow, Strickland, p.118-123

    DayCounter dc = Actual360();
    Date today = Date::todaysDate();

    boost::shared_ptr<SimpleQuote> spot(new SimpleQuote(100.0));
    boost::shared_ptr<SimpleQuote> qRate(new SimpleQuote(0.03));
    boost::shared_ptr<YieldTermStructure> qTS = flatRate(today, qRate, dc);
    boost::shared_ptr<SimpleQuote> rRate(new SimpleQuote(0.06));
    boost::shared_ptr<YieldTermStructure> rTS = flatRate(today, rRate, dc);
    boost::shared_ptr<SimpleQuote> vol(new SimpleQuote(0.20));
    boost::shared_ptr<BlackVolTermStructure> volTS = flatVol(today, vol, dc);

    boost::shared_ptr<BlackScholesMertonProcess> stochProcess(new
        BlackScholesMertonProcess(Handle<Quote>(spot),
                                  Handle<YieldTermStructure>(qTS),
                                  Handle<YieldTermStructure>(rTS),
                                  Handle<BlackVolTermStructure>(volTS)));

    Real tolerance = 4.0e-3;

    boost::shared_ptr<PricingEngine> engine =
        MakeMCDiscreteGeometricAPEngine<LowDiscrepancy>(stochProcess)
        .withStepsPerYear(1)
        .withSamples(8191);

    Average::Type averageType = Average::Geometric;
    Real runningAccumulator = 1.0;
    Size pastFixings = 0;
    Size futureFixings = 10;
    Option::Type type = Option::Call;
    Real strike = 100.0;
    boost::shared_ptr<StrikedTypePayoff> payoff(
                                        new PlainVanillaPayoff(type, strike));

    Date exerciseDate = today + 360;
    boost::shared_ptr<Exercise> exercise(new EuropeanExercise(exerciseDate));

    std::vector<Date> fixingDates(futureFixings);
    Integer dt = Integer(360/futureFixings+0.5);
    fixingDates[0] = today + dt;
    for (Size j=1; j<futureFixings; j++)
        fixingDates[j] = fixingDates[j-1] + dt;

    DiscreteAveragingAsianOption option(averageType, runningAccumulator,
                                        pastFixings, fixingDates,
                                        payoff, exercise);
    option.setPricingEngine(engine);

    Real calculated = option.NPV();

    boost::shared_ptr<PricingEngine> engine2(
          new AnalyticDiscreteGeometricAveragePriceAsianEngine(stochProcess));
    option.setPricingEngine(engine2);
    Real expected = option.NPV();

    if (std::fabs(calculated-expected) > tolerance) {
        REPORT_FAILURE("value", averageType, runningAccumulator, pastFixings,
                       fixingDates, payoff, exercise, spot->value(),
                       qRate->value(), rRate->value(), today,
                       vol->value(), expected, calculated, tolerance);
    }
}


namespace {

    struct DiscreteAverageData {
        Option::Type type;
        Real underlying;
        Real strike;
        Rate dividendYield;
        Rate riskFreeRate;
        Time first;
        Time length;
        Size fixings;
        Volatility volatility;
        bool controlVariate;
        Real result;
    };

}


void AsianOptionTest::testMCDiscreteArithmeticAveragePrice() {

    BOOST_MESSAGE("Testing Monte Carlo discrete arithmetic average-price "
                  "Asians...");

    QL_TEST_START_TIMING

    // data from "Asian Option", Levy, 1997
    // in "Exotic Options: The State of the Art",
    // edited by Clewlow, Strickland
    DiscreteAverageData cases4[] = {
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 0.0, 11.0/12.0, 2,
          0.13, true, 1.3942835683 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 0.0, 11.0/12.0, 4,
          0.13, true, 1.5852442983 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 0.0, 11.0/12.0, 8,
          0.13, true, 1.66970673 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 0.0, 11.0/12.0, 12,
          0.13, true, 1.6980019214 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 0.0, 11.0/12.0, 26,
          0.13, true, 1.7255070456 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 0.0, 11.0/12.0, 52,
          0.13, true, 1.7401553533 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 0.0, 11.0/12.0, 100,
          0.13, true, 1.7478303712 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 0.0, 11.0/12.0, 250,
          0.13, true, 1.7490291943 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 0.0, 11.0/12.0, 500,
          0.13, true, 1.7515113291 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 0.0, 11.0/12.0, 1000,
          0.13, true, 1.7537344885 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0, 11.0/12.0, 2,
          0.13, true, 1.8496053697 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0, 11.0/12.0, 4,
          0.13, true, 2.0111495205 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0, 11.0/12.0, 8,
          0.13, true, 2.0852138818 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0, 11.0/12.0, 12,
          0.13, true, 2.1105094397 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0, 11.0/12.0, 26,
          0.13, true, 2.1346526695 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0, 11.0/12.0, 52,
          0.13, true, 2.147489651 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0, 11.0/12.0, 100,
          0.13, true, 2.154728109 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0, 11.0/12.0, 250,
          0.13, true, 2.1564276565 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0, 11.0/12.0, 500,
          0.13, true, 2.1594238588 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0, 11.0/12.0, 1000,
          0.13, true, 2.1595367326 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0, 11.0/12.0, 2,
          0.13, true, 2.63315092584 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0, 11.0/12.0, 4,
          0.13, true, 2.76723962361 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0, 11.0/12.0, 8,
          0.13, true, 2.83124836881 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0, 11.0/12.0, 12,
          0.13, true, 2.84290301412 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0, 11.0/12.0, 26,
          0.13, true, 2.88179560417 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0, 11.0/12.0, 52,
          0.13, true, 2.88447044543 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0, 11.0/12.0, 100,
          0.13, true, 2.89985329603 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0, 11.0/12.0, 250,
          0.13, true, 2.90047296063 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0, 11.0/12.0, 500,
          0.13, true, 2.89813412160 },
        { Option::Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0, 11.0/12.0, 1000,
          0.13, true, 2.89703362437 }
   };

    DayCounter dc = Actual360();
    Date today = Date::todaysDate();

    boost::shared_ptr<SimpleQuote> spot(new SimpleQuote(100.0));
    boost::shared_ptr<SimpleQuote> qRate(new SimpleQuote(0.03));
    boost::shared_ptr<YieldTermStructure> qTS = flatRate(today, qRate, dc);
    boost::shared_ptr<SimpleQuote> rRate(new SimpleQuote(0.06));
    boost::shared_ptr<YieldTermStructure> rTS = flatRate(today, rRate, dc);
    boost::shared_ptr<SimpleQuote> vol(new SimpleQuote(0.20));
    boost::shared_ptr<BlackVolTermStructure> volTS = flatVol(today, vol, dc);



    Average::Type averageType = Average::Arithmetic;
    Real runningSum = 0.0;
    Size pastFixings = 0;
    for (Size l=0; l<LENGTH(cases4); l++) {

        boost::shared_ptr<StrikedTypePayoff> payoff(new
            PlainVanillaPayoff(cases4[l].type, cases4[l].strike));

        Time dt = cases4[l].length/(cases4[l].fixings-1);
        std::vector<Time> timeIncrements(cases4[l].fixings);
        std::vector<Date> fixingDates(cases4[l].fixings);
        timeIncrements[0] = cases4[l].first;
        fixingDates[0] = today + Integer(timeIncrements[0]*360+0.5);
        for (Size i=1; i<cases4[l].fixings; i++) {
            timeIncrements[i] = i*dt + cases4[l].first;
            fixingDates[i] = today + Integer(timeIncrements[i]*360+0.5);
        }
        boost::shared_ptr<Exercise> exercise(new
            EuropeanExercise(fixingDates[cases4[l].fixings-1]));

        spot ->setValue(cases4[l].underlying);
        qRate->setValue(cases4[l].dividendYield);
        rRate->setValue(cases4[l].riskFreeRate);
        vol  ->setValue(cases4[l].volatility);

        boost::shared_ptr<BlackScholesMertonProcess> stochProcess(new
            BlackScholesMertonProcess(Handle<Quote>(spot),
                                      Handle<YieldTermStructure>(qTS),
                                      Handle<YieldTermStructure>(rTS),
                                      Handle<BlackVolTermStructure>(volTS)));


        boost::shared_ptr<PricingEngine> engine =
            MakeMCDiscreteArithmeticAPEngine<LowDiscrepancy>(stochProcess)
            .withStepsPerYear(1)
            .withSamples(2047)
            .withControlVariate();

        DiscreteAveragingAsianOption option(averageType, runningSum,
                                            pastFixings, fixingDates,
                                            payoff, exercise);
        option.setPricingEngine(engine);

        Real calculated = option.NPV();
        Real expected = cases4[l].result;
        Real tolerance = 2.0e-2;
        if (std::fabs(calculated-expected) > tolerance) {
            REPORT_FAILURE("value", averageType, runningSum, pastFixings,
                        fixingDates, payoff, exercise, spot->value(),
                        qRate->value(), rRate->value(), today,
                        vol->value(), expected, calculated, tolerance);
        }
    }

}

void AsianOptionTest::testAnalyticDiscreteGeometricAveragePriceGreeks() {

    BOOST_MESSAGE("Testing discrete-averaging geometric Asian greeks...");

    SavedSettings backup;

    std::map<std::string,Real> calculated, expected, tolerance;
    tolerance["delta"]  = 1.0e-5;
    tolerance["gamma"]  = 1.0e-5;
    tolerance["theta"]  = 1.0e-5;
    tolerance["rho"]    = 1.0e-5;
    tolerance["divRho"] = 1.0e-5;
    tolerance["vega"]   = 1.0e-5;

    Option::Type types[] = { Option::Call, Option::Put };
    Real underlyings[] = { 100.0 };
    Real strikes[] = { 90.0, 100.0, 110.0 };
    Rate qRates[] = { 0.04, 0.05, 0.06 };
    Rate rRates[] = { 0.01, 0.05, 0.15 };
    Integer lengths[] = { 1, 2 };
    Volatility vols[] = { 0.11, 0.50, 1.20 };

    DayCounter dc = Actual360();
    Date today = Date::todaysDate();
    Settings::instance().evaluationDate() = today;

    boost::shared_ptr<SimpleQuote> spot(new SimpleQuote(0.0));
    boost::shared_ptr<SimpleQuote> qRate(new SimpleQuote(0.0));
    Handle<YieldTermStructure> qTS(flatRate(qRate, dc));
    boost::shared_ptr<SimpleQuote> rRate(new SimpleQuote(0.0));
    Handle<YieldTermStructure> rTS(flatRate(rRate, dc));
    boost::shared_ptr<SimpleQuote> vol(new SimpleQuote(0.0));
    Handle<BlackVolTermStructure> volTS(flatVol(vol, dc));

    boost::shared_ptr<BlackScholesMertonProcess> process(
         new BlackScholesMertonProcess(Handle<Quote>(spot), qTS, rTS, volTS));

    for (Size i=0; i<LENGTH(types); i++) {
      for (Size j=0; j<LENGTH(strikes); j++) {
        for (Size k=0; k<LENGTH(lengths); k++) {

            boost::shared_ptr<EuropeanExercise> maturity(
                              new EuropeanExercise(today + lengths[k]*Years));

            boost::shared_ptr<PlainVanillaPayoff> payoff(
                                new PlainVanillaPayoff(types[i], strikes[j]));

            Real runningAverage = 120;
            Size pastFixings = 1;

            std::vector<Date> fixingDates;
            for (Date d = today + 3*Months;
                      d <= maturity->lastDate();
                      d += 3*Months)
                fixingDates.push_back(d);


            boost::shared_ptr<PricingEngine> engine(
               new AnalyticDiscreteGeometricAveragePriceAsianEngine(process));

            DiscreteAveragingAsianOption option(Average::Geometric,
                                                runningAverage, pastFixings,
                                                fixingDates, payoff, maturity);
            option.setPricingEngine(engine);

            for (Size l=0; l<LENGTH(underlyings); l++) {
              for (Size m=0; m<LENGTH(qRates); m++) {
                for (Size n=0; n<LENGTH(rRates); n++) {
                  for (Size p=0; p<LENGTH(vols); p++) {

                      Real u = underlyings[l];
                      Rate q = qRates[m],
                           r = rRates[n];
                      Volatility v = vols[p];
                      spot->setValue(u);
                      qRate->setValue(q);
                      rRate->setValue(r);
                      vol->setValue(v);

                      Real value = option.NPV();
                      calculated["delta"]  = option.delta();
                      calculated["gamma"]  = option.gamma();
                      calculated["theta"]  = option.theta();
                      calculated["rho"]    = option.rho();
                      calculated["divRho"] = option.dividendRho();
                      calculated["vega"]   = option.vega();

                      if (value > spot->value()*1.0e-5) {
                          // perturb spot and get delta and gamma
                          Real du = u*1.0e-4;
                          spot->setValue(u+du);
                          Real value_p = option.NPV(),
                               delta_p = option.delta();
                          spot->setValue(u-du);
                          Real value_m = option.NPV(),
                               delta_m = option.delta();
                          spot->setValue(u);
                          expected["delta"] = (value_p - value_m)/(2*du);
                          expected["gamma"] = (delta_p - delta_m)/(2*du);

                          // perturb rates and get rho and dividend rho
                          Spread dr = r*1.0e-4;
                          rRate->setValue(r+dr);
                          value_p = option.NPV();
                          rRate->setValue(r-dr);
                          value_m = option.NPV();
                          rRate->setValue(r);
                          expected["rho"] = (value_p - value_m)/(2*dr);

                          Spread dq = q*1.0e-4;
                          qRate->setValue(q+dq);
                          value_p = option.NPV();
                          qRate->setValue(q-dq);
                          value_m = option.NPV();
                          qRate->setValue(q);
                          expected["divRho"] = (value_p - value_m)/(2*dq);

                          // perturb volatility and get vega
                          Volatility dv = v*1.0e-4;
                          vol->setValue(v+dv);
                          value_p = option.NPV();
                          vol->setValue(v-dv);
                          value_m = option.NPV();
                          vol->setValue(v);
                          expected["vega"] = (value_p - value_m)/(2*dv);

                          // perturb date and get theta
                          Time dT = dc.yearFraction(today-1, today+1);
                          Settings::instance().evaluationDate() = today-1;
                          value_m = option.NPV();
                          Settings::instance().evaluationDate() = today+1;
                          value_p = option.NPV();
                          Settings::instance().evaluationDate() = today;
                          expected["theta"] = (value_p - value_m)/dT;

                          // compare
                          std::map<std::string,Real>::iterator it;
                          for (it = calculated.begin();
                               it != calculated.end(); ++it) {
                              std::string greek = it->first;
                              Real expct = expected  [greek],
                                   calcl = calculated[greek],
                                   tol   = tolerance [greek];
                              Real error = relativeError(expct,calcl,u);
                              if (error>tol) {
                                  REPORT_FAILURE(greek, Average::Geometric,
                                                 runningAverage, pastFixings,
                                                 std::vector<Date>(),
                                                 payoff, maturity,
                                                 u, q, r, today, v,
                                                 expct, calcl, tol);
                              }
                          }
                      }
                  }
                }
              }
            }
        }
      }
    }
}

test_suite* AsianOptionTest::suite() {
    test_suite* suite = BOOST_TEST_SUITE("Asian option tests");

    suite->add(BOOST_TEST_CASE(
        &AsianOptionTest::testAnalyticContinuousGeometricAveragePrice));
    suite->add(BOOST_TEST_CASE(
        &AsianOptionTest::testAnalyticContinuousGeometricAveragePriceGreeks));
    suite->add(BOOST_TEST_CASE(
        &AsianOptionTest::testAnalyticDiscreteGeometricAveragePrice));
    suite->add(BOOST_TEST_CASE(
        &AsianOptionTest::testMCDiscreteGeometricAveragePrice));
    suite->add(BOOST_TEST_CASE(
        &AsianOptionTest::testMCDiscreteArithmeticAveragePrice));
    suite->add(BOOST_TEST_CASE(
        &AsianOptionTest::testAnalyticDiscreteGeometricAveragePriceGreeks));

    return suite;
}
