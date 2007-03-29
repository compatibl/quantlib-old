/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Ferdinando Ametrano
 Copyright (C) 2006 Marco Bianchetti
 Copyright (C) 2006 Cristina Duminuco
 Copyright (C) 2006 StatPro Italia srl

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/reference/license.html>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include "curvestates.hpp"
#include "utilities.hpp"

#include <ql/MarketModels/CurveStates/lmmcurvestate.hpp>
#include <ql/MarketModels/CurveStates/coterminalswapcurvestate.hpp>
#include <ql/MarketModels/CurveStates/cmswapcurvestate.hpp>
#include <ql/MarketModels/evolutiondescription.hpp>
#include <ql/MarketModels/DriftComputation/lmmdriftcalculator.hpp>
#include <ql/Utilities/dataformatters.hpp>

#include <ql/Math/matrix.hpp>
#include <ql/schedule.hpp>
#include <ql/DayCounters/simpledaycounter.hpp>
//#include <iostream>
#include <sstream>

#if defined(BOOST_MSVC)
#include <float.h>
//namespace { unsigned int u = _controlfp(_EM_INEXACT, _MCW_EM); }
#endif

using namespace QuantLib;
using namespace boost::unit_test_framework;

QL_BEGIN_TEST_LOCALS(CurveStatesTest)


#define BEGIN(x) (x+0)
#define END(x) (x+LENGTH(x))


Date todaysDate, startDate, endDate;
std::vector<Time> rateTimes, paymentTimes;
std::vector<Real> accruals;
Calendar calendar;
DayCounter dayCounter;
std::vector<Rate> todaysForwards, todaysCoterminalSwapRates;
std::vector<Real> coterminalAnnuity;
Spread displacement;
std::vector<DiscountFactor> todaysDiscounts;

void setup() {
    // Times
    calendar = NullCalendar();
    todaysDate = Settings::instance().evaluationDate();
    //startDate = todaysDate + 5*Years;
    endDate = todaysDate + 10*Years;
    Schedule dates(todaysDate, endDate, Period(Semiannual),
                   calendar, Following, Following, true, false);
    rateTimes = std::vector<Time>(dates.size()-1);
    paymentTimes = std::vector<Time>(rateTimes.size()-1);
    accruals = std::vector<Real>(rateTimes.size()-1);
    dayCounter = SimpleDayCounter();
    for (Size i=1; i<dates.size(); ++i)
        rateTimes[i-1] = dayCounter.yearFraction(todaysDate, dates[i]);
    std::copy(rateTimes.begin()+1, rateTimes.end(), paymentTimes.begin());
    for (Size i=1; i<rateTimes.size(); ++i)
        accruals[i-1] = rateTimes[i] - rateTimes[i-1];

    // Rates & displacement
    todaysForwards = std::vector<Rate>(paymentTimes.size());
    displacement = 0.0;
    for (Size i=0; i<todaysForwards.size(); ++i)
        todaysForwards[i] = 0.03 + 0.0010*i;

    // Discounts
    todaysDiscounts = std::vector<DiscountFactor>(rateTimes.size());
    todaysDiscounts[0] = 0.95;
    for (Size i=1; i<rateTimes.size(); ++i)
        todaysDiscounts[i] = todaysDiscounts[i-1] /
            (1.0+todaysForwards[i-1]*accruals[i-1]);

    // Coterminal swap rates & annuities
    Size N = todaysForwards.size();
    todaysCoterminalSwapRates = std::vector<Rate>(N);
    coterminalAnnuity = std::vector<Real>(N);
    Real floatingLeg = 0.0;
    for (Size i=1; i<=N; ++i) {
        if (i==1) {
            coterminalAnnuity[N-1] = accruals[N-1]*todaysDiscounts[N];
        } else {
            coterminalAnnuity[N-i] = coterminalAnnuity[N-i+1] +
                                     accruals[N-i]*todaysDiscounts[N-i+1];
        }
        floatingLeg = todaysDiscounts[N-i]-todaysDiscounts[N];
        todaysCoterminalSwapRates[N-i] = floatingLeg/coterminalAnnuity[N-i];
    }

    std::vector<Time> evolutionTimes(rateTimes.size()-1);
    std::copy(rateTimes.begin(), rateTimes.end()-1, evolutionTimes.begin());
    EvolutionDescription evolution(rateTimes,evolutionTimes);
    std::vector<Real> rateTaus = evolution.rateTaus();
    std::vector<Size> alive = evolution.firstAliveRate();
}

QL_END_TEST_LOCALS(CurveStatesTest)

void CurveStatesTest::testLMMCurveState() {

    BOOST_MESSAGE("Testing LMMCurveState class"
                  "in a LIBOR market model...");
    QL_TEST_SETUP
}

void CurveStatesTest::testCoterminalSwapCurveState() {

    BOOST_MESSAGE("Testing CoterminalSwapCurveState class"
                  "in a Swap market model...");
    QL_TEST_SETUP
}


void CurveStatesTest::testCMSwapCurveState() {

    BOOST_MESSAGE("Testing CoterminalSwapCurveState class"
                  "in a Swap market model...");
    QL_TEST_SETUP
    Size nbRates = todaysForwards.size();
    Size factors = nbRates;
    Matrix pseudo(nbRates, factors, 0.1);
    std::vector<Spread> displacements(nbRates, .0);
    std::vector<Time> rateTimes(nbRates+1);
    std::vector<Time> taus(nbRates, .5);
    std::vector<Rate> forwards(nbRates, 0.0);

    //std::cout << "rate value:"<< std::endl;

    for (Size i = 0; i < forwards.size(); ++i){
        forwards[i] = static_cast<Time>(i)*.001+.04;

    }

    for (Size i = 0; i < rateTimes.size(); ++i){
        rateTimes[i] = static_cast<Time>(i)*.5;

    }

    //BOOST_MESSAGE( << "Rates\nTime\tValue:"<< std::endl;)
    //for (Size i = 0; i < rateTimes.size()-1; ++i){
    //    std::cout << rateTimes[i+1] << "\t"<<io::rate(forwards[i]) << std::endl;
    //}

    Size numeraire = nbRates;
    Size alive = 0;

    Size spanningFwds = 1;

    CMSMMDriftCalculator cmsDriftcalulator(pseudo, displacements, taus, numeraire,
                                        alive, spanningFwds);

    CMSwapCurveState cmsCs(rateTimes, spanningFwds);
    cmsCs.setOnCMSwapRates(forwards);
    std::vector<Real> cmsDrifts(nbRates);
    cmsDriftcalulator.compute(cmsCs,cmsDrifts);

    LMMDriftCalculator  lmmDriftcalulator(pseudo, displacements, taus, numeraire,
                                        alive);
    LMMCurveState lmmCs(rateTimes);
    lmmCs.setOnForwardRates(forwards);
    std::vector<Real> lmmDrifts(nbRates);


 /*   std::cout << "drifts:"<< std::endl;
    std::cout << "LMM\t\tCMS"<< std::endl;
    for (Size i = 0; i<nbRates; ++i){
         std::cout << lmmDrifts[i] << "\t\t"<< cmsDrifts[i] << std::endl;
    }*/

//    const std::vector<Rate>& dfs = cs.discountRatios();
    //std::cout << "discounts ratios:"<< std::endl;
    //std::cout << "LMM\tCMS"<< std::endl;
    /*for (Size i = 0; i <nbRates; ++i){
        std::cout << lmmCs.discountRatio(i, nbRates) << "\t"<< cmsCs.discountRatio(i, nbRates) << std::endl;
    }*/
}

// --- Call the desired tests
test_suite* CurveStatesTest::suite() {
    test_suite* suite = BOOST_TEST_SUITE("Curve States tests");
    suite->add(BOOST_TEST_CASE(&CurveStatesTest::testLMMCurveState));
    suite->add(BOOST_TEST_CASE(&CurveStatesTest::testCoterminalSwapCurveState));
    suite->add(BOOST_TEST_CASE(&CurveStatesTest::testCMSwapCurveState));
    return suite;
}
