/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2004, 2005, 2006, 2007 Ferdinando Ametrano
 Copyright (C) 2004, 2005, 2006, 2007, 2008 StatPro Italia srl
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

#include <ql/types.hpp>
#include <ql/settings.hpp>
#include <ql/version.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/timer.hpp>

/* Use BOOST_MSVC instead of _MSC_VER since some other vendors (Metrowerks,
   for example) also #define _MSC_VER
*/
#ifdef BOOST_MSVC
#  include <ql/auto_link.hpp>
#  define BOOST_LIB_NAME boost_unit_test_framework
#  include <boost/config/auto_link.hpp>
#  undef BOOST_LIB_NAME

/* uncomment the following lines to unmask floating-point exceptions.
   See http://www.wilmott.com/messageview.cfm?catid=10&threadid=9481
*/
//#  include <float.h>
//   namespace { unsigned int u = _controlfp(_EM_INEXACT, _MCW_EM); }

#endif

#include <iostream>
#include <iomanip>
#include "utilities.hpp"

#include "adjointcreditdefaultswaptest.hpp"
#include "adjointzerospreadedtermstructuretest.hpp"
#include "adjointtermstructuretest.hpp"
#include "adjointbermudanswaptiontest.hpp"
#include "adjointeuropeanoptionportfoliotest.hpp"
#include "adjointmatricestest.hpp"
#include "adjointswaptionvolatilitycubetest.hpp"
#include "adjointpathgeneratortest.hpp"
#include "adjointpiecewiseyieldcurvetest.hpp"
#include "adjointmarketmodelcalibrationtest.hpp"
#include "adjointswaptiontest.hpp"
#include "adjointgjrgarchmodeltest.hpp"
#include "adjointgreekstest.hpp"
#include "adjointpiecewiseyieldcurvetest.hpp"
#include "adjointshortratemodelstest.hpp"
#include "adjointblackformulatest.hpp"
#include "adjointbondportfoliotest.hpp"
#include "adjointdistributiontest.hpp"
#include "adjointrealmathtest.hpp"
#include "adjointdefaultprobabilitycurvetest.hpp"
#include "adjointswaptest.hpp"
#include "adjointvariategeneratorstest.hpp"
#include "adjointfraportfoliotest.hpp"

# if defined FIXED_
#include "adjointarraytest.hpp"
# endif

#include "adjointcomplextest.hpp"
#include "adjointfastfouriertransformtest.hpp"
#include "adjointhestonprocesstest.hpp"
#include "adjointbatesmodeltest.hpp"
#include "adjointspecialfunctionstest.hpp"

using namespace boost::unit_test_framework;

namespace {

    boost::timer t;

    void startTimer() { t.restart(); }
    void stopTimer() {
        double seconds = t.elapsed();
        int hours = int(seconds/3600);
        seconds -= hours * 3600;
        int minutes = int(seconds/60);
        seconds -= minutes * 60;
        std::cout << " \nTests completed in ";
        if (hours > 0)
            std::cout << hours << " h ";
        if (hours > 0 || minutes > 0)
            std::cout << minutes << " m ";
        std::cout << std::fixed << std::setprecision(0)
                  << seconds << " s\n" << std::endl;
    }

    void configure() {
        /* if needed, either or both the lines below can be
           uncommented and/or changed to run the test suite with a
           different configuration. In the future, we'll need a
           mechanism that doesn't force us to recompile (possibly a
           couple of command-line flags for the test suite?)
        */

        //QuantLib::Settings::instance().includeReferenceDateCashFlows() = true;
        //QuantLib::Settings::instance().includeTodaysCashFlows() = boost::none;
    }

}

#if defined(QL_ENABLE_SESSIONS)
namespace QuantLib {

    Integer sessionId() { return 0; }

}
#endif

test_suite* init_unit_test_suite(int, char* []) {

    std::string header =
        " Testing "
            #ifdef BOOST_MSVC
            QL_LIB_NAME
            #else
            "QuantLib " QL_VERSION
            #endif
        "\n  QL_NEGATIVE_RATES "
            #ifdef QL_NEGATIVE_RATES
            "       defined"
            #else
            "     undefined"
            #endif
        "\n  QL_EXTRA_SAFETY_CHECKS "
            #ifdef QL_EXTRA_SAFETY_CHECKS
            "  defined"
            #else
            "undefined"
            #endif
        "\n  QL_USE_INDEXED_COUPON "
            #ifdef QL_USE_INDEXED_COUPON
            "   defined"
            #else
            " undefined"
            #endif
         ;
# pragma message ("AdjointArrayTest::suite() was exclude from project in 2015 because it does not buildable.")
    std::string rule = std::string(35, '=');

    BOOST_TEST_MESSAGE(rule);
    BOOST_TEST_MESSAGE(header);
    BOOST_TEST_MESSAGE(rule);
    test_suite* test = BOOST_TEST_SUITE("QuantLib test suite");

    test->add(QUANTLIB_TEST_CASE(startTimer));
    test->add(QUANTLIB_TEST_CASE(configure));

//#   define CL_CERTAIN_TEST AdjointArrayTest

#if defined CL_CERTAIN_TEST
    test->add(CL_CERTAIN_TEST::suite());
#endif

#if !defined CL_ENABLE_BOOST_TEST_ADAPTER && !defined CL_CERTAIN_TEST
    test->add(AdjointBondPortfolioTest::suite());
    test->add(AdjointCreditDefaultSwapTest::suite());
    test->add(AdjointEuropeanOptionPortfolioTest::suite());
    test->add(AdjointGjrgarchModelTest::suite());
    test->add(AdjointFRAPortfolioTest::suite());
    test->add(AdjointMarketModelCalibrationTest::suite());
    test->add(AdjointMatricesTest::suite());
    test->add(AdjointPathGeneratorTest::suite());
    test->add(AdjointPiecewiseYieldCurveTest::suite());
    test->add(AdjointZeroSpreadedTermStructureTest::suite());
    test->add(AdjointShortRateModelsTest::suite());
    test->add(AdjointSwaptionTest::suite());
    test->add(AdjointSwaptionVolatilityCubeTest::suite());
    test->add(AdjointDistributionTest::suite());
    test->add(AdjointRealMathTest::suite());
    test->add(AdjointDefaultProbabilityCurveTest::suite());
    test->add(AdjointSwapTest::suite());
    test->add(AdjointVariateGeneratorsTest::suite());


# if defined FIXED_
    test->add(AdjointArrayTest::suite());
# endif

    test->add(AdjointGreeksTest::suite());

    test->add(AdjointHestonProcessTest::suite());

    // Complex Differentiation test
#if defined CL_TAPE_COMPLEX_ENABLED
    test->add(AdjointComplexTest::suite());
    test->add(AdjointFastFourierTransformTest::suite());
    test->add(AdjointBatesModelTest::suite());
    test->add(AdjointSpecialFunctionsTest::suite());
#endif

    // very slow
    test->add(AdjointBermudanSwaptionTest::suite());
    test->add(AdjointBlackFormulaTest::suite());
    test->add(AdjointTermStructureTest::suite());
#endif

    test->add(QUANTLIB_TEST_CASE(stopTimer));

    return test;
}
