/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2006 Cristina Duminuco
Copyright (C) 2006, 2008 Ferdinando Ametrano
Copyright (C) 2006 Katiuscia Manzoni
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

#include "adjointswaptionvolatilitycubetest.hpp"
#include <test-suite/swaptionvolstructuresutilities.hpp>
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include <ql/indexes/swap/euriborswap.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/volatility/swaption/swaptionvolcube2.hpp>
#include <ql/termstructures/volatility/swaption/swaptionvolcube1.hpp>
#include <ql/termstructures/volatility/swaption/spreadedswaptionvol.hpp>
#include <ql/utilities/dataformatters.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;

namespace {

    struct CommonVars
    {
        CommonVars()
            : vegaWeighedSmileFit(false)
        {
            conventions.setConventions();

            // ATM swaptionvolmatrix
            atm.setMarketData();

            atmVolMatrix = RelinkableHandle<SwaptionVolatilityStructure>(
                boost::shared_ptr<SwaptionVolatilityStructure>(new
                SwaptionVolatilityMatrix(conventions.calendar,
                conventions.optionBdc,
                atm.tenors.options,
                atm.tenors.swaps,
                atm.volsHandle,
                conventions.dayCounter)));
            // Swaptionvolcube
            cube.setMarketData();

            termStructure.linkTo(flatRate(0.05, Actual365Fixed()));

            swapIndexBase = boost::shared_ptr<SwapIndex>(new
                EuriborSwapIsdaFixA(2 * Years, termStructure));
            shortSwapIndexBase = boost::shared_ptr<SwapIndex>(new
                EuriborSwapIsdaFixA(1 * Years, termStructure));
        }

        bool makeAtmVolTest(const SwaptionVolatilityCube& volCube, Real tolerance)
        {
            bool result = true;;
            for (Size i = 0; i < atm.tenors.options.size(); i++) {
                for (Size j = 0; j < atm.tenors.swaps.size(); j++) {
                    Rate strike = volCube.atmStrike(atm.tenors.options[i],
                        atm.tenors.swaps[j]);
                    std::vector<cl::TapeDouble> X(1, strike);

                    //Beginning of tape recording
                    Independent(X);
                    Volatility actVol = volCube.volatility(atm.tenors.options[i],
                        atm.tenors.swaps[j],
                        X[0], true);
                    std::vector<cl::TapeDouble> Y(1);
                    Y[0] = actVol;

                    //End of tape recording. Differentiaion wiil be held with respect to  the independent variables vector
                    cl::TapeFunction<double> f(X, Y);

                    //Differentiation in Forward mode
                    std::vector<double> dY = f.Forward(1, std::vector<double>(1, 1));

                    //Finite differences (forward difference)
                    Real h = 1e-7 * strike;
                    Rate strikeShifted = strike + h;
                    Volatility actVolShifted = volCube.volatility(atm.tenors.options[i],
                        atm.tenors.swaps[j],
                        strikeShifted, true);
                    std::vector<Real> finDiff;
                    finDiff.push_back((actVolShifted - actVol) / h);

                    if (std::abs(dY[0] - finDiff[0]) > tolerance) {
                        result = false;
                        BOOST_ERROR("\ndifferettiation of recovered atm vols failed:"
                            "\nexpiry time     = " << atm.tenors.options[i] <<
                            "\nswap length     = " << atm.tenors.swaps[j] <<
                            "\n atm strike     = " << io::rate(strike) <<
                            "\n actual vol     = " << io::volatility(actVol) <<
                            "\n adjoint res    = " << io::volatility(actVol) <<
                            "\n fin. diff. res = " << io::volatility(actVol) <<
                            "\n  tolerance     = " << tolerance);
                    }
                }
            }
            return result;
        }

        // global data
        SwaptionMarketConventions conventions;
        AtmVolatility atm;
        RelinkableHandle<SwaptionVolatilityStructure> atmVolMatrix;
        VolatilityCube cube;
        RelinkableHandle<YieldTermStructure> termStructure;
        boost::shared_ptr<SwapIndex> swapIndexBase;
        boost::shared_ptr<SwapIndex> shortSwapIndexBase;
        bool vegaWeighedSmileFit;

        // cleanup
        SavedSettings backup;

    };

}


bool AdjointSwaptionVolatilityCubeTest::testAtmVols() {

    BOOST_TEST_MESSAGE("Testing swaption volatility cube (implied vols differentiation)...");

#ifdef CL_TAPE_CPPAD

    CommonVars vars;

    SwaptionVolCube2 volCube(vars.atmVolMatrix,
        vars.cube.tenors.options,
        vars.cube.tenors.swaps,
        vars.cube.strikeSpreads,
        vars.cube.volSpreadsHandle,
        vars.swapIndexBase,
        vars.shortSwapIndexBase,
        vars.vegaWeighedSmileFit);

    Real tolerance = 1.0e-6;
    return vars.makeAtmVolTest(volCube, tolerance);
#endif
    return true;
}


test_suite* AdjointSwaptionVolatilityCubeTest::suite() {
    test_suite* suite = BOOST_TEST_SUITE("Swaption Volatility Cube tests");

    suite->add(QUANTLIB_TEST_CASE(&AdjointSwaptionVolatilityCubeTest::testAtmVols));

    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_swaption_volatility_cube)

BOOST_AUTO_TEST_CASE(testSwaptionVolatilityCubeAtmVols)
{
    BOOST_CHECK(AdjointSwaptionVolatilityCubeTest::testAtmVols());
}

BOOST_AUTO_TEST_SUITE_END()

#endif