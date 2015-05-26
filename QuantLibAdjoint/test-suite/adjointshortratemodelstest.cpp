/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2007 Marco Bianchetti
Copyright (C) 2007 Giorgio Facchinetti
Copyright (C) 2006 Chiara Fornarola
Copyright (C) 2005 StatPro Italia srl
Copyright (C) 2013 Peter Caspers
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

// Based on shortratemodels.cpp file from Quantlib/test-suite.

#include "adjointshortratemodelstest.hpp"
#include <test-suite/shortratemodels.hpp>
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include <ql/quantlib.hpp>
#include <ql/models/shortrate/onefactormodels/hullwhite.hpp>
#include <ql/models/shortrate/calibrationhelpers/swaptionhelper.hpp>
#include <ql/pricingengines/swaption/jamshidianswaptionengine.hpp>
#include <ql/pricingengines/swap/treeswapengine.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>
#include <ql/indexes/ibor/euribor.hpp>
#include <ql/indexes/indexmanager.hpp>
#include <ql/math/optimization/simplex.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
#include <ql/termstructures/yield/discountcurve.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/thirty360.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/schedule.hpp>
#include <ql/quotes/simplequote.hpp>
#include <boost/timer.hpp>


using namespace QuantLib;
using namespace boost::unit_test_framework;

// Hull and White in 1990 proposed the following model for the short rate:
// $dr(t)=(\theta(t)-a(t)*r(t))dt +\sigma(t) dW_t$, where
// a(t) is a mean reversion and \sigma(t) is a volatility
// (a, r and $\sigma$ are assumed to be constants in the testing  model);
// The goal of approximating algorithm is to  find the values of  a and  $\sigma$ in
// the  Hull-White model that best fit the market prices of the swaptions 
// or find the minimum: $\sum_{i=1}^N (MarketPrice_i-PriceHW_i)^2$, where:
// N is a number of swaptions, $PriceHW= PriceHW(a,\sigma))$ is a swaption prices vector,
// calculated using Hull-White model, MarketPrice is market swaption prices vector.
// Black-Sholes implied volatilities for swaptions are provided in
// std::vector<cl::TapeDouble> vol_ vector of independent variables. 
// Calculated using the Hull-White calibration  model   $\sigma$
// (see $model\righarrow \sigma()$ in the code below) is differentiated 
// with respect to input volatilities vector.
// Model $\sigma$ dependence on the first element of input volatilties vector is plotted.

namespace {

    struct CalibrationData {
        Integer start_;
        Integer length_;
        Real volatility_;
    };

    struct CommonVars
    {
        // setup
        CommonVars(Size size) : today_(15, February, 2002)
            , settlement_(19, February, 2002)
            , termStructure_(flatRate(settlement_, 0.04875825, Actual365Fixed()))
            , model_(new HullWhite(termStructure_))
            , index_(new Euribor6M(termStructure_))
            , engine_(new JamshidianSwaptionEngine(model_))
            // Size of volatilities vector
            , vol_(size)
            // Step in finite-differences method.
            , h_(1.0e-3)
            , relTolerance_(1e-1)
            , absTolerance_(1e-4)
        {
            Settings::instance().evaluationDate() = today_;
            initValues();
        }

        void initValues()
        {
            std::vector<cl::TapeDouble>::iterator it = vol_.begin();
            for (Size i = 0; it != vol_.end(); it++, i++)
                *it = 0.1148 - 0.004*i;
        }

        // common data
        Date today_;
        Date settlement_;
        Handle<YieldTermStructure> termStructure_;
        boost::shared_ptr<HullWhite> model_;
        boost::shared_ptr<IborIndex> index_;
        boost::shared_ptr<PricingEngine> engine_;
        std::vector<boost::shared_ptr<CalibrationHelper> > swaptions_;
        std::vector<CalibrationData> data_;

        //Adjoint routine vectors
        //Vector of independent variables
        std::vector<cl::TapeDouble> vol_;

        //finite-differences scheme
        double h_;
        double relTolerance_;
        double absTolerance_;

        //results
        std::vector<PerformanceTime> performanceTime_;
        std::vector<AdjointTime> adjointTime_;
    };

    void calibrationData(std::vector<cl::TapeDouble>& vol, std::vector<CalibrationData>& data)
    {
        Size pos = 0;
        Size size = vol.size();
        data.clear();
        for (std::vector<cl::TapeDouble>::iterator it = vol.begin(); it != vol.end(); it++, pos++)
        {
            data.push_back(CalibrationData{ pos + 1, size - pos, *it });
        }
    }
    
    // Calibration type 1: Testing Hull-White calibration against cached values using swaptions with start delay
    // type 2: Testing Hull-White calibration with fixed reversion against cached values
 
    Real calibrate(Handle<YieldTermStructure>& termStructure
        , boost::shared_ptr<HullWhite>& model
        , boost::shared_ptr<IborIndex>& index
        , boost::shared_ptr<PricingEngine>& engine
        , std::vector<boost::shared_ptr<CalibrationHelper> >& swaptions
        , std::vector<CalibrationData>& data
        , Size calibration_type
        , Size sizeof_indep)
    {
        swaptions.clear();
        for (Size i = 0; i < sizeof_indep; i++) {
            boost::shared_ptr<Quote> volatil(new SimpleQuote(data[i].volatility_));
            boost::shared_ptr<CalibrationHelper> helper(
                new SwaptionHelper(Period(data[i].start_, Years),
                Period(data[i].length_, Years),
                Handle<Quote>(volatil),
                index,
                Period(1, Years), Thirty360(),
                Actual360(), termStructure));
            helper->setPricingEngine(engine);
            swaptions.push_back(helper);
        }

        // Set up the optimization problem
        // Real simplexLambda = 0.1;
        // Simplex optimizationMethod(simplexLambda);
        LevenbergMarquardt optimizationMethod(1.0e-8, 1.0e-8, 1.0e-8);
        EndCriteria endCriteria(1000, 100, 1e-6, 1e-8, 1e-8);

        //Optimize
        switch (calibration_type)
        {
        case 1:
            model->calibrate(swaptions, optimizationMethod, endCriteria);
            break;
        case 2:
            model->calibrate(swaptions, optimizationMethod, endCriteria, Constraint(), std::vector<Real>(),
                HullWhite::FixedReversion());
            // The difference is in the choice of  HullWhite::FixedReversion() to calibrate the model below
        }

        return  model->sigma();
    }

    double finiteDiff(Handle<YieldTermStructure>& termStructure
        , boost::shared_ptr<HullWhite>& model
        , boost::shared_ptr<IborIndex>& index
        , boost::shared_ptr<PricingEngine>& engine
        , std::vector<boost::shared_ptr<CalibrationHelper> >& swaptions
        , std::vector<CalibrationData>& data
        , Size calibration_type
        , Size  sizeIndep
        , double h
        , std::vector<Real>& sfFinite
        , cl::AdjointTestOutput& out)
    {
        sfFinite.resize(sizeIndep);
        out.log() << "Start differentiation using finite-differences:\t " << currentTime() << std::endl;

        //Finite differences (central formula)
        boost::timer timerDiff;

        std::vector<CalibrationData>::iterator it;
        std::vector<Real>::iterator it_fin;
        for (it = data.begin(), it_fin = sfFinite.begin(); it != data.end(); it++, it_fin++)
        {
            (*it).volatility_ -= h;
            calibrate(termStructure, model, index, engine, swaptions, data, calibration_type, sizeIndep);
            cl::TapeDouble stepbackward = model->sigma();
            (*it).volatility_ += 2 * h;
            calibrate(termStructure, model, index, engine, swaptions, data,  calibration_type, sizeIndep);
            cl::TapeDouble stepforward = model->sigma();
            *it_fin = (stepforward - stepbackward) / (2 * h);
            (*it).volatility_ -= h;
        }
        double timeCalculated_ = timerDiff.elapsed();
        out.log() << "Time for differentiation using finite-differences:\t " << timeCalculated_ << std::endl;
        // Return time  to calculate all derivatives in sfFinite vector.
        return timeCalculated_;
    }
}

struct VolatilDependence
{
    static std::deque<std::string > get_columns()
    {
        static std::deque<std::string > columns =
        {
            "Volatility", "ModelSigma"
        };

        return columns;
    }

    template <typename stream_type>
    friend inline stream_type&
        operator << (stream_type& stm, VolatilDependence& v)
    {
            stm << v.inputVolatil_
                << ";" << v.modelSigma_ << std::endl;

            return stm;
        }

    Real inputVolatil_;
    Real modelSigma_;
};

void rateDependencePlot(CommonVars& vars, Size size, Size calibration_type, cl::AdjointTestOutput& output)
{
    std::vector<VolatilDependence> volatilDependence_;

    //Init range for the  rate and output portfolio price
    for (Size i = 0; i < size; i++)
    {
        Real volatil_ = 0.06 + 0.0003*i;
        vars.vol_[0] = volatil_;
        calibrationData(vars.vol_, vars.data_);
        Real sigma = calibrate(vars.termStructure_
            , vars.model_
            , vars.index_
            , vars.engine_
            , vars.swaptions_
            , vars.data_
            , calibration_type
            , vars.vol_.size());
        volatilDependence_.push_back(VolatilDependence{ volatil_, sigma });
    }

    //Reset vol-[0] to initial value
    vars.vol_[0] = 0.1148;
    output << volatilDependence_;
    output.log() << "Plot Sigma on volatility dependence successfully generated" << std::endl;
}

bool testHullWhite(CommonVars& vars, Size size, Size type)
{
    // Vector to store derivatives calculated in Forward, Reverse modes and using finite-differences method.
    std::vector<double>sfForward;
    std::vector<double> sfReverse;
    std::vector<Real> sfFinite;

    // TapeSize on Number of Inpependent Variables dependence.
    std::vector<TapeSize> tapeSize;

    std::string path,pathOut;
    switch (type)
    {
    case 1:
        path = "AdjointShortRateModels\\TestCachedHullWhite";
        pathOut = "AdjointShortRateModels\\TestCachedHullWhite\\output";
        break;
    case 2:
        path = "AdjointShortRateModels\\TestCachedHullWhiteFixedReversion";
        pathOut = "AdjointShortRateModels\\TestCachedHullWhiteFixedReversion\\output";
    }
   
    // Plots streams.
    cl::AdjointTestOutput out(pathOut, {
         { "filename", "SigmaonVolatil" }
        ,{ "not_clear", "Not" }
        ,{ "title", "Model sigma on volatilty dependence" }
        ,{ "ylabel", "Model Sigma" }
    });
    cl::AdjointTestOutput outPerform(path, {
         { "filename", "ModelSigma" }
        ,{ "not_clear", "Not" }
        ,{ "title", "Performance Time" }
        ,{ "ylabel", "Time (s)" }
    });
    cl::AdjointTestOutput outAdjoint(path, {
         { "filename", "Adjoint" }
        ,{ "not_clear", "Not" }
        ,{ "title", "Adjoint time dependence on number of volatilities in an input vector" }
        ,{ "ylabel", "Adjoint Time (s)" }
    });
    cl::AdjointTestOutput outSize(path, {
         { "filename", "TapeSize" }
        ,{ "not_clear", "Not" }
        ,{ "title", "Tape size dependence on  number of independent variables" }
        ,{ "ylabel", "Memory (MB)" }
    });

    // Starting number of votilities
    Size volatilStart = size;

#ifdef CL_GRAPH_GEN
    // Plot output dependence of model sigma on the first component of the  volatilities vector
    Size pointsNumber_ = 100;
    rateDependencePlot(vars, pointsNumber_, type, out);

    // Starting number of volatilities
    volatilStart = 5;
#endif

    // Running Adjoint Differentiation for various numbers of bonds in a portfolio
    for (Size pos = volatilStart; pos <= size; pos++)
    {
        boost::timer timer;
        std::vector<cl::TapeDouble> indepVar(vars.vol_.begin(), vars.vol_.begin() + pos);

        // Beginning of tape recording
        outPerform.log() << "Start of tape recording: " <<currentTime()<< std::endl;
        Independent(indepVar);
        calibrationData(indepVar, vars.data_);
        std::vector<cl::TapeDouble> output(1);
        output[0] = calibrate(vars.termStructure_, vars.model_, vars.index_, vars.engine_, vars.swaptions_, vars.data_, type, pos);
        cl::TapeFunction<double > f(indepVar, output);

        // Memory for tape calculation
        tapeSize.push_back(TapeSize{ pos, f.Memory() });

        // End of tape recording. Differentiaion wiil be held with respect to  the independent variables vector
        double timeTapeRecording = timer.elapsed();
        outPerform.log() << "Time for tape recording:\t " << timeTapeRecording <<" s"<< std::endl;

        // Start differentiation in Forward mode
        gradForward(f, sfForward, outPerform,true, false);

        // Start differentiation in Reverse mode
        double timeAdjoint = gradReverse(f, sfReverse, outPerform, true, false);

        // Finite differences (central formula)
        double timeAnalytical = finiteDiff(vars.termStructure_, vars.model_, vars.index_, vars.engine_, vars.swaptions_, vars.data_, type, pos, vars.h_, sfFinite, outPerform);

        // Adding new data to the performace result vector
        vars.performanceTime_.push_back(PerformanceTime{ timeTapeRecording, timeAdjoint, timeAnalytical, pos });
        vars.adjointTime_.push_back({ timeAdjoint, pos });

        outPerform.log() << std::endl;
    }

    // Write performance results to a plt file
    outSize << tapeSize;
    outPerform << vars.performanceTime_;
    vars.performanceTime_.clear();
    outAdjoint << vars.adjointTime_;
    vars.adjointTime_.clear();

    Settings::instance().resetEvaluationDate();

    return checkWithFiniteDiff(sfForward, sfReverse, sfFinite, outPerform, vars.relTolerance_, vars.absTolerance_);
}
bool AdjointShortRateModelsTest::testCachedHullWhite() {
    bool result = false;
    std::cout << "Testing Hull-White calibration against cached values using swaptions with start delay..." << std::endl;
    std::cout << "Testing Hull-White model dependency on volatilities" << std::endl;
#ifdef CL_TAPE_CPPAD
    Size size = 24;
    CommonVars vars(size);
    result = testHullWhite(vars, size, 1);
#endif
    return result;
}

bool  AdjointShortRateModelsTest::testCachedHullWhiteFixedReversion() {
    bool result = false;
    std::cout << "Testing Hull-White calibration with fixed reversion against cached values..." << std::endl;
#ifdef CL_TAPE_CPPAD
    Size size = 24;
    CommonVars vars(size);
    result = testHullWhite(vars, size, 2);
#endif
    return result;
}

bool AdjointShortRateModelsTest::testFuturesConvexityBias() {
    std::cout<<"Testing Hull-White futures convexity bias...";
#ifdef CL_TAPE_CPPAD

    // G. Kirikos, D. Novak, "Convexity Conundrums", Risk Magazine, March 1997
    Real futureQuote = 94.0;
    Real a = 0.03;
    Real sigma = 0.015;
    Time t = 5.0;
    Time T = 5.25;
    std::vector<cl::TapeDouble> parameters = {a, sigma, futureQuote};
    size_t sizeof_indep = parameters.size();

    // Beginning of tape recording
    boost::timer timer;

    Independent(parameters);
    std::vector<cl::TapeDouble> output(1);
    output[0] = (100.0 - parameters[sizeof_indep - 1]) / 100.0 -
        HullWhite::convexityBias(parameters[sizeof_indep - 1], t, T, parameters[1], parameters[0]);
    cl::TapeFunction<double > f(parameters, output);

    // End of tape recording. Differentiaion wiil be held with respect to  the independent variables vector
    std::cout << "Time for tape recording : " << timer.elapsed() << std::endl;

    std::vector<double> sf_Forward(sizeof_indep);
    std::vector<double> sf_Reverse(sizeof_indep);

    // Start of differentiation in Forward mode
    gradForward(f, sf_Forward, true, false);

    // Start of  differentiation in Reverse mode
    gradReverse(f, sf_Reverse, true, false);

    // Finite differences (central formula)
    double h = 1.0e-3;
    std::vector<Real> sf_Finite(sizeof_indep);
    for (Size i = 0; i < sizeof_indep; i++)
    {
        parameters[i]-=h;
        Real stepbackward = (100.0 - parameters[sizeof_indep - 1]) / 100.0 -
            HullWhite::convexityBias(parameters[sizeof_indep - 1], t, T, parameters[1], parameters[0]);
        parameters[i]+= 2*h;
        Real stepforward = (100.0 - parameters[sizeof_indep - 1]) / 100.0 -
            HullWhite::convexityBias(parameters[sizeof_indep - 1], t, T, parameters[1], parameters[0]);
        sf_Finite[i] = (stepforward - stepbackward) / (2 * h);
        parameters[i]-=h;
    }
    double tolerance = 1e-2;
    double abs_tolerance = 1e-4;
    return checkWithFiniteDiff(sf_Forward, sf_Reverse, sf_Finite, tolerance, abs_tolerance);
#endif
}

test_suite*  AdjointShortRateModelsTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("CppAD Hull-White model calibration  tests");
    suite->add(QUANTLIB_TEST_CASE(&AdjointShortRateModelsTest::testCachedHullWhite));
    suite->add(QUANTLIB_TEST_CASE(&AdjointShortRateModelsTest::testCachedHullWhiteFixedReversion));
    suite->add(QUANTLIB_TEST_CASE(&AdjointShortRateModelsTest::testFuturesConvexityBias));
    return suite;
}
#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_ShortRateModelsTest)

BOOST_AUTO_TEST_CASE(testCachedHullWhite)
{
    BOOST_CHECK(AdjointShortRateModelsTest::testCachedHullWhite());
}

BOOST_AUTO_TEST_CASE(testCachedHullWhiteFixedReversion)
{
    BOOST_CHECK(AdjointShortRateModelsTest::testCachedHullWhiteFixedReversion());
}

BOOST_AUTO_TEST_CASE(testFuturesConvexityBias)
{
    BOOST_CHECK(AdjointShortRateModelsTest::testFuturesConvexityBias());
}

BOOST_AUTO_TEST_SUITE_END()

#endif