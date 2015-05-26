/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2008 Yee Man Chan
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
#include "adjointtestutilities.hpp"
#include "adjointgjrgarchmodeltest.hpp"

//based on gjrgarchmodel.cpp file from test-suite
using namespace QuantLib;
using namespace boost::unit_test_framework;

namespace
{
    struct CommonVars
    {
        // Cleanup.
        SavedSettings backup_;

        // Global data.
        Date settlement_;

        DayCounter dayCounter_;
        Calendar calendar_;

        std::vector<Integer> t_;
        std::vector<Real> r_;
        std::vector<Real> v_;

        Real strike;

        const Real omega_ = 2.0e-6;
        const Real alpha_ = 0.024;
        const Real beta_ = 0.93;
        const Real gamma_ = 0.059;
        const Real lambda_ = 0.1;
        const Real daysPerYear_ = 365.0;

        CommonVars()
        {
            settlement_ = Date(5, July, 2002);
            Settings::instance().evaluationDate() = settlement_;

            dayCounter_ = Actual365Fixed();
            calendar_ = TARGET();

            t_ = { 13, 41, 75, 165, 235 };
            r_ = { 0.0357, 0.0349, 0.0341, 0.0338, 0.0329 };
            v_ = { 0.6645, 0.4875, 0.4204, 0.3667, 0.3432 };

            strike = 3400;
        }
        
        // Get model calibration error with given volatilities.
        Real getGJRGARCHmodelCalibrationError(std::vector<Real> vol)
        {
            std::vector<Date> dates;
            std::vector<Rate> rates;
            dates.push_back(settlement_);
            rates.push_back(0.0357);
            for (Size i = 0; i < vol.size(); i++)
            {
                dates.push_back(settlement_ + t_[i]);
                rates.push_back(r_[i]);
            }

            // Create handle for YieldTermStructure. 
            Handle<YieldTermStructure> riskFreeTS(boost::shared_ptr<YieldTermStructure>(new ZeroCurve(dates, rates, dayCounter_)));
            Handle<YieldTermStructure> dividendTS(boost::shared_ptr<YieldTermStructure>(new FlatForward(settlement_, Handle<Quote>(boost::shared_ptr<Quote>(new SimpleQuote(0.0))), dayCounter_)));

            Handle<Quote> s0(boost::shared_ptr<Quote>(new SimpleQuote(4468.17)));

            // Calculate coef.
            const Real m1 = beta_ + (alpha_ + gamma_*CumulativeNormalDistribution()(lambda_))
                *(1.0 + lambda_*lambda_) + gamma_*lambda_*std::exp(-lambda_*lambda_ / 2.0)
                / std::sqrt(2.0*M_PI);
            const Real v0 = omega_ / (1.0 - m1);

            // Create gjrgarch process for model.
            boost::shared_ptr<GJRGARCHProcess> process(new GJRGARCHProcess(
                riskFreeTS, dividendTS, s0, v0,
                omega_, alpha_, beta_, gamma_, lambda_, daysPerYear_));

            // Create gjrgarch model using gjrgarch process.
            boost::shared_ptr<GJRGARCHModel> model(new GJRGARCHModel(process));

            // Create engine.
            boost::shared_ptr<PricingEngine> engine(new AnalyticGJRGARCHEngine(boost::shared_ptr<GJRGARCHModel>(model)));

            std::vector<boost::shared_ptr<CalibrationHelper>> options;

            for (Size i = 0; i < vol.size(); i++)
            {
                Handle<Quote> vol(boost::shared_ptr<Quote>(new SimpleQuote(vol[i])));
                Period maturity((int)(t_[i] / 7.), Weeks);

                // Create calibration helper.
                boost::shared_ptr<CalibrationHelper> helper(
                    new HestonModelHelper(maturity, calendar_,
                    s0->value(), strike, vol,
                    riskFreeTS, dividendTS,
                    CalibrationHelper::ImpliedVolError));

                // Set engine in helper.
                helper->setPricingEngine(engine);
                options.push_back(helper);
            }

            Real error = 0;
            // Calculate calibration error before calibration.
            for (Size i = 0; i < vol.size(); ++i)
            {
                const Real diff = options[i]->calibrationError()*100.0;
                error += diff*diff;
            }
            //std::cout << "Error before calibration: " << error << std::endl;

            // Set optimization method.
            Simplex om(0.01);

            // Calibrate model.
            model->calibrate(options, om, EndCriteria(50, 10, 1.0e-2, 1.0e-2, 1.0e-2));

            error = 0;
            // Calculate calibration error after calibration.
            for (Size i = 0; i < options.size(); ++i)
            {
                const Real diff = options[i]->calibrationError()*100.0;
                error += diff*diff;
            }

            //std::cout << "Error after calibration: " << error << std::endl;
            return error;
        }
    };

    struct Variation
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "Volatility", ""
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, Variation& v)
        {
                stm << v.volatility_
                    << ";" << v.calibrationError_
                    << std::endl;
                return stm;
            }

        Real volatility_;
        Real calibrationError_;
    };

}

// Test GJRGARCH model calibration using DAX volatility data.
bool AdjointGjrgarchModelTest::testGJRGARCHmodel()
{
    BOOST_TEST_MESSAGE("Testing adjoint differentiation for GJRGARCH model calibration using DAX volatility data...\n");
    bool result = false;
#ifdef CL_TAPE_CPPAD
    CommonVars vars;

    // Create output streams for plots.
    cl::AdjointTestOutput outPerform("AdjointGjrgarchModel//GJRGARCHmodel"
                                                , { { "filename", "AdjointPerformance" }
                                                  , { "not_clear", "Not" }
                                                  , { "line_box_width", "-5" }
                                                  , { "title", "Calibration error differentiation performance with respect to volatility" }
                                                  , { "ylabel", "Time (s)" }
                                                  , { "xlabel", "Number of volatilities" } });

    cl::AdjointTestOutput outAdjoint("AdjointGjrgarchModel//GJRGARCHmodel"
                                     , { { "filename", "Adjoint" }
                                       , { "not_clear", "Not" }
                                       , { "title", "Calibration error adjoint differentiation with respect to volatility" }
                                       , { "cleanlog", "false" }
                                       , { "ylabel", "Time (s)" }
                                       , { "xlabel", "Number of volatilities" } });

    cl::AdjointTestOutput outSize("AdjointGjrgarchModel//GJRGARCHmodel"
                                      , { { "filename", "TapeSize" }
                                        , { "not_clear", "Not" }
                                        , { "title", "Tape size dependence on number of volatilities" }
                                        , { "cleanlog", "false" }
                                        , { "ylabel", "Size (MB)" }
                                        , { "xlabel", "Number of volatilities" } });

    cl::AdjointTestOutput out("AdjointGjrgarchModel//GJRGARCHmodel//output"
                                       , { { "filename", "CalibrErronVol" }
                                         , { "ylabel", "Calibration Error" }
                                         , { "not_clear", "Not" }
                                         , { "title", "Calibration error on volatility" } 
                                         , { "cleanlog", "false" }
                                         , { "xlabel", "Volatility" } });


    // Create vectors for store performance results.
    std::vector<PerformanceTime> performResults;
    std::vector<AdjointTime> adjointResults;
    std::vector<TapeSize> tapeSizeResults;

#if defined CL_GRAPH_GEN
    Size startSize = 0;
    Size maxSize = vars.v_.size();
#else
    Size startSize = 1;
    Size maxSize = 2;
#endif

    outPerform.log() << "Testing adjoint differentiation for GJRGARCH model calibration using DAX volatility data." << std::endl << std::endl;

    // Variate number of independent variables.
    for (Size s = startSize; s < maxSize; s++)
    {
        // Create temporary structures for store performance results.
        PerformanceTime performTime;
        AdjointTime adTime;
        TapeSize tSize;

        Size sizeof_indep = s + 1;

        // Store number of independent variables
        performTime.indepVarNumber_ = sizeof_indep;
        adTime.indepVarNumber_ = sizeof_indep;
        tSize.indepVarNumber_ = sizeof_indep;

        outPerform.log() << "Number of independent variables: " << sizeof_indep << std::endl;
        outPerform.log() << "Number of dependent variables: 1" << std::endl;

        outPerform.log() << "Start of tape recording: " << currentTime() << std::endl;

        // Start timing of tape recording.
        boost::timer timer;

        // Create vector of independent variables.
        std::vector<cl::TapeDouble> vol(sizeof_indep);

        // Initialize vector of independent variables.
        for (Size i = 0; i < sizeof_indep; i++)
            vol[i] = vars.v_[i];

        // Start taping. Declare strikes as independent variables.
        Independent(vol);

        // Create vector of dependent variables.
        std::vector<cl::TapeDouble> Y(1);

        // Calculate calibration error.
        Y[0] = vars.getGJRGARCHmodelCalibrationError(vol);

        // End of tape recording. Declare calibration error as dependent variable.
        // Differentiaion will be held with respect to the independent variables vector.
        cl::TapeFunction<double> f(vol, Y);

        // Store time of tape recording.
        performTime.timeTapeRecording_ = timer.elapsed();

        outPerform.log() << "End of tape recording. " << std::endl;
        outPerform.log() << "Time for tape recording: " << performTime.timeTapeRecording_ << std::endl;

        // Store size of tape.
        tSize.memory_ = f.Memory();

        std::vector<double> sf_Forward(sizeof_indep), sf_Reverse(sizeof_indep);

        // Start differentiation in Forward mode.
        double tf = gradForward(f, sf_Forward, outPerform, false, false);

        // Start differentiation in Reverse mode.
        double tr = gradReverse(f, sf_Reverse, outPerform, false, false);

        // Store time of adjoint differntiation.
        performTime.timeAdjoint_ = std::min(tf, tr);
        adTime.timeAdjoint_ = performTime.timeAdjoint_;

        outPerform.log() << "Start of  differentiation using finite differences method: " << currentTime() << std::endl;

        // Start timing for calculating derivatives by finite difference.
        timer.restart();

        //Start differentiation using finite differences method.
        double h = 1.0e-10;

        //  Create vector for derivatives.
        std::vector<Real> sf_Finite(sizeof_indep);

        // Calculate calibration error.
        Real YF = vars.getGJRGARCHmodelCalibrationError(vol);

        for (Size i = 0; i < sizeof_indep; i++)
        {
            vol[i] += h;
            //Evaluate derivative using finite difference.
            sf_Finite[i] = (vars.getGJRGARCHmodelCalibrationError(vol) - YF) / h;
            vol[i] -= h;
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
        result = checkWithFiniteDiff(sf_Forward, sf_Reverse, sf_Finite, outPerform, 1e-2, 1e-2);

        outPerform.log() << std::endl << std::endl;
    }

    std::vector<Variation> variationResults;
#if defined CL_GRAPH_GEN
    std::vector<Real> volatility;

    // Create input volatilities vector.
    for (Size i = 0; i < vars.v_.size(); i++)
    {
        volatility.push_back(vars.v_[i]);
    }

    Real step = (volatility[1] - volatility[0]) / 5;

    // Calculate calibration error on variation of volatility.
    for (Size i = 0; i < 5; i++)
    {
        volatility[0] += i * step;

        Variation var;
        var.volatility_ = volatility[0];
        var.calibrationError_ = vars.getGJRGARCHmodelCalibrationError(volatility);
        variationResults.push_back(var);

        volatility[0] -= i * step;
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

test_suite* AdjointGjrgarchModelTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("CppAD GJRGARCH model test");
    suite->add(QUANTLIB_TEST_CASE(&AdjointGjrgarchModelTest::testGJRGARCHmodel));
    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_gjrgarchmodel)

BOOST_AUTO_TEST_CASE(testGJRGarchModelCalibration)
{
    BOOST_CHECK(AdjointGjrgarchModelTest::testGJRGARCHmodel());
}

BOOST_AUTO_TEST_SUITE_END()

#endif
