/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2007 Ferdinando Ametrano
Copyright (C) 2007 Marco Bianchetti
Copyright (C) 2007 Cristina Duminuco
Copyright (C) 2007 Mark Joshi
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

//based on marketmodel_smmcapletalphacalibration.cpp from test-suite

#include "adjointmarketmodelcalibrationtest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include <ql/models/marketmodels/products/multiproductcomposite.hpp>
#include <ql/quantlib.hpp>
#include <boost/timer.hpp>
#include <iostream>

using namespace QuantLib;
using namespace std;
using namespace boost::unit_test_framework;

namespace
{
    bool checkWithFiniteDiff(
        std::vector<double>& forward_mode,
        std::vector<double>& reverse_mode,
        std::vector<Real>& finite_diff, Size n,
        double relativeTol, double absTol)
    {
        bool result = true;
        for (size_t i = 0; i < n; i++)
        {
            Real maxabs = std::max(std::abs(finite_diff[i]), std::max(std::abs(forward_mode[i]), std::abs(reverse_mode[i])));
            Real tol = relativeTol * maxabs + absTol;
            if (std::abs(forward_mode[i] - finite_diff[i]) > tol)
            {
                result = false;
                BOOST_ERROR("\nForward mode and Finite difference derivative[" << i << "] mismatch."
                    << "\n    forward mode: " << forward_mode[i]
                    << "\n    finite difference:   " << finite_diff[i]
                    << "\n    tolerance:  " << tol);
            }
            if (std::abs(reverse_mode[i] - finite_diff[i]) > tol)
            {
                result = false;
                BOOST_ERROR("\nReverse mode and Finite difference derivative[" << i << "] mismatch."
                    << "\n    reverse mode: " << reverse_mode[i]
                    << "\n    finite difference:   " << finite_diff[i]
                    << "\n    tolerance:  " << tol);
            }
            if (std::abs(reverse_mode[i] - forward_mode[i]) > tol)
            {
                result = false;
                BOOST_ERROR("\nForward mode and Reverse mode derivative[" << i << "] mismatch."
                    << "\n    forward mode: " << forward_mode[i]
                    << "\n    reverse mode:   " << reverse_mode[i]
                    << "\n    tolerance:  " << tol);
            }
        }
        return result;
    }

    // Function std::vector<Real> jacobian(cl::TapeFunction<double>& f, std::vector<Real> const& x, size_t iterNum, double& time)
    // Adapter with a timer for Jacobian.
    std::vector<Real> jacobian(cl::TapeFunction<double>& f, std::vector<Real> const& x, size_t iterNum, double& time)
    {
        std::vector<double> xd, yd;
        std::for_each(x.cbegin(), x.cend(),
            [&](Real const& a){ xd.push_back(double(a)); }
        );

        boost::timer timer;
        for (size_t i = 0; i < iterNum; ++i)
            yd = f.Jacobian(xd);
        time = timer.elapsed() / iterNum;

        std::vector<Real> y;
        std::for_each(yd.cbegin(), yd.cend(),
            [&](double const& a){ y.push_back(Real(a)); }
        );
        return y;
    }

    struct RealPair
    {
        Real x, y;

        static std::deque<std::string> get_columns()
        {
            static std::deque<std::string> columns = { " ", "" };
            return columns;
        }

        template <typename stream_type>
        friend inline stream_type& operator << (stream_type& stm, RealPair& p)
        {
            stm << p.x << ";" << p.y << std::endl;
            return stm;
        }
    };
}
//example based on marketmodel_smmcapletalphacalibration.cpp
//algorithm of calibration may be found here:
//http://fbe.unimelb.edu.au/__data/assets/pdf_file/0019/806302/172.pdf

namespace {

    Date todaysDate_, startDate_, endDate_;
    std::vector<Time> rateTimes_;
    std::vector<Real> accruals_;
    Calendar calendar_;
    DayCounter dayCounter_;
    std::vector<Rate> todaysForwards_, todaysSwaps_;
    std::vector<Real> coterminalAnnuity_;
    Size numberOfFactors_;
    Real alpha_, alphaMax_, alphaMin_;
    Spread displacement_;
    std::vector<DiscountFactor> todaysDiscounts_;
    std::vector<Volatility> swaptionDisplacedVols_, swaptionVols_;
    std::vector<Volatility> capletDisplacedVols_, capletVols_;
    Real a_, b_, c_, d_;
    Real longTermCorrelation_, beta_;
    Size measureOffset_;
    unsigned long seed_;
    Size paths_, trainingPaths_;

    void setup(Size size) {
        // Times
        calendar_ = NullCalendar();
        todaysDate_ = Settings::instance().evaluationDate();
        //startDate = todaysDate + 5*Years;
        endDate_ = todaysDate_ + (size + 1) * Months;
        Schedule dates(todaysDate_, endDate_, Period(Monthly),
            calendar_, Following, Following, DateGeneration::Backward, false);
        rateTimes_ = std::vector<Time>(dates.size() - 1);
        accruals_ = std::vector<Real>(rateTimes_.size() - 1);
        dayCounter_ = SimpleDayCounter();
        for (Size i = 1; i<dates.size(); ++i)
            rateTimes_[i - 1] = dayCounter_.yearFraction(todaysDate_, dates[i]);
        for (Size i = 1; i<rateTimes_.size(); ++i)
            accruals_[i - 1] = rateTimes_[i] - rateTimes_[i - 1];

        // Rates & displacement
        todaysForwards_ = std::vector<Rate>(accruals_.size());
        numberOfFactors_ = 3;
        alpha_ = 0.0;
        alphaMax_ = 1.0;
        alphaMin_ = -1.0;
        displacement_ = 0.0;
        for (Size i = 0; i<todaysForwards_.size(); ++i) {
            todaysForwards_[i] = 0.03 + 0.0025 * i / size;
            //    todaysForwards_[i] = 0.03;
        }
        LMMCurveState curveState_lmm(rateTimes_);
        curveState_lmm.setOnForwardRates(todaysForwards_);
        todaysSwaps_ = curveState_lmm.coterminalSwapRates();

        // Discounts
        todaysDiscounts_ = std::vector<DiscountFactor>(rateTimes_.size());
        todaysDiscounts_[0] = 0.95;
        for (Size i = 1; i<rateTimes_.size(); ++i)
            todaysDiscounts_[i] = todaysDiscounts_[i - 1] /
            (1.0 + todaysForwards_[i - 1] * accruals_[i - 1]);

        a_ = 0.0;
        b_ = 0.17;
        c_ = 1.0;
        d_ = 0.10;

        // Cap/Floor Correlation
        longTermCorrelation_ = 0.5;
        beta_ = 0.2;
        measureOffset_ = 5;
    }
}


std::vector<Real> AdjointMarketModelCalibrationTest::calibrationErrors(const std::vector<Volatility>& capletVols)
{
    Size numberOfRates = todaysForwards_.size();

    EvolutionDescription evolution(rateTimes_);

    boost::shared_ptr<PiecewiseConstantCorrelation> fwdCorr(new
        ExponentialForwardCorrelation(rateTimes_,
        longTermCorrelation_,
        beta_));

    boost::shared_ptr<LMMCurveState> cs(new LMMCurveState(rateTimes_));
    cs->setOnForwardRates(todaysForwards_);

    boost::shared_ptr<PiecewiseConstantCorrelation> corr(new
        CotSwapFromFwdCorrelation(fwdCorr, *cs, displacement_));

    std::vector<boost::shared_ptr<PiecewiseConstantVariance> >
        swapVariances(numberOfRates);
    for (Size i = 0; i<numberOfRates; ++i) {
        swapVariances[i] = boost::shared_ptr<PiecewiseConstantVariance>(new
            PiecewiseConstantAbcdVariance(a_, b_, c_, d_,
            i, rateTimes_));
    }

    // create calibrator
    std::vector<Real> alphaInitial(numberOfRates, alpha_);
    std::vector<Real> alphaMax(numberOfRates, 1.0);
    std::vector<Real> alphaMin(numberOfRates, -1.0);
    bool maximizeHomogeneity = false;

    CTSMMCapletAlphaFormCalibration calibrator(evolution,
        corr,
        swapVariances,
        capletVols,
        cs,
        displacement_,
        alphaInitial,
        alphaMax,
        alphaMin,
        maximizeHomogeneity);

    // calibrate
    Natural maxIterations = 10;
    Real capletTolerance = 1e-4; // i.e. 1 bp
    Natural innerMaxIterations = 100;
    Real innerTolerance = 1e-8;

    bool result = calibrator.calibrate(numberOfFactors_,
        maxIterations,
        capletTolerance,
        innerMaxIterations,
        innerTolerance);

    if (!result)
        BOOST_ERROR("calibration failed");

    std::vector<Real> calibrErrors(4);
    calibrErrors[0] = calibrator.capletRmsError();
    calibrErrors[1] = calibrator.capletMaxError();
    calibrErrors[2] = calibrator.swaptionRmsError();
    calibrErrors[3] = calibrator.swaptionMaxError();

    return calibrErrors;
}


bool AdjointMarketModelCalibrationTest::testFunction()
{
    const Real relativeError = 1.0e-4,
        absoluteError = 1.0e-5,
        delta = 1.0e-6;

    // Sizes of independent vectors
    const size_t minSize = 12,
        maxSize = 100,
        step = 12;

    bool result = true;

    std::vector<Volatility> capletVolsAll(maxSize);
    Real coef = 0.1 / maxSize;
    for (size_t i = 0; i < maxSize; ++i) {
        capletVolsAll[i] = 0.2 - coef * i;
    }

    for (size_t size = minSize; size <= maxSize; size += step)
    {
        setup(size);
        std::vector<Volatility> capletVols(capletVolsAll.begin(), capletVolsAll.begin() + size);
        
        boost::timer timer;
        Independent(capletVols);
        std::vector<Real> calibrErrors = calibrationErrors(capletVols);
        cl::TapeFunction<double> f(capletVols, calibrErrors);
        double timeTapeRecording = timer.elapsed();
        
        // Compute derivatives using Jacobian
        double timeAdjoint;
        std::vector<Real> jacob = jacobian(f, capletVols, 1000, timeAdjoint);

        // Compute finite difference derivatives using step-forward formula
        // Function derivative = (step-forward function value - current function value) / step
        std::vector<Real> jacobFinDiff(4 * size);
        timer.restart();
        const size_t iterNum = 10;
        for (size_t k = 0; k < iterNum; ++k)
        {
            for (size_t i = 0; i < size; ++i)
            {
                capletVols[i] += delta;
                std::vector<Real> calibrErrorsNext = calibrationErrors(capletVols);
                for (size_t j = 0; j < 4; ++j)
                    jacobFinDiff[i + j * size] = (calibrErrorsNext[j] - calibrErrors[j]) / delta;
                capletVols[i] -= delta;
            }
        }
        double timeAnalytical = timer.elapsed() / iterNum;

        // Check results
        result &= checkWithFiniteDiff(jacob, jacobFinDiff, relativeError, absoluteError);
    }
    return result;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////


Real AdjointMarketModelCalibrationTest::findCapletRmsError(
    std::vector<Rate>& todaysForwards_,
    std::vector<cl::TapeDouble>& mktCapletVols,
    std::vector<Real>& h)
{
    capletVols_.resize(todaysSwaps_.size());
    for (Size i = 0; i < todaysSwaps_.size(); i++)
        capletVols_[i] = mktCapletVols[i] + h[i];

    // Cap/Floor Correlation

    longTermCorrelation_ = 0.5;
    beta_ = 0.2;
    measureOffset_ = 5;

    Size numberOfRates = todaysForwards_.size();

    EvolutionDescription evolution(rateTimes_);

    boost::shared_ptr<PiecewiseConstantCorrelation> fwdCorr(new
                                                            ExponentialForwardCorrelation(rateTimes_,
                                                            longTermCorrelation_,
                                                            beta_));

    boost::shared_ptr<LMMCurveState> cs(new LMMCurveState(rateTimes_));
    cs->setOnForwardRates(todaysForwards_);

    boost::shared_ptr<PiecewiseConstantCorrelation> corr(new
                                                         CotSwapFromFwdCorrelation(fwdCorr, *cs, displacement_));

    std::vector<boost::shared_ptr<PiecewiseConstantVariance> >
        swapVariances(numberOfRates);
    for (Size i = 0; i < numberOfRates; ++i)
        swapVariances[i] = boost::shared_ptr<PiecewiseConstantVariance>(new
        PiecewiseConstantAbcdVariance(a_, b_, c_, d_,
        i, rateTimes_));

    // create calibrator

    std::vector<Real> alpha(numberOfRates, alpha_);
    bool lowestRoot = true;
    bool useFullApprox = false;
    CTSMMCapletOriginalCalibration calibrator(evolution,
                                              corr,
                                              swapVariances,
                                              capletVols_,
                                              cs,
                                              displacement_,
                                              alpha,
                                              lowestRoot,
                                              useFullApprox);
    // calibrate

    Natural maxIterations = 2;
    Real capletTolerance = (maxIterations == 1 ? 0.0032 : 0.0001);
    Natural innerMaxIterations = 50;
    Real innerTolerance = 1e-9;
    bool result = calibrator.calibrate(numberOfFactors_,
                                       maxIterations,
                                       capletTolerance / 10,
                                       innerMaxIterations,
                                       innerTolerance);
    if (!result)
        BOOST_ERROR("calibration failed");

    return calibrator.capletRmsError().value();
}

bool AdjointMarketModelCalibrationTest::testCapletCalibration()
{
    BOOST_TEST_MESSAGE("Testing GHLS caplet calibration "
                       "in a lognormal coterminal swap market model...");
    bool result = false;

#ifdef CL_TAPE_CPPAD

    Size n = 10;
    double tolerance = 1e-3;
    boost::timer timer;
    // Times
    calendar_ = NullCalendar();
    todaysDate_ = Settings::instance().evaluationDate();

    endDate_ = todaysDate_ + 66 * Months;
    Schedule dates(todaysDate_, endDate_, Period(Semiannual),
                   calendar_, Following, Following, DateGeneration::Backward, false);

    rateTimes_ = std::vector<Time>(dates.size() - 1);
    accruals_ = std::vector<Real>(rateTimes_.size() - 1);
    dayCounter_ = SimpleDayCounter();
    for (Size i = 1; i < dates.size(); ++i)
        rateTimes_[i - 1] = dayCounter_.yearFraction(todaysDate_, dates[i]);

    for (Size i = 1; i < rateTimes_.size(); ++i)
        accruals_[i - 1] = rateTimes_[i] - rateTimes_[i - 1];

    // Rates & displacement
    std::vector<Rate> todaysForwards_(accruals_.size());
    numberOfFactors_ = 3;
    alpha_ = -0.05;
    alphaMax_ = 1.0;
    alphaMin_ = -1.0;
    displacement_ = 0.0;
    for (Size i = 0; i < todaysForwards_.size(); ++i)
        todaysForwards_[i] = 0.03 + 0.0025*i;

    LMMCurveState curveState_lmm(rateTimes_);
    curveState_lmm.setOnForwardRates(todaysForwards_);
    todaysSwaps_ = curveState_lmm.coterminalSwapRates();

    a_ = 0.0;
    b_ = 0.17;
    c_ = 1.0;
    d_ = 0.10;

    std::vector<cl::TapeDouble> mktCapletVols(n);
    mktCapletVols[0] = 0.1640;
    mktCapletVols[1] = 0.1740;
    mktCapletVols[2] = 0.1840;
    mktCapletVols[3] = 0.1940;
    mktCapletVols[4] = 0.1840;
    mktCapletVols[5] = 0.1740;
    mktCapletVols[6] = 0.1640;
    mktCapletVols[7] = 0.1540;
    mktCapletVols[8] = 0.1440;
    mktCapletVols[9] = 0.1340376439125532;

    std::vector<cl::TapeDouble> capletRmsError(1);
    std::vector<Real> H(n);

    cl::Independent(mktCapletVols);

    capletRmsError[0] = findCapletRmsError(todaysForwards_, mktCapletVols, H);
    cl::TapeFunction<double> f(mktCapletVols, capletRmsError);

    std::cout << "Time for taping : " << timer.elapsed() << std::endl;

    //finite differences
    //d(capletRmsError)/d(mktCapletVols[0]) = (3.31037e-005 - 3.31025e-005)/ 0.1e-005 = 0.0012;

    std::vector<double> sx(n);
    std::vector<double> sy(1);
    std::vector<double> sf, sf_Forward(n), sf_Reverse(n);

    //Forward mode
    timer.restart();
    for (Size k = 0; k < 1000; k++)
    for (size_t i = 0; i < n; i++)
    {
        sx[i] = 1.;
        sf = f.Forward(1, sx);
        sx[i] = 0;
        sf_Forward[i] = sf[0];
    }

    std::cout << "Time for Forward mode : " << timer.elapsed() / 1000 << std::endl;

    //Reverse  mode
    timer.restart();
    for (Size k = 0; k < 1000; k++)
    {
        sy[0] = 1;
        sf_Reverse = f.Reverse(1, sy);
    }
    std::cout << "Time for Reverse mode : " << timer.elapsed() / 1000 << std::endl;

    //Finite differences
    double h = 1.0e-10;
    std::vector<Real> sf_Finite(n);
    for (Size i = 0; i < n / 2; i++)
    {
        H[i] = h;
        sf_Finite[i] = (findCapletRmsError(todaysForwards_, mktCapletVols, H) - capletRmsError[0]) / h;
        H[i] = 0;
    }

    result = checkWithFiniteDiff(sf_Forward, sf_Reverse, sf_Finite, n / 2, tolerance, tolerance);
#endif
    return result;
}

test_suite* AdjointMarketModelCalibrationTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("ADMarketModelCalibration test");

    //suite->add(QUANTLIB_TEST_CASE(&AdjointMarketModelCalibrationTest::testFunction));
    suite->add(QUANTLIB_TEST_CASE(&AdjointMarketModelCalibrationTest::testCapletCalibration));

    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_market_model_calibration)

BOOST_AUTO_TEST_CASE(testCapletCalibration)
{
    BOOST_CHECK(AdjointMarketModelCalibrationTest::testCapletCalibration());
}

BOOST_AUTO_TEST_SUITE_END()

#endif
