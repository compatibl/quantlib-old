/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2013 Gary Kennedy
Copyright (C) 2015 Peter Caspers
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

//based on blackformula.cpp file from test-suite

#include "adjointblackformulatest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include <ql/pricingengines/blackformula.hpp>

using namespace QuantLib;
using namespace std;
using namespace boost::unit_test_framework;


struct CommonVarBS
{
    CommonVarBS() :
    forwardPrice_(1.0)
    , bpvol_(0.01)
    , tte_(10.0)
    , stdDev_(bpvol_*std::sqrt(tte_))
    , optionType_(Option::Call)
    , discount_(0.95)
    , n_(10)
    , h_(1e-4)
    , tol_(1e-8)
    {
        d_.resize(n_);
        for (Size i = 0; i < n_; i++)
            d_[i] = -3.0 + 6.0 * i / n_;
    }
    double forwardPrice_ = 1.0;
    double bpvol_ = 0.01;
    double tte_ = 10.0;
    double stdDev_ = bpvol_*std::sqrt(tte_);
    Option::Type optionType_ = Option::Call;
    double discount_ = 0.95;
    vector<double> d_;
    Size n_;
    double h_;
    double tol_;
    double timeTapeRecording_;
    double timeAdjoint_;
    double timeAnalytical_;
    vector<PerformanceTime> performanceTime_;
    vector<TapeSize> tapeMemory_;
    boost::timer timer;

    Real bachelierImpliedVol(Real forward, Real tte, Real stdDev, Real discount, Real strike)
    {
        Real callPrem = bachelierBlackFormula(optionType_, strike, forward, stdDev, discount);
        Real impliedBpVol = bachelierBlackFormulaImpliedVol(optionType_, strike, forward, tte, callPrem, discount);
        return impliedBpVol;
    }
};

bool AdjointBlackFormulaTest::testBachelierImpliedVolStrike()
{
    BOOST_TEST_MESSAGE("Testing Bachelier implied vol...");
    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVarBS vars;
    std::vector<cl::TapeDouble> strike;
    vector<double> strikeDouble;
    std::vector<cl::TapeDouble> impliedBpVol;
    std::vector<double> Jacobian;
    std::vector<Real> sf_Finite;
    vars.performanceTime_.clear();
    cl::AdjointTestOutput outPerform("AdjointBlackFormula//BachelierImpliedVolStrike");

    //Running Adjoint differentiation with respect to various numbers of strike prices
    for (Size k = 0; k < vars.n_; k++)
    {
        outPerform.log() << "\nNumber of strike prices : k = " << k + 1 << std::endl;
        // Add new strike price.
        strikeDouble.push_back(vars.forwardPrice_ - (vars.d_[k]) * vars.bpvol_ * std::sqrt(vars.tte_));
        strike.push_back(strikeDouble[k]);

        // Start timing of tape recording.
        outPerform.log() << "Start taping : " << currentTime() << std::endl;
        vars.timer.restart();

        // Start taping. Declare forward rates as independent variables.
        Independent(strike);
        impliedBpVol.resize(k + 1, 0);

        // Calculate implied bachelier volatility vector.
        outPerform.log() << "Calculate implied bachelier volatility vector." << std::endl;
        for (Size i = 0; i < k + 1; i++)
        {
            impliedBpVol[i] = vars.bachelierImpliedVol(vars.forwardPrice_, vars.tte_, vars.stdDev_, vars.discount_, strike[i]);
            if (std::fabs(vars.bpvol_ - impliedBpVol[i])>1.0e-10)
            {
                BOOST_ERROR("Failed, expected " << vars.bpvol_ << " realised " << impliedBpVol[i]);
            }
        }

        cl::TapeFunction<double> f(strike, impliedBpVol);

        //End of tape recording. Differentiaion wiil be held with respect to  the independent variables vector
        vars.timeTapeRecording_ = vars.timer.elapsed();
        outPerform.log() << "End of tape recording. Time for tape recording : " << vars.timeTapeRecording_ << std::endl;
        //Store size of tape.
        vars.tapeMemory_.push_back(TapeSize { k + 1, f.Memory() });

        //Calculate derivatives using Jacobian
        outPerform.log() << "Start differentiation using Jacobian : " << currentTime() << std::endl;
        vars.timer.restart();
        for (Size s = 0; s < 10; s++)
            Jacobian = f.Jacobian(strikeDouble);
        vars.timeAdjoint_ = vars.timer.elapsed() / 10;
        outPerform.log() << "Time for differentiation using Jacobian : " << vars.timeAdjoint_ << std::endl;

        //Central-difference scheme
        outPerform.log() << "Start differentiation using central finite differences : " << currentTime() << std::endl;
        vars.timer.restart();
        sf_Finite.resize((k + 1)*(k + 1), 0.0);
        Real rightValue, leftValue;

        for (Size l = 0; l < k + 1; l++)
        {
            for (Size m = 0; m < k + 1; m++)
            {
                if (l == m)
                {
                    rightValue = vars.bachelierImpliedVol(vars.forwardPrice_, vars.tte_, vars.stdDev_, vars.discount_, strike[m] + vars.h_);
                    leftValue = vars.bachelierImpliedVol(vars.forwardPrice_, vars.tte_, vars.stdDev_, vars.discount_, strike[m] - vars.h_);
                    sf_Finite[l*(k + 1) + m] = (rightValue - leftValue) / (2 * vars.h_);
                }
                else
                    sf_Finite[l*(k + 1) + m] = 0;
            }
        }


        vars.timeAnalytical_ = vars.timer.elapsed();
        outPerform.log() << "Time for differentiation using central finite differences : " << vars.timeAnalytical_ << " s" << std::endl;

        //Check results
        outPerform.log() << "Check derivatives calculated by Jacobian and finite differences methods." << std::endl;
        result = checkWithFiniteDiff(Jacobian, sf_Finite, vars.tol_, vars.tol_);

        //Adding new data to the performace result vector
        vars.performanceTime_.push_back(PerformanceTime { vars.timeTapeRecording_, vars.timeAdjoint_, vars.timeAnalytical_, k + 1 });
    }
    outPerform << vars.performanceTime_;
#endif
    return result;
}


bool AdjointBlackFormulaTest::testBachelierImpliedVolD()
{
    BOOST_TEST_MESSAGE("Testing Bachelier implied vol...");
    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVarBS vars;
    std::vector<cl::TapeDouble> strike(vars.n_);
    std::vector<cl::TapeDouble> dVector;
    vector<double> dVectorDouble;
    std::vector<cl::TapeDouble> impliedBpVol;
    vars.performanceTime_.clear();
    cl::AdjointTestOutput outPerform("AdjointBlackFormula//BachelierImpliedVolD");
    //Running Adjoint differentiation with respect to various numbers of parameters d = (F-K)/(vol*sqrt(tte))
    for (Size k = 0; k < vars.n_; k++)
    {
        outPerform.log() << "\nNumber of parameters d = (F-K)/(vol*sqrt(tte)) : k = " << k + 1 << std::endl;
        // Add new parameter d.
        dVectorDouble.push_back(vars.d_[k]);
        dVector.push_back(dVectorDouble[k]);

        // Start timing of tape recording.
        outPerform.log() << "Start taping : " << currentTime() << std::endl;
        vars.timer.restart();

        // Start taping. Declare forward rates as independent variables.
        Independent(dVector);
        impliedBpVol.resize(k + 1, 0);
        outPerform.log() << "Calculate implied bachelier volatility vector." << std::endl;
        for (Size i = 0; i < k + 1; i++)
        {
            strike[i] = vars.forwardPrice_ - dVector[i] * vars.bpvol_ * std::sqrt(vars.tte_);
            impliedBpVol[i] = vars.bachelierImpliedVol(vars.forwardPrice_, vars.tte_, vars.stdDev_, vars.discount_, strike[i]);
            if (std::fabs(vars.bpvol_ - impliedBpVol[i])>1.0e-10)
            {
                BOOST_ERROR("Failed, expected " << vars.bpvol_ << " realised " << impliedBpVol[i]);
            }
        }

        cl::TapeFunction<double> f(dVector, impliedBpVol);

        //End of tape recording. Differentiaion wiil be held with respect to  the independent variables vector
        vars.timeTapeRecording_ = vars.timer.elapsed();
        outPerform.log() << "End of tape recording. Time for tape recording : " << vars.timeTapeRecording_ << std::endl;
        //Store size of tape.
        vars.tapeMemory_.push_back(TapeSize { k + 1, f.Memory() });

        //Calculate derivatives using Jacobian
        outPerform.log() << "Start differentiation using Jacobian : " << currentTime() << std::endl;
        std::vector<double> Jacobian;
        vars.timer.restart();
        for (Size s = 0; s < 10; s++)
            Jacobian = f.Jacobian(dVectorDouble);
        vars.timeAdjoint_ = vars.timer.elapsed() / 10;
        outPerform.log() << "Time for differentiation using Jacobian : " << vars.timeAdjoint_ << std::endl;

        //Central-difference scheme
        outPerform.log() << "Start differentiation using central finite differences : " << currentTime() << std::endl;
        vars.timer.restart();
        vector<Real> sf_Finite((k + 1)*(k + 1), 0.0);
        Real rightValue, leftValue;
        for (Size l = 0; l < k + 1; l++)
        {
            for (Size m = 0; m < k + 1; m++)
            {
                if (l == m)
                {
                    strike[m] = vars.forwardPrice_ - (dVector[m] + vars.h_)* vars.bpvol_ * std::sqrt(vars.tte_);
                    rightValue = vars.bachelierImpliedVol(vars.forwardPrice_, vars.tte_, vars.stdDev_, vars.discount_, strike[m]);

                    strike[m] = vars.forwardPrice_ - (dVector[m] - vars.h_)* vars.bpvol_ * std::sqrt(vars.tte_);
                    leftValue = vars.bachelierImpliedVol(vars.forwardPrice_, vars.tte_, vars.stdDev_, vars.discount_, strike[m]);

                    sf_Finite[l*(k + 1) + m] = (rightValue - leftValue) / (2 * vars.h_);
                }
                else
                    sf_Finite[l*(k + 1) + m] = 0;
            }
        }

        vars.timeAnalytical_ = vars.timer.elapsed();
        outPerform.log() << "Time for differentiation using central finite differences : " << vars.timeAnalytical_ << " s" << std::endl;

        //Check results
        outPerform.log() << "Check derivatives calculated by Jacobian and finite differences methods." << std::endl;
        result = checkWithFiniteDiff(Jacobian, sf_Finite, vars.tol_, vars.tol_);

        //Adding new data to the performace result vector
        vars.performanceTime_.push_back(PerformanceTime { vars.timeTapeRecording_, vars.timeAdjoint_, vars.timeAnalytical_, k + 1 });
    }
    outPerform << vars.performanceTime_;
#endif;
    return result;
}

bool AdjointBlackFormulaTest::testBachelierImpliedVolForward()
{
    BOOST_TEST_MESSAGE("Testing Bachelier implied vol... Derivatives with respect to forward = 0");
    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVarBS vars;
    std::vector<cl::TapeDouble> strike(vars.n_);
    std::vector<cl::TapeDouble> forwardPrice(1);
    forwardPrice[0] = vars.forwardPrice_;
    std::vector<cl::TapeDouble> impliedBpVol(vars.n_);
    cl::AdjointTestOutput outPerform("AdjointBlackFormula//BachelierImpliedVolForward");

    //Running Adjoint differentiation with respect to forward price
    outPerform.log() << "\nNumber of forward prices : k = " << 1 << std::endl;

    // Start timing of tape recording.
    outPerform.log() << "Start taping : " << currentTime() << std::endl;
    vars.timer.restart();

    // Start taping. Declare forward rates as independent variables.
    Independent(forwardPrice);

    // Calculate implied bachelier volatility vector.
    outPerform.log() << "Calculate implied bachelier volatility vector." << std::endl;
    for (Size i = 0; i < vars.n_; i++)
    {
        strike[i] = forwardPrice[0] - vars.d_[i] * vars.bpvol_ * std::sqrt(vars.tte_);
        impliedBpVol[i] = vars.bachelierImpliedVol(forwardPrice[0], vars.tte_, vars.stdDev_, vars.discount_, strike[i]);
        if (std::fabs(vars.bpvol_ - impliedBpVol[i])>1.0e-10)
        {
            BOOST_ERROR("Failed, expected " << vars.bpvol_ << " realised " << impliedBpVol[i]);
        }
    }

    cl::TapeFunction<double> f(forwardPrice, impliedBpVol);

    //End of tape recording. Differentiaion wiil be held with respect to  the independent variables vector
    vars.timeTapeRecording_ = vars.timer.elapsed();
    outPerform.log() << "End of tape recording. Time for tape recording : " << vars.timeTapeRecording_ << std::endl;
    //Store size of tape.
    vars.tapeMemory_.push_back(TapeSize { 1, f.Memory() });

    //Start differentiation in Forward mode
    vector<double> sf_Forward;
    outPerform.log() << "Start differentiation in Forward mode : " << currentTime() << std::endl;
    vars.timeAdjoint_ = gradForward(f, sf_Forward, false, false);
    outPerform.log() << "Time for differentiation in Forward mode : " << vars.timeAdjoint_ << " s" << std::endl;

    //Start differentiation using central finite differences.
    outPerform.log() << "Start differentiation using central finite differences : " << currentTime() << std::endl;
    vars.timer.restart();
    std::vector<Real> sf_Finite(vars.n_, 0.0);
    Real rightValue, leftValue;
    for (Size i = 0; i < vars.n_; i++)
    {
        forwardPrice[0] += vars.h_;
        strike[i] = forwardPrice[0] - vars.d_[i] * vars.bpvol_ * std::sqrt(vars.tte_);
        rightValue = vars.bachelierImpliedVol(forwardPrice[0], vars.tte_, vars.stdDev_, vars.discount_, strike[i]);

        forwardPrice[0] -= 2 * vars.h_;
        strike[i] = forwardPrice[0] - vars.d_[i] * vars.bpvol_ * std::sqrt(vars.tte_);
        leftValue = vars.bachelierImpliedVol(forwardPrice[0], vars.tte_, vars.stdDev_, vars.discount_, strike[i]);

        sf_Finite[i] = (rightValue - leftValue) / (2 * vars.h_);
    }
    vars.timeAnalytical_ = vars.timer.elapsed();
    outPerform.log() << "Time for differentiation using central finite differences : " << vars.timeAnalytical_ << " s" << std::endl;

    // Check derivatives calculated by forward, reverse and finite differences methods.
    outPerform.log() << "Check derivatives calculated by forward and finite differences methods." << std::endl;
    result = checkWithFiniteDiff(sf_Forward, sf_Finite, vars.tol_, vars.tol_);
    //Adding new data to the performace result vector
    vars.performanceTime_.push_back(PerformanceTime { vars.timeTapeRecording_, vars.timeAdjoint_, vars.timeAnalytical_, 1 });
#endif
    return result;
}

bool AdjointBlackFormulaTest::testBachelierImpliedVolTte()
{
    BOOST_TEST_MESSAGE("Testing Bachelier implied vol...");
    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVarBS vars;
    std::vector<cl::TapeDouble> strike(vars.n_);
    std::vector<cl::TapeDouble> tteVector(1);
    tteVector[0] = vars.tte_;
    std::vector<cl::TapeDouble> impliedBpVol(vars.n_);
    cl::AdjointTestOutput outPerform("AdjointBlackFormula//BachelierImpliedVolTte");

    //Running Adjoint differentiation with respect to tte
    outPerform.log() << "\nNumber of ttes : k = " << 1 << std::endl;

    // Start timing of tape recording.
    outPerform.log() << "Start taping : " << currentTime() << std::endl;
    vars.timer.restart();

    // Start taping. Declare forward rates as independent variables.
    Independent(tteVector);

    // Calculate implied bachelier volatility vector.
    outPerform.log() << "Calculate implied bachelier volatility vector." << std::endl;
    Real stdDev_tte = vars.bpvol_*std::sqrt(vars.tte_);
    for (Size i = 0; i < vars.n_; i++)
    {
        strike[i] = vars.forwardPrice_ - vars.d_[i] * vars.bpvol_ * std::sqrt(tteVector[0]);
        impliedBpVol[i] = vars.bachelierImpliedVol(vars.forwardPrice_, tteVector[0], stdDev_tte, vars.discount_, strike[i]);
        if (std::fabs(vars.bpvol_ - impliedBpVol[i])>1.0e-10)
        {
            BOOST_ERROR("Failed, expected " << vars.bpvol_ << " realised " << impliedBpVol[i]);
        }
    }

    cl::TapeFunction<double> f(tteVector, impliedBpVol);

    vars.timeTapeRecording_ = vars.timer.elapsed();
    //End of tape recording. Differentiaion wiil be held with respect to  the independent variables vector
    vars.timeTapeRecording_ = vars.timer.elapsed();
    outPerform.log() << "End of tape recording. Time for tape recording : " << vars.timeTapeRecording_ << std::endl;
    //Store size of tape.
    vars.tapeMemory_.push_back(TapeSize { 1, f.Memory() });

    //Start differentiation in Forward mode
    vector<double> sf_Forward;
    outPerform.log() << "Start differentiation in Forward mode : " << currentTime() << std::endl;
    vars.timeAdjoint_ = gradForward(f, sf_Forward, false, false);
    outPerform.log() << "Time for differentiation in Forward mode : " << vars.timeAdjoint_ << " s" << std::endl;

    //Start differentiation using central finite differences.
    outPerform.log() << "Start differentiation using central finite differences : " << currentTime() << std::endl;
    vars.timer.restart();
    vector<Real> sf_Finite(vars.n_, 0.0);
    Real rightValue, leftValue;
    for (Size i = 0; i < vars.n_; i++)
    {
        tteVector[0] += vars.h_;
        strike[i] = vars.forwardPrice_ - vars.d_[i] * vars.bpvol_ * std::sqrt(tteVector[0]);
        rightValue = vars.bachelierImpliedVol(vars.forwardPrice_, tteVector[0], vars.stdDev_, vars.discount_, strike[i]);

        tteVector[0] -= 2 * vars.h_;
        strike[i] = vars.forwardPrice_ - vars.d_[i] * vars.bpvol_ * std::sqrt(tteVector[0]);
        leftValue = vars.bachelierImpliedVol(vars.forwardPrice_, tteVector[0], vars.stdDev_, vars.discount_, strike[i]);

        sf_Finite[i] = (rightValue - leftValue) / (2 * vars.h_);
    }
    vars.timeAnalytical_ = vars.timer.elapsed();
    outPerform.log() << "Time for differentiation using central finite differences : " << vars.timeAnalytical_ << " s" << std::endl;

    // Check derivatives calculated by forward, reverse and finite differences methods.
    outPerform.log() << "Check derivatives calculated by forward and finite differences methods." << std::endl;
    result = checkWithFiniteDiff(sf_Forward, sf_Finite, 1e-5, 1e-5);

    //Adding new data to the performace result vector
    vars.performanceTime_.push_back(PerformanceTime { vars.timeTapeRecording_, vars.timeAdjoint_, vars.timeAnalytical_, 1 });
#endif
    return result;
}

bool AdjointBlackFormulaTest::testBachelierImpliedVolStdDev()
{
    BOOST_TEST_MESSAGE("Testing Bachelier implied vol...");
    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVarBS vars;
    std::vector<cl::TapeDouble> strike(vars.n_);
    std::vector<cl::TapeDouble> stdDevVector(1);
    stdDevVector[0] = vars.bpvol_*std::sqrt(vars.tte_);
    std::vector<cl::TapeDouble> impliedBpVol(vars.n_);
    cl::AdjointTestOutput outPerform("AdjointBlackFormula//BachelierImpliedVolStdDev");

    //Running Adjoint differentiation with respect to std
    outPerform.log() << "\nNumber of std : k = " << 1 << std::endl;

    // Start timing of tape recording.
    outPerform.log() << "Start taping : " << currentTime() << std::endl;
    vars.timer.restart();

    // Start taping. Declare forward rates as independent variables.
    Independent(stdDevVector);

    // Calculate implied bachelier volatility vector.
    outPerform.log() << "Calculate implied bachelier volatility vector." << std::endl;
    for (Size i = 0; i < vars.n_; i++)
    {
        strike[i] = vars.forwardPrice_ - vars.d_[i] * vars.bpvol_ * std::sqrt(vars.tte_);
        impliedBpVol[i] = vars.bachelierImpliedVol(vars.forwardPrice_, vars.tte_, stdDevVector[0], vars.discount_, strike[i]);
        if (std::fabs(vars.bpvol_ - impliedBpVol[i])>1.0e-10)
        {
            BOOST_ERROR("Failed, expected " << vars.bpvol_ << " realised " << impliedBpVol[i]);
        }
    }

    cl::TapeFunction<double> f(stdDevVector, impliedBpVol);

    vars.timeTapeRecording_ = vars.timer.elapsed();
    outPerform.log() << "End of tape recording. Time for tape recording : " << vars.timeTapeRecording_ << std::endl;
    //Store size of tape.
    vars.tapeMemory_.push_back(TapeSize { 1, f.Memory() });

    //Start differentiation in Forward mode
    vector<double> sf_Forward(vars.n_);
    outPerform.log() << "Start differentiation in Forward mode : " << currentTime() << std::endl;
    vars.timeAdjoint_ = gradForward(f, sf_Forward, false, false);
    outPerform.log() << "Time for differentiation in Forward mode : " << vars.timeAdjoint_ << " s" << std::endl;

    //Start differentiation using central finite differences.
    outPerform.log() << "Start differentiation using central finite differences : " << currentTime() << std::endl;
    vars.timer.restart();
    vector<Real> sf_Finite(vars.n_);
    Real rightValue, leftValue;
    for (Size i = 0; i < vars.n_; i++)
    {
        rightValue = vars.bachelierImpliedVol(vars.forwardPrice_, vars.tte_, stdDevVector[0] + vars.h_, vars.discount_, strike[i]);
        leftValue = vars.bachelierImpliedVol(vars.forwardPrice_, vars.tte_, stdDevVector[0] - vars.h_, vars.discount_, strike[i]);
        sf_Finite[i] = (rightValue - leftValue) / (2 * vars.h_);
    }
    vars.timeAnalytical_ = vars.timer.elapsed();
    outPerform.log() << "Time for differentiation using central finite differences : " << vars.timeAnalytical_ << " s" << std::endl;

    // Check derivatives calculated by forward, reverse and finite differences methods.
    outPerform.log() << "Check derivatives calculated by forward and finite differences methods." << std::endl;
    result = checkWithFiniteDiff(sf_Forward, sf_Finite, vars.tol_, vars.tol_);

    //Adding new data to the performace result vector
    vars.performanceTime_.push_back(PerformanceTime { vars.timeTapeRecording_, vars.timeAdjoint_, vars.timeAnalytical_, 1 });
#endif
    return result;
}

bool AdjointBlackFormulaTest::testBachelierImpliedVolDiscount()
{
    BOOST_TEST_MESSAGE("Testing Bachelier implied vol...Derivatives with respect to discount = 0");
    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVarBS vars;
    std::vector<cl::TapeDouble> strike(vars.n_);
    std::vector<cl::TapeDouble> discountVector(1);
    discountVector[0] = vars.discount_;
    std::vector<cl::TapeDouble> impliedBpVol(vars.n_);
    cl::AdjointTestOutput outPerform("AdjointBlackFormula//BachelierImpliedVolDiscount");

    //Running Adjoint differentiation with respect to discount
    outPerform.log() << "\nNumber of discounts : k = " << 1 << std::endl;

    // Start timing of tape recording.
    outPerform.log() << "Start taping : " << currentTime() << std::endl;
    vars.timer.restart();

    // Start taping. Declare forward rates as independent variables.
    Independent(discountVector);

    // Calculate implied bachelier volatility vector.
    outPerform.log() << "Calculate implied bachelier volatility vector." << std::endl;
    for (Size i = 0; i < vars.n_; i++)
    {
        strike[i] = vars.forwardPrice_ - vars.d_[i] * vars.bpvol_ * std::sqrt(vars.tte_);
        impliedBpVol[i] = vars.bachelierImpliedVol(vars.forwardPrice_, vars.tte_, vars.stdDev_, discountVector[0], strike[i]);
        if (std::fabs(vars.bpvol_ - impliedBpVol[i])>1.0e-10)
        {
            BOOST_ERROR("Failed, expected " << vars.bpvol_ << " realised " << impliedBpVol[i]);
        }
    }

    cl::TapeFunction<double> f(discountVector, impliedBpVol);
    //End of tape recording. Differentiaion wiil be held with respect to  the independent variables vector
    vars.timeTapeRecording_ = vars.timer.elapsed();
    outPerform.log() << "End of tape recording. Time for tape recording : " << vars.timeTapeRecording_ << std::endl;
    //Store size of tape.
    vars.tapeMemory_.push_back(TapeSize { 1, f.Memory() });

    //Start differentiation in Forward mode
    vector<double> sf_Forward(vars.n_);
    outPerform.log() << "Start differentiation in Forward mode : " << currentTime() << std::endl;
    vars.timeAdjoint_ = gradForward(f, sf_Forward, false, false);
    outPerform.log() << "Time for differentiation in Forward mode : " << vars.timeAdjoint_ << " s" << std::endl;

    //Start differentiation using central finite differences.
    outPerform.log() << "Start differentiation using central finite differences : " << currentTime() << std::endl;
    vars.timer.restart();
    vector<Real> sf_Finite(vars.n_);
    Real rightValue, leftValue;
    for (Size i = 0; i < vars.n_; i++)
    {
        rightValue = vars.bachelierImpliedVol(vars.forwardPrice_, vars.tte_, vars.stdDev_, discountVector[0] + vars.h_, strike[i]);
        leftValue = vars.bachelierImpliedVol(vars.forwardPrice_, vars.tte_, vars.stdDev_, discountVector[0] - vars.h_, strike[i]);
        sf_Finite[i] = (rightValue - leftValue) / (2 * vars.h_);
    }
    vars.timeAnalytical_ = vars.timer.elapsed();
    outPerform.log() << "Time for differentiation using central finite differences : " << vars.timeAnalytical_ << " s" << std::endl;

    // Check derivatives calculated by forward, reverse and finite differences methods.
    outPerform.log() << "Check derivatives calculated by forward and finite differences methods." << std::endl;
    result = checkWithFiniteDiff(sf_Forward, sf_Finite, vars.tol_, vars.tol_);

    //Adding new data to the performace result vector
    vars.performanceTime_.push_back(PerformanceTime { vars.timeTapeRecording_, vars.timeAdjoint_, vars.timeAnalytical_, 1 });
#endif
    return result;
}

bool AdjointBlackFormulaTest::testBachelierImpliedVolBpVol()
{
    BOOST_TEST_MESSAGE("Testing Bachelier implied vol...");
    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVarBS vars;
    std::vector<cl::TapeDouble> strike(vars.n_);
    std::vector<cl::TapeDouble> bpvolVector(1);
    bpvolVector[0] = vars.bpvol_;
    cl::AdjointTestOutput outPerform("AdjointBlackFormula//BachelierImpliedVolBpVol");
    std::vector<cl::TapeDouble> impliedBpVol(vars.n_);

    //Running Adjoint differentiation with respect to bpvol (or normal volatility)
    outPerform.log() << "\nNumber of bpvol (or normal volatility) : k = " << 1 << std::endl;

    // Start timing of tape recording.
    outPerform.log() << "Start taping : " << currentTime() << std::endl;
    vars.timer.restart();

    // Start taping. Declare forward rates as independent variables.
    Independent(bpvolVector);

    // Calculate implied bachelier volatility vector.
    outPerform.log() << "Calculate implied bachelier volatility vector." << std::endl;
    Real stdDev_ = bpvolVector[0] * std::sqrt(vars.tte_);
    for (Size i = 0; i < vars.n_; i++)
    {
        strike[i] = vars.forwardPrice_ - vars.d_[i] * bpvolVector[0] * std::sqrt(vars.tte_);
        impliedBpVol[i] = vars.bachelierImpliedVol(vars.forwardPrice_, vars.tte_, stdDev_, vars.discount_, strike[i]);

        if (std::fabs(vars.bpvol_ - impliedBpVol[i])>1.0e-10)
        {
            BOOST_ERROR("Failed, expected " << vars.bpvol_ << " realised " << impliedBpVol[i]);
        }
    }

    cl::TapeFunction<double> f(bpvolVector, impliedBpVol);
    //End of tape recording. Differentiaion wiil be held with respect to  the independent variables vector
    vars.timeTapeRecording_ = vars.timer.elapsed();
    outPerform.log() << "End of tape recording. Time for tape recording : " << vars.timeTapeRecording_ << std::endl;
    //Store size of tape.
    vars.tapeMemory_.push_back(TapeSize { 1, f.Memory() });

    //Start differentiation in Forward mode
    vector<double> sf_Forward(vars.n_);
    outPerform.log() << "Start differentiation in Forward mode : " << currentTime() << std::endl;
    vars.timeAdjoint_ = gradForward(f, sf_Forward, false, false);
    outPerform.log() << "Time for differentiation in Forward mode : " << vars.timeAdjoint_ << " s" << std::endl;

    //Start differentiation using central finite differences.
    outPerform.log() << "Start differentiation using central finite differences : " << currentTime() << std::endl;
    vars.timer.restart();
    vector<Real> sf_Finite(vars.n_);
    Real rightValue, leftValue;
    for (Size i = 0; i < vars.n_; i++)
    {
        stdDev_ = (bpvolVector[0] + vars.h_)* std::sqrt(vars.tte_);
        strike[i] = vars.forwardPrice_ - vars.d_[i] * (bpvolVector[0] + vars.h_) * std::sqrt(vars.tte_);
        rightValue = vars.bachelierImpliedVol(vars.forwardPrice_, vars.tte_, stdDev_, vars.discount_, strike[i]);

        stdDev_ = (bpvolVector[0] - vars.h_)* std::sqrt(vars.tte_);
        strike[i] = vars.forwardPrice_ - vars.d_[i] * (bpvolVector[0] - vars.h_) * std::sqrt(vars.tte_);
        leftValue = vars.bachelierImpliedVol(vars.forwardPrice_, vars.tte_, stdDev_, vars.discount_, strike[i]);

        sf_Finite[i] = (rightValue - leftValue) / (2 * vars.h_);
    }
    vars.timeAnalytical_ = vars.timer.elapsed();
    outPerform.log() << "Time for differentiation using central finite differences : " << vars.timeAnalytical_ << " s" << std::endl;

    // Check derivatives calculated by forward, reverse and finite differences methods.
    outPerform.log() << "Check derivatives calculated by forward and finite differences methods." << std::endl;
    result = checkWithFiniteDiff(sf_Forward, sf_Finite, vars.tol_, vars.tol_);

    //Adding new data to the performace result vector
    vars.performanceTime_.push_back(PerformanceTime { vars.timeTapeRecording_, vars.timeAdjoint_, vars.timeAnalytical_, 1 });
#endif
    return result;
}


test_suite* AdjointBlackFormulaTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Black formula tests");

    suite->add(QUANTLIB_TEST_CASE(&AdjointBlackFormulaTest::testBachelierImpliedVolStrike));
    suite->add(QUANTLIB_TEST_CASE(&AdjointBlackFormulaTest::testBachelierImpliedVolForward));
    suite->add(QUANTLIB_TEST_CASE(&AdjointBlackFormulaTest::testBachelierImpliedVolTte));
    suite->add(QUANTLIB_TEST_CASE(&AdjointBlackFormulaTest::testBachelierImpliedVolStdDev));
    suite->add(QUANTLIB_TEST_CASE(&AdjointBlackFormulaTest::testBachelierImpliedVolDiscount));
    suite->add(QUANTLIB_TEST_CASE(&AdjointBlackFormulaTest::testBachelierImpliedVolD));
    suite->add(QUANTLIB_TEST_CASE(&AdjointBlackFormulaTest::testBachelierImpliedVolBpVol));

    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_blackformula)

BOOST_AUTO_TEST_CASE(testBachelierImpliedVolStrike)
{
    BOOST_CHECK(AdjointBlackFormulaTest::testBachelierImpliedVolStrike());
}

BOOST_AUTO_TEST_CASE(testBachelierImpliedVolForward)
{
    BOOST_CHECK(AdjointBlackFormulaTest::testBachelierImpliedVolForward());
}
BOOST_AUTO_TEST_CASE(testBachelierImpliedVolMaturityTime)
{
    BOOST_CHECK(AdjointBlackFormulaTest::testBachelierImpliedVolTte());
}

BOOST_AUTO_TEST_CASE(testBachelierImpliedVolStdDev)
{
    BOOST_CHECK(AdjointBlackFormulaTest::testBachelierImpliedVolStdDev());
}

BOOST_AUTO_TEST_CASE(testBachelierImpliedVolDiscount)
{
    BOOST_CHECK(AdjointBlackFormulaTest::testBachelierImpliedVolDiscount());
}

BOOST_AUTO_TEST_CASE(testBachelierImpliedVolTermD)
{
    BOOST_CHECK(AdjointBlackFormulaTest::testBachelierImpliedVolD());
}

BOOST_AUTO_TEST_CASE(testBachelierImpliedVolBpVol)
{
    BOOST_CHECK(AdjointBlackFormulaTest::testBachelierImpliedVolBpVol());
}

BOOST_AUTO_TEST_SUITE_END()

#endif