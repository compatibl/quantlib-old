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

#ifndef cl_adjoint_shortratemodels_impl_hpp
#define cl_adjoint_shortratemodels_impl_hpp
#pragma once

#include <ql/quantlib.hpp>
#include "utilities.hpp"
#include "adjointshortratemodelstest.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include <boost/shared_ptr.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;

#define OUTPUT_FOLDER_NAME "AdjointShortRateModels"

namespace
{
    struct CalibrationData
    {
        Integer start_;
        Integer length_;
        Real volatility_;
    };

    struct VolatilityDependence
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
            operator << (stream_type& stm, VolatilityDependence& v)
        {
                stm << v.inputVolatility_
                    << ";" << v.modelSigma_ << std::endl;

                return stm;
            }

        Real inputVolatility_;
        Real modelSigma_;
    };

    struct ModelData
    {
        ModelData()
        : today_(15, February, 2002)
        , settlement_(19, February, 2002)
        , termStructure_(flatRate(settlement_, 0.04875825, Actual365Fixed()))
        , model_(new HullWhite(termStructure_))
        , index_(new Euribor6M(termStructure_))
        , engine_(new JamshidianSwaptionEngine(model_))
        , swaptions_()
        , data_()
        {
            Settings::instance().evaluationDate() = today_;
        }

        void calibrationData(std::vector<cl::tape_double>& vol, std::vector<CalibrationData>& data)
        {
            Integer pos = 0;
            Integer size = vol.size();
            data.clear();
            for (std::vector<cl::tape_double>::iterator it = vol.begin(); it != vol.end(); it++, pos++)
            {
                data.push_back(CalibrationData { pos + 1, size - pos, *it });
            }
        }

        std::string getPath(Size calibrationType)
        {
            if (calibrationType == 1)
                return "\\TestCachedHullWhite";
            else
                return "\\TestCachedHullWhiteFixedReversion";
        }

        // Calibration type 1: Testing Hull-White calibration against cached values using swaptions with start delay
        // type 2: Testing Hull-White calibration with fixed reversion against cached values.
        Real calibrate(std::vector<Real>& volatility
                       , Size calibrationType)
        {
            std::vector<CalibrationData> data;
            calibrationData(volatility, data);
            Size sizeof_indep = data.size();
            swaptions_.clear();
            for (Size i = 0; i < sizeof_indep; i++)
            {
                boost::shared_ptr<Quote> volatil(new SimpleQuote(data[i].volatility_));
                boost::shared_ptr<CalibrationHelper> helper(
                    new SwaptionHelper(Period(data[i].start_, Years),
                    Period(data[i].length_, Years),
                    Handle<Quote>(volatil),
                    index_,
                    Period(1, Years), Thirty360(),
                    Actual360(), termStructure_));
                helper->setPricingEngine(engine_);
                swaptions_.push_back(helper);
            }

            // Set up the optimization problem
            // Real simplexLambda = 0.1;
            // Simplex optimizationMethod(simplexLambda);
            LevenbergMarquardt optimizationMethod(1.0e-8, 1.0e-8, 1.0e-8);
            EndCriteria endCriteria(1000, 100, 1e-6, 1e-8, 1e-8);

            // Optimize.
            switch (calibrationType)
            {

                case 1:
                    model_->calibrate(swaptions_, optimizationMethod, endCriteria);
                    break;
                case 2:
                    model_->calibrate(swaptions_, optimizationMethod, endCriteria, Constraint(), std::vector<Real>(),
                                     HullWhite::FixedReversion());
                    // The difference is in the choice of  HullWhite::FixedReversion() to calibrate the model below.
            }

            return  model_->sigma();
        }

        void finiteDiff(std::vector<Real>& vol
                          , Size calibrationType
                          , double h
                          , std::vector<Real>& sfFinite)
        {
            std::vector<Real>::iterator it;
            std::vector<Real>::iterator it_fin;
            for (it = vol.begin(), it_fin = sfFinite.begin(); it != vol.end(); it++, it_fin++)
            {
                *it -= h;
                calibrate(vol, calibrationType);
                cl::tape_double stepbackward = model_->sigma();
                *it += 2 * h;
                calibrate(vol, calibrationType);
                cl::tape_double stepforward = model_->sigma();
                *it_fin = (stepforward - stepbackward) / (2 * h);
                *it -= h;
            }
        }

        // Common data.
        Date today_;
        Date settlement_;
        Handle<YieldTermStructure> termStructure_;
        boost::shared_ptr<HullWhite> model_;
        boost::shared_ptr<IborIndex> index_;
        boost::shared_ptr<PricingEngine> engine_;
        std::vector<boost::shared_ptr<CalibrationHelper>> swaptions_;
        std::vector<CalibrationData> data_;

    };

    struct CachedHullWhiteTestData
        : public ModelData
    {
        struct Test
        : public cl::AdjointTest<Test>
        {
            Test(Size size, CachedHullWhiteTestData* data)
            : size_(size)
            , data_(data)
            , volatility_(size)
            , iterNumFactor_(1)
            {
                setLogger(&data_->outPerform_);

                for (Size i = 0; i < size; i++)
                    volatility_[i] = 0.1148 - 0.004*i;
            }

            Size indepVarNumber() { return size_;}

            Size depVarNumber() { return 1; }

            Size minPerfIteration() { return iterNumFactor_; }

            void recordTape()
            {
                cl::Independent(volatility_);
                calculateModelSigma();
                f_ = std::make_unique<cl::tape_function<double>>(volatility_, modelSigma_);
            }

            // Calculates total calibration error.
            void  calculateModelSigma()
            {
                modelSigma_.push_back(data_->calibrate(volatility_, data_->calibrationType_));
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                // Shift for finite diff. method.
                double h = 1.0e-3;
                analyticalResults_.resize(size_);
                data_->finiteDiff(volatility_, data_->calibrationType_, h, analyticalResults_);
            }

            double relativeTol() const { return 1e-2; }

            double absTol() const { return 1e-2; }

            Size size_;
            Size iterNumFactor_;
            CachedHullWhiteTestData* data_;
            std::vector<cl::tape_double> volatility_;
            std::vector<cl::tape_double> modelSigma_;
        };

        CachedHullWhiteTestData(Size calibrationType)
            : ModelData()
            , calibrationType_(calibrationType)

            , outPerform_(OUTPUT_FOLDER_NAME + getPath(calibrationType)
            , { { "filename", "AdjointPerformance" }
              , { "not_clear", "Not" }
              , { "line_box_width", "-5" }
              , { "cleanlog", "true" }
              , { "title", "Model sigma differentiation performance with respect to calibration implied swaption volatilities" }
              , { "xlabel", "Size of  calibration  volatilities vector" }
              , { "ylabel", "Time (s)" }
              , { "smooth", "12" } })

            , outAdjoint_(OUTPUT_FOLDER_NAME + getPath(calibrationType)
            , { { "filename", "Adjoint" }
              , { "not_clear", "Not" }
              , { "cleanlog", "false" }
              , { "title", "Adjoint differentiation performance with respect to calibration implied swaption volatilities" }
              , { "xlabel", "Size of  calibration  volatilities vector" }
              , { "ylabel", "Time (s)" }
              , { "smooth", "12" } })

            , outSize_(OUTPUT_FOLDER_NAME + getPath(calibrationType)
            , { { "filename", "TapeSize" }
              , { "not_clear", "Not" }
              , { "cleanlog", "false" }
              , { "title", "Tape size dependence on  size of  calibration  volatilities vector" }
              , { "xlabel", "Size of  calibration  volatilities vector" }
              , { "ylabel", "Memory (MB)" }
              , { "smooth", "12" } })

            , out_(OUTPUT_FOLDER_NAME + getPath(calibrationType) + "//output"
            , { { "filename", "SigmaonVolatil" }
              , { "not_clear", "Not" }
              , { "title", "Model sigma on calibration swaption implied volatilty dependence" }
              , { "xlabel", "Volatility" }
              , { "ylabel", "Model Sigma" }
              , { "cleanlog", "false" } })

#if defined CL_GRAPH_GEN
            , pointNo_(20)
            , iterNo_(20)
            , step_(1)
#else
            , pointNo_(1)
            , iterNo_(1)
            , step_(1)
#endif
        {
        }

        bool makeOutput()
        {
            bool ok = true;
            if (pointNo_ > 0)
            {
                ok &= recordDependencePlot();
            }
            ok &= cl::recordPerformance(*this, iterNo_, step_);
            return ok;
        }

        std::shared_ptr<Test> getTest(size_t size)
        {
            return std::make_shared<Test>(size + 4, this);
        }

        // Makes plots for strike sensitivity dependence.
        bool recordDependencePlot()
        {
            std::vector<VolatilityDependence> outData(pointNo_);
            auto test = getTest(pointNo_);
            for (Size i = 0; i < pointNo_; i++)
            {
                Real vol = 0.06 + 0.0003*i*5;
                test->volatility_[0] = vol;
                outData[i] = { vol, test->data_->calibrate(test->volatility_, test->data_->calibrationType_) };
            }
            out_ << outData;

            // Reset volatility_[0] to initial value.
            test->volatility_[0] = 0.1148;

            return true;
        }

        Size pointNo_;
        Size iterNo_;
        Size step_;

        Size calibrationType_;

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output out_;
    };
    typedef CachedHullWhiteTestData::Test CachedHullWhiteTest;

    struct FuturesConvexityBiasTest
        : public cl::AdjointTest<FuturesConvexityBiasTest>
    {
        FuturesConvexityBiasTest()
        : iterNumFactor_(1)
        , futureQuote_(94.0)
        , a_(0.03)
        , sigma_(0.015)
        , t_(5.0)
        , T_(5.25)
        , parameters_(3)
        , calculatedFunction_(1)
        {
            parameters_ = { a_, sigma_, futureQuote_ };
        }

        Size indepVarNumber() { return 3; }

        Size depVarNumber() { return 1; }

        Size minPerfIteration() { return iterNumFactor_; }

        void recordTape()
        {
            cl::Independent(parameters_);
            calculateFunction();
            f_ = std::make_unique<cl::tape_function<double>>(parameters_, calculatedFunction_);
        }

        void calculateFunction()
        {
            calculatedFunction_[0] = (100.0 - parameters_[2]) / 100.0 -
                HullWhite::convexityBias(parameters_[2], t_, T_, parameters_[1], parameters_[0]);
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            double h = 1e-4;
            Size size = indepVarNumber();
            analyticalResults_.resize(size);
            for (Size i = 0; i < size; i++)
            {
                parameters_[i] -= h;
                Real stepbackward = (100.0 - parameters_[2]) / 100.0 -
                    HullWhite::convexityBias(parameters_[2], t_, T_, parameters_[1], parameters_[0]);
                parameters_[i] += 2 * h;
                Real stepforward = (100.0 - parameters_[2]) / 100.0 -
                    HullWhite::convexityBias(parameters_[2], t_, T_, parameters_[1], parameters_[0]);
                analyticalResults_[i] = (stepforward - stepbackward) / (2 * h);
                parameters_[i] -= h;
            }

        }

        double relativeTol() const { return 1e-2; }

        double absTol() const { return 1e-4; }

        Size iterNumFactor_;

        Real futureQuote_;
        Real a_;
        Real sigma_;
        Time t_;
        Time T_;

        std::vector<cl::tape_double> parameters_;
        std::vector<cl::tape_double> calculatedFunction_;
    };
}

#endif