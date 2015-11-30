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

// based on gjrgarchmodel.cpp file from test-suite

#ifndef cl_adjoint_gjrgarch_model_impl_hpp
#define cl_adjoint_gjrgarch_model_impl_hpp
#pragma once

#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include "adjointgjrgarchmodeltest.hpp"
#include "adjointtestbase.hpp"

#include <boost/make_shared.hpp>
#include <ql/quantlib.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;


#define OUTPUT_FOLDER_NAME "AdjointGjrgarchmodel"


namespace
{
    enum
    {
#if defined CL_GRAPH_GEN
        // Number of points for dependency plots.
        pointNo = 15,
        // Number of points for performance plot.
        iterNo = 15,
        // Step for portfolio size for performance testing .
        step = 1,
#else
        // Number of points for dependency plots.
        pointNo = 1,
        // Number of points for performance plot.
        iterNo = 1,
        // Step for portfolio size for performance testing .
        step = 1,
#endif
        // Defines performance accuracy. Its value is a minimum number
        // of calling of O(1) complexity methods per one performance test.
        iterNumFactor = 0,
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

    struct GJRGARCHmodelData
    {

        GJRGARCHmodelData()
            : settlement_(Date(5, July, 2002))
            , dayCounter_(Actual365Fixed())
            , calendar_(TARGET())
            , strike_(3400)
        {
            Settings::instance().evaluationDate() = settlement_;

            t_ = { 7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 98, 105 };

            r_ = { 0.0341
                , 0.0334
                , 0.0323
                , 0.0378
                , 0.0362
                , 0.0305
                , 0.0381
                , 0.0361
                , 0.0395
                , 0.0327
                , 0.0391
                , 0.0302
                , 0.0392
                , 0.0321
                , 0.0318 };

            v_ = { 0.3145
                , 0.3153
                , 0.3382
                , 0.3464
                , 0.3491
                , 0.3827
                , 0.4358
                , 0.4436
                , 0.5467
                , 0.55
                , 0.5604
                , 0.5716
                , 0.6724
                , 0.6895
                , 0.6942};
        }

        // Calculate model calibration error with given volatilities.
        Real calculateGJRGARCHmodelCalibrationError(std::vector<Real> vol)
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
                    s0->value(), strike_, vol,
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

            // Set optimization method.
            Simplex om(0.01);

            // Calibrate model.
            model->calibrate(options, om, EndCriteria(20, 10, 1.0e-2, 1.0e-2, 1.0e-2));

            error = 0;
            // Calculate calibration error after calibration.
            for (Size i = 0; i < options.size(); ++i)
            {
                const Real diff = options[i]->calibrationError()*100.0;
                error += diff*diff;
            }

            return error;
        }

        void calculateFinDiff(std::vector<Real>& vol, double h, std::vector<Real>& sf_Finite)
        {
            Size size = vol.size();

            Real totalCalibrationErrorValue = calculateGJRGARCHmodelCalibrationError(vol);

            for (Size i = 0; i < size; i++)
            {
                vol[i] += h;
                sf_Finite[i] = (calculateGJRGARCHmodelCalibrationError(vol) - totalCalibrationErrorValue) / h;
                vol[i] -= h;
            }
        }

        // Cleanup.
        SavedSettings backup_;

        // Global data.
        Date settlement_;

        DayCounter dayCounter_;
        Calendar calendar_;

        std::vector<Integer> t_;
        std::vector<Real> r_;
        std::vector<Real> v_;

        Real strike_;

        const Real omega_ = 2.0e-6;
        const Real alpha_ = 0.024;
        const Real beta_ = 0.93;
        const Real gamma_ = 0.059;
        const Real lambda_ = 0.1;
        const Real daysPerYear_ = 365.0;
    };


    struct TestData
        : public GJRGARCHmodelData
    {
        struct Test
        : public cl::AdjointTest<Test>
        {
            Test(Size size, TestData* data)
            : size_(size)
            , data_(data)
            , volatility_(size)
            , totalCalibrationErrorValue_()
            {
                setLogger(&data_->outPerform_);

                for (Size i = 0; i < size_; i++)
                {
                    volatility_[i] = data->v_[i];
                }

            }

            Size indepVarNumber() { return size_; }

            Size depVarNumber() { return 1; }

            Size minPerfIteration() { return iterNumFactor; }

            void recordTape()
            {
                cl::Independent(volatility_);
                calculateTotalError();
                f_ = std::make_unique<cl::tape_function<double>>(volatility_, totalCalibrationErrorValue_);
            }

            // Calculates total calibration error.
            void calculateTotalError()
            {
                totalCalibrationErrorValue_.push_back(data_->calculateGJRGARCHmodelCalibrationError(volatility_));
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                double h = 1.0e-10;  // shift for finite diff. method
                analyticalResults_.resize(size_);
                data_->calculateFinDiff(volatility_, h, analyticalResults_);
            }

            double relativeTol() const { return 1e-2; }

            double absTol() const { return 1e-10; }

            Size size_;
            TestData* data_;
            std::vector<cl::tape_double> volatility_;
            std::vector<cl::tape_double> totalCalibrationErrorValue_;
        };

        TestData()
            : GJRGARCHmodelData()

            , outPerform_(OUTPUT_FOLDER_NAME "//"
              , { { "filename", "AdjointPerformance" }
                , { "not_clear", "Not" }
                , { "line_box_width", "-5" }
                , { "title", "Calibration error differentiation performance with respect to volatility" }
                , { "ylabel", "Time (s)" }
                , { "xlabel", "Number of volatilities" }
                , { "smooth", "default" }
                , { "cleanlog", "true" } })

            , outAdjoint_(OUTPUT_FOLDER_NAME "//"
              , { { "filename", "Adjoint" }
                , { "not_clear", "Not" }
                , { "smooth", "2" }
                , { "title", "Calibration error adjoint differentiation performance with respect to volatility" }
                , { "cleanlog", "false" }
                , { "ylabel", "Time (s)" }
                , { "xlabel", "Number of volatilities" } })

            , outSize_(OUTPUT_FOLDER_NAME "//"
              , { { "filename", "TapeSize" }
                , { "not_clear", "Not" }
                , { "title", "Tape size dependence on number of volatilities" }
                , { "cleanlog", "false" }
                , { "smooth", "default" }
                , { "ylabel", "Size (MB)" }
                , { "xlabel", "Number of volatilities" } })

            , out_(OUTPUT_FOLDER_NAME "//output"
              , { { "filename", "CalibrErronVol" }
                , { "ylabel", "Calibration Error" }
                , { "not_clear", "Not" }
                , { "title", "Calibration error on volatility" }
                , { "cleanlog", "false" }
                , { "smooth", "10" }
                , { "xlabel", "Volatility" } })
        {
        }

        bool makeOutput()
        {
            bool ok = true;
            if (pointNo > 0)
            {
                ok &= recordDependencePlot();
            }
            ok &= cl::recordPerformance(*this, iterNo, step);
            return ok;
        }

        std::shared_ptr<Test> getTest(size_t size)
        {
            return std::make_shared<Test>(size, this);
        }

        // Makes plots for strike sensitivity dependence.
        bool recordDependencePlot()
        {
            std::vector<Variation> outData(pointNo);
            auto test = getTest(pointNo);
            for (Size i = 0; i < pointNo; i++)
            {
                outData[i] = { test->volatility_[i], test->data_->calculateGJRGARCHmodelCalibrationError(std::vector<Real>(1, test->volatility_[i])) };
            }
            out_ << outData;
            return true;
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output out_;
    };


    typedef TestData::Test GJRGARCHmodelTest;
}

#endif