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

// Based on marketmodel_smmcapletalphacalibration.cpp from test-suite

// Algorithm of calibration can be found here:
// http://fbe.unimelb.edu.au/__data/assets/pdf_file/0019/806302/172.pdf


#ifndef cl_adjoint_market_model_calibration_impl_hpp
#define cl_adjoint_market_model_calibration_impl_hpp
#pragma once

#include "adjointmarketmodelcalibrationtest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include <ql/models/marketmodels/products/multiproductcomposite.hpp>
#include <ql/quantlib.hpp>
#include <boost/make_shared.hpp>
#include <iostream>

using namespace QuantLib;
using namespace std;
using namespace boost::unit_test_framework;


#define OUTPUT_FOLDER_NAME "AdjointMarketModelCalibration"


namespace
{
    enum
    {
#if defined CL_GRAPH_GEN
        // Number of points for dependency plots.
        pointNo = 10,
        // Start size.
        startSize = 4,
        // Number of points for performance plot.
        iterNo = 14,
        // Step for portfolio size for performance testing .
        step = 1,
        // Defines performance accuracy. Its value is a minimum number
        // of calling of O(1) complexity methods per one performance test.
        iterNumFactor = 1000,
#else
        // Number of points for dependency plots.
        pointNo = 1,
        // Start size.
        startSize = 12,
        // Number of points for performance plot.
        iterNo = 1,
        // Step for portfolio size for performance testing .
        step = 1,
        // Defines performance accuracy. Its value is a minimum number
        // of calling of O(1) complexity methods per one performance test.
        iterNumFactor = 1,
#endif
    };

    struct RealPair
    {
        static std::deque<std::string> get_columns()
        {
            static std::deque<std::string> columns = { " ", "" };
            return columns;
        }

        template <typename stream_type>
        friend inline stream_type& operator << (stream_type& stm, RealPair& p)
        {
            stm << p.x_ << ";" << p.y_ << std::endl;
            return stm;
        }

        Real x_, y_;

    };

    struct CommonVars
    {
        explicit CommonVars(Size size):
          size_(size)
        , displacement_(0.0)
        , rateTimes_(initRateTimes())
        , evolution_(EvolutionDescription(rateTimes_))
        , cs_(createCS())
        , corr_(calcCorrelation())
        , swapVariances_(initSwapVariance())
        , alphaInitial_(size, 0.0)
        , alphaMax_(size, 1.0)
        , alphaMin_(size, -1.0)
        {
        }

        std::vector<Real> initRateTimes()
        {
            Date todaysDate = Date(3, June, 2015);
            Date endDate = todaysDate + (size_ + 1) * Weeks;

            Schedule dates(
                todaysDate,
                endDate,
                Period(Weekly),
                NullCalendar(),
                Following,
                Following,
                DateGeneration::Backward,
                false);

            std::vector<Real> rateTimes(dates.size() - 1);
            DayCounter dayCounter = SimpleDayCounter();
            for (size_t i = 1; i < dates.size(); ++i)
                rateTimes[i - 1] = dayCounter.yearFraction(todaysDate, dates[i]);
            return rateTimes;
        }

        boost::shared_ptr<LMMCurveState> createCS()
        {
            std::vector<Rate> todaysForwards(rateTimes_.size() - 1);
            for (Size i = 0; i < todaysForwards.size(); ++i)
            {
                todaysForwards[i] = 0.03 + 0.0025 * i / size_;
            }
            boost::shared_ptr<LMMCurveState> cs = boost::make_shared<LMMCurveState>(rateTimes_);
            cs->setOnForwardRates(todaysForwards);
            return cs;
        }

        boost::shared_ptr<PiecewiseConstantCorrelation> calcCorrelation()
        {
            Real longTermCorrelation = 0.5;
            Real beta = 0.2;

            boost::shared_ptr<PiecewiseConstantCorrelation> fwdCorr =
                boost::make_shared<ExponentialForwardCorrelation>(
                rateTimes_, longTermCorrelation, beta);

            return boost::make_shared<CotSwapFromFwdCorrelation>(fwdCorr, *cs_, displacement_);
        }

        std::vector<boost::shared_ptr<PiecewiseConstantVariance>> initSwapVariance()
        {
            Real a_ = 0.0, b_ = 0.17, c_ = 1.0, d_ = 0.10;
            std::vector<boost::shared_ptr<PiecewiseConstantVariance>> swapVariances(size_);
            for (Size i = 0; i < size_; ++i)
            {
                swapVariances[i] = boost::make_shared< PiecewiseConstantAbcdVariance>(a_, b_, c_, d_, i, rateTimes_);
            }
            return swapVariances;
        }

        std::vector<Volatility> generateCapletVolatilities()
        {
            std::vector<Real> volsDelta(13, 0.01);
            volsDelta.resize(19);
            volsDelta[13] = 0.0090;
            volsDelta[14] = 0.0057;
            volsDelta[15] = 0.0045;
            volsDelta[16] = 0.0044;
            volsDelta[17] = 0.0015;
            volsDelta[18] = 0.0004;

            std::vector<Volatility> capletVols(size_);
            capletVols[size_ - 1] = swapVariances_[size_ - 1]->totalVolatility(size_ - 1);
            Real delta = (size_ <= 18) ? volsDelta[size_] : 0;
            for (int i = size_ - 2; i >= 0; --i)
            {
                if (i >= size_ / 3)
                    capletVols[i] = capletVols[i + 1] + delta;
                else
                    capletVols[i] = capletVols[i + 1] - delta;
            }
            return capletVols;
        }

        std::vector<Real> calibrationVolatilitiesAndErrors(const std::vector<Volatility>& capletVols)
        {
            // Create calibrator
            CTSMMCapletAlphaFormCalibration calibrator(evolution_,
                                                       corr_,
                                                       swapVariances_,
                                                       capletVols,
                                                       cs_,
                                                       displacement_,
                                                       alphaInitial_,
                                                       alphaMax_,
                                                       alphaMin_,
                                                       false);

            // Calibrate
            Natural numberOfFactors = 3;
            Natural maxIterations = 10;
            Real capletTolerance = 1e-4;
            Natural innerMaxIterations = 100;
            Real innerTolerance = 1e-8;

            bool result = calibrator.calibrate(numberOfFactors,
                                               maxIterations,
                                               capletTolerance,
                                               innerMaxIterations,
                                               innerTolerance);
            if (!result)
                BOOST_ERROR("Calibration failed");

            const std::vector<Matrix>& swapPseudoRoots = calibrator.swapPseudoRoots();
            boost::shared_ptr<MarketModel> smm = boost::make_shared<PseudoRootFacade>(
                swapPseudoRoots,
                rateTimes_,
                cs_->coterminalSwapRates(),
                std::vector<Spread>(size_, displacement_));
            CotSwapToFwdAdapter flmm(smm);
            Matrix capletTotCovariance = flmm.totalCovariance(size_ - 1);

            std::vector<Volatility> capletVolsAndErrors(size_ + 2);
            for (Size i = 0; i < size_; ++i)
            {
                capletVolsAndErrors[i] = std::sqrt(capletTotCovariance[i][i] / rateTimes_[i]);
            }
            capletVolsAndErrors[size_] = calibrator.capletRmsError();
            capletVolsAndErrors[size_ + 1] = calibrator.capletMaxError();

            return capletVolsAndErrors;
        }

        Size size_;
        Spread displacement_;
        std::vector<Time> rateTimes_;
        EvolutionDescription evolution_;
        boost::shared_ptr<LMMCurveState> cs_;
        boost::shared_ptr<PiecewiseConstantCorrelation> corr_;
        std::vector<boost::shared_ptr<PiecewiseConstantVariance>> swapVariances_;
        std::vector<Real> alphaInitial_, alphaMax_, alphaMin_;
    };


    struct TestData
    {
        struct MarketModelTest
        : public cl::AdjointTest<MarketModelTest>
        {
            static const CalcMethod default_method = other;

            explicit MarketModelTest(Size size, TestData* data)
            : size_(size)
            , data_(data)
            , capletVols_(size)
            , doubleCapletVols_(size)
            , volsAndErrs_()
            , vars(size)
            {
                setLogger(&data_->outPerform_);
                capletVols_ = vars.generateCapletVolatilities();
                for (Size i = 0; i < size_; i++)
                    doubleCapletVols_[i] = double(capletVols_[i]);
            }

            Size indepVarNumber() { return size_; }

            Size depVarNumber() { return size_ + 2; }

            Size minPerfIteration() { return iterNumFactor; }

            void recordTape()
            {
                cl::Independent(capletVols_);
                calibrate();
                f_ = std::make_unique<cl::tape_function<double>>(capletVols_, volsAndErrs_);
            }

            // Calculates price of portfolio and each option.
            void calibrate()
            {
                volsAndErrs_ = vars.calibrationVolatilitiesAndErrors(capletVols_);
            }

            void calcAdjoint()
            {
                adjointResults_ = f_->Jacobian(doubleCapletVols_);
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                double h = 1.0e-7;  // shift for finite diff. method
                analyticalResults_.resize(size_*(size_+2));
                for (size_t i = 0; i < size_ - 1; ++i)
                {
                    capletVols_[i] += h;
                    std::vector<Real> volsAndErrsRight = vars.calibrationVolatilitiesAndErrors(capletVols_);
                    capletVols_[i] -= 2*h;
                    std::vector<Real> volsAndErrsLeft = vars.calibrationVolatilitiesAndErrors(capletVols_);
                    for (size_t j = 0; j < volsAndErrsLeft.size(); ++j)
                    {
                        analyticalResults_[i + j * size_] = (volsAndErrsRight[j] - volsAndErrsLeft[j]) / (2*h);
                    }
                    capletVols_[i] += h;
                }
            }

            double relativeTol() const { return 1.0e-2; }

            double absTol() const { return 1.0e-3; }

            Size size_;
            TestData* data_;
            std::vector<cl::tape_double> capletVols_;
            std::vector<double> doubleCapletVols_;
            std::vector<cl::tape_double> volsAndErrs_;
            CommonVars vars;
        };

        TestData()
            : outPerform_(OUTPUT_FOLDER_NAME "//",
            { { "filename", "AdjointPerformance" }
            , { "not_clear", "Not" }
            , { "title", "Caplet volatilities from CTSMMCapletAlphaFormCalibration \\ndifferentiation performance" }
            , { "xlabel", "Number of caplet volatilities" }
            , { "ylabel", "Time (s)" }
            , { "line_box_width", "-5" }
            , { "smooth", "default" }
            , { "cleanlog", "true" }
            })
            , outAdjoint_(OUTPUT_FOLDER_NAME "//",
            {{ "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "title", "Caplet volatilities from CTSMMCapletAlphaFormCalibration \\nadjoint differentiation performance" }
            , { "xlabel", "Number of caplet volatilities" }
            , { "ylabel", "Time (s)" }
            , { "smooth", "default" }
            , { "cleanlog", "false" }
            })
            , outSize_(OUTPUT_FOLDER_NAME "//",
            { { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "title", "Tape size dependence on number of caplet volatilities" }
            , { "xlabel", "Number of caplet volatilities" }
            , { "ylabel", "Size (MB)" }
            , { "cleanlog", "false" }
            })
            , out_(OUTPUT_FOLDER_NAME "//output//",
            { { "filename", "VolOnCapletVol" }
            , { "not_clear", "Not" }
            , { "title", "Volatility on caplet volatilities from CTSMMCapletAlphaFormCalibration" }
            , { "xlabel", "Caplet volatility" }
            , { "ylabel", "Volatility" }
            , { "smooth", "default" }
            , { "cleanlog", "false" }
            })
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

        std::shared_ptr<MarketModelTest> getTest(size_t size)
        {
            return std::make_shared<MarketModelTest>(startSize + size, this);
        }

        bool recordDependencePlot()
        {
            std::vector<RealPair> outData(pointNo);
            auto test = getTest(pointNo);
            for (Size i = 0; i < pointNo; i++)
            {
                test->capletVols_[0] += 0.00001;
                test->calibrate();
                outData[i] = { test->capletVols_[0], test->volsAndErrs_[0] };

            }
            out_ << outData;
            return true;
        }
        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output out_;
    };

    typedef TestData::MarketModelTest MarketModelTest;
}
#endif