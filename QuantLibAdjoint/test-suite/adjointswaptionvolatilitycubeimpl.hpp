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

// Based on swaptionvolatilitycube.cpp file from test-suite.
#ifndef cl_adjoint_swaption_volatility_cube_impl_hpp
#define cl_adjoint_swaption_volatility_cube_impl_hpp
#pragma once


#include "adjointswaptionvolatilitycubetest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include <test-suite/swaptionvolstructuresutilities.hpp>
#include <ql/termstructures/volatility/swaption/swaptionvolcube2.hpp>
#include <ql/indexes/swap/euriborswap.hpp>
#include <boost/make_shared.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;


#define OUTPUT_FOLDER_NAME "SwaptionVolatilityCube"

/*! The swaption vol cube is made up of ordered swaption vol surface
layers, each layer refers to a swap index of a given length
(in years), all indexes belong to the same family. In order
to identify the family (and its market conventions) an index of
whatever length from that family must be passed in as
swapIndexBase.

Tested functions are volatility and blackvariance.

Output of Volatility function is implied volatility.
Input components are strike rate, option time, and swaption period.
Implied volatility was obtained using interpolation from volatility cube values.
Volatility is a piecewise linear function of the strike rate. And the Volatility function is differentiated
with respect to strike.
Also Volatility is differentiated with respect to volatilities from volatility cube obtained using interpolation.

Output of BlackVariance function is Black variance.
Input components are strike rate, option time, and swaption period.
Black variance was derived from the formula:
volatility * volatility * optionTime.
The BlackVariance function is differentiated with respect to strike.
Derivative of blackVariance with respect to strike:
2 * volatility * optionTime * (volatility)'
Also blackVariance is differentiated with respect to volatilities from volatility cube
because volatility is differentiated  with respect to them.
*/

namespace
{
    enum
    {
#if defined CL_GRAPH_GEN
        // Number of points for dependency plots.
        pointNo = 100,
        // Number of points for performance plot.
        iterNo = 20,
        // Step for portfolio size for performance testing .
        step = 50,
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
        iterNumFactor = 1000,
    };


    enum TestType { VolAtmVol, StrikeAtmVol, VolBlackVariance, StrikeBlackVariance };

    std::string toString(TestType type)
    {
        switch (type)
        {
        case VolAtmVol:
            return "ImpliedVolatilitiesDiffVols";
        case StrikeAtmVol:
            return "ImpliedVolatilitiesDiffStrike";
        case VolBlackVariance:
            return "BlackVariancesDiffVols";
        case StrikeBlackVariance:
            return "BlackVariancesDiffStrike";
        default:
            return "";
        }
    }

    // Output variable which is differentiated.
    std::string variableType(TestType type)
    {
        switch (type)
        {
        case VolAtmVol:
            return "Volatility";
        case StrikeAtmVol:
            return "Volatility";
        case VolBlackVariance:
            return "BlackVariance";
        case StrikeBlackVariance:
            return "BlackVariance";
        default:
            return "";
        }
    }

    struct CommonVars1
    {
        CommonVars1()
        : backup_()
        , conventions_()
        , vegaWeighedSmileFit_(false)
        , atm_()
        , cube_()
        , termStructure_(flatRate(0.05, Actual365Fixed()))
        {
            conventions_.setConventions();

            // ATM swaptionvolmatrix.
            atm_.setMarketData();

            // Swaptionvolcube.
            cube_.setMarketData();
        }

        // Cleanup.
        SavedSettings backup_;
        // Global data.
        SwaptionMarketConventions conventions_;
        bool vegaWeighedSmileFit_;
        AtmVolatility atm_;
        VolatilityCube cube_;
        RelinkableHandle<YieldTermStructure> termStructure_;
    };

    struct CommonVars : public CommonVars1
    {
        CommonVars()
        :CommonVars1()
        , atmVolMatrix_(boost::make_shared<SwaptionVolatilityMatrix>(conventions_.calendar,
        conventions_.optionBdc,
        atm_.tenors.options,
        atm_.tenors.swaps,
        atm_.volsHandle,
        conventions_.dayCounter))
        , swapIndexBase_(boost::make_shared<EuriborSwapIsdaFixA>(2 * Years, termStructure_))
        , shortSwapIndexBase_(boost::make_shared<EuriborSwapIsdaFixA>(1 * Years, termStructure_))
        {
        }

        RelinkableHandle<SwaptionVolatilityStructure> atmVolMatrix_;
        boost::shared_ptr<SwapIndex> swapIndexBase_;
        boost::shared_ptr<SwapIndex> shortSwapIndexBase_;
    };

    // Struct for graphics recording.
    struct StrikeVariance
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "Strike rate", ""
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, StrikeVariance& v)
        {
                stm << v.strike_
                    << ";" << v.blackVariance_
                    << std::endl;
                return stm;
            }

        Real strike_;
        Real blackVariance_;
    };


    struct SwaptionVolatilityCubeTestData : public CommonVars
    {
        explicit SwaptionVolatilityCubeTestData(TestType type, Real shift = 1e-10)
        : CommonVars()
        , test_type_(type)
        , generalVar_()
        , shift_(shift)
        , optionSize_(atm_.tenors.options.size())
        , swapsSize_(atm_.tenors.swaps.size())
        {}

        // Calculate volatility in case of test of function Volatility
        // and blackVariance in case of test of function blackVariance.
        Volatility calcVol(const std::vector<std::vector<Handle<Quote> > >& vols,
            const Period& optionTenor, const Period& swapTenor)
        {
            RelinkableHandle<SwaptionVolatilityStructure> atmVolMatrix = RelinkableHandle<SwaptionVolatilityStructure>(
                boost::make_shared<SwaptionVolatilityMatrix>(conventions_.calendar,
                conventions_.optionBdc,
                atm_.tenors.options,
                atm_.tenors.swaps,
                vols,
                conventions_.dayCounter));
            return calcVol(atmVolMatrix, optionTenor, swapTenor);
        }

        // Calculate volatility in case of test of function Volatility
        // and blackVariance in case of test of function blackVariance.
        Volatility calcVol(const RelinkableHandle<SwaptionVolatilityStructure>& atmVolMatrix,
            const Period& optionTenor, const Period& swapTenor)
        {
            SwaptionVolCube2 volCube(atmVolMatrix,
                cube_.tenors.options,
                cube_.tenors.swaps,
                cube_.strikeSpreads,
                cube_.volSpreadsHandle,
                swapIndexBase_,
                shortSwapIndexBase_,
                vegaWeighedSmileFit_);

            return calcVol(volCube, optionTenor, swapTenor);
        }

        // Calculate volatility in case of test of function Volatility
        // and blackVariance in case of test of function blackVariance.
        Volatility calcVol(const SwaptionVolCube2& volCube, const Period& optionTenor, const Period& swapTenor)
        {
            Rate strike = volCube.atmStrike(optionTenor,
                swapTenor);

            if (test_type_ == VolAtmVol)
                return volCube.volatility(optionTenor, swapTenor, strike, true);
            if (test_type_ == VolBlackVariance)
                return volCube.blackVariance(optionTenor, swapTenor, strike, true);

            return Volatility(0.0);
        }


        TestType test_type_;
        std::vector<Real> generalVar_;
        Real shift_;
        unsigned optionSize_, swapsSize_;
    };

    struct SwaptionVolatilityCubeTestDiffVols
        : public cl::AdjointTest<SwaptionVolatilityCubeTestDiffVols>
    {
        static const CalcMethod default_method = other;

        // Constructor parameters:
        // data to construct swaption volatility cube
        // logger - pointer to output stream that contains log stream.
        explicit SwaptionVolatilityCubeTestDiffVols(SwaptionVolatilityCubeTestData* data, cl::tape_empty_test_output* logger = nullptr)
        : AdjointTest()
        , data_(data)
        , vols_()
        , indepDouble_()
        , BlackVariance_()
        , swapSize_(data_->atm_.tenors.swaps.size())
        , optionSize_(data_->atm_.tenors.options.size())
        , optionTenorNumber_(optionSize_ / 2)
        , optionTenor_(data_->atm_.tenors.options[optionTenorNumber_])
        {
            setUniformData();
        }

        Size indepVarNumber() { return swapSize_; }

        Size depVarNumber() { return swapSize_; }

        Size minPerfIteration() { return 0; }

        double absTol() const { return 1e-4; }

        void recordTape()
        {
            cl::Independent(vols_[optionTenorNumber_]);
            calculateBlackVar();
            f_ = std::make_unique<cl::tape_function<double>>(vols_[optionTenorNumber_], BlackVariance_);
        }

        // Calculate volatility or blackVariance depending on test type
        // for different volatility cubes based on different input volatilities.
        void calculateBlackVar()
        {
            std::vector<std::vector<Handle<Quote> > > vols;
            vols.resize(optionSize_);
            for (Size i = 0; i < optionSize_; i++) {
                vols[i].resize(swapSize_);
                for (Size j = 0; j < swapSize_; j++)
                {
                    vols[i][j] = Handle<Quote>(boost::make_shared<SimpleQuote>(vols_[i][j]));
                }
            }

            BlackVariance_.resize(swapSize_);
            for (int i = 0; i < swapSize_; i++)
            {
                BlackVariance_[i] = data_->calcVol(vols, optionTenor_, data_->atm_.tenors.swaps[i]);
            }
        }



        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            analyticalResults_.resize(swapSize_ * swapSize_, 0.0);
            Real h = data_->shift_;
            std::vector<std::vector<cl::tape_double>> volsLeft, volsRight;
            volsLeft = vols_;
            volsRight = vols_;
            for (Size k = 0; k < swapSize_; k++)
            {
                volsRight[optionTenorNumber_][k] += h;
                volsLeft[optionTenorNumber_][k] -= h;

                std::vector<std::vector<Handle<Quote> > > volsShiftedRight;
                volsShiftedRight.resize(optionSize_);
                std::vector<std::vector<Handle<Quote> > > volsShiftedLeft;
                volsShiftedLeft.resize(optionSize_);
                for (Size i = 0; i < optionSize_; i++) {
                    volsShiftedRight[i].resize(swapSize_);
                    volsShiftedLeft[i].resize(swapSize_);
                    for (Size j = 0; j < swapSize_; j++)
                    {
                        volsShiftedRight[i][j] = Handle<Quote>(boost::make_shared<SimpleQuote>(volsRight[i][j]));
                        volsShiftedLeft[i][j] = Handle<Quote>(boost::make_shared<SimpleQuote>(volsLeft[i][j]));
                    }
                }

                analyticalResults_[k*swapSize_ + k] = (data_->calcVol(volsShiftedRight, optionTenor_, data_->atm_.tenors.swaps[k]) -
                    data_->calcVol(volsShiftedLeft, optionTenor_, data_->atm_.tenors.swaps[k])) / (2 * h);

                volsRight[optionTenorNumber_][k] -= h;
                volsLeft[optionTenorNumber_][k] += h;
            }
        }

        // Calculates derivatives using adjoint.
        void calcAdjoint()
        {
            adjointResults_ = f_->Jacobian(indepDouble_);
        }

        void setUniformData()
        {
            vols_.resize(optionSize_);
            for (Size i = 0; i < data_->optionSize_; i++) {
                vols_[i].resize(swapSize_);
                for (Size j = 0; j < swapSize_; j++)
                {
                    vols_[i][j] = data_->atm_.vols[i][j];
                }
            }

            indepDouble_.resize(swapSize_);
            for (Size i = 0; i < swapSize_; i++)
            {
                indepDouble_[i] = double(data_->atm_.vols[optionTenorNumber_][i]);
            }
        }


        SwaptionVolatilityCubeTestData* data_;
        std::vector<std::vector<cl::tape_double>> vols_;
        std::vector<double> indepDouble_;
        std::vector<Real> BlackVariance_;
        Size swapSize_;
        Size optionSize_;
        Size optionTenorNumber_;
        Period optionTenor_;
    };

    struct SwaptionVolatilityCubeTest
        : public cl::AdjointTest<SwaptionVolatilityCubeTest>
    {
        // Constructor parameters:
        // size - size of strike rates vector for perfomance plot
        // data to construct swaption volatility cube
        // logger - pointer to output stream which contains log stream.
        explicit SwaptionVolatilityCubeTest(Size size, SwaptionVolatilityCubeTestData* data, cl::tape_empty_test_output* logger = nullptr)
        : AdjointTest()
        , size_(size)
        , strikes_(size)
        , data_(data)
        , BlackVariance_()
        , SumBlackVariance_()
        , optPer_(data_->optionSize_ / 2)
        , swPer_(data_->swapsSize_ / 2)
        {
            assert(size_);
            setLogger(logger);
            if (size_ != 1)
            {
                for (Size i = 0; i < size_; i++)
                {
                    strikes_[i] = (0.4 + (0.2 * i) / (size_ - 1));
                }
            }
            else
            {
                strikes_[0] = 0.5;
            }
        }

        Size indepVarNumber() { return size_; }

        Size depVarNumber() { return 1; }

        Size minPerfIteration() { return iterNumFactor; }

        double absTol() const { return 1e-4; }

        void recordTape()
        {
            Independent(strikes_);
            calculateBlackVariance();
            f_ = std::make_unique<cl::tape_function<double>>(strikes_, SumBlackVariance_);
        }

        // Calculate volatility in case of test of function Volatility
        // and blackVariance in case of test of function blackVariance.
        void calculateBlackVariance()
        {
            SwaptionVolCube2 volCube(data_->atmVolMatrix_,
                data_->cube_.tenors.options,
                data_->cube_.tenors.swaps,
                data_->cube_.strikeSpreads,
                data_->cube_.volSpreadsHandle,
                data_->swapIndexBase_,
                data_->shortSwapIndexBase_,
                data_->vegaWeighedSmileFit_);

            SumBlackVariance_.resize(1, 0);
            BlackVariance_.resize(size_);
            for (int i = 0; i < size_; i++)
            {
                if (data_->test_type_ == StrikeAtmVol)
                {
                    // Calculate volatility.
                    BlackVariance_[i] = volCube.volatility(data_->atm_.tenors.options[optPer_],
                        data_->atm_.tenors.swaps[swPer_],
                        strikes_[i], true);
                    SumBlackVariance_[0] += BlackVariance_[i];
                }

                if (data_->test_type_ == StrikeBlackVariance)
                {
                    // Calculate black variance.
                    BlackVariance_[i] += volCube.blackVariance(data_->atm_.tenors.options[optPer_],
                        data_->atm_.tenors.swaps[swPer_],
                        strikes_[i], true);
                    SumBlackVariance_[0] += BlackVariance_[i];
                }
            }
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            SwaptionVolCube2 volCube(data_->atmVolMatrix_,
                data_->cube_.tenors.options,
                data_->cube_.tenors.swaps,
                data_->cube_.strikeSpreads,
                data_->cube_.volSpreadsHandle,
                data_->swapIndexBase_,
                data_->shortSwapIndexBase_,
                data_->vegaWeighedSmileFit_);

            analyticalResults_.resize(size_);
            for (int i = 0; i < size_; i++)
            {
                if (data_->test_type_ == StrikeAtmVol)
                {
                    // Calculate volatility.
                    Real blackDown = volCube.volatility(data_->atm_.tenors.options[optPer_],
                        data_->atm_.tenors.swaps[swPer_],
                        strikes_[i] - data_->shift_, true);

                    Real blackShifted = volCube.volatility(data_->atm_.tenors.options[optPer_],
                        data_->atm_.tenors.swaps[swPer_],
                        strikes_[i] + data_->shift_, true);

                    analyticalResults_[i] = (blackShifted - blackDown) / (2 * data_->shift_);
                }

                if (data_->test_type_ == StrikeBlackVariance)
                {
                    // Calculate black variance.
                    Real black = volCube.blackVariance(data_->atm_.tenors.options[optPer_],
                        data_->atm_.tenors.swaps[swPer_],
                        strikes_[i], true);

                    Real blackShifted = volCube.blackVariance(data_->atm_.tenors.options[optPer_],
                        data_->atm_.tenors.swaps[swPer_],
                        strikes_[i] + data_->shift_, true);

                    analyticalResults_[i] = (blackShifted - black) / data_->shift_;
                }
            }
        }


        SwaptionVolatilityCubeTestData* data_;
        Size size_;
        std::vector<cl::tape_double> strikes_;
        std::vector<Real> BlackVariance_;
        std::vector<Real> SumBlackVariance_;
        int optPer_;
        int swPer_;
    };

    struct TestData
        : public SwaptionVolatilityCubeTestData
    {
        explicit TestData(TestType type, Real shift = 1e-10)
        : SwaptionVolatilityCubeTestData(type, shift)
        , outPerform_title_("Swaption  " + variableType(type) + " differentiation performance with respect to strike rate")
        , outAdjoint_title_("Swaption " + variableType(type) + " adjoint differentiation with respect to fixed strike rate")
        , out_title_(variableType(type) + " dependence on strike rate")
        , out_ylabel_(variableType(type))
        , outPerform_(OUTPUT_FOLDER_NAME "//" + toString(type), {
            { "title", outPerform_title_ }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of strike rates" }
            , { "line_box_width", "-5" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        })
            , outAdjoint_(OUTPUT_FOLDER_NAME "//" + toString(type), {
                { "title", outAdjoint_title_ }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of strikes" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        })
            , outSize_(OUTPUT_FOLDER_NAME "//" + toString(type), {
                { "title", "Tape size dependence on number of strike rates" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
            , { "cleanlog", "false" }
        })


            , out_(OUTPUT_FOLDER_NAME "//" + toString(type) + "//output", {
                { "filename", "StrikeDependence" }
            , { "not_clear", "Not" }
            , { "title", out_title_ }
            , { "ylabel", out_ylabel_ }
            , { "xlabel", "Strike rate" }
            , { "cleanlog", "false" }
        })
        { }

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



        std::shared_ptr<SwaptionVolatilityCubeTest> getTest(size_t size)
        {
            return std::make_shared<SwaptionVolatilityCubeTest>(size, this, &outPerform_);
        }

        // Makes plots for dependence BlackVariance on strikes rates.
        bool recordDependencePlot()
        {
            std::vector<StrikeVariance> outData(pointNo);
            auto test = getTest(pointNo);
            bool ok = test->testReverse();
            for (Size i = 0; i < pointNo; i++)
            {
                outData[i] = { test->strikes_[i], test->BlackVariance_[i] };
            }

            out_ << outData;
            return ok;

            return true;
        }

        std::string outPerform_title_;
        std::string outAdjoint_title_;
        std::string out_title_;
        std::string out_ylabel_;
        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output out_;
    };

}

#endif