/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2006, 2008 Ferdinando Ametrano
Copyright (C) 2006 François du Vignaud
Copyright (C) 2007 Cristina Duminuco
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

// Based on swaptionvolatilitymatrix.cpp project from Quantlib.

#ifndef cl_adjoint_swaption_matrix_impl_hpp
#define cl_adjoint_swaption_matrix_impl_hpp
#pragma once

#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include "utilities.hpp"
#include <ql/indexes/swap/euriborswap.hpp>
#include <ql/instruments/makeswaption.hpp>
#include <ql/pricingengines/swaption/blackswaptionengine.hpp>
#include <string>
#include <boost/make_shared.hpp>
#include <ql/quantlib.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;


#define OUTPUT_FOLDER_NAME "AdjointSwaptionVolatilityMatrix"


namespace
{
    enum
    {
#if defined CL_GRAPH_GEN
        // Number of points for dependency plots.
        pointNo = 50,
        // Number of points for performance plot.
        iterNo = 50,
        // Step for portfolio size for performance testing .
        step = 1,
#else
        // Number of points for dependency plots.
        pointNo = 5,
        // Number of points for performance plot.
        iterNo = 5,
        // Step for portfolio size for performance testing .
        step = 1,
#endif
        // Defines performance accuracy. Its value is a minimum number
        // of calling of O(1) complexity methods per one performance test.
        iterNumFactor = 0,
        maxOptionSize = 100,
        maxSwapSize = 2,
    };

    enum DataType { floatDateFloatData, floatDateFixedData, fixedDateFloatData, fixedDateFixedData };

    std::string toString(DataType description)
    {
        switch (description)
        {
            case floatDateFloatData:
                return "FloatDateFloatData";
            case floatDateFixedData:
                return "FloatDateFixedData";
            case fixedDateFloatData:
                return "FixedDateFloatData";
            case fixedDateFixedData:
                return "FixedDateFixedData";
            default:
                return "";
        }
    }

    struct ExpVolVariation
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "Expected Vol", ""
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, ExpVolVariation& v)
        {
                stm << v.expVol_
                    << ";" << v.forwardValue_
                    << std::endl;
                return stm;
            }

        Real expVol_;
        Real forwardValue_;
    };

    struct SwaptionTenor
    {
        std::vector<Period> options_;
        std::vector<Period> swaps_;
    };

    struct SwaptionMarketConvention
    {
        SwaptionMarketConvention():
          calendar_(TARGET())
        , optionBdc_ (ModifiedFollowing)
        , dayCounter_(Actual365Fixed())
        {
        }
        Calendar calendar_;
        BusinessDayConvention optionBdc_;
        DayCounter dayCounter_;
    };

    struct AtmVolatilityData
    {
        AtmVolatilityData()
        {
            tenors_.options_.resize(maxOptionSize);
            for (Size i = 0; i < maxOptionSize; i++)
                tenors_.options_[i] = Period(5 * i + 1, Months);

            tenors_.swaps_.resize(maxSwapSize);
            for (Size i = 0; i < maxSwapSize; i++)
                tenors_.swaps_[i] = Period(5 * i + 1, Years);

            vols_ = Matrix(maxOptionSize, maxSwapSize);
            for (Size i = 0; i < maxOptionSize; i++)
            for (Size j = 0; j < maxSwapSize; j++)
                vols_[i][j] = 0.14 + i*0.002 - j*0.001;

            volsHandle_.resize(tenors_.options_.size());
            for (Size i = 0; i < maxOptionSize; i++)
            {
                volsHandle_[i].resize(maxSwapSize);
                for (Size j = 0; j < maxSwapSize; j++)
                    volsHandle_[i][j] = Handle<Quote>(boost::make_shared<SimpleQuote>(vols_[i][j]));
            }
        };
        SwaptionTenor tenors_;
        Matrix vols_;
        std::vector<std::vector<Handle<Quote> > > volsHandle_;
    };

    struct CommonVars
    {
        explicit CommonVars(DataType description) :
          conventions_()
        , atm_()
        , termStructure_(boost::make_shared<FlatForward>(0, conventions_.calendar_, 0.05, Actual365Fixed()))
        {
            Settings::instance().evaluationDate() =
                conventions_.calendar_.adjust(Date::todaysDate());
            createVolatilityMatrix(description, vol_);
        }

        void createVolatilityMatrix(DataType description, boost::shared_ptr<SwaptionVolatilityMatrix>& vol_)
        {
            switch (description)
            {
                case floatDateFloatData:
                    vol_ = boost::make_shared<SwaptionVolatilityMatrix>(conventions_.calendar_,
                                                                      conventions_.optionBdc_,
                                                                      atm_.tenors_.options_,
                                                                      atm_.tenors_.swaps_,
                                                                      atm_.volsHandle_,
                                                                      conventions_.dayCounter_);
                    break;
                case floatDateFixedData:
                    vol_ = boost::make_shared<SwaptionVolatilityMatrix>(conventions_.calendar_,
                                                                      conventions_.optionBdc_,
                                                                      atm_.tenors_.options_,
                                                                      atm_.tenors_.swaps_,
                                                                      atm_.volsHandle_,
                                                                      conventions_.dayCounter_);
                    break;
                case fixedDateFloatData:
                    vol_ = boost::make_shared<SwaptionVolatilityMatrix>(Settings::instance().evaluationDate(),
                                                                      conventions_.calendar_,
                                                                      conventions_.optionBdc_,
                                                                      atm_.tenors_.options_,
                                                                      atm_.tenors_.swaps_,
                                                                      atm_.volsHandle_,
                                                                      conventions_.dayCounter_);
                    break;
                case fixedDateFixedData:
                    vol_ = boost::make_shared<SwaptionVolatilityMatrix>(Settings::instance().evaluationDate(),
                                                                      conventions_.calendar_,
                                                                      conventions_.optionBdc_,
                                                                      atm_.tenors_.options_,
                                                                      atm_.tenors_.swaps_,
                                                                      atm_.volsHandle_,
                                                                      conventions_.dayCounter_);
                    break;
                default:
                    BOOST_ERROR("Define type of initial data.");
              }
        }

        // Find implied volatility based on expected volatility.
        // Input parameters:
        // - expVol - vector of expected volatilities;
        // - implVol - output vector of implied volatilities.
        void findImpliedVol(std::vector<Volatility>& expVol, std::vector<Volatility>& implVol)
        {
            Size m = expVol.size();
            Size n = m / maxSwapSize;
            implVol.resize(m);
            boost::shared_ptr<BlackSwaptionEngine> engine =
                                    boost::make_shared<BlackSwaptionEngine>(termStructure_,
                                    Handle<SwaptionVolatilityStructure>(vol_));
            for (Size j = 0; j < maxSwapSize; j++)
            {
                boost::shared_ptr<SwapIndex> swapIndex = boost::make_shared<EuriborSwapIsdaFixA>(Period(j + 1, Years), termStructure_);

                for (Size i = 0; i < n; ++i)
                {
                    // Create swaption.
                    Swaption swaption =
                        MakeSwaption(swapIndex, Period(i + 1, Months))
                        .withPricingEngine(engine);
                    // Calculate NPV.
                    Real npv = swaption.NPV();
                    // Calculate implied volatility and put it to the vector.
                    implVol[j*n + i] = swaption.impliedVolatility(npv, termStructure_,
                        expVol[j*n + i] * 0.98, 10e-6,
                        100, 10.0e-7, 4.0, 0.0);
                }
            }
        }

        // Use central finite difference for approximation of the derivatives:
        // dy/dx = (y(x+h) - y(x-h))/(2*h)
        // Input parameters:
        //  - expVol - vector of expected volatilities;
        //  - sf_Analytical - vector for calculated derivatives.
        //  - h - step size for finite difference method;
        void calculateCentralFinDiff(std::vector<Volatility>& expVol, std::vector<Real>& sf_Analytical, double h)
        {
            Size n = expVol.size();
            sf_Analytical.resize(n*n);
            std::vector<Volatility> implVolRight(n, 0.0), implVolLeft(n, 0.0);
            for (Size i = 0; i < n; i++)
            {
                // Find implied volatility with shifting expected volatility with a step size of +h.
                expVol[i] += h;
                findImpliedVol(expVol, implVolRight);
                // Find implied volatility with shifting expected volatility with a step size of -h.
                expVol[i] -= 2 * h;
                findImpliedVol(expVol, implVolLeft);
                //Evaluate derivatives using central finite difference
                sf_Analytical[i*n + i] = (implVolRight[i] - implVolLeft[i]) / (2 * h);
                expVol[i] += h;
            }
        }

        SwaptionMarketConvention conventions_;
        AtmVolatilityData atm_;
        RelinkableHandle<YieldTermStructure> termStructure_;
        SavedSettings backup_;
        boost::shared_ptr<SwaptionVolatilityMatrix> vol_;
    };

    struct TestData
        : public CommonVars
    {
        struct Test
        : public cl::AdjointTest<Test>
        {
            Test(Size optionSize, TestData* data)
            : optionSize_(optionSize)
            , data_(data)
            , expVol_(optionSize_*maxSwapSize)
            , implVol_(optionSize_*maxSwapSize)
            {
                setLogger(&data_->outPerform_);

                for (Size i = 0; i < maxSwapSize; i++)
                {
                    for (Size j = 0; j < optionSize_; j++)
                    {
                        expVol_[i*optionSize_ + j] = 0.14 + i*0.002 - j*0.001;
                    }
                }
            }

            Size indepVarNumber() { return optionSize_*maxSwapSize; }

            Size depVarNumber() { return optionSize_*maxSwapSize; }

            Size minPerfIteration() { return iterNumFactor; }

            void recordTape()
            {
                cl::Independent(expVol_);
                calculateImplVol();
                f_ = std::make_unique<cl::tape_function<double>>(expVol_, implVol_);
            }

            void calcReverse()
            {
                using namespace std;
                size_t n = indepVarNumber();
                size_t m = depVarNumber();
                reverseResults_.resize(m*n);
                vector<double> dw(m), dy(n);

                for (size_t i = 0; i < m; i++)
                {
                    dw[i] = 1;
                    dy = f_->Reverse(1, dw);
                    dw[i] = 0;
                    for (size_t j = 0; j < n; j++)
                        reverseResults_[i*n + j] = dy[j];
                }
            }

            // Calculates price of portfolio and each option.
            void calculateImplVol()
            {
                implVol_.resize(depVarNumber(), 0);
                data_->findImpliedVol(expVol_, implVol_);
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                double h = 1e-3;  // shift for finite diff. method
                analyticalResults_.resize(depVarNumber()*indepVarNumber());
                data_->calculateCentralFinDiff(expVol_, analyticalResults_, h);
            }

            double relativeTol() const { return 1e-5; }

            double absTol() const { return 1e-8; }

            TestData* data_;
            Size optionSize_;
            std::vector<cl::tape_double> expVol_;
            std::vector<cl::tape_double> implVol_;
        };

        explicit TestData(DataType description)
            : CommonVars(description)
            , outPerform_(OUTPUT_FOLDER_NAME "//"+ toString(description)
            , { { "filename", "AdjointPerformance" }
            , { "not_clear", "Not" }
            , { "title", "Implied volatility differentiation performance with respect to expected volatility"}
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of expected volatilities" }
            , { "line_box_width", "-5" }
            , { "smooth", "default" }
            , { "cleanlog", "true" }
        })
            , outAdjoint_(OUTPUT_FOLDER_NAME "//" + toString(description)
            , { { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "title", "Implied volatility adjoint performance with respect to expected volatility"}
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of expected volatilities" }
            , { "smooth", "default" }
            , { "cleanlog", "false" }
        })
            , outSize_(OUTPUT_FOLDER_NAME "//" + toString(description)
            , { { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "title", "Tape size dependence on  number of expected volatilities" }
            , { "ylabel", "Memory (MB)" }
            , { "xlabel", "Number of expected volatilities" }
            , { "smooth", "default" }
            , { "cleanlog", "false" }
        })
            , out_(OUTPUT_FOLDER_NAME "//" + toString(description) + "//output",
            { { "filename", "ImpliedVolOnExpectedVol" }
            , { "not_clear", "Not" }
            , { "title", "Implied volatility dependence on expected volatility" }
            , { "ylabel", "Implied volatility" }
            , { "xlabel", "Expected volatility" }
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

        std::shared_ptr<Test> getTest(size_t size)
        {
            return std::make_shared<Test>(size, this);
        }

        // Makes plots for strike sensitivity dependence.
        bool recordDependencePlot()
        {
            std::vector<ExpVolVariation> outData(pointNo);
            auto test = getTest(pointNo);
            test->doTape();
            for (Size i = 0; i < pointNo; i++)
            {
                outData[i] = { test->expVol_[i], test->implVol_[i] };
            }
            out_ << outData;
            return true;
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output out_;
    };

    typedef TestData::Test SwaptionMatrixTest;
}
#endif