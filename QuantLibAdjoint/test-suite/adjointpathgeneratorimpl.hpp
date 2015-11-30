/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2005 StatPro Italia srl
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

#ifndef cl_adjoint_path_generator_impl_hpp
#define cl_adjoint_path_generator_impl_hpp
#pragma once

#include <boost/test/unit_test.hpp>
#include "utilities.hpp"
#include "adjointpathgeneratortest.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include <ql/methods/montecarlo/mctraits.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/processes/ornsteinuhlenbeckprocess.hpp>
#include <ql/processes/squarerootprocess.hpp>
#include <ql/processes/stochasticprocessarray.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/utilities/dataformatters.hpp>
#include <ql/stochasticprocess.hpp>
#include <ql/processes/eulerdiscretization.hpp>


using namespace QuantLib;
using namespace boost::unit_test_framework;

#define OUTPUT_FOLDER_NAME "AdjointPathGenerator"

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
        pointNo = 1,
        // Number of points for performance plot.
        iterNo = 1,
        // Step for portfolio size for performance testing .
        step = 1,
#endif
        // Defines performance accuracy. Its value is a minimum number
        // of calling of O(1) complexity methods per one performance test.
        iterNumFactor = 10,
    };

    inline Matrix generateCorrelation(Size size)
    {
        Matrix M(size, size);
        for (Size i = 0; i < size; i++)
        {
            M[i][i] = 1;
            for (Size j = i + 1; j < size; j++)
            {
                M[i][j] = 1 - j*0.6 / size;
                M[j][i] = M[i][j];
            }
        }
        return M;
    }


    struct BlackScholes
    {
        BlackScholes()
        {
        }

        std::string getName() { return "BlackScholes"; }

        static boost::shared_ptr<StochasticProcess> createProcess(std::vector<Real> sigmaValue)
        {
            std::vector<boost::shared_ptr<StochasticProcess1D> > processes;
            Size size = sigmaValue.size();

            Handle<Quote> x0(boost::shared_ptr<Quote>(new SimpleQuote(100.0)));
            Handle<YieldTermStructure> r(flatRate(0.05, Actual360()));
            Handle<YieldTermStructure> q(flatRate(0.02, Actual360()));

            Matrix correlation = generateCorrelation(size);

            for (Size i = 0; i < size; i++)
            {
                Handle<BlackVolTermStructure> sigma(flatVol(sigmaValue[i], Actual360()));
                processes.push_back(boost::shared_ptr<StochasticProcess1D>(
                    new BlackScholesMertonProcess(x0, q, r, sigma)));
            }

            return boost::shared_ptr<StochasticProcess>(
                new StochasticProcessArray(processes, correlation));
        }
    };

    struct OrnstUhlenbeck
    {
        OrnstUhlenbeck()
        {
        }

        std::string getName() { return "OrnsteinUhlenbeck"; }

        static boost::shared_ptr<StochasticProcess> createProcess(std::vector<Real> sigmaValue)
        {
            std::vector<boost::shared_ptr<StochasticProcess1D> > processes;
            Size size = sigmaValue.size();

            Matrix correlation = generateCorrelation(size);
            Real speed = 0.1;

            for (Size i = 0; i < size; i++)
            {
                Real sigma = sigmaValue[i];
                processes.push_back(boost::shared_ptr<StochasticProcess1D>(
                    new QuantLib::OrnsteinUhlenbeckProcess(0.1, sigma)));
            }

            return boost::shared_ptr<StochasticProcess>(
                new StochasticProcessArray(processes, correlation));
        }
    };

    struct SqRoot
    {
        SqRoot()
        {
        }

        std::string getName() { return "SquareRoot"; }

        static boost::shared_ptr<StochasticProcess> createProcess(std::vector<Real> sigmaValue)
        {
            std::vector<boost::shared_ptr<StochasticProcess1D>> processes;
            Size size = sigmaValue.size();

            Matrix correlation = generateCorrelation(size);
            Real mean = 0.5;
            Real speed = 0.1;
            Real x0 = 1.0;

            for (Size i = 0; i < size; i++)
            {
                Real sigma = sigmaValue[i];
                processes.push_back(boost::shared_ptr<StochasticProcess1D>(new QuantLib::SquareRootProcess(mean, speed, sigma, 1.0)));
            }

            return boost::shared_ptr<StochasticProcess>(
                new StochasticProcessArray(processes, correlation));
        }
    };

    template <class Process>
    struct PathGenerators : public Process
    {
        PathGenerators()
        {
        };

        // Generates a multipath from a random number generator.
        // Input parameter:
        //  - sigmaValue - sigma for stochastic process for generating random paths.
        static std::vector<Real> testMultiple(std::vector<Real> sigmaValue)
        {
            typedef PseudoRandom::rsg_type rsg_type;
            typedef MultiPathGenerator<rsg_type>::sample_type sample_type;

            boost::shared_ptr<StochasticProcess> process = createProcess(sigmaValue);

            BigNatural seed = 42;
            Time length = 10;
            Size timeSteps = 12;
            Size assets = process->size();
            std::vector<Real> calculated(2 * assets);
            rsg_type rsg = PseudoRandom::make_sequence_generator(timeSteps*assets,
                                                                 seed);
            MultiPathGenerator<rsg_type> generator(process,
                                                   TimeGrid(length, timeSteps),
                                                   rsg, false);
            Size i, j;
            for (i = 0; i < 100; i++)
                generator.next();

            sample_type sample = generator.next();

            Real error, tolerance = 2.0e-7;
            for (j = 0; j < assets; j++)
                calculated[j] = sample.value[j].back();

            sample = generator.antithetic();
            for (j = 0; j < assets; j++)
                calculated[assets + j] = sample.value[j].back();

            return calculated;
        }

        // Use forward finite difference for approximation of the derivatives.
        // dy/dx = (y(x+h) - y(x))/h
        // Input parameters:
        //  - sigma - vector of input sigma values;
        //  - sf_Finite - vector of calculated derivatives;
        //  - h - step size for finite difference method.
        static void calculateFiniteDifference(std::vector<Real>& sigma, std::vector<Real>& sf_Finite, double h)
        {
            Size numOfIndep = sigma.size();
            Size n = 2 * numOfIndep;
            sf_Finite.resize(n*numOfIndep, 0.0);

            std::vector<Real> value = testMultiple(sigma);
            std::vector<Real> rightValue;

            for (Size i = 0; i < numOfIndep; i++)
            {
                sigma[i] += h;
                rightValue = testMultiple(sigma);
                for (Size j = 0; j < 2; j++)
                    sf_Finite[(i + numOfIndep * j)*numOfIndep + i] = (rightValue[i + j*numOfIndep] - value[i + j*numOfIndep]) / h;
                sigma[i] -= h;
            }
        }
    };

    // Struct for plots recording.
    struct SigmaSens
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "Sigma", ""
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, SigmaSens& v)
        {
                stm << v.sigma_
                    << ";" << v.sensitivity_
                    << std::endl;
                return stm;
            }

        Real sigma_;
        Real sensitivity_;
    };

    template <class Process>
    struct TestData
        : public PathGenerators<Process>
    {
        struct Test
        : public cl::AdjointTest<Test>
        {
            Test(Size size, TestData<Process>* data)
            : size_(size)
            , data_(data)
            , sigma_(size)
            , assets_()
            {
                setLogger(&data_->outPerform_);

                if (size_ != 1)
                {
                    for (Size i = 0; i < size_; i++)
                    {
                        sigma_[i] = 0.20 + i*0.001;
                    }
                }
                else
                {
                    sigma_[0] = 0.20;
                }
            }

            Size indepVarNumber() { return size_; }

            Size depVarNumber() { return 2*size_; }

            Size minPerfIteration() { return iterNumFactor; }

            void recordTape()
            {
                cl::Independent(sigma_);
                generateAssets();
                f_ = std::make_unique<cl::tape_function<double>>(sigma_, assets_);
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
            void generateAssets()
            {
                assets_ = data_->testMultiple(sigma_);
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                double h = 1.0e-10;  // shift for finite diff. method
                analyticalResults_.resize(size_);
                data_->calculateFiniteDifference(sigma_, analyticalResults_, h);
            }

            double relativeTol() const { return 1e-2; }

            double absTol() const { return 1e-8; }

            Size size_;
            TestData<Process>* data_;
            std::vector<cl::tape_double> sigma_;
            std::vector<cl::tape_double> assets_;
        };

        TestData()
            : PathGenerators<Process>()
        , outPerform_(OUTPUT_FOLDER_NAME"//" + getName() + "//"
        , { { "filename", "AdjointPerformance" }
        , { "not_clear", "Not" }
        , { "title", "Last-time step asset value differentiation performance with respect to sigma" }
        , { "ylabel", "Time (s)" }
        , { "xlabel", "Number of sigmas" }
        , { "line_box_width", "-5" }
        , { "smooth", "default" }
        , { "cleanlog", "true" } })
            , outAdjoint_(OUTPUT_FOLDER_NAME"//" +  getName() + "//"
        , { { "filename", "Adjoint" }
        , { "not_clear", "Not" }
        , { "title", "Last-time step asset value adjoint differentiation performance with respect to sigma" }
        , { "ylabel", "Time (s)" }
        , { "xlabel", "Number of sigmas" }
        , { "smooth", "default" }
        , { "cleanlog", "false" } })
            , outSize_(OUTPUT_FOLDER_NAME"//" +  getName() + "//"
        , { { "filename", "TapeSize" }
        , { "not_clear", "Not" }
        , { "title", "Tape size dependence on number of sigma" }
        , { "ylabel", "Memory (MB)" }
        , { "xlabel", "Number of sigmas" }
        , { "smooth", "default" }
        , { "cleanlog", "false" } })
            , out_(OUTPUT_FOLDER_NAME"//" +  getName() + "//output"
        , { { "filename", "SigmaDependence" }
        , { "not_clear", "Not" }
        , { "title", "Last-time step asset value sensitivity dependence on sigma" }
        , { "ylabel", "Last-time step asset value sensitivity" }
        , { "xlabel", "Sigma" }
        , { "smooth", "default" }
        , { "cleanlog", "false" } })
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
            std::vector<SigmaSens> outData(pointNo);
            auto test = getTest(pointNo);
            for (Size i = 0; i < pointNo; i++)
            {
                test->sigma_[0] += 0.05;
                test->generateAssets();
                outData[i] = { test->sigma_[0], test->assets_[i] };
            }
            out_ << outData;
            return true;
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output out_;

    };

}

#endif
