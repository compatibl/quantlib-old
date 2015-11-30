/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2005, 2007, 2009, 2014 Klaus Spanderen
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

#ifndef cl_adjoint_heston_process_impl_hpp
#define cl_adjoint_heston_process_impl_hpp
#pragma once


#include "adjointhestonprocesstest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include <ql/quotes/simplequote.hpp>
#include <ql/math/modifiedbessel.hpp>
#include <ql/processes/hestonprocess.hpp>
#include <ql/time/daycounters/actualactual.hpp>
#include <ql/time/schedule.hpp>
#include <boost/make_shared.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;


#define OUTPUT_FOLDER_NAME "AdjointHestonProcess"

/*!
AdjointHestonProcessTest is tested using adjoint differentiation for the class hestonprocess.

Heston model describes the evolution of the volatility of an underlying asset.
It is a stochastic volatility model: such a model assumes that the
volatility of the asset is neither constant, nor even deterministic,
but is a random process.

This class describes the square root stochastic volatility
process governed by
\f[
\begin{array}{rcl}
dS(t, S)  &=& \mu S dt + \sqrt{v} S dW_1 \\
dv(t, S)  &=& \kappa (\theta - v) dt + \sigma \sqrt{v} dW_2 \\
dW_1 dW_2 &=& \rho dt
\end{array}
\f]

Tested functions are expectation, evolve, phi, drift, stdDeviation, covariance.

Phi is the continuous version of a characteristic function
for the exact sampling of the Heston process, s. page 8, formula 13,
M. Broadie, O. Kaya, Exact Simulation of Stochastic Volatility and
other Affine Jump Diffusion Processes
http://finmath.stanford.edu/seminars/documents/Broadie.pdf
Output of function Phi is value of haracteristic function.
Inputs are process, value of variable, start and finish variance distibution, and delta time.
Phi is obtained from fourier transform, so Phi has complex derivative with respect to input variable.

Output of fuction expectation is Expectation
\f$ E(\mathrm{x}_{t_0 + \Delta t}
| \mathrm{x}_{t_0} = \mathrm{x}_0) \f$
of the process after a time interval \f$ \Delta t \f$
according to the given discretization.
Inputs are time, start value, and delta time.
Expectation is differentiated with respect to delta time.

Output of fuction evolve is Evolve,
which returns the asset value after a time interval \f$ \Delta t
\f$ according to the given discretization. By default, it
returns
\f[
E(x_0,t_0,\Delta t) + S(x_0,t_0,\Delta t) \cdot \Delta w
\f]
where \f$ E \f$ is the expectation and \f$ S \f$ the
standard deviation.
Inputs are time, start value, and delta time.
Evolve is differentiated with respect to delta time.

Output of fuction drift is Drift (the change of the average value)
Inputs are time, array of average values.
Drift is differentiated with respect to average values.

Output of fuction stdDeviation is Standart Deviation.
Inputs are time, array of average values, and delta time.
StdDeviation is differentiated with respect to delta time.

Output of fuction covariance is covariances matrix
Inputs are time, array of average values, and delta time.
Covariance is differentiated with respect to delta time.
*/


namespace
{
    enum
    {
#if defined CL_GRAPH_GEN
        // Number of points for dependency plots.
        pointNo = 100,
        // Number of points for performance plot.
        iterNo = 30,
        // Step for portfolio size for performance testing .
        step = 1,
#else
        // Number of points for dependency plots.
        pointNo = 0,
        // Number of points for performance plot.
        iterNo = 0,
        // Step for portfolio size for performance testing .
        step = 1,
#endif
        // Defines performance accuracy. Its value is a minimum number
        // of calling of O(1) complexity methods per one performance test.
        iterNumFactor = 1000,
    };

    // copy of function Phi from hestonprocess.hpp
    // This is the continuous version of a characteristic function
    // for the exact sampling of the Heston process, s. page 8, formula 13,
    // M. Broadie, O. Kaya, Exact Simulation of Stochastic Volatility and
    // other Affine Jump Diffusion Processes
    // http://finmath.stanford.edu/seminars/documents/Broadie.pdf
    std::complex<Real> Phi1(const HestonProcess& process,
        const std::complex<Real>& a,
        Real nu_0, Real nu_t, Time dt) {
        const Real theta = process.theta();
        const Real kappa = process.kappa();
        const Real sigma = process.sigma();

        const Volatility sigma2 = sigma*sigma;
        const std::complex<Real> ga = std::sqrt(
            kappa*kappa - 2 * sigma2*a*std::complex<Real>(0.0, 1.0));
        const Real d = 4 * theta*kappa / sigma2;

        const Real nu = 0.5*d - 1;
        const std::complex<Real> z
            = ga*std::exp(-0.5*ga*dt) / (1.0 - std::exp(-ga*dt));
        const std::complex<Real> log_z
            = -0.5*ga*dt + std::log(ga / (1.0 - std::exp(-ga*dt)));

        const std::complex<Real> alpha
            = 4.0*ga*std::exp(-0.5*ga*dt) / (sigma2*(1.0 - std::exp(-ga*dt)));

        const std::complex<Real> beta = 4.0*kappa*std::exp(-0.5*kappa*dt)
            / (sigma2*(1.0 - std::exp(-kappa*dt)));

        return ga*std::exp(-0.5*(ga - kappa)*dt)*(1 - std::exp(-kappa*dt))
            / (kappa*(1.0 - std::exp(-ga*dt)))
            *std::exp((nu_0 + nu_t) / sigma2 * (
            kappa*(1.0 + std::exp(-kappa*dt)) / (1.0 - std::exp(-kappa*dt))
            - ga*(1.0 + std::exp(-ga*dt)) / (1.0 - std::exp(-ga*dt))))
            *std::exp(nu*log_z) / std::pow(z, nu)
            *((nu_t > 1e-8)
            ? modifiedBesselFunction_i(
            nu, std::sqrt(nu_0*nu_t)*alpha)
            / modifiedBesselFunction_i(
            nu, std::sqrt(nu_0*nu_t)*beta)
            : std::pow(alpha / beta, nu)
            );
    }

    struct CommonVars0
    {
        CommonVars0()
        : dayCounter_()
        {
            dayCounter_ = ActualActual();
        }
        DayCounter dayCounter_;

    };

    struct CommonVars : public CommonVars0
    {
        explicit CommonVars()
        : CommonVars0()
        , s0_(boost::make_shared<SimpleQuote>(SimpleQuote(1.05)))
        , process_(Handle<YieldTermStructure>(flatRate(0.04, dayCounter_)), Handle<YieldTermStructure>(flatRate(0.5, dayCounter_)), s0_, 0.3, 1.16, 0.2, 0.8, 0.8,
        HestonProcess::QuadraticExponentialMartingale)
        {
        }

        Handle<Quote> s0_;
        HestonProcess process_;
    };

    enum VariableName { expectation, evolve, phi, drift, stdDeviation, covariance };

    std::string toString(VariableName name)
    {
        switch (name)
        {
        case expectation:
            return "Expectation";
        case evolve:
            return "Evolve";
        case phi:
            return "Phi";
        case drift:
            return "Drift";
        case stdDeviation:
            return "stdDeviation";
        case covariance:
            return "covariance";
        default:
            return "";
        }
    }

    struct HestonProcessTestData : public CommonVars
    {
        explicit HestonProcessTestData(VariableName type)
        : CommonVars()
        , type_(type)
        , t(0.5)
        , x(2, 0.3, 0.2)
        , dw(2, 0.1, 0.4)
        {}

        // Calculate the an expectation in case of test of function Expectation
        // and an evolve in case of test of function Evolve.
        Real calculateRealFromDeltaTime(Real deltaTimes)
        {
            if (type_ == expectation)
            {
                return process_.expectation(t, x, deltaTimes)[0];
            }

            if (type_ == evolve)
            {
                return process_.evolve(t, x, deltaTimes, dw)[0];
            }

            return Real(0.0);
        }

        Matrix calculateMatrixFromDeltaTime(Real deltaTimes)
        {
            if (type_ == stdDeviation)
            {
                return process_.stdDeviation(t, x, deltaTimes);
            }

            if (type_ == covariance)
            {
                return process_.covariance(t, x, deltaTimes);
            }

            return Matrix();
        }

        VariableName type_;
        Real t;
        Array x;
        Array dw;
    };

    // Test of functions expectation and evolve.
    // Testing derivative with respect to delta time.
    struct HestonProcessRespectToDeltaTimeTest
        : public cl::AdjointTest<HestonProcessRespectToDeltaTimeTest>
    {
        explicit HestonProcessRespectToDeltaTimeTest(HestonProcessTestData* data, Size size, cl::tape_empty_test_output* logger = nullptr, double shift = 1e-10)
        : AdjointTest()
        , data_(data)
        , size_(size)
        , shift_(shift)
        , deltaTimes_(size)
        , Variables_(size)
        , SumVariables_(1)
        {
            setLogger(logger);
            for (int i = 0; i < size_; i++)
            {
                deltaTimes_[i] = (i + 1) * 0.2 / double(size);
            }
        }

        Size indepVarNumber() { return size_; }

        Size depVarNumber() { return 1; }

        Size minPerfIteration() { return iterNumFactor; }

        double absTol() const { return 1e-5; }

        void calcAnalytical()
        {
            analyticalResults_.resize(size_);
            for (int i = 0; i < size_; i++)
            {
                double deltaTimeUp = double(deltaTimes_[i]) + shift_;
                double deltaTimeDown = double(deltaTimes_[i]) - shift_;
                double varUp = double(data_->calculateRealFromDeltaTime(deltaTimeUp));
                double varDown = double(data_->calculateRealFromDeltaTime(deltaTimeDown));
                analyticalResults_[i] = (varUp - varDown) / 2 / shift_;
            }
        }

        void calculateVariables()
        {
            SumVariables_[0] = 0.0;
            for (int i = 0; i < size_; i++)
            {
                Variables_[i] = data_->calculateRealFromDeltaTime(deltaTimes_[i]);
                SumVariables_[0] += Variables_[i];
            }
        }

        void recordTape()
        {
            Independent(deltaTimes_);
            calculateVariables();
            f_ = std::make_unique<cl::tape_function<double>>(deltaTimes_, SumVariables_);
        }

        HestonProcessTestData* data_;
        Size size_;
        double shift_;
        std::vector<Real> deltaTimes_;
        std::vector<Real> Variables_;
        std::vector<Real> SumVariables_;
    };

    // Test of matrix functions stdDeviation and covariance.
    struct HestonProcessMatrixFunctionsTest :
        public cl::AdjointTest<HestonProcessMatrixFunctionsTest>
    {
        static const CalcMethod default_method = other;

        explicit HestonProcessMatrixFunctionsTest(HestonProcessTestData* data, cl::tape_empty_test_output* logger = nullptr, double shift = 1e-10)
        : AdjointTest()
        , data_(data)
        , shift_(shift)
        , deltaTimes_(1, 0.1)
        , size_(4)
        {}

        Size indepVarNumber() { return 1; }

        Size depVarNumber() { return size_; }

        Size minPerfIteration() { return iterNumFactor; }

        double absTol() const { return 1e-6; }

        void calculateMatrixRes()
        {
            Matrix StdDeviation = data_->calculateMatrixFromDeltaTime(deltaTimes_[0]);
            Results_.assign(StdDeviation.begin(), StdDeviation.end());
        }

        void recordTape()
        {
            Independent(deltaTimes_);
            calculateMatrixRes();
            f_ = std::make_unique<cl::tape_function<double>>(deltaTimes_, Results_);
        }

        void calcAdjoint()
        {
            int num = Results_.size();
            std::vector<double> w(1, 1);
            adjointResults_ = f_->Forward(1, w);
        }

        void calcAnalitical()
        {
            Time dtUp = deltaTimes_[0] + shift_;
            Time dtDown = deltaTimes_[0] - shift_;

            Matrix ResUp = data_->calculateMatrixFromDeltaTime(dtUp);
            Matrix ResDown = data_->calculateMatrixFromDeltaTime(dtDown);

            for (int i = 0; i < ResUp.columns(); i++)
            {
                for (int j = 0; j < ResUp.rows(); j++)
                {
                    double analyticalRes = double((ResUp[i][j] - ResDown[i][j]) / 2 / shift_);
                    analyticalResults_.push_back(analyticalRes);
                }
            }
        }

        HestonProcessTestData* data_;
        double shift_;
        std::vector<Real> deltaTimes_;
        std::vector<Real> Results_;
        Size size_;
    };

    struct HestonProcessDriftTest : public cl::AdjointTest<HestonProcessDriftTest>
    {
        static const CalcMethod default_method = other;

        explicit HestonProcessDriftTest(CommonVars* data, cl::tape_empty_test_output* logger = nullptr, double shift = 1e-10)
            : AdjointTest()
            , data_(data)
            , shift_(shift)
            , size_(2)
            , values_(size_)
            , x_(size_)
            , drift_(2)
            , t_(0.5)
            , indepDouble_(2)
        {
            for (int i = 0; i < size_; i++)
            {
                indepDouble_[i] = 0.3 + 0.2 * i;
                values_[i] = indepDouble_[i];
            }
        }

        Size indepVarNumber() { return size_; }

        Size depVarNumber() { return size_; }

        Size minPerfIteration() { return iterNumFactor; }

        double absTol() const { return 1e-6; }

        void calculateDrift()
        {
            for (int i = 0; i < size_; i++)
            {
                x_[i] = values_[i];
            }
            Array drift = data_->process_.drift(t_, x_);
            for (int i = 0; i < size_; i++)
            {
                drift_[i] = drift[i];
            }
        }

        void recordTape()
        {
            cl::Independent(values_);
            calculateDrift();
            f_ = std::make_unique<cl::tape_function<double>>(values_, drift_);
        }

        void calcAdjoint()
        {
            adjointResults_ = f_->Jacobian(indepDouble_);
        }

        void calcAnalitical()
        {
            analyticalResults_.resize(size_ * size_);
            Array xUp = x_;
            Array xDown = x_;
            for (int i = 0; i < size_; i++)
            {
                xUp[i] += shift_;
                Array DriftUp = data_->process_.drift(t_, xUp);
                xUp[i] -= shift_;

                xDown[i] -= shift_;
                Array DriftDown = data_->process_.drift(t_, xDown);
                xDown[i] += shift_;

                for (int j = 0; j < size_; j++)
                {
                    analyticalResults_[i + size_ * j] = (DriftUp[j] - DriftDown[j]) / 2 / shift_;
                }
            }
        }

        CommonVars* data_;
        double shift_;
        Size size_;
        std::vector<Real> values_;
        Array x_;
        std::vector<Real> drift_;
        Time t_;
        std::vector<double> indepDouble_;
    };


    // Struct for plots recording.
    struct DeltaTimeResult
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "Delta time", ""
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, DeltaTimeResult& v)
        {
                stm << v.deltaTime_
                    << ";" << v.result_
                    << std::endl;
                return stm;
            }

        Real deltaTime_;
        Real result_;
    };


    struct TestDataRespectDeltaTime
        : public HestonProcessTestData
    {
        explicit TestDataRespectDeltaTime(VariableName name)
        : HestonProcessTestData(name)
        , outPerform_title_("Swaption " + toString(name) + " differentiation performance with respect to delta time")
        , outAdjoint_title_("Swaption " + toString(name) + " adjoint differentiation with respect to fixed delta time")
        , out_title_(toString(name) + " dependence on delta time")
        , out_ylabel_(toString(name))
        , outPerform_(OUTPUT_FOLDER_NAME "//" + toString(name), {
            { "title", outPerform_title_ }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of delta times" }
            , { "line_box_width", "-5" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        })
            , outAdjoint_(OUTPUT_FOLDER_NAME "//" + toString(name), {
                { "title", outAdjoint_title_ }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of delta times" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        })
            , outSize_(OUTPUT_FOLDER_NAME "//" + toString(name), {
                { "title", "Tape size dependence on number of delta time" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
            , { "cleanlog", "false" }
        })

            , out_(OUTPUT_FOLDER_NAME "//" + toString(name) + "//output", {
                { "filename", "DeltaTimeDependence" }
            , { "not_clear", "Not" }
            , { "title", out_title_ }
            , { "ylabel", out_ylabel_ }
            , { "xlabel", "delta time" }
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



        std::shared_ptr<HestonProcessRespectToDeltaTimeTest> getTest(size_t size)
        {
            return std::make_shared<HestonProcessRespectToDeltaTimeTest>(this, size, &outPerform_);
        }

        // Makes plots for derivative of blackVariance dependence.
        bool recordDependencePlot()
        {
            std::vector<DeltaTimeResult> outData(pointNo);
            auto test = getTest(pointNo);
            bool ok = test->testReverse();
            for (Size i = 0; i < pointNo; i++)
            {
                outData[i] = { test->deltaTimes_[i], test->Variables_[i] };
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