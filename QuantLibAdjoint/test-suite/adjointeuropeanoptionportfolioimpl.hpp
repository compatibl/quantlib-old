/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
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

#ifndef cl_adjoint_european_option_portfolio_impl_hpp
#define cl_adjoint_european_option_portfolio_impl_hpp
#pragma once


#include <ql/quantlib.hpp>
#include <iostream>

#include "adjointeuropeanoptionportfoliotest.hpp"
#include "adjointtestbase.hpp"
#include "utilities.hpp"

using namespace QuantLib;
using namespace boost::unit_test_framework;


#define OUTPUT_FOLDER_NAME "AdjointEuropeanOptionPortfolio"

namespace
{
    // Test settings.
    // Finite difference numerical method is used to check calculation of some derivatives.
    // If the following settings are changed the test might fail
    enum
    {
#if defined CL_GRAPH_GEN
#if defined NDEBUG
        // Number of points for performance plot.
        iterNo = 100,
        // Step for portfolio size for performance testing.
        step = 1,
        // Defines performance accuracy. Its value is a minimum number
        // of calling of O(1) complexity methods per one performance test.
        // With iterNumFactor = 100000 one AdjointTestBase::recordPerformance()
        // call runs about 1 second in release mode.
        iterNumFactor = 10000,
        // Number of points for dependency plots.
        pointNo = 1000,
#else
        iterNo = 20,
        step = 1,
        iterNumFactor = 100,
        pointNo = 100,
#endif
#else
        iterNo = 1,
        step = 1,
        iterNumFactor = 1,
        pointNo = 1,
#endif
        testPortfolioSize = 100
    };
}


// Greek calculation and testing logic.
namespace
{
    // Parameters of single option.
    enum DataType
    {
        stock
        , strike
        , sigma
        , time
        , rate
    };

    std::string toString(DataType d, bool fullName = false, bool xlabel = false)
    {
        static const std::string name[] = {
            "Stock",
            "Strike",
            "Sigma",
            "Time",
            "Rate",
        };
        static const std::string full[] = {
            "underlying stock price",
            "stock strike",
            "volatility",
            "maturity time",
            "interest rate",
        };
        static const std::string label[] = {
            "Stock price",
            "Strike",
            "Volatility",
            "Maturity time",
            "Interest rate"
        };
        if (fullName)
        {
            return full[d];
        }
        if (xlabel)
        {
            return label[d];
        }
        return name[d];
    }

    enum Greek
    {
        Delta
        , Vega
        , Theta
        , Rho
        , Gamma
        , Vomma
        , Vanna
        , Charm
        , Veta
        , Vera
        , Speed
        , Ultima
        , Zomma
        , Color
    };


    template <size_t order = 1, bool is_mixed = false, int sign_factor = 1>
    struct DerivCalculator;

    // Calculates first derivatives.
    template <int signFactor>
    struct DerivCalculator<1, false, signFactor>
    {
        static void calcForward(cl::tape_function<double>& f, std::vector<double>& forwardResults)
        {
            forwardResults.resize(f.Domain());
            // Direction for directional derivative.
            std::vector<double> dX(f.Domain(), 0);
            auto px = dX.begin();
            std::generate(forwardResults.begin(), forwardResults.end(), [&]()
            {
                *px = 1;
                double temp = signFactor * f.Forward(1, dX)[0];
                *px++ = 0;
                return temp;
            });
        }

        static void calcReverse(cl::tape_function<double>& f, std::vector<double>& reverseResults)
        {
            reverseResults = f.Reverse(1, std::vector<double>(1, 1));
            for (double& res : reverseResults)
            {
                res *= signFactor;
            }
        }
    };


    // Calculates second derivatives, assuming that all mixed derivatives are zero.
    template <int signFactor>
    struct DerivCalculator<2, false, signFactor>
    {
        static void calcForward(cl::tape_function<double>& f, std::vector<double>& forwardResults)
        {
            forwardResults.resize(f.Domain());
            // Direction for directional derivative.
            std::vector<double> dX(f.Domain(), 0);
            std::vector<double>::iterator px = dX.begin();
            std::generate(forwardResults.begin(), forwardResults.end(), [&]()
            {
                *px = 1;
                f.Forward(1, dX);
                *px++ = 0;
                return signFactor * 2 * f.Forward(2, dX)[0];
            });
        }

        static void calcReverse(cl::tape_function<double>& f, std::vector<double>& reverseResults)
        {
            const size_t n = f.Domain();
            f.Forward(1, std::vector<double>(n, 1));
            std::vector<double> dw = f.Reverse(2, std::vector<double>(1, 1));
            reverseResults.resize(n);
            for (size_t i = 0; i < n; i++)
            {
                reverseResults[i] = signFactor * dw[2 * i + 1];
            }
        }
    };

    // Calculates third derivatives, assuming that all mixed derivatives are zero.
    template <int signFactor>
    struct DerivCalculator<3, false, signFactor>
    {
        static void calcForward(cl::tape_function<double>& f, std::vector<double>& forwardResults)
        {
            forwardResults.resize(f.Domain());
            // Direction for directional derivative.
            std::vector<double> dX(f.Domain(), 0);
            std::vector<double>::iterator px = dX.begin();
            std::generate(forwardResults.begin(), forwardResults.end(), [&]()
            {
                *px = 1;
                f.Forward(1, dX);
                *px++ = 0;
                f.Forward(2, dX);
                return signFactor * 6 * f.Forward(3, dX)[0];
            });
        }

        static void calcReverse(cl::tape_function<double>& f, std::vector<double>& reverseResults)
        {
            const size_t n = f.Domain();
            f.Forward(1, std::vector<double>(n, 1));
            f.Forward(2, std::vector<double>(n, 0));
            std::vector<double> dw = f.Reverse(3, std::vector<double>(1, 1));
            reverseResults.resize(n);
            for (size_t i = 0; i < n; i++)
            {
                reverseResults[i] = signFactor * 2 * dw[3 * i + 2];
            }
        }
    };

    // Calculates second order mixed derivatives with respect to i-th and (i+n)-th
    // variables from 2*n dimensional domain vector. For correct work all other
    // mixed derivatives should be zero.
    template <int signFactor>
    struct DerivCalculator<2, true, signFactor>
    {
        static void calcForward(cl::tape_function<double>& f, std::vector<double>& forwardResults)
        {
            size_t n = f.Domain() / 2;
            forwardResults.resize(n);
            // Direction for directional derivative.
            std::vector<double> dX(f.Domain(), 0);
            std::vector<double>::iterator px = dX.begin();
            std::vector<double>::iterator pxn = dX.begin() + n;
            std::generate(forwardResults.begin(), forwardResults.end(), [&]()
            {
                // Mixed derivative calculated from second order derivatives like this:
                // z = x + y
                // w = x - y
                // ddF/dzdz = ddF/dxdx + 2 * ddF/dxdy + ddF/dydy
                // ddF/dwdw = ddF/dxdx - 2 * ddF/dxdy + ddF/dydy
                // ddF/dxdy = (ddF/dzdz - ddF/dwdw) / 4
                *px = 1;
                *pxn = 1;
                f.Forward(1, dX);
                *px = 0;
                *pxn = 0;
                double dzdz = 2 * f.Forward(2, dX)[0];
                *px = 1;
                *pxn = -1;
                f.Forward(1, dX);
                *px++ = 0;
                *pxn++ = 0;
                double dwdw = 2 * f.Forward(2, dX)[0];
                return signFactor * (dzdz - dwdw) / 4;
            });
        }

        static void calcReverse(cl::tape_function<double>& f, std::vector<double>& reverseResults)
        {
            const size_t n = f.Domain() / 2;
            std::vector<double> dX(f.Domain(), 0);
            std::fill(dX.begin(), dX.begin() + n, 1);
            f.Forward(1, dX);
            std::vector<double> dw = f.Reverse(2, std::vector<double>(1, 1));
            reverseResults.resize(n);
            for (size_t i = 0; i < n; i++)
            {
                reverseResults[i] = signFactor * dw[2 * (n + i) + 1];
            }
        }
    };

    // Calculates third order mixed derivatives with respect to i-th variable
    // one time and (i+n)-th variable two times. Where domain vector is 2*n
    // dimensional. For correct work all other mixed derivatives should be zero.
    template <int signFactor>
    struct DerivCalculator<3, true, signFactor>
    {
        static void calcForward(cl::tape_function<double>& f, std::vector<double>& forwardResults)
        {
            size_t n = f.Domain() / 2;
            forwardResults.resize(n);
            // Direction for directional derivative.
            std::vector<double> dX(f.Domain(), 0);
            std::vector<double>::iterator px = dX.begin();
            std::vector<double>::iterator pxn = dX.begin() + n;
            std::generate(forwardResults.begin(), forwardResults.end(), [&]()
            {
                // Mixed derivative calculated from third order derivatives like this:
                // z = x + y
                // w = x - y
                // d3F/dz3 = d3F/dx3 + 3 * d3F/dx2dy + 3 * d3F/dxdy2 + d3F/dy3
                // d3F/dw3 = d3F/dx3 - 3 * d3F/dx2dy + 3 * d3F/dxdy2 - d3F/dy3
                // d3F/dx2dy = (d3F/dz3 - d3F/dw3 - 2 * d3F/dy3) / 6
                *px = 1;
                *pxn = 1;
                f.Forward(1, dX);
                *px = 0;
                *pxn = 0;
                f.Forward(2, dX);
                double dz3 = 6 * f.Forward(3, dX)[0];
                *px = -1;
                *pxn = 1;
                f.Forward(1, dX);
                *px = 0;
                *pxn = 0;
                f.Forward(2, dX);
                double dw3 = 6 * f.Forward(3, dX)[0];
                *px = 1;
                f.Forward(1, dX);
                *px = 0;
                f.Forward(2, dX);
                ++px;
                ++pxn;
                double dy3 = 6 * f.Forward(3, dX)[0];
                return signFactor * (dz3 - dw3 - 2 * dy3) / 6;
            });
        }

        static void calcReverse(cl::tape_function<double>& f, std::vector<double>& reverseResults)
        {
            const size_t n = f.Domain() / 2;
            std::vector<double> dX(f.Domain(), 0);
            std::fill(dX.begin() + n, dX.end(), 1);
            f.Forward(1, dX);
            std::fill(dX.begin() + n, dX.end(), 0);
            f.Forward(2, dX);
            std::vector<double> dw = f.Reverse(3, std::vector<double>(1, 1));
            reverseResults.resize(n);
            for (size_t i = 0; i < n; i++)
            {
                reverseResults[i] = signFactor * 2 * dw[3 * i + 2];
            }
        }
    };


    template <Greek>
    struct GreekTraits;

#define REGISTER_GREEK_TRAIT(Gr, ind1, ord, sf)                 \
    template <>                                                 \
    struct GreekTraits<(Gr)>                                    \
    {                                                           \
        static const Greek greek = (Gr);                        \
        static const DataType independent1 = (ind1);            \
        static const bool is_mixed = false;                     \
        static const size_t order = (ord);                      \
        static const int sign_factor = (sf);                    \
        static std::string name() { return (#Gr); }             \
    };

#define REGISTER_MIXED_GREEK_TRAIT(Gr, ind1, ind2, ord, sf)     \
    template <>                                                 \
    struct GreekTraits<(Gr)>                                    \
    {                                                           \
        static const Greek greek = (Gr);                        \
        static const DataType independent1 = (ind1);            \
        static const DataType independent2 = (ind2);            \
        static const bool is_mixed = true;                      \
        static const size_t order = (ord);                      \
        static const int sign_factor = (sf);                    \
        static std::string name() { return (#Gr); }             \
    };

    REGISTER_GREEK_TRAIT(Delta, stock, 1, 1);
    REGISTER_GREEK_TRAIT(Vega, sigma, 1, 1);
    REGISTER_GREEK_TRAIT(Theta, time, 1, -1);
    REGISTER_GREEK_TRAIT(Rho, rate, 1, 1);
    REGISTER_GREEK_TRAIT(Gamma, stock, 2, 1);
    REGISTER_GREEK_TRAIT(Vomma, sigma, 2, 1);
    REGISTER_GREEK_TRAIT(Speed, stock, 3, 1);
    REGISTER_GREEK_TRAIT(Ultima, sigma, 3, 1);
    REGISTER_MIXED_GREEK_TRAIT(Vanna, sigma, stock, 2, 1);
    REGISTER_MIXED_GREEK_TRAIT(Charm, stock, time, 2, -1);
    REGISTER_MIXED_GREEK_TRAIT(Veta, sigma, time, 2, 1);
    REGISTER_MIXED_GREEK_TRAIT(Vera, sigma, rate, 2, 1);
    REGISTER_MIXED_GREEK_TRAIT(Zomma, sigma, stock, 3, 1);
    REGISTER_MIXED_GREEK_TRAIT(Color, time, stock, 3, 1);


    // Returns BlackCalculator for option with parameters determined by
    // optionData vector in format { stock, strike, sigma, time, rate}.
    inline BlackCalculator getBC(const std::vector<Real>& optionData)
    {
        return BlackCalculator(Option::Call,
            optionData[strike],
            optionData[stock] * std::exp(optionData[rate] * optionData[time]),
            optionData[sigma] * std::sqrt(optionData[time]),
            std::exp(-optionData[rate] * optionData[time]));
    }


    // Abstract class that determine method of adjoint result checking.
    class CheckingModel
    {
    public:
        // Relative tolerance and absolute tolerance settings.
        CheckingModel(double relativeTol, double absTol, std::string const& description = "BlackCalculator")
            : description_(description)
            , relativeTol_(relativeTol)
            , absTol_(absTol)
        { }

        virtual ~CheckingModel() {}

        // Should calculate Greek without adjoint.
        virtual Real operator()(const std::vector<Real>& optionData) const = 0;

        std::string description_;
        double relativeTol_;
        double absTol_;
    };

    // Provides Greek value from given function.
    class SimpleCheckingModel
        : public CheckingModel
    {
    public:
        typedef Real(*func_type)(const std::vector<Real>&);
        // Constructor parameters:
        // function - function which calculates Greek.
        SimpleCheckingModel(func_type function,
            double relativeTol = 1e-8, double absTol = 1e-10)
            : CheckingModel(relativeTol, absTol)
            , function_(function)
        { }

        // Returns Greek value.
        inline Real operator()(const std::vector<Real>& optionData) const
        {
            return (*function_)(optionData);
        }

    protected:
        func_type function_;
    };

    // Provides Greek value with finite difference derivative calculation.
    class FdCheckingModel
        : public CheckingModel
    {
    public:
        typedef Real(*func_type)(const std::vector<Real>&);
        // Constructor parameters:
        // function - function which is differentiated,
        // param - function parameter to differentiate,
        // shift - shift for finite difference calculation.
        FdCheckingModel(func_type function, DataType param,
            double shift = 1e-5,
            double relativeTol = 1e-8, double absTol = 1e-10,
            std::string const& description = "Finite diff.")
            : CheckingModel(relativeTol, absTol, description)
            , function_(function)
            , param_(param)
            , shift_(shift)
        { }

        // Returns finite difference derivative of function_ with respect to param_.
        inline Real operator()(const std::vector<Real>& optionData) const
        {
            std::vector<Real> optionData_left(optionData);
            optionData_left[param_] -= shift_;
            std::vector<Real> optionData_right(optionData);
            optionData_right[param_] += shift_;
            return ((*function_)(optionData_right)-(*function_)(optionData_left)) / (2 * shift_);
        }

    protected:
        func_type function_;
        DataType param_;
        double shift_;
    };

    // Provides Greek value with finite difference second order derivative calculation.
    class Fd2CheckingModel
        : public CheckingModel
    {
    public:
        typedef Real(*func_type)(const std::vector<Real>&);
        // Constructor parameters:
        // function - function which is differentiated,
        // param - function parameter to differentiate,
        // shift - shift for finite difference calculation.
        Fd2CheckingModel(func_type function, DataType param,
            double shift = 1e-5,
            double relativeTol = 1e-2, double absTol = 1e-7,
            std::string const& description = "Finite diff.")
            : CheckingModel(relativeTol, absTol, description)
            , function_(function)
            , param_(param)
            , shift_(shift)
        { }

        // Returns second order finite difference derivative of function_ with respect to param_.
        inline Real operator()(const std::vector<Real>& optionData) const
        {
            std::vector<Real> optionData_left(optionData);
            optionData_left[param_] -= shift_;
            std::vector<Real> optionData_right(optionData);
            optionData_right[param_] += shift_;
            return ((*function_)(optionData_right)+(*function_)(optionData_left)-2 * (*function_)(optionData)) / (shift_ * shift_);
        }

    protected:
        func_type function_;
        DataType param_;
        double shift_;
    };


    inline Real delta(const std::vector<Real>& optionData)
    {
        return getBC(optionData).delta(optionData[stock]);
    }
    inline Real gamma(const std::vector<Real>& optionData)
    {
        return getBC(optionData).gamma(optionData[stock]);
    }
    inline Real vega(const std::vector<Real>& optionData)
    {
        return getBC(optionData).vega(optionData[DataType::time]);
    }
    inline Real theta(const std::vector<Real>& optionData)
    {
        return getBC(optionData).theta(optionData[stock], optionData[DataType::time]);
    }
    inline Real rho(const std::vector<Real>& optionData)
    {
        return getBC(optionData).rho(optionData[DataType::time]);
    }
    inline Real bcValue(const std::vector<Real>& optionData)
    {
        return getBC(optionData).value();
    }
    inline Real minusBcValue(const std::vector<Real>& optionData)
    {
        return -getBC(optionData).value();
    }

    // Struct for plots recording.
    template <DataType var, Greek greek>
    struct OutputData
    {
        static std::deque<std::string > get_columns()
        {
            return{ toString(var, false, true), GreekTraits<greek>::name() };
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, OutputData& v)
        {
                for (Size i = 0; i < v.data_.size() - 1; i++)
                {
                    stm << v.data_[i] << ";";
                }
                stm << v.data_.back() << std::endl;
                return stm;
            }

        std::deque<Real> data_;
    };


    struct GreekTestData
    {
        GreekTestData()
            : size_()
            , data_(5)
            , totalPrice_()
            , model_()
            , prices_()
        { }

        // parameters:
        // size - size of portfolio,
        // var - uniformly distributed option parameter,
        // left, rigth - bounds of distribution.
        void setUniformData(size_t size, DataType var, Real left, Real right)
        {
            size_ = size;
            data_[rate].resize(size, 0.05);
            data_[stock].resize(size, 100.0);
            data_[strike].resize(size, 120.0);
            data_[sigma].resize(size, 0.3);
            data_[time].resize(size, 1.0);
            if (size != 1)
            {
                for (Size i = 0; i < size; i++)
                {
                    data_[var][i] = left + ((right - left) * i) / (size - 1);
                }
            }
            else
            {
                data_[var][0] = (left + right) / 2;
            }
        }

        void setPseudorandomData(size_t size)
        {
            size_ = size;
            data_[rate].resize(size, 0.05);
            data_[stock].resize(size);
            data_[strike].resize(size);
            data_[sigma].resize(size);
            data_[time].resize(size);
            // Pseudorandom values.
            for (Size i = 0; i < size; i++)
            {
                data_[stock][i] = 100.0 + std::log(i + 1);
                data_[strike][i] = data_[stock][i] * (1 + 0.1 * (0.87 * std::cos(1.26 * i) + 1.23 * std::sin(2.48 * i)));
                data_[sigma][i] = 0.3 + 0.2 * std::sin(i);
                data_[time][i] = 2.0 + std::cos(2 * i);
            }
        }

        void setAnalyticalModel(std::shared_ptr<CheckingModel> model)
        {
            model_ = model;
        }

        // Calculates price of portfolio and each option.
        void calculatePrices()
        {
            totalPrice_.resize(1, 0);
            prices_.resize(size_);
            for (Size i = 0; i < size_; i++)
            {
                prices_[i] = getBC(optionData(i)).value();
                totalPrice_.front() += prices_[i];
            }
        }

        // Calculates Greeks using analytical method.
        std::vector<Real> analytical()
        {
            std::vector<Real> res(size_);
            for (Size i = 0; i < size_; i++)
            {
                res[i] = (*model_)(optionData(i));
            }
            return res;
        }

        // Returns parameters for one option in format { stock, strike, sigma, time, rate}.
        std::vector<Real> optionData(Size index)
        {
            std::vector<Real> optionData(5);
            for (Size dt = stock; dt <= rate; dt++)
            {
                optionData[dt] = data_[dt][index];
            }
            return optionData;
        }

        Size size_;
        std::vector<std::vector<cl::tape_double>> data_;
        std::vector<cl::tape_double> totalPrice_;
        std::shared_ptr<CheckingModel> model_;
        std::vector<Real> prices_;
    };

    template <class GreekTraits, class Meta = void>
    struct GreekTestBase
        : GreekTestData
    {
        typedef GreekTraits greek_traits;
        typedef DerivCalculator<
            greek_traits::order
            , greek_traits::is_mixed
            , greek_traits::sign_factor
        > calculator;
        static const DataType independent1 = greek_traits::independent1;

        void recordTape(std::unique_ptr<cl::tape_function<double>>& f)
        {
            cl::Independent(data_[independent1]);
            calculatePrices();
            f = std::make_unique<cl::tape_function<double>>(data_[independent1], totalPrice_);
        }

        Size indepVarNumber()
        {
            return size_;
        }
    };

    // Specialization for mixed Greeks.
    template <class GreekTraits>
    struct GreekTestBase<GreekTraits
        , typename std::enable_if<
            GreekTraits::is_mixed
        >::type>
        : GreekTestData
    {
        typedef GreekTraits greek_traits;
        typedef DerivCalculator<
            greek_traits::order
            , greek_traits::is_mixed
            , greek_traits::sign_factor
        > calculator;
        static const DataType independent1 = greek_traits::independent1;
        static const DataType independent2 = greek_traits::independent2;

        void recordTape(std::unique_ptr<cl::tape_function<double>>& f)
        {
            std::vector<cl::tape_double> indVector(2 * size_);
            std::copy(data_[independent1].begin(), data_[independent1].end(), indVector.begin());
            std::copy(data_[independent2].begin(), data_[independent2].end(), indVector.begin() + size_);
            cl::Independent(indVector);
            std::copy(indVector.begin(), indVector.begin() + size_, data_[independent1].begin());
            std::copy(indVector.begin() + size_, indVector.end(), data_[independent2].begin());

            calculatePrices();

            f = std::make_unique<cl::tape_function<double>>(indVector, totalPrice_);
        }

        Size indepVarNumber()
        {
            return 2 * size_;
        }
    };

    template <class GreekTraits>
    struct GreekTest
        : GreekTestBase<GreekTraits>
        , cl::AdjointTest<GreekTest<GreekTraits>>
    {
        // Constructor parameters:
        // logger - pointer to output stream which contains log stream.
        explicit GreekTest(cl::tape_empty_test_output* logger = nullptr)
            : GreekTestBase()
            , AdjointTest()
        {
            setLogger(logger);
            setAnalyticalModel();
        }

        // Constructor parameters:
        // size - size of portfolio.
        // logger - pointer to output stream that contains log stream.
        explicit GreekTest(Size size, cl::tape_empty_test_output* logger = nullptr)
            : GreekTestBase()
            , AdjointTest()
        {
            setLogger(logger);
            setPseudorandomData(size);
            setAnalyticalModel();
        }

        void setAnalyticalModel()
        {
            switch (GreekTraits::greek)
            {
            case Delta:
                GreekTestData::setAnalyticalModel(std::make_shared<FdCheckingModel>(&bcValue, stock, 1e-6, 1e-3 ));
                break;
            case Vega:
                GreekTestData::setAnalyticalModel(std::make_shared<FdCheckingModel>(&bcValue, sigma, 1e-6, 1e-4));
                break;
            case Theta:
                GreekTestData::setAnalyticalModel(std::make_shared<FdCheckingModel>(&minusBcValue, time, 1e-6, 1e-4));
                break;
            case Rho:
                GreekTestData::setAnalyticalModel(std::make_shared<FdCheckingModel>(&bcValue, rate, 1e-6, 1e-4));
                break;
            case Gamma:
                GreekTestData::setAnalyticalModel(std::make_shared<FdCheckingModel>(&delta, stock, 1e-6, 1e-4));
                break;
            case Vomma:
                GreekTestData::setAnalyticalModel(std::make_shared<FdCheckingModel>(&vega, sigma, 3e-7, 1e-4));
                break;
            case Vanna:
                GreekTestData::setAnalyticalModel(std::make_shared<FdCheckingModel>(&delta, sigma, 3e-7, 1e-4, 1e-10));
                break;
            case Charm:
                GreekTestData::setAnalyticalModel(std::make_shared<FdCheckingModel>(&theta, stock, 1e-6, 1e-4));
                break;
            case Veta:
                GreekTestData::setAnalyticalModel(std::make_shared<FdCheckingModel>(&vega, time, 1e-6, 1e-6, 1e-8));
                break;
            case Vera:
                GreekTestData::setAnalyticalModel(std::make_shared<FdCheckingModel>(&rho, sigma, 3e-7, 1e-4));
                break;
            case Speed:
                GreekTestData::setAnalyticalModel(std::make_shared<FdCheckingModel>(&gamma, stock, 1e-5, 1e-7));
                break;
            case Ultima:
                GreekTestData::setAnalyticalModel(std::make_shared<Fd2CheckingModel>(&vega, sigma, 1e-4));
                break;
            case Zomma:
                GreekTestData::setAnalyticalModel(std::make_shared<FdCheckingModel>(&gamma, sigma, 1e-6, 1e-4));
                break;
            case Color:
                GreekTestData::setAnalyticalModel(std::make_shared<FdCheckingModel>(&gamma, time, 1e-6, 1e-6, 1e-8));
                break;
            }
        }

        void recordTape()
        {
            assert(!data_[independent1].empty());

            GreekTestBase::recordTape(f_);
        }

        void calcReverse()
        {
            calculator::calcReverse(*f_, reverseResults_);
        }

        void calcForward()
        {
            calculator::calcForward(*f_, forwardResults_);
        }

        void calcAnalytical()
        {
            analyticalResults_ = GreekTestData::analytical();
        }

        Size indepVarNumber()
        {
            return GreekTestBase::indepVarNumber();
        }

        Size minPerfIteration() { return iterNumFactor; }

        double relativeTol() const { return model_->relativeTol_; }

        double absTol() const { return model_->absTol_; }

        std::string analyticalName() const { return model_->description_; }
    };

    template <Greek greek>
    using Test = GreekTest<GreekTraits<greek>>;


    template <Greek greek>
    struct TestData
    {
        typedef Test<greek> test_type;

        TestData()
        {}

        bool makeOutput()
        {
            bool ok = true;
            if (pointNo > 0)
            {
                ok &= recordDependencePlots();
            }
            return ok;
        }

        std::shared_ptr<test_type> makeTest()
        {
            return std::make_shared<test_type>();
        }

        std::shared_ptr<test_type> getTest(size_t size)
        {
            auto test = makeTest();
            test->setPseudorandomData(size);
            return test;
        }

        std::shared_ptr<test_type> getUniformTest(size_t size, DataType var, Real left, Real right)
        {
            auto test = makeTest();
            test->setUniformData(size, var, left, right);
            return test;
        }

        bool recordDependencePlots()
        {
            bool ok = true;
            switch (greek)
            {
            case Delta:
                ok &= recordDependencePlot<stock>(30, 250);
                break;
            case Vega:
                ok &= recordDependencePlot<sigma>(0.015, 0.5);
                break;
            case Theta:
                ok &= recordDependencePlot<time >(0.001, 2);
                break;
            case Rho:
                ok &= recordDependencePlot<rate >(0.0001, 0.5);
                break;
            case Gamma:
                ok &= recordDependencePlot<stock>(10, 200);
                break;
            case Vomma:
                ok &= recordDependencePlot<sigma>(0.015, 0.25);
                break;
            case Vanna:
                ok &= recordDependencePlot<stock>(10, 300);
                ok &= recordDependencePlot<sigma>(0.015, 0.9);
                break;
            case Charm:
                ok &= recordDependencePlot<stock>(50, 200);
                ok &= recordDependencePlot<time >(0.001, 0.5);
                break;
            case Veta:
                ok &= recordDependencePlot<sigma>(0.02, 0.9);
                ok &= recordDependencePlot<time >(0.01, 3);
                break;
            case Vera:
                ok &= recordDependencePlot<rate >(0.0001, 1);
                ok &= recordDependencePlot<sigma>(0.0003, 0.5);
                break;
            case Speed:
                ok &= recordDependencePlot<stock>(40, 180);
                break;
            case Ultima:
                ok &= recordDependencePlot<sigma>(0.015, 0.2);
                break;
            case Zomma:
                ok &= recordDependencePlot<stock>(10, 200);
                ok &= recordDependencePlot<sigma>(0.03, 0.15);
                break;
            case Color:
                ok &= recordDependencePlot<stock>(40, 250);
                ok &= recordDependencePlot<time >(0.01, 0.5);
                break;
            default:
                ok = false;
            }
            return ok;
        }

        // Makes plots for greek dependence.
        // var - option parameter on which we record dependence.
        // Function parameters:
        // left, right - bounds of variation of var parameter.
        template <DataType var>
        bool recordDependencePlot(Real left, Real right)
        {
            std::string filename = toString(var) + "Dependence";
            std::string title = GreekTraits<greek>::name() + " dependence on " + toString(var, true);
            cl::tape_empty_test_output out(OUTPUT_FOLDER_NAME "//" + GreekTraits<greek>::name() + "CallPortfolio//output", {
                { "filename", filename }
                , { "not_clear", "Not" }
                , { "title", title }
                , { "ylabel", "Greek value" }
                , { "line_box_width", "2" }
                , { "cleanlog", "false" }
            });
            std::vector<OutputData<var, greek>> outData(pointNo);
            auto test = getUniformTest(pointNo, var, left, right);
            test->setLogger(&out);
            bool ok = test->testReverse();
            for (Size i = 0; i < pointNo; i++)
            {
                outData[i].data_.push_back(test->data_[var][i]);
                outData[i].data_.push_back(test->reverseResults_[i]);
            }
            out << outData;
            return ok;
        }
    };
}

#endif