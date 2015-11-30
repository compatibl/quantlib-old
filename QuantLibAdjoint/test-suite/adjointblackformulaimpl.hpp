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

#ifndef cl_adjoint_black_formula_impl_hpp
#define cl_adjoint_black_formula_impl_hpp
#pragma once

#include "adjointblackformulatest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include <ql/pricingengines/blackformula.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;

#define OUTPUT_FOLDER_NAME "AdjointBlackFormula"

namespace
{
    enum DataType
    {
          tte
        , stdDev
        , moneyness
    };

    inline Real bachelierImpliedVol(const std::vector<Real>& data)
    {
        Real strike = 1.0 - 0.01 * data[moneyness] * std::sqrt(data[tte]);
        Real callPrem = bachelierBlackFormula(Option::Call
                                              , strike
                                              , 1.0
                                              , data[stdDev]
                                              , 0.95);
        return bachelierBlackFormulaImpliedVol(Option::Call
                                               , strike
                                               , 1.0
                                               , data[tte]
                                               , callPrem
                                               , 0.95);

    }

    struct Test :
        cl::AdjointTest<Test>
    {
        static const CalcMethod default_method = other;

        explicit Test(size_t size, DataType var)
            : AdjointTest()
            , size_(size)
            , data_(3)
            , indepDouble_(size)
            , implBpVol_(size)
            , var_(var)
        {
            switch (var_)
            {
                case tte:
                    setUniformData(size_, DataType::tte, 10.0, 20.0);
                    break;
                case stdDev:
                    setUniformData(size_, DataType::stdDev, 0.03, 0.08);
                    break;
                case moneyness:
                    setUniformData(size_, DataType::moneyness, -3.0, 3.0);
                    break;
            }
        }

        Size indepVarNumber() { return size_; }

        Size depVarNumber() { return size_; }

        Size minPerfIteration() { return 0; }

        void recordTape()
        {
            cl::Independent(data_[var_]);
            calculateImplBpVol();
            f_ = std::make_unique<cl::tape_function<double>>(data_[var_], implBpVol_);
        }

        void calculateImplBpVol()
        {
            implBpVol_.resize(size_, 0);
            for (Size i = 0; i < size_; i++)
            {
                implBpVol_[i] = bachelierImpliedVol(optionData(i));
            }
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            analyticalResults_.resize(size_*size_, 0.0);
            std::vector<Real> shiftRight, shiftLeft;
            Real h = 1e-4;
            for (Size i = 0; i < size_; i++)
            {
                shiftRight = shiftLeft = optionData(i);
                shiftRight[var_] += h;
                shiftLeft[var_] -= h;
                analyticalResults_[i*size_ + i] = (bachelierImpliedVol(shiftRight) - bachelierImpliedVol(shiftLeft)) / (2 * h);
            }
        }

        // Calculates derivatives using adjoint.
        void calcAdjoint()
        {
            adjointResults_ = f_->Jacobian(indepDouble_);
        }

        void setUniformData(size_t size, DataType var, double left, double right)
        {
            size_ = size;
            data_[tte].resize(size, 10.0);
            data_[stdDev].resize(size, 0.03);
            data_[moneyness].resize(size, -3);

            if (size != 1)
            {
                for (Size i = 0; i < size; i++)
                {
                    data_[var][i] = left + ((right - left) * i) / (size - 1);
                    indepDouble_[i] = left + ((right - left) * i) / (size - 1);
                }
            }
            else
            {
                data_[var][0] = (left + right) / 2;
                indepDouble_[0] = (left + right) / 2;
            }
        }

        // Returns parameters for one option in format { tte, stdDev, moneyness}.
        std::vector<Real> optionData(Size index)
        {
            std::vector<Real> optionData(3);
            for (Size dt = tte; dt <= moneyness; dt++)
            {
                optionData[dt] = data_[dt][index];
            }
            return optionData;
        }

        double relativeTol() const { return 1e-3; }

        double absTol() const { return 1e-5; }

        Size size_;
        std::vector<std::vector<cl::tape_double>> data_;
        std::vector<double> indepDouble_;
        std::vector<cl::tape_double> implBpVol_;
        DataType var_;
    };
}

#endif