/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2003 Ferdinando Ametrano
Copyright (C) 2007, 2008 Klaus Spanderen
Copyright (C) 2007 Neil Firth
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

//based on matrices.cpp from test-suite

#ifndef cl_adjoint_matrices_impl_hpp
#define cl_adjoint_matrices_impl_hpp
#pragma once

#include "adjointmatricestest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include <ql/math/matrix.hpp>
#include <ql/math/matrixutilities/pseudosqrt.hpp>
#include <ql/math/matrixutilities/symmetricschurdecomposition.hpp>
#include <ql/math/randomnumbers/mt19937uniformrng.hpp>
#include <ql/math/matrixutilities/basisincompleteordered.hpp>


using namespace QuantLib;
using namespace boost::unit_test_framework;


namespace
{
    enum
    {
#if defined CL_GRAPH_GEN
        // Number of points for performance plot.
        iterNo = 30,
        // Step for portfolio size for performance testing .
        step = 1,
#else
        // Number of points for performance plot.
        iterNo = 0,
        // Step for portfolio size for performance testing .
        step = 1,
#endif
    };


    struct CommonVars
    {
        double timeTapeRecording_;
        double timeAdjoint_;
        double timeAnalytical_;
        std::vector<PerformanceTime> performanceTime_;
        std::vector<AdjointTime> performanceAdjointTime_;
        std::vector<TapeSize> tapeSize_;
    };

    template <typename T = cl::tape_double>
    inline std::vector<T> generateVector(Size dim, double seed = 1.2345)
    {
        seed += 0.412176013753151;
        std::vector<T> vec;
        vec.reserve(dim * dim);
        int size = int(dim);
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                vec.push_back((1.4 + sin(i + seed) + cos(seed * j)) * exp(-(i - j) * (i - j)));
            }
        }
        return vec;
    }

    inline std::vector<double> generateSymVector(Size size)
    {
        std::vector<double> vec(size * size, 0);
        for (Size i = 0; i < size; i++)
        {
            vec[i*size + i] = 1;
            for (Size j = i + 1; j < size; j++)
            {
                vec[i*size + j] = 1 - j*0.7 / size;
                vec[j*size + i] = 1 - j*0.7 / size;
            }
        }
        return vec;
    }

    template <class Cont>
    typename std::iterator_traits<typename Cont::iterator>::value_type norm(const Cont& c)
    {
        // Matrix does not define value_type directly, but defines iterator.
        typedef std::iterator_traits<typename Cont::iterator>::value_type value_type;
        value_type sum = 0;
        for (auto iter = c.begin(); iter != c.end(); ++iter)
        {
            sum += (*iter) * (*iter);
        }
        return std::sqrt(sum);
    }

    // return $\sum_k a_k M^k$
    Matrix polynom(const Matrix& M, std::vector<Real> a)
    {
        Matrix result(M.rows(), M.columns(), 0);
        Matrix power(M.rows(), M.columns(), 0);
        for (Size i = 0; i < result.rows(); i++)
        {
            power[i][i] = 1;
        }
        for (Size k = 0; k < a.size(); k++)
        {
            result += a[k] * power;
            power = M * power;
        }
        return result;
    }

    Real cofactor(Matrix& A, Size i, Size j)
    {
        std::vector<Real> row(A.row_begin(i), A.row_end(i));
        std::fill(A.row_begin(i), A.row_end(i), 0);
        A[i][j] = 1;
        Real result = determinant(A);
        std::copy(row.begin(), row.end(), A.row_begin(i));
        return result;
    }

    std::vector<Real> calculateEigen(Size n, Matrix& M1)
    {
        std::vector<Real> det(n);
        Array eigenValues;
        SymmetricSchurDecomposition dec(M1);
        eigenValues = dec.eigenvalues();
        for (Size i = 0; i < n; i++)
        {
            det[i] = eigenValues[i];
        }
        return det;
    }


    struct InverseTest
        : cl::AdjointTest<InverseTest>
    {
        static const CalcMethod default_method = other;

        explicit InverseTest(Size dim, cl::tape_empty_test_output* logger = nullptr)
            : AdjointTest()
            , dim_(dim)
            , size_(dim_ * dim_)
            , doubleX_()
            , X_(size_)
            , Y_()
            , A_()
            , invA_()
        {
            setLogger(logger);

            doubleX_ = generateVector<double>(dim_);
            std::copy(doubleX_.begin(), doubleX_.end(), X_.begin());
        }

        Size indepVarNumber() { return size_; }

        Size depVarNumber() { return size_; }

        Size minPerfIteration() { return 0; }

        void recordTape()
        {
            cl::Independent(X_);
            A_ = vector2matrix(X_, dim_, dim_);
            invA_ = inverse(A_);
            Y_ = matrix2vector(invA_);
            f_ = std::make_unique<cl::tape_function<double>>(X_, Y_);
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            analyticalResults_.resize(size_ * size_);
            Matrix invA = inverse(vector2matrix(X_, dim_, dim_));
            Matrix dA(dim_, dim_, 0);
            auto pda = dA.begin();
            for (Size k = 0; k < size_; k++)
            {
                *pda = 1;
                Matrix dInvA = -1 * invA * dA * invA;
                *pda++ = 0;
                for (Size i = 0; i < size_; i++)
                {
                    analyticalResults_[i * size_ + k] = *(dInvA.begin() + i);
                }
            }
        }

        // Calculates derivatives using adjoint.
        void calcAdjoint()
        {
            adjointResults_ = f_->Jacobian(doubleX_);
        }

        Size dim_;
        Size size_;
        std::vector<double> doubleX_;
        std::vector<cl::tape_double> X_;
        std::vector<cl::tape_double> Y_;
        Matrix A_;
        Matrix invA_;
    };

    struct InverseTestData
    {
        explicit InverseTestData()
        : outPerform_("AdjointMatrices//Inverse", {
            { "title", "Inverse matrices differentiation performance" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of matrix elements" }
            , { "line_box_width", "-5" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        })
            , outAdjoint_("AdjointMatrices//Inverse", {
                { "title", "Inverse matrices adjoint differentiation performance" }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of matrix elements" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        })
            , outSize_("AdjointMatrices//Inverse", {
                { "title", "Tape size dependence on number of matrix elements" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
            , { "cleanlog", "false" }
        })
        {
        }

        std::shared_ptr<InverseTest> getTest(size_t dim)
        {
            return std::make_shared<InverseTest>(dim, &outPerform_);
        }

        bool makeOutput()
        {
            return cl::recordPerformance(*this, iterNo, step);
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
    };

    struct DeterminantTest
        : cl::AdjointTest<DeterminantTest>
    {
        explicit DeterminantTest(Size dim, cl::tape_empty_test_output* logger = nullptr)
        : AdjointTest()
        , dim_(dim)
        , size_(dim_ * dim_)
        , X_()
        , det_()
        {
            setLogger(logger);

            X_ = generateVector(dim_);
        }

        Size indepVarNumber() { return size_; }

        Size minPerfIteration() { return 0; }

        void recordTape()
        {
            cl::Independent(X_);
            Matrix A = vector2matrix(X_, dim_, dim_);
            det_.resize(1, determinant(A));
            f_ = std::make_unique<cl::tape_function<double>>(X_, det_);
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            analyticalResults_.resize(size_);
            Matrix A = vector2matrix(X_, dim_, dim_);
            for (Size i = 0; i < dim_; i++)
            {
                for (Size j = 0; j < dim_; j++)
                {
                    analyticalResults_[i * dim_ + j] = cofactor(A, i, j);
                }
            }
        }

        Size dim_;
        Size size_;
        std::vector<cl::tape_double> X_;
        std::vector<cl::tape_double> det_;
    };

    struct DeterminantTestData
    {
        explicit DeterminantTestData()
        : outPerform_("AdjointMatrices//Determinant", {
            { "title", "Matrix determinant differentiation performance" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of matrix elements" }
            , { "line_box_width", "-5" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        })
            , outAdjoint_("AdjointMatrices//Determinant", {
                { "title", "Matrix determinant adjoint differentiation performance" }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of matrix elements" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        })
            , outSize_("AdjointMatrices//Determinant", {
                { "title", "Tape size dependence on number of matrix elements" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
            , { "cleanlog", "false" }
        })
        {
        }

        std::shared_ptr<DeterminantTest> getTest(size_t dim)
        {
            return std::make_shared<DeterminantTest>(dim, &outPerform_);
        }

        bool makeOutput()
        {
            return cl::recordPerformance(*this, iterNo, step);
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
    };

    struct SqrtTest
        : cl::AdjointTest<SqrtTest>
    {
        explicit SqrtTest(Size dim)
        : AdjointTest()
        , dim_(dim)
        , size_(dim_ * dim_)
        , doubleX_()
        , X_(size_)
        , Y_(1)
        , A_()
        , sqrtA_()
        {
            doubleX_ = { 5.0, 0.9, 0.8, 0.7,
                0.9, 2.0, 0.6, 0.5,
                0.8, 0.6, 3.0, 0.4,
                0.7, 0.5, 0.4, 4.0 };

            std::copy(doubleX_.begin(), doubleX_.end(), X_.begin());
            A_ = vector2matrix(X_, dim_, dim_);
        }

        Size indepVarNumber() { return size_; }

        Size depVarNumber() { return 1; }

        Size minPerfIteration() { return 0; }

        void recordTape()
        {
            cl::Independent(X_);
            sqrtA_ = pseudoSqrt(A_, SalvagingAlgorithm::None);
            Y_[0] = determinant(sqrtA_);
            f_ = std::make_unique<cl::tape_function<double>>(X_, Y_);
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            double h = 1.0e-6;
            analyticalResults_.resize(size_, 0.0);
            Size ii, jj;
            Matrix fm;
            for (Size i = 0; i < size_; i++)
            {
                ii = i / dim_;
                jj = i %  dim_;

                if (ii <= jj)
                {
                    A_[ii][jj] = A_[jj][ii] += h;
                    fm = pseudoSqrt(A_, SalvagingAlgorithm::None);
                    analyticalResults_[i] = (determinant(fm) - Y_[0]) / h;
                    A_[ii][jj] = A_[jj][ii] -= h;
                }
            }
        }

        double relativeTol() const { return 1e-4; }

        double absTol() const { return 1e-10; }

        Size dim_;
        Size size_;
        std::vector<Real> doubleX_;
        std::vector<cl::tape_double> X_;
        std::vector<cl::tape_double> Y_;
        Matrix A_;
        Matrix sqrtA_;
    };

    struct PolynomTest
        : cl::AdjointTest<PolynomTest>
    {
        explicit PolynomTest(Size dim)
        : AdjointTest()
        , dim_(dim)
        , size_(dim_ * dim_)
        , X_(1, 1.0)
        , P_({ 5.4, 1.2, 0.3, 0.1 })
        , Y_()
        , A_()
        , sqrtA_()
        {
            A_ = X_[0] * vector2matrix(generateVector<double>(dim_), dim_, dim_);
        }

        Size indepVarNumber() { return 1; }

        Size depVarNumber() { return size_; }

        Size minPerfIteration() { return 0; }

        void recordTape()
        {
            cl::Independent(X_);
            Matrix PA = polynom(A_, P_);
            Y_ = matrix2vector(PA);
            f_ = std::make_unique<cl::tape_function<double>>(X_, Y_);
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            std::vector<Real> dP(P_.size());
            for (Size i = 1; i < P_.size(); i++)
            {
                dP[i - 1] = (i * P_[i]);
            }
            Matrix expected = polynom(A_, dP) * A_;
            analyticalResults_ = matrix2vector(expected);
        }

        double relativeTol() const { return 1e-4; }

        double absTol() const { return 1e-10; }

        Size dim_;
        Size size_;
        std::vector<cl::tape_double> X_;
        std::vector<Real> P_;
        std::vector<cl::tape_double> Y_;
        Matrix A_;
        Matrix sqrtA_;
    };

    struct EigenvectorsTest
        : cl::AdjointTest<EigenvectorsTest>
    {
        static const CalcMethod default_method = other;

        explicit EigenvectorsTest(Size dim, cl::tape_empty_test_output* logger = nullptr)
            : AdjointTest()
            , dim_(dim)
            , size_(dim_ * dim_)
            , doubleX_()
            , X_(size_)
            , Y_(dim_)
            , A_()
        {
            setLogger(logger);

            doubleX_ = generateSymVector(dim_);
            std::copy(doubleX_.begin(), doubleX_.end(), X_.begin());
        }

        Size indepVarNumber() { return size_; }

        Size depVarNumber() { return dim_; }

        Size minPerfIteration() { return 0; }

        void recordTape()
        {
            cl::Independent(X_);
            A_ = vector2matrix(X_, dim_, dim_);
            Y_ = calculateEigen(dim_, A_);
            f_ = std::make_unique<cl::tape_function<double>>(X_, Y_);
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            double h = 1.0e-10;
            analyticalResults_.resize(size_ * dim_, 0);
            std::vector<Real> eigen_values_right, eigen_values_left;
            A_ = vector2matrix(X_, dim_, dim_);
            Y_ = calculateEigen(dim_, A_);
            for (Size k = 0; k < dim_; k++)
            {
                for (Size l = k; l < dim_; l++)
                {
                    A_[k][l] += h;
                    eigen_values_right = calculateEigen(dim_, A_);
                    A_[k][l] -= 2 * h;
                    eigen_values_left = calculateEigen(dim_, A_);
                    for (Size s = 0; s < dim_; s++)
                    {
                        analyticalResults_[s*dim_*dim_ + k*dim_ + l] = (eigen_values_right[s] - eigen_values_left[s]) / (2 * h);
                    }
                    A_[k][l] += h;
                }
            }
        }

        // Calculates derivatives using adjoint.
        void calcAdjoint()
        {
            adjointResults_ = f_->Jacobian(doubleX_);
        }

        template <CalcMethod method = test_type::default_method>
        bool checkAdjoint()
        {
            for (Size k = 0; k < dim_; k++)
            {
                for (Size l = k + 1; l < dim_; l++)
                {
                    for (Size s = 0; s < dim_; s++)
                    {
                        if (fabs(adjointResults_[s*dim_*dim_ + k*dim_ + l] - analyticalResults_[s*dim_*dim_ + k*dim_ + l].value()) > 1e-4)
                        {
                            BOOST_FAIL("Derivatives of eigenvalues is not satisfied: "
                                       << "\nDerivatives of eigenvalues using Quantlib function  = " << adjointResults_[s*dim_*dim_ + k*dim_ + l]
                                       << "\nDerivatives of eigenvalues using finite differences = " << analyticalResults_[s*dim_*dim_ + k*dim_ + l]);
                            return false;
                        }
                    }
                }
            }
            return true;
        }

        double relativeTol() const { return 1e-2; }

        double absTol() const { return 1e-4; }

        Size dim_;
        Size size_;
        std::vector<double> doubleX_;
        std::vector<cl::tape_double> X_;
        std::vector<cl::tape_double> Y_;
        Matrix A_;
    };

    struct EigenvectorsTestData
    {
        explicit EigenvectorsTestData()
        : outPerform_("AdjointMatrices//Eigenvectors", {
            { "title", "Eigenvectors differentiation performance" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of matrix elements" }
            , { "line_box_width", "-5" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        })
            , outAdjoint_("AdjointMatrices//Eigenvectors", {
                { "title", "Eigenvectors adjoint differentiation performance" }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of matrix elements" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        })
            , outSize_("AdjointMatrices//Eigenvectors", {
                { "title", "Tape size dependence on number of matrix elements" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
            , { "cleanlog", "false" }
        })
        {
        }

        std::shared_ptr<EigenvectorsTest> getTest(size_t dim)
        {
            return std::make_shared<EigenvectorsTest>(dim, &outPerform_);
        }

        bool makeOutput()
        {
            return cl::recordPerformance(*this, iterNo, step);
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
    };

    struct OrthogonalProjectionTest
        : cl::AdjointTest<OrthogonalProjectionTest>
    {

        explicit OrthogonalProjectionTest(Size dim, cl::tape_empty_test_output* logger = nullptr)
        : AdjointTest()
        , dim_(dim)
        , numberVectors_(8)
        , X_(numberVectors_* dim_)
        , Y_(1)
        , A_()
        , multiplier_(100)
        {
            setLogger(logger);
            randomInitVector(X_);
        }

        void randomInitVector(std::vector<cl::tape_double>& A)
        {
            MersenneTwisterUniformRng rng(1);

            for (Size i = 0; i < numberVectors_; ++i)
            for (Size j = 0; j < dim_; ++j)
                A[i*dim_ + j] = rng.next().value;
        }

        Size indepVarNumber() { return numberVectors_* dim_; }

        Size depVarNumber() { return 1; }

        Size minPerfIteration() { return 10; }

        void recordTape()
        {
            cl::Independent(X_);
            A_ = vector2matrix(X_, numberVectors_, dim_);
            Matrix projA_ = calcProjection(A_);
            Y_[0] = norm(projA_);
            f_ = std::make_unique<cl::tape_function<double>>(X_, Y_);
        }

        Matrix calcProjection(Matrix& A)
        {
            Matrix projA(numberVectors_, dim_);
            OrthogonalProjections projector(A, multiplier_, 1e-4);
            for (Size i = 0; i < numberVectors_; ++i)
            for (Size j = 0; j < dim_; ++j)
                projA[i][j] = projector.GetVector(i)[j];
            return projA;
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            double h = 1e-8;
            analyticalResults_.resize(numberVectors_* dim_);
            Matrix projRight, projLeft;
            A_ = vector2matrix(X_, numberVectors_, dim_);
            for (Size i = 0; i < numberVectors_; ++i)
            {
                for (Size j = 0; j < dim_; ++j)
                {
                    A_[i][j] += h;
                    projRight = calcProjection(A_);
                    A_[i][j] -= 2 * h;
                    projLeft = calcProjection(A_);
                    analyticalResults_[i*dim_ + j] = (norm(projRight) - norm(projLeft)) / (2 * h);
                    A_[i][j] += h;
                }
            }
        }

        double relativeTol() const { return 1e-2; }

        double absTol() const { return 1e-4; }

        Size dim_;
        Size numberVectors_;
        std::vector<cl::tape_double> X_;
        std::vector<cl::tape_double> Y_;
        Matrix A_;
        Real multiplier_;
    };

    struct OrthogonalProjectionTestData
    {
        explicit OrthogonalProjectionTestData()
        : outPerform_("AdjointMatrices//OrthogonalProjection", {
            { "title", "Orthogonal projection differentiation performance" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of matrix elements" }
            , { "line_box_width", "-5" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        })
            , outAdjoint_("AdjointMatrices//OrthogonalProjection", {
                { "title", "Orthogonal projection adjoint differentiation performance" }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of matrix elements" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        })
            , outSize_("AdjointMatrices//OrthogonalProjection", {
                { "title", "Tape size dependence on number of matrix elements" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
            , { "cleanlog", "false" }
        })
        {
        }

        std::shared_ptr<OrthogonalProjectionTest> getTest(size_t dim)
        {
            return std::make_shared<OrthogonalProjectionTest>(dim, &outPerform_);
        }

        bool makeOutput()
        {
            return cl::recordPerformance(*this, iterNo, step);
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
    };
}
#endif
