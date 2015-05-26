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

#include "adjointmatricestest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include <ql/math/matrix.hpp>
#include <ql/math/matrixutilities/pseudosqrt.hpp>
//#include <ql/math/matrixutilities/svd.hpp>
#include <ql/math/matrixutilities/symmetricschurdecomposition.hpp>
#include <ql/math/randomnumbers/mt19937uniformrng.hpp>
//#include <ql/math/matrixutilities/qrdecomposition.hpp>
#include <ql/math/matrixutilities/basisincompleteordered.hpp>


using namespace QuantLib;
using namespace boost::unit_test_framework;


#ifdef CL_TAPE_CPPAD
namespace
{
    struct CommonVars
    {
        double timeTapeRecording_;
        double timeAdjoint_;
        double timeAnalytical_;
        std::vector<PerformanceTime> performanceTime_;
        std::vector<AdjointTime> performanceAdjointTime_;
        std::vector<TapeSize> tapeSize_;
    };

    inline std::vector<cl::TapeDouble> generateVector(Size sqrt_size, double seed = 1.2345)
    {
        std::vector<cl::TapeDouble> vec;
        vec.reserve(sqrt_size * sqrt_size);
        int size = (int)sqrt_size;
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                vec.push_back((1.4 + sin(i + seed) + cos(seed * j)) * exp(-(i - j) * (i - j)));
            }
        }
        return vec;
    }

    inline std::vector<cl::TapeDouble> generateSymVector(Size size)
    {
        std::vector<cl::TapeDouble> vec(size * size, 0);
        for (Size i = 0; i < size; i++)
        {
            vec[i*size + i] = 1;
            for (Size j = i+1; j < size; j++)
            {
                vec[i*size + j] = 1 - j*0.7 / size;
                vec[j*size + i] = 1 - j*0.7 / size;
            }
        }
        return vec;
    }

    inline Real norm(const Matrix& m)
    {
        Real sum = 0;
        for (Matrix::const_iterator iter = m.begin(); iter != m.end(); ++iter)
        {
            sum += (*iter) * (*iter);
        }
        return std::sqrt(sum);
    }

    template <class Cont>
    typename Cont::value_type norm(const Cont& c)
    {
        Cont::value_type sum = 0;
        for (Cont::const_iterator iter = c.begin(); iter != c.end(); ++iter)
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

    std::vector<cl::TapeDouble> calculateEigen(Size n, Matrix& M1)
    {
        std::vector<cl::TapeDouble> det(n);
        Array eigenValues;
        SymmetricSchurDecomposition dec(M1);
        eigenValues = dec.eigenvalues();
        for (Size i = 0; i < n; i++)
        {
            det[i] = eigenValues[i];
        }
        return det;
    }
}
#endif


bool AdjointMatricesTest::testInverse()
{
#ifdef CL_TAPE_CPPAD
    BOOST_TEST_MESSAGE("Testing Adjoint differentiation for inverse matrices...");
    Real tol = 1.0e-10;
    bool ok = true;
    boost::timer timer;
    CommonVars varsInv;
    for (Size dim = 1; dim <= 50; dim++)
    {
#ifndef CL_GRAPH_GEN
        dim = 10;
#endif
        std::vector<cl::TapeDouble> X = generateVector(dim);
        cl::Independent(X);
        timer.restart();
        Matrix A = vector2matrix(X, dim, dim);
        Matrix invA = inverse(A);
        std::vector<cl::TapeDouble> Y = matrix2vector(invA);
        cl::TapeFunction<double> f(X, Y);
        varsInv.timeTapeRecording_ = timer.elapsed();

        //Store size of tape
        varsInv.tapeSize_.push_back(TapeSize { dim, f.Memory() });

        // checking all partial derivatives
        std::vector<double> dX(X.size(), 0);
        varsInv.timeAdjoint_ = 0;
        varsInv.timeAnalytical_ = 0;
        for (Size k = 0; k < dX.size(); k++)
        {
            dX[k] = 1;
            timer.restart();
            std::vector<double> dY = f.Forward(1, dX);
            Matrix d_invA = vector2matrix(dY, dim, dim);
            varsInv.timeAdjoint_ += timer.elapsed();

            timer.restart();
            invA = inverse(A);
            Matrix dA = vector2matrix(dX, dim, dim);
            Matrix expected = -1 * invA * dA * invA;
            varsInv.timeAnalytical_ += timer.elapsed();


            Real error = norm(expected - d_invA);
            if (error > tol)
            {
                BOOST_ERROR("Failed to reproduce expected derivative of inverse martix"
                            << "\n    error:      " << error
                            << "\n    tolerance:  " << tol);
                ok = false;
            }
            dX[k] = 0;
        }
        varsInv.performanceTime_.push_back(PerformanceTime { varsInv.timeTapeRecording_, varsInv.timeAdjoint_, varsInv.timeAnalytical_, dim });
        varsInv.performanceAdjointTime_.push_back(AdjointTime { varsInv.timeAdjoint_, dim});

        cl::AdjointTestOutput outPerform("AdjointMatrices//Invert", { { "filename", "AdjointPerformance" }, { "not_clear", "Not" }, { "title", "Performance for Adjoint Inverting Matrices" }, { "ylabel", "Time" } });
        outPerform << varsInv.performanceTime_;


        cl::AdjointTestOutput outAdjoint("AdjointMatrices//Invert", { { "filename", "Adjoint" }, { "not_clear", "Not" }, { "title", "Adjoint" }, { "ylabel", "Time" } });
        outAdjoint << varsInv.performanceAdjointTime_;

        cl::AdjointTestOutput outSize("AdjointMatrices//Invert", { { "filename", "TapeSize" }, { "not_clear", "Not" }, { "title", "Size of Tape" }, { "ylabel", "Memory (bytes)" } });
        outSize << varsInv.tapeSize_;
        X = generateVector(dim, 3.3);
        cl::Independent(X);
        A = vector2matrix(X, dim, dim);
        invA = inverse(A);
        Matrix I1 = invA * A;
        Matrix I2 = A * invA;
        Real detA = determinant(A);
        Y.clear();
        Y.push_back(norm(I1));
        Y.push_back(norm(I2));
        Y.push_back(determinant(I1));
        Y.push_back(determinant(I2));
        Y.push_back(detA);
        Y.push_back(determinant(invA));
        cl::TapeFunction<double> f2(X, Y);

        std::vector<std::vector<double> > w;
        // checking that following gradients are zero:
        // I1 and I2 are constant so their norms and determinants are constant
        w.push_back({ 1, 0, 0, 0, 0, 0 });
        w.push_back({ 0, 1, 0, 0, 0, 0 });
        w.push_back({ 0, 0, 1, 0, 0, 0 });
        w.push_back({ 0, 0, 0, 1, 0, 0 });
        // checking that $\frac{\partial}{\partial A_i^j}detA^{-1} = -\frac{1}{(detA)^2}\frac{\partial}{\partial A_i^j}detA$
        w.push_back({ 0, 0, 0, 0, 1, static_cast<double>(detA * detA) });
        for (Size i = 0; i < w.size(); i++)
        {
            double gradient_norm = norm(f2.Reverse(1, w[i]));
            if (gradient_norm > tol)
            {
                BOOST_ERROR("Failed to reproduce expected 0 gradient (norm = "
                            << gradient_norm << ")"
                            << "\n    tolerance:  " << tol
                            << "\n    test case #" << i);
                ok = false;
            }
        }
#ifndef CL_GRAPH_GEN
        break;
#endif
    }
    return ok;
#endif
    return true;
}


bool AdjointMatricesTest::testDeterminant()
{
#ifdef CL_TAPE_CPPAD
    BOOST_TEST_MESSAGE("Testing Adjoint differentiation for matrix determinant...");
    double tol = 1e-9;
    bool ok = true;
    boost::timer timer;

    CommonVars vars;
    for (Size dim = 1; dim <= 50; dim ++)
    {
#ifndef CL_GRAPH_GEN
        dim = 10;
#endif
        std::vector<cl::TapeDouble> X = generateVector(dim);
        timer.restart();
        Independent(X);
        Matrix A = vector2matrix(X, dim, dim);
        Real detA = determinant(A);
        Matrix B = vector2matrix(generateVector(dim, 4.321), dim, dim);
        Real detB = determinant(B);
        Real detAB = determinant(A * B);
        Real detA3 = determinant(A * A * A);
        std::vector<cl::TapeDouble> Y;
        Y.push_back(detA);
        Y.push_back(detAB);
        Y.push_back(detA3);
        cl::TapeFunction<double> f(X, Y);
        vars.timeTapeRecording_ = timer.elapsed();

        //Store size of tape
        vars.tapeSize_.push_back(TapeSize { dim, f.Memory() });

        timer.restart();
        //Calculate derivatives in Reverse mode
        std::vector<double> grad_detA = f.Reverse(1, std::vector<double>({ 1, 0, 0 }));
        vars.timeAdjoint_ = timer.elapsed();
        timer.restart();
        for (Size i = 0; i < dim; i++)
        {
            for (Size j = 0; j < dim; j++)
            {
                //Calculate derivatives in analytical mode
                double expected = static_cast<double>(cofactor(A, i, j));
                if (abs(expected - grad_detA[dim * i + j]) > tol)
                {
                    BOOST_ERROR("Failed to reproduce expected derivative of determinant"
                                << "\n    i = " << i << ";  j = " << j
                                << "\n    expected:    " << expected
                                << "\n    calculated:  " << grad_detA[dim * i + j]
                                << "\n    tolerance:  " << tol);
                    ok = false;
                }
            }
        }
        vars.timeAnalytical_ = timer.elapsed();

        vars.performanceTime_.push_back(PerformanceTime { vars.timeTapeRecording_, vars.timeAdjoint_, vars.timeAnalytical_, dim });
        vars.performanceAdjointTime_.push_back(AdjointTime { vars.timeAdjoint_, dim});

        std::vector<std::vector<double> > w;
        // checking some equalities
        // $\frac{\partial}{\partial A_i^j}detAB = detB\frac{\partial}{\partial A_i^j}detA$
        w.push_back({ static_cast<double>(detB), -1, 0 });
        // $\frac{\partial}{\partial A_i^j}detA^3 = 3(detA)^2\frac{\partial}{\partial A_i^j}detA$
        w.push_back({ static_cast<double>(3 * detA * detA), 0, -1 });
        for (Size i = 0; i < w.size(); i++)
        {
            double gradient_norm = norm(f.Reverse(1, w[i]));
            if (gradient_norm > tol)
            {
                BOOST_ERROR("Failed to reproduce expected 0 gradient (norm = "
                            << gradient_norm << ")"
                            << "\n    tolerance:  " << tol
                            << "\n    test case #" << i);
                ok = false;
            }
        }
#ifndef CL_GRAPH_GEN
        break;
#endif
    }
    cl::AdjointTestOutput outPerform("AdjointMatrices//Determinant", { { "filename", "AdjointPerformance" }, { "not_clear", "Not" }, { "title", "Performance for Adjoint Determinant" }, { "ylabel", "Time" } });
    outPerform << vars.performanceTime_;


    cl::AdjointTestOutput outAdjoint("AdjointMatrices//Determinant", { { "filename", "Adjoint" }, { "not_clear", "Not" }, { "title", "Adjoint" }, { "ylabel", "Time" } });
    outAdjoint << vars.performanceAdjointTime_;

    cl::AdjointTestOutput outSize("AdjointMatrices//Determinant", { { "filename", "TapeSize" }, { "not_clear", "Not" }, { "title", "Size of Tape" }, { "ylabel", "Memory (bytes)" } });
    outSize << vars.tapeSize_;
    return ok;
#endif
    return true;
}


bool AdjointMatricesTest::testSqrt()
{
    bool result = false;
#ifdef CL_TAPE_CPPAD

    BOOST_TEST_MESSAGE("Testing matricial pseudo square root...");

    std::vector<Real> v = { 5.0, 0.9, 0.8, 0.7,
        0.9, 2.0, 0.6, 0.5,
        0.8, 0.6, 3.0, 0.4,
        0.7, 0.5, 0.4, 4.0 };

    Size sizeof_indep = v.size();
    Size sizeof_matrix = 4;
    SalvagingAlgorithm::Type salvagingAlgorithm = SalvagingAlgorithm::None;

    std::vector<cl::TapeDouble> val(sizeof_indep);
    for (Size i = 0; i < sizeof_indep; i++)
        val[i] = v[i];

    Independent(val);

    Matrix M = vector2matrix(val, sizeof_matrix, sizeof_matrix);
    Matrix m = pseudoSqrt(M, salvagingAlgorithm);

    std::vector<cl::TapeDouble> Y(1);
    Y[0] = determinant(m);

    cl::TapeFunction<double> f(val, Y);

    std::vector<double> sf_Forward, sf_Reverse;

    //Start differentiation in Forward mode
    gradForward(f, sf_Forward, true, false);

    //Start differentiation in Reverse mode
    gradReverse(f, sf_Reverse, true, false);

    //Finite differences
    double h = 1.0e-6;
    std::vector<Real> sf_Finite(sizeof_indep);

    for (Size i = 0; i < sizeof_indep; i++)
    {
        std::vector<Real> val(sizeof_indep);
        for (Size j = 0; j < sizeof_indep; j++)
            val[j] = v[j];

        Matrix FM = vector2matrix(val, sizeof_matrix, sizeof_matrix);
        Size ii = i / sizeof_matrix;
        Size jj = i % sizeof_matrix;

        if (ii > jj)
        {
            sf_Finite[i] = 0;
        }
        else
        {
            FM[ii][jj] = FM[jj][ii] += h;
            Matrix fm = pseudoSqrt(FM, salvagingAlgorithm);
            sf_Finite[i] = (determinant(fm) - Y[0]) / h;
        }
    }
    result = checkWithFiniteDiff(sf_Forward, sf_Reverse, sf_Finite, 1e-4, 1e-10);

#endif
    return result;
}


bool AdjointMatricesTest::testPolynom()
{
#ifdef CL_TAPE_CPPAD
    BOOST_TEST_MESSAGE("Testing matricial polynom differentiation...");
    double tol = 1e-10;
    Size dim = 10;
    std::vector<Real> P = { 5.4, 1.2, 0.3, 0.1 };
    std::vector<Real> dP(P.size());
    for (Size i = 1; i < P.size(); i++)
    {
        dP[i - 1] = (i * P[i]);
    }
    std::vector<cl::TapeDouble> t = { 1 };

    Independent(t);

    Matrix A = t[0] * vector2matrix(generateVector(dim), dim, dim);
    Matrix PA = polynom(A, P);
    std::vector<cl::TapeDouble> Y = matrix2vector(PA);
    cl::TapeFunction<double> f(t, Y);

    // $dP(tA) = dtP'(tA)A$
    bool ok = true;
    //Start differentiation in Forward mode
    std::vector<double> dY = f.Forward(1, std::vector<double>(1, 1));
    Matrix dPA = vector2matrix(dY, dim, dim);

    //Calculate derivatives in analytical mode
    Matrix expected = polynom(A, dP) * A;
    std::vector<cl::TapeDouble> exp = matrix2vector(expected);
    Real error = norm(dPA - expected);
    if (error > tol)
    {
        BOOST_ERROR("Failed to reproduce expected derivative"
                    << "\n    error:      " << error
                    << "\n    tolerance:  " << tol);
        ok = false;
    }
    return ok;
#endif
    return true;
}


bool AdjointMatricesTest::testEigenvectors()
{
#ifdef CL_TAPE_CPPAD
    BOOST_TEST_MESSAGE("Testing eigenvalues and eigenvectors calculation with AD...");
    bool result = true;
    //Size n = 3;
    double tol = 1.0e-5;
    boost::timer timer;
    CommonVars vars;
    for (Size dim = 2; dim <= 50; dim++)
    {
#ifndef CL_GRAPH_GEN
        dim = 50;
#endif
        std::vector<cl::TapeDouble> X = generateSymVector(dim);
        std::vector<cl::TapeDouble> eigen_values(dim, 0.0);

        std::vector<Real> sf_Rev(dim*dim*dim);
        std::vector<Real> sf_Finite(dim*dim*dim);

        timer.restart();
        cl::Independent(X);

        Matrix M = vector2matrix(X, dim, dim);
        eigen_values = calculateEigen(dim, M);

        cl::TapeFunction<double> f(X, eigen_values);
        vars.timeTapeRecording_ = timer.elapsed();

        //Store size of tape
        vars.tapeSize_.push_back(TapeSize { dim, f.Memory() });
        std::vector<double> sy(dim);
        std::vector<double> sf(dim*dim);

        //Start differentiation in Reverse mode
        timer.restart();
        for (Size i = 0; i < dim; i++)
        {
            sy[i] = 1;
            sf = f.Reverse(1, sy);
            sy[i] = 0;
            for (Size j = 0; j < dim*dim; j++)
                sf_Rev[i*dim*dim + j] = sf[j];
        }
        vars.timeAdjoint_ = timer.elapsed();
        //Finite differences
        timer.restart();
        double h = 1.0e-5;
        std::vector<cl::TapeDouble> eigen_values_right, eigen_values_left;
        for (Size k = 0; k < dim; k++)
        {
            for (Size l = k + 1; l < dim; l++)
            {
                M[k][l] += h;
                eigen_values_right = calculateEigen(dim, M);
                M[k][l] -= 2*h;
                eigen_values_left = calculateEigen(dim, M);
                for (Size s = 0; s < dim; s++)
                {
                    sf_Finite[s*dim*dim + k*dim + l] = (eigen_values_right[s] - eigen_values_left[s]) / (2*h);
                    if (fabs(sf_Rev[s*dim*dim + k*dim + l].value() - sf_Finite[s*dim*dim + k*dim + l].value()) > tol)
                    {
                        BOOST_FAIL("Derivatives of eigenvalues is not satisfied: "
                                   << "\nDerivatives of eigenvalues using Quantlib function  = " << sf_Rev[s*dim*dim + k*dim + l]
                                   << "\nDerivatives of eigenvalues using finite differences = " << sf_Finite[s*dim*dim + k*dim + l]);
                        result = false;
                    }
                }
                M[k][l] += h;
            }
        }
        vars.timeAnalytical_ = timer.elapsed();

        vars.performanceTime_.push_back(PerformanceTime { vars.timeTapeRecording_, vars.timeAdjoint_, vars.timeAnalytical_, dim });
        vars.performanceAdjointTime_.push_back(AdjointTime { vars.timeAdjoint_ , dim });
    }

    cl::AdjointTestOutput outPerform("AdjointMatrices//Eigenvectors", { { "filename", "AdjointPerformance" }, { "not_clear", "Not" }, { "title", "Performance for Adjoint Eigenvectors" }, { "ylabel", "Time" } });
    outPerform << vars.performanceTime_;


    cl::AdjointTestOutput outAdjoint("AdjointMatrices//Eigenvectors", { { "filename", "Adjoint" }, { "not_clear", "Not" }, { "title", "Adjoint" }, { "ylabel", "Time" } });
    outAdjoint << vars.performanceAdjointTime_;

    cl::AdjointTestOutput outSize("AdjointMatrices//Eigenvectors", { { "filename", "TapeSize" }, { "not_clear", "Not" }, { "title", "Size of Tape" }, { "ylabel", "Memory (bytes)" } });
    outSize << vars.tapeSize_;
    return result;
#endif
    return true;
}

bool AdjointMatricesTest::testOrthogonalProjection()
{
    BOOST_TEST_MESSAGE("Testing orthogonal projections...");
    bool result = true;
#ifdef CL_TAPE_CPPAD
    //Size dimension = 8;
    Size numberVectors = 8;
    Real multiplier = 100;
    Real tolerance = 1e-4;
    unsigned long seed = 1;
    boost::timer timer;
    CommonVars vars;
    for (Size dimension = 1; dimension < 200; dimension += 10)
    {
#ifndef CL_GRAPH_GEN
        dimension = 200;
#endif
        std::vector<cl::TapeDouble> A(numberVectors* dimension);
        MersenneTwisterUniformRng rng(seed);

        for (Size i = 0; i < numberVectors; ++i)
        for (Size j = 0; j < dimension; ++j)
            A[i*dimension + j] = rng.next().value;
        timer.restart();
        Independent(A);
        Matrix AM(numberVectors, dimension), AMP(numberVectors, dimension);
        AM = vector2matrix(A, numberVectors, dimension);
        OrthogonalProjections projector(AM, multiplier, tolerance);

        for (Size i = 0; i < numberVectors; ++i)
        for (Size j = 0; j < dimension; ++j)
            AMP[i][j] = projector.GetVector(i)[j];
        std::vector<cl::TapeDouble> Y(1);
        Y[0] = norm(AMP);
        cl::TapeFunction<double> f(A, Y);
        vars.timeTapeRecording_ = timer.elapsed();

        //Store size of tape
        vars.tapeSize_.push_back(TapeSize { dimension, f.Memory() });

        std::vector<double> sf_Reverse;

        //Start differentiation in Reverse mode
        vars.timeAdjoint_ = gradReverse(f, sf_Reverse, false, false);

        //Finite differences
        double h = 1e-8;
        timer.restart();
        std::vector<Real> sf_Finite(numberVectors* dimension);
        Matrix AMR(numberVectors, dimension), AML(numberVectors, dimension);
        for (Size i = 0; i < numberVectors; ++i)
        {

            for (Size j = 0; j < dimension; ++j)
            {
                AM[i][j] += h;
                OrthogonalProjections projector_h(AM, multiplier, tolerance);
                for (Size is = 0; is < numberVectors; ++is)
                for (Size js = 0; js < dimension; ++js)
                    AMR[is][js] = projector_h.GetVector(is)[js];
                AM[i][j] -= 2*h;
                projector_h = OrthogonalProjections(AM, multiplier, tolerance);
                for (Size is = 0; is < numberVectors; ++is)
                for (Size js = 0; js < dimension; ++js)
                    AML[is][js] = projector_h.GetVector(is)[js];
                sf_Finite[i*dimension + j] = (norm(AMR) - norm(AML)) / (2*h);
                AM[i][j] += h;
            }
        }
        vars.timeAnalytical_ = timer.elapsed();

        vars.performanceTime_.push_back(PerformanceTime { vars.timeTapeRecording_, vars.timeAdjoint_, vars.timeAnalytical_, dimension });
        vars.performanceAdjointTime_.push_back(AdjointTime { vars.timeAdjoint_, dimension });
        //Check
        for (Size i = 0; i < numberVectors* dimension; ++i)
        {
            if (std::fabs(sf_Reverse[i] - sf_Finite[i])>tolerance)
            {
                BOOST_FAIL("Derivatives of norm is not satisfied: "
                           << "\nusing Quantlib function  = " << sf_Reverse[i]
                           << "\nusing finite differences = " << sf_Finite[i]);
                result = false;
            }
        }
    }
    cl::AdjointTestOutput outPerform("AdjointMatrices//OrthogonalProjection", { { "filename", "AdjointPerformance" }, { "not_clear", "Not" }, { "title", "Performance for Adjoint OrthogonalProjection" }, { "ylabel", "Time" } });
    outPerform << vars.performanceTime_;


    cl::AdjointTestOutput outAdjoint("AdjointMatrices//OrthogonalProjection", { { "filename", "Adjoint" }, { "not_clear", "Not" }, { "title", "Adjoint" }, { "ylabel", "Time" } });
    outAdjoint << vars.performanceAdjointTime_;

    cl::AdjointTestOutput outSize("AdjointMatrices//OrthogonalProjection", { { "filename", "TapeSize" }, { "not_clear", "Not" }, { "title", "Size of Tape" }, { "ylabel", "Memory (bytes)" } });
    outSize << vars.tapeSize_;
#endif
    return result;
}




test_suite* AdjointMatricesTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Adjoint with Matrices tests");
    suite->add(QUANTLIB_TEST_CASE(&AdjointMatricesTest::testInverse));
    suite->add(QUANTLIB_TEST_CASE(&AdjointMatricesTest::testDeterminant));
    suite->add(QUANTLIB_TEST_CASE(&AdjointMatricesTest::testSqrt));
    suite->add(QUANTLIB_TEST_CASE(&AdjointMatricesTest::testPolynom));
    suite->add(QUANTLIB_TEST_CASE(&AdjointMatricesTest::testEigenvectors));
    suite->add(QUANTLIB_TEST_CASE(&AdjointMatricesTest::testOrthogonalProjection));
    return suite;
}


#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_matrices)

BOOST_AUTO_TEST_CASE(testMatricesInverse)
{
    BOOST_CHECK(AdjointMatricesTest::testInverse());
}
BOOST_AUTO_TEST_CASE(testMatricesDeterminant)
{
    BOOST_CHECK(AdjointMatricesTest::testDeterminant());
}
BOOST_AUTO_TEST_CASE(testMatricesSqrt)
{
    BOOST_CHECK(AdjointMatricesTest::testSqrt());
}
BOOST_AUTO_TEST_CASE(testMatricesPolynom)
{
    BOOST_CHECK(AdjointMatricesTest::testPolynom());
}
BOOST_AUTO_TEST_CASE(testMatricesEigenvectors)
{
    BOOST_CHECK(AdjointMatricesTest::testEigenvectors());
}
BOOST_AUTO_TEST_CASE(testMatricesOrthogonalProjection)
{
    BOOST_CHECK(AdjointMatricesTest::testOrthogonalProjection());
}

BOOST_AUTO_TEST_SUITE_END()

#endif




