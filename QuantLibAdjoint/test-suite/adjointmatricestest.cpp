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


#include "adjointmatricesimpl.hpp"


bool AdjointMatricesTest::testInverse()
{
    BOOST_TEST_MESSAGE("Testing Adjoint differentiation for inverse matrices...");

    Size n = 10;
    InverseTest test(n);

    std::vector<double> doubleX = generateVector<double>(n);
    // Vector representation of matrix.
    std::vector<cl::tape_double> X(doubleX.begin(), doubleX.end());

    // Tape recording.
    cl::Independent(X);

    Matrix A = vector2matrix(X, n, n);
    Matrix invA = inverse(A);
    std::vector<cl::tape_double> Y = matrix2vector(invA);

    // End of tape recording.
    cl::tape_function<double> f(X, Y);

    // Calculating derivatives using adjoint.
    test.adjointResults_ = f.Jacobian(doubleX);

    // Calculate derivatives using analytical formula.
    test.calcAnalytical();

    bool ok = test.checkAdjoint();

    InverseTestData testData;

    return ok && testData.makeOutput();
}


bool AdjointMatricesTest::testDeterminant()
{
    BOOST_TEST_MESSAGE("Testing Adjoint differentiation for matrix determinant...");

    Size n = 10;
    DeterminantTest test(n);

    // X_ is vector representation of input Matrix.
    cl::Independent(test.X_);
    Matrix A = vector2matrix(test.X_, n, n);
    test.det_ = { determinant(A) };
    cl::tape_function<double> f(test.X_, test.det_);

    // Forward mode derivatives calculation.
    test.forwardResults_.resize(n * n);
    std::vector<double> dX(n * n, 0);
    for (size_t i = 0; i < n * n; i++)
    {
        dX[i] = 1;
        test.forwardResults_[i] = f.Forward(1, dX)[0];
        dX[i] = 0;
    }

    // Reverse mode derivatives calulation.
    test.reverseResults_ = f.Reverse(1, std::vector<double>(1, 1));

    // Calculation of analytical results.
    test.calcAnalytical();

    bool ok = test.check();

    DeterminantTestData testData;

    return ok && testData.makeOutput();
}


bool AdjointMatricesTest::testSqrt()
{
    BOOST_TEST_MESSAGE("Testing matricial pseudo square root...");

    Size n = 4;
    SqrtTest test(n);

    // X_ is vector representation of input Matrix.
    // Tape recording.
    cl::Independent(test.X_);

    Matrix A = vector2matrix(test.X_, n, n);
    Matrix sqrtA = pseudoSqrt(A, SalvagingAlgorithm::None);
    test.Y_[0] = determinant(sqrtA);

    // End of tape recording.
    cl::tape_function<double> f(test.X_, test.Y_);

    // Forward mode derivatives calculation.
    test.forwardResults_.resize(n * n);
    std::vector<double> dX(n * n, 0);
    for (size_t i = 0; i < n * n; i++)
    {
        dX[i] = 1;
        test.forwardResults_[i] = f.Forward(1, dX)[0];
        dX[i] = 0;
    }

    // Reverse mode derivatives calulation.
    test.reverseResults_ = f.Reverse(1, std::vector<double>(1, 1));

    // Calculation of analytical results.
    test.calcAnalytical();

    return test.check();
}


bool AdjointMatricesTest::testPolynom()
{

    BOOST_TEST_MESSAGE("Testing matricial polynom differentiation...");
    Size n = 10;
    PolynomTest test(n);
    std::vector<Real> P = { 5.4, 1.2, 0.3, 0.1 };

    // Tape recording.
    cl::Independent(test.X_);
    Matrix A = test.X_[0] * vector2matrix(generateVector<double>(n), n, n);
    Matrix PA = polynom(A, P);
    test.Y_ = matrix2vector(PA);

    // End of tape recording.
    cl::tape_function<double> f(test.X_, test.Y_);

    // Forward mode derivatives calculation.
    test.forwardResults_ = f.Forward(1, std::vector<double>(1, 1));

    // Reverse mode derivatives calulation.
    test.reverseResults_.resize(n * n);
    std::vector<double> dw(n * n, 0);
    for (size_t i = 0; i < n * n; i++)
    {
        dw[i] = 1;
        test.reverseResults_[i] = f.Reverse(1, dw)[0];
        dw[i] = 0;
    }
    // Calculation of analytical results.
    test.calcAnalytical();

    return test.check();
}


bool AdjointMatricesTest::testEigenvectors()
{

   BOOST_TEST_MESSAGE("Testing eigenvalues and eigenvectors calculation with AD...");

   Size n = 10;
   EigenvectorsTest test(n);

   // Vector representation of matrix.
   std::vector<double> doubleX = generateSymVector(n);
   std::vector<cl::tape_double> X(doubleX.begin(), doubleX.end());

   // Tape recording.
   cl::Independent(X);

   Matrix A = vector2matrix(X, n, n);
   std::vector<cl::tape_double> Y = calculateEigen(n, A);

   // End of tape recording.
   cl::tape_function<double> f(X, Y);

   // Calculating derivatives with adjoint.
   test.adjointResults_ = f.Jacobian(doubleX);

   // Calculate derivatives from analytical formula.
   test.calcAnalytical();

   EigenvectorsTestData testData;

   return test.checkAdjoint() && testData.makeOutput();
}


bool AdjointMatricesTest::testOrthogonalProjection()
{
    BOOST_TEST_MESSAGE("Testing orthogonal projections...");
    Size n = 10;
    OrthogonalProjectionTest test(n);

    // X_ is vector representation of input Matrix.
    // Tape recording.
    cl::Independent(test.X_);

    Matrix A = vector2matrix(test.X_, test.numberVectors_, n);
    Matrix projA_ = test.calcProjection(A);
    test.Y_[0] = norm(projA_);

    // End of tape recording.
    cl::tape_function<double> f(test.X_, test.Y_);

    // Forward mode derivatives calculation.
    test.forwardResults_.resize(test.numberVectors_*n);
    std::vector<double> dX(test.numberVectors_*n, 0);
    for (size_t i = 0; i < test.numberVectors_*n; i++)
    {
        dX[i] = 1;
        test.forwardResults_[i] = f.Forward(1, dX)[0];
        dX[i] = 0;
    }

    // Reverse mode derivatives calulation.
    test.reverseResults_ = f.Reverse(1, std::vector<double>(1, 1));

    // Calculation of analytical results.
    test.calcAnalytical();
    OrthogonalProjectionTestData testData;

    return test.check() && testData.makeOutput();
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
