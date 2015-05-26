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

#include <ql/quantlib.hpp>
#include <boost/timer.hpp>
#include "utilities.hpp"
#include "adjointcurvefittingtest.hpp"
#include "adjointtestutilities.hpp"

using namespace QuantLib;
using namespace boost::unit_test_framework;

namespace
{
	class BraceVolatilityFunction
	{
		Real a_, b_, c_, d_;
	public:
		inline BraceVolatilityFunction(Real a, Real b, Real c, Real d)
			: a_(a), b_(b), c_(c), d_(d)
		{
		}
		inline ~BraceVolatilityFunction() {}

		inline Real operator() (Real t)
		{
			return (a_ + b_*t)*std::exp(-c_*t) + d_;
		}
	};


	struct LSCurvesValue
	{
		static std::deque<std::string > get_columns()
		{
			static std::deque<std::string > columns =
			{
				"x", "Observed", "Estimated", "True"
			};

			return columns;
		}

		template <typename stream_type>
		friend inline stream_type&
			operator << (stream_type& stm, LSCurvesValue& v)
		{
				stm << v.xvalue_
					<< ";" << v.observed_
					<< ";" << v.estimated_
					<< ";" << v.true_ << std::endl;

				return stm;
			}

		Real observed_;
		Real estimated_;
		Real true_;
		Real xvalue_;
	};

	/*!
	Curve Fitting Problem
	*/
	class CurveFittingProblem : public LeastSquareProblem
	{
		const Array &ttm_;
		const Array &target_;

	public:
		/*!
		Default constructor : set time to maturity vector
		and target value
		*/

		static Size counter;

		CurveFittingProblem(const Array &ttm, const Array &target) : ttm_(ttm), target_(target)
		{
		}

		//! Destructor
		virtual ~CurveFittingProblem() {}

		//! Size of the least square problem
		virtual size_t size()
		{
			return ttm_.size();
		}

		//! return function and target values
		virtual void targetAndValue(const Array& coefficients, Array& target, Array& fct2fit)
		{
			BraceVolatilityFunction bvf(coefficients[0], coefficients[1], coefficients[2], coefficients[3]);

			target = target_;// target values
			for (Size i = 0; i < ttm_.size(); ++i)
				fct2fit[i] = bvf(ttm_[i]);
		}

		//! return function, target and first derivatives values
		virtual void targetValueAndGradient(const Array& coefficients, Matrix& grad_fct2fit, Array& target, Array& fct2fit)
		{
			Size n = coefficients.size();
			std::vector<cl::TapeDouble> coef(n);
			std::vector<double> coefD(n);
			Size m = ttm_.size();
			for (size_t i = 0; i < n; i++)
			{
				coef[i] = coefficients[i].value();
				coefD[i] = (double)coefficients[i];
			}

			cl::Independent(coef);

			BraceVolatilityFunction bvf(coef[0], coef[1], coef[2], coef[3]);

			std::vector<cl::TapeDouble> y(m);
			for (Size i = 0; i < m; ++i)
			{
				y[i] = bvf(ttm_[i]).value();
			}
			cl::TapeFunction<double> f(coef, y);
			//Compute derivatives using Jacobian
			std::vector<double> Jacobian(m*n);
			Jacobian = f.Jacobian(coefD);

			target = target_;// target values
			for (Size i = 0; i < ttm_.size(); ++i)
			{
				// function value at current point coefficients
				fct2fit[i] = bvf(ttm_[i]);
			}

			/*
			matrix of first derivatives :
			the derivatives with respect to the parameter a,b,c,d
			are stored by row.
			*/
			counter++;
			for (size_t i = 0; i < m; i++)
			for (size_t j = 0; j < n; j++)
				grad_fct2fit[i][j] = Jacobian[i*n + j];

		}
	};

	Size CurveFittingProblem::counter = 0;
}

/*
We define here an inverse problem to show how to fit
parametric function to data.
*/
bool AdjointCurveFittingTest::testCurveFitting()
{
	BOOST_TEST_MESSAGE("Testing curve fitting with AAD...\n");
	bool result = false;
#ifdef CL_TAPE_CPPAD
	boost::mt19937 rng;
	boost::normal_distribution<> nd(0.0, 0.01);
	boost::variate_generator<boost::mt19937&, boost::normal_distribution<>> noise(rng, nd);

	/*
	Parameter values that produce the volatility hump.
	Consider it as optimal values of the curve fitting
	problem.
	*/
	std::vector<Real> coefficients_ =
	{
		0.147014,
		0.057302,
		0.249964,
		0.148556
	};


	Array coefficients(coefficients_.begin(), coefficients_.end());

	// Define the target volatility function
	BraceVolatilityFunction bvf(coefficients[0], coefficients[1], coefficients[2], coefficients[3]);

	// start date of volatilty
	const Real startDate = 0.0;
	// end date of volatility
	const Real endDate = 20.;
	// period length between values (in year fraction : quarterly)
	const Real period = 0.1;

	// number of period
	size_t periodNumber = (size_t)(endDate / period);

	Array targetValue(periodNumber);
	Array timeToMaturity(periodNumber);



	// Fill target and time to maturity arrays
	for (size_t i = 0; i < periodNumber; ++i)
	{
		const Real t = startDate + i * period;
		timeToMaturity[i] = t;
		targetValue[i] = bvf(t) + noise();
	}

	// Accuracy of the optimization method
	const Real accuracy = 1e-5;// It is the square of the accuracy
	// Maximum number of iterations
	Size maxiter = 10000;

	Array initialValue(4, 0.1);

	// Least square optimizer
	NoConstraint nc;
	NonLinearLeastSquare lsqnonlin(nc, accuracy, maxiter, boost::shared_ptr<OptimizationMethod>(new ConjugateGradient()));

	// Define the least square problem
	CurveFittingProblem cfp(timeToMaturity, targetValue);

	// Set initial values
	lsqnonlin.setInitialValue(initialValue);
	// perform fitting
	Array solution = lsqnonlin.perform(cfp);
	BraceVolatilityFunction bvf_est(solution[0], solution[1], solution[2], solution[3]);

	// Plot stream
	cl::AdjointTestOutput out("AdjointCurveFitting"
							  , { { "filename", "CurveFitting" }
	, { "ylabel", "F(x)" }
	, { "not_clear", "Not" }
	, { "title", "Least-squares curve fitting" }
	, { "cleanlog", "false" }
	, { "xlabel", "x" } });

	cl::AdjointTestOutput outPerform("AdjointCurveFitting/"
									 , { { "filename", "LSCurves" }
	, { "not_clear", "Not" }
	, { "line_box_width", "-5" }
	, { "title", "Swaption NPV differentiation performance with respect to volatility" }
	, { "ylabel", "Time (s)" }
	, { "xlabel", "Number of volatilities" } });

#if defined CL_GRAPH_GEN

	// The observed volatility function 
	std::vector<LSCurvesValue> Volatility;

	LSCurvesValue temp;
	// Fill target and time to maturity arrays
	for (size_t i = 0; i < periodNumber; ++i)
	{
		const Real t = startDate + i * period;
		temp.xvalue_ = t;
		timeToMaturity[i] = t;
		temp.observed_ = targetValue[i];
		temp.estimated_ = bvf_est(t);
		temp.true_ = bvf(t);
		Volatility.push_back(temp);
	}
	std::cout << Volatility.size() << std::endl;
	outPerform << Volatility;
#endif
	

	// check the result with defined tolerance
	const Real tollerance = 1e-2;
	result = true;
	for (Size i = 0; i < solution.size(); i++)
	{
		if (solution[i] - coefficients[i] > tollerance)
		{
			result = false;
			BOOST_ERROR("Optimal and solution values [" << i << "] mismatch."
						<< "\nOptimal values  : " << coefficients[i]
						<< "\nSolution values : " << solution[i]
						<< "\nError  : " << (solution[i] - coefficients[i])
						<< "\nTolerance:  " << tollerance);
		}
	}

#endif
	return result;
}

test_suite* AdjointCurveFittingTest::suite()
{
	test_suite* suite = BOOST_TEST_SUITE("AAD curve fitting test");
	suite->add(QUANTLIB_TEST_CASE(&AdjointCurveFittingTest::testCurveFitting));
	return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_curve_fitting)

BOOST_AUTO_TEST_CASE(testCurveFitting)
{
	BOOST_CHECK(AdjointCurveFittingTest::testCurveFitting());
}

BOOST_AUTO_TEST_SUITE_END()

#endif
