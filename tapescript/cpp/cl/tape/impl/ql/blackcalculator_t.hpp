/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003, 2004, 2005, 2006 Ferdinando Ametrano
 Copyright (C) 2006 StatPro Italia srl
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

/*! \file blackcalculator_t.hpp
    \brief Black-formula calculator class
*/

#ifndef quantlib_blackcalculator_t_hpp
#define quantlib_blackcalculator_t_hpp

#include <cl/tape/impl/ql/payoffs_t.hpp>
#include <cl/tape/impl/ql/normaldistribution_t.hpp>
#include <ql/math/comparison.hpp>

namespace QuantLib
{
    //! Black 1976 calculator class
    /*! \bug When the variance is null, division by zero occur during
             the calculation of delta, delta forward, gamma, gamma
             forward, rho, dividend rho, vega, and strike sensitivity.
    */
    template <class T = Real>
    class BlackCalculator_T
    {
    private:
        class Calculator;
        friend class Calculator;
    public:
        typedef T value_type;

        BlackCalculator_T(const boost::shared_ptr<StrikedTypePayoff_T<T>>& payoff,
            T const& forward,
            T const& stdDev,
            T const& discount = 1.0);
        BlackCalculator_T(Option::Type optionType,
            T const& strike,
            T const& forward,
            T const& stdDev,
            T const& discount = 1.0);
        virtual ~BlackCalculator_T() {}

        T value() const;

        /*! Sensitivity to change in the underlying forward price. */
        T deltaForward() const;
        /*! Sensitivity to change in the underlying spot price. */
        virtual T delta(T const& spot) const;

        /*! Sensitivity in percent to a percent change in the
            underlying forward price. */
        T elasticityForward() const;
        /*! Sensitivity in percent to a percent change in the
            underlying spot price. */
        virtual T elasticity(T const& spot) const;

        /*! Second order derivative with respect to change in the
            underlying forward price. */
        T gammaForward() const;
        /*! Second order derivative with respect to change in the
            underlying spot price. */
        virtual T gamma(T const& spot) const;

        /*! Sensitivity to time to maturity. */
        virtual T theta(T const& spot,
            T const& maturity) const;
        /*! Sensitivity to time to maturity per day,
            assuming 365 day per year. */
        virtual T thetaPerDay(T const& spot,
            T const& maturity) const;

        /*! Sensitivity to volatility. */
        T vega(T const& maturity) const;

        /*! Sensitivity to discounting rate. */
        T rho(T const& maturity) const;

        /*! Sensitivity to dividend/growth rate. */
        T dividendRho(T const& maturity) const;

        /*! Probability of being in the money in the bond martingale
            measure, i.e. N(d2).
            It is a risk-neutral probability, not the real world one.
            */
        T itmCashProbability() const;

        /*! Probability of being in the money in the asset martingale
            measure, i.e. N(d1).
            It is a risk-neutral probability, not the real world one.
        */
        T itmAssetProbability() const;

        /*! Sensitivity to strike. */
        T strikeSensitivity() const;

        T alpha() const;
        T beta() const;
    protected:
        void initialize(const boost::shared_ptr<StrikedTypePayoff_T<T>>& p);
        T strike_, forward_, stdDev_, discount_, variance_;
        T d1_, d2_;
        T alpha_, beta_, DalphaDd1_, DbetaDd2_;
        T n_d1_, cum_d1_, n_d2_, cum_d2_;
        T x_, DxDs_, DxDstrike_;
    };

    // inline
    template <class T>
    inline T BlackCalculator_T<T>::thetaPerDay(T const& spot,
        T const& maturity) const
    {
        return theta(spot, maturity) / 365.0;
    }

    template <class T>
    inline T BlackCalculator_T<T>::itmCashProbability() const {
        return cum_d2_;
    }

    template <class T>
    inline T BlackCalculator_T<T>::itmAssetProbability() const {
        return cum_d1_;
    }

    template <class T>
    inline T BlackCalculator_T<T>::alpha() const
    {
        return alpha_;
    }

    template <class T>
    inline T BlackCalculator_T<T>::beta() const
    {
        return beta_;
    }

    template <class T>
    class BlackCalculator_T<T>::Calculator
        : public AcyclicVisitor
        //, public Visitor<Payoff_T<T>>
        , public Visitor<PlainVanillaPayoff_T<T>>
        //, public Visitor<CashOrNothingPayoff_T<T>>
        //, public Visitor<AssetOrNothingPayoff_T<T>>
        //, public Visitor<GapPayoff_T<T>>
    {
    private:
        BlackCalculator_T<T>& black_;
    public:
        Calculator(BlackCalculator_T<T>& black) : black_(black) {}
        //void visit(Payoff_T<T>&);
        void visit(PlainVanillaPayoff_T<T>&) {}
        //void visit(CashOrNothingPayoff_T<T>&);
        //void visit(AssetOrNothingPayoff_T<T>&);
        //void visit(GapPayoff_T<T>&);
    };


    template <class T>
    BlackCalculator_T<T>::BlackCalculator_T(const boost::shared_ptr<StrikedTypePayoff_T<T>>& p,
        T const& forward,
        T const& stdDev,
        T const& discount)
        : strike_(p->strike()), forward_(forward), stdDev_(stdDev),
        discount_(discount), variance_(stdDev*stdDev)
    {
        initialize(p);
    }

    template <class T>
    BlackCalculator_T<T>::BlackCalculator_T(Option::Type optionType,
        T const& strike,
        T const& forward,
        T const& stdDev,
        T const& discount)
        : strike_(strike), forward_(forward), stdDev_(stdDev),
        discount_(discount), variance_(stdDev*stdDev)
    {
        initialize(boost::make_shared<PlainVanillaPayoff_T<T>>(optionType, strike));
    }

    template <class T>
    void BlackCalculator_T<T>::initialize(const boost::shared_ptr<StrikedTypePayoff_T<T>>& p)
    {
        QL_REQUIRE(strike_ >= 0.0,
            "strike (" << strike_ << ") must be non-negative");
        QL_REQUIRE(forward_ > 0.0,
            "forward (" << forward_ << ") must be positive");
        //QL_REQUIRE(displacement_>=0.0,
        //           "displacement (" << displacement_ << ") must be non-negative");
        QL_REQUIRE(stdDev_ >= 0.0,
            "stdDev (" << stdDev_ << ") must be non-negative");
        QL_REQUIRE(discount_ > 0.0,
            "discount (" << discount_ << ") must be positive");

        if (stdDev_ >= QL_EPSILON)
        {
            d1_ = std::log(forward_ / strike_) / stdDev_ + 0.5*stdDev_;
            d2_ = d1_ - stdDev_;
            CumulativeNormalDistribution_T<T> f;
            cum_d1_ = f(d1_);
            cum_d2_ = f(d2_);
            n_d1_ = f.derivative(d1_);
            n_d2_ = f.derivative(d2_);
        }
        else
        {
            if (forward_ > strike_)
            {
                d1_ = QL_MAX_REAL;
                d2_ = QL_MAX_REAL;
                cum_d1_ = 1.0;
                cum_d2_ = 1.0;
                n_d1_ = 0.0;
                n_d2_ = 0.0;
            }
            else
            {
                d1_ = QL_MIN_REAL;
                d2_ = QL_MIN_REAL;
                cum_d1_ = 0.0;
                cum_d2_ = 0.0;
                n_d1_ = 0.0;
                n_d2_ = 0.0;
            }
        }

        x_ = strike_;
        DxDstrike_ = 1.0;

        // the following one will probably disappear as soon as
        // super-share will be properly handled
        DxDs_ = 0.0;

        // this part is always executed.
        // in case of plain-vanilla payoffs, it is also the only part
        // which is executed.
        switch (p->optionType())
        {
        case Option::Call:
            alpha_ = cum_d1_;//  N(d1)
            DalphaDd1_ = n_d1_;//  n(d1)
            beta_ = -cum_d2_;// -N(d2)
            DbetaDd2_ = -n_d2_;// -n(d2)
            break;
        case Option::Put:
            alpha_ = -1.0 + cum_d1_;// -N(-d1)
            DalphaDd1_ = n_d1_;//  n( d1)
            beta_ = 1.0 - cum_d2_;//  N(-d2)
            DbetaDd2_ = -n_d2_;// -n( d2)
            break;
        default:
            QL_FAIL("invalid option type");
        }

        // now dispatch on type.

        Calculator calc(*this);
        p->accept(calc);
    }


    template <class T>
    T BlackCalculator_T<T>::value() const {
        T result = discount_ * (forward_ * alpha_ + x_ * beta_);
        return result;
    }

    template <class T>
    T BlackCalculator_T<T>::delta(T const& spot) const {

        QL_REQUIRE(spot > 0.0, "positive spot value required: " <<
            spot << " not allowed");

        T DforwardDs = forward_ / spot;

        T temp = stdDev_*spot;
        T DalphaDs = DalphaDd1_ / temp;
        T DbetaDs = DbetaDd2_ / temp;
        T temp2 = DalphaDs * forward_ + alpha_ * DforwardDs
            + DbetaDs  * x_ + beta_  * DxDs_;

        return discount_ * temp2;
    }

    template <class T>
    T BlackCalculator_T<T>::deltaForward() const {

        T temp = stdDev_*forward_;
        T DalphaDforward = DalphaDd1_ / temp;
        T DbetaDforward = DbetaDd2_ / temp;
        T temp2 = DalphaDforward * forward_ + alpha_
            + DbetaDforward  * x_; // DXDforward = 0.0

        return discount_ * temp2;
    }

    template <class T>
    T BlackCalculator_T<T>::elasticity(T const& spot) const
    {
        T val = value();
        T del = delta(spot);
        if (val > QL_EPSILON)
            return del / val*spot;
        else if (std::fabs(del)<QL_EPSILON)
            return 0.0;
        else if (del>0.0)
            return QL_MAX_REAL;
        else
            return QL_MIN_REAL;
    }

    template <class T>
    T BlackCalculator_T<T>::elasticityForward() const {
        T val = value();
        T del = deltaForward();
        if (val > QL_EPSILON)
            return del / val*forward_;
        else if (std::fabs(del)<QL_EPSILON)
            return 0.0;
        else if (del>0.0)
            return QL_MAX_REAL;
        else
            return QL_MIN_REAL;
    }

    template <class T>
    T BlackCalculator_T<T>::gamma(T const& spot) const {

        QL_REQUIRE(spot > 0.0, "positive spot value required: " <<
            spot << " not allowed");

        T DforwardDs = forward_ / spot;

        T temp = stdDev_*spot;
        T DalphaDs = DalphaDd1_ / temp;
        T DbetaDs = DbetaDd2_ / temp;

        T D2alphaDs2 = -DalphaDs / spot*(1 + d1_ / stdDev_);
        T D2betaDs2 = -DbetaDs / spot*(1 + d2_ / stdDev_);

        T temp2 = D2alphaDs2 * forward_ + 2.0 * DalphaDs * DforwardDs
            + D2betaDs2  * x_ + 2.0 * DbetaDs  * DxDs_;

        return  discount_ * temp2;
    }

    template <class T>
    T BlackCalculator_T<T>::gammaForward() const {

        T temp = stdDev_*forward_;
        T DalphaDforward = DalphaDd1_ / temp;
        T DbetaDforward = DbetaDd2_ / temp;

        T D2alphaDforward2 = -DalphaDforward / forward_*(1 + d1_ / stdDev_);
        T D2betaDforward2 = -DbetaDforward / forward_*(1 + d2_ / stdDev_);

        T temp2 = D2alphaDforward2 * forward_ + 2.0 * DalphaDforward
            + D2betaDforward2  * x_; // DXDforward = 0.0

        return discount_ * temp2;
    }

    template <class T>
    T BlackCalculator_T<T>::theta(T const& spot,
        T const& maturity) const {

        QL_REQUIRE(maturity >= 0.0,
            "maturity (" << maturity << ") must be non-negative");
        //if (close(maturity, 0.0)) return 0.0;
        return -(std::log(discount_)            * value()
            + std::log(forward_ / spot) * spot * delta(spot)
            + 0.5*variance_ * spot  * spot * gamma(spot)) / maturity;
    }

    template <class T>
    T BlackCalculator_T<T>::vega(T const& maturity) const {
        QL_REQUIRE(maturity >= 0.0,
            "negative maturity not allowed");

        T temp = std::log(strike_ / forward_) / variance_;
        // actually DalphaDsigma / SQRT(T)
        T DalphaDsigma = DalphaDd1_*(temp + 0.5);
        T DbetaDsigma = DbetaDd2_ *(temp - 0.5);

        T temp2 = DalphaDsigma * forward_ + DbetaDsigma * x_;

        return discount_ * std::sqrt(maturity) * temp2;

    }

    template <class T>
    T BlackCalculator_T<T>::rho(T const& maturity) const {
        QL_REQUIRE(maturity >= 0.0,
            "negative maturity not allowed");

        // actually DalphaDr / T
        T DalphaDr = DalphaDd1_ / stdDev_;
        T DbetaDr = DbetaDd2_ / stdDev_;
        T temp = DalphaDr * forward_ + alpha_ * forward_ + DbetaDr * x_;

        return maturity * (discount_ * temp - value());
    }

    template <class T>
    T BlackCalculator_T<T>::dividendRho(T const& maturity) const {
        QL_REQUIRE(maturity >= 0.0,
            "negative maturity not allowed");

        // actually DalphaDq / T
        T DalphaDq = -DalphaDd1_ / stdDev_;
        T DbetaDq = -DbetaDd2_ / stdDev_;

        T temp = DalphaDq * forward_ - alpha_ * forward_ + DbetaDq * x_;

        return maturity * discount_ * temp;
    }

    template <class T>
    T BlackCalculator_T<T>::strikeSensitivity() const {

        T temp = stdDev_*strike_;
        T DalphaDstrike = -DalphaDd1_ / temp;
        T DbetaDstrike = -DbetaDd2_ / temp;

        T temp2 =
            DalphaDstrike * forward_ + DbetaDstrike * x_ + beta_ * DxDstrike_;

        return discount_ * temp2;
    }

}

#endif
