/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003, 2006 Ferdinando Ametrano
 Copyright (C) 2006 Warren Chou
 Copyright (C) 2006, 2008 StatPro Italia srl
 Copyright (C) 2006 Chiara Fornarola

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

/*! \file payoffs.hpp
    \brief Payoffs for various options
    */

#ifndef quantlib_payoffs_t_hpp
#define quantlib_payoffs_t_hpp

#include <ql/option.hpp>
#include <cl/tape/impl/ql/payoff_t.hpp>

namespace QuantLib {

    //! Intermediate class for put/call payoffs
    template <class T>
    class TypePayoff_T : public Payoff_T<T>
    {
    public:
        Option::Type optionType() const { return type_; };
    protected:
        TypePayoff_T(Option::Type type) : type_(type) {}
        Option::Type type_;
    };

    //! Intermediate class for payoffs based on a fixed strike
    template <class T>
    class StrikedTypePayoff_T : public TypePayoff_T<T>
    {
    public:
        T strike() const { return strike_; };
    protected:
        StrikedTypePayoff_T(Option::Type type, T strike)
            : TypePayoff_T(type), strike_(strike) {}
        T strike_;
    };

    //! Plain-vanilla payoff
    template <class T>
    class PlainVanillaPayoff_T : public StrikedTypePayoff_T<T> {
    public:
        PlainVanillaPayoff_T(Option::Type type, T strike)
            : StrikedTypePayoff_T<T>(type, strike) {}

        std::string name() const { return "Vanilla"; }
        T operator()(const T& price) const;
        virtual void accept(AcyclicVisitor&);
    };

    template <class T>
    T PlainVanillaPayoff_T<T>::operator()(const T& price) const
    {
        switch (type_)
        {
        case Option::Call:
            return std::max(price - strike_, 0.0);
        case Option::Put:
            return std::max(strike_ - price, 0.0);
        default:
            QL_FAIL("unknown/illegal option type");
        }
    }

    template <class T>
    void PlainVanillaPayoff_T<T>::accept(AcyclicVisitor& v)
    {
        Visitor<PlainVanillaPayoff_T>* v1 =
            dynamic_cast<Visitor<PlainVanillaPayoff_T>*>(&v);
        if (v1 != 0)
            v1->visit(*this);
        else
            Payoff_T::accept(v);
    }
}

#endif
