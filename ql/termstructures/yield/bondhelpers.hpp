/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2008 Ferdinando Ametrano
 Copyright (C) 2005 Toyin Akin
 Copyright (C) 2007 StatPro Italia srl

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

/*! \file bondhelpers.hpp
    \brief bond rate helpers
*/

#ifndef quantlib_bond_helpers_hpp
#define quantlib_bond_helpers_hpp

#include <ql/termstructures/yield/ratehelpers.hpp>
#include <ql/instruments/bonds/fixedratebond.hpp>

namespace QuantLib {

    //! fixed-coupon bond helper
    /*! \warning This class assumes that the reference date
                 does not change between calls of setTermStructure().
    */
    class BondHelper : public RelativeDateRateHelper {
      public:
        /*! \warning Setting a pricing engine to the passed bond from
                     external code will cause the bootstrap to fail or
                     to give wrong results. It is advised to discard
                     the bond after creating the helper, so that the
                     helper has sole ownership of it.
        */
        BondHelper(const Handle<Quote>& cleanPrice,
                   const boost::shared_ptr<Bond>& bond);
        //! \name BootstrapHelper interface
        //@{
        Real impliedQuote() const;
        void setTermStructure(YieldTermStructure*);
        //@}
        //! \name additional inspectors
        //@{
        boost::shared_ptr<Bond> bond() const;
        //@}
        //! \name Visitability
        //@{
        void accept(AcyclicVisitor&);
        //@}
      protected:
        boost::shared_ptr<Bond> bond_;
        RelinkableHandle<YieldTermStructure> termStructureHandle_;
      private:
        void initializeDates();
    };

    class FixedRateBondHelper : public BondHelper {
      public:
        FixedRateBondHelper(const Handle<Quote>& cleanPrice,
                            Natural settlementDays,
                            Real faceAmount,
                            const Schedule& schedule,
                            const std::vector<Rate>& coupons,
                            const DayCounter& dayCounter,
                            BusinessDayConvention paymentConv = Following,
                            Real redemption = 100.0,
                            const Date& issueDate = Date());
        //! \name additional inspectors
        //@{
        boost::shared_ptr<FixedRateBond> fixedRateBond() const;
        //@}
        //! \name Visitability
        //@{
        void accept(AcyclicVisitor&);
        //@}
      protected:
        boost::shared_ptr<FixedRateBond> fixedRateBond_;
    };

}

#endif