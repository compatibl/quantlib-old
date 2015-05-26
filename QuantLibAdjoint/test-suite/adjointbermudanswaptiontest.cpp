/* based on test-suite\bermudanswaption.cpp */

/*
Copyright (C) 2005, 2007 StatPro Italia srl
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

#include "adjointbermudanswaptiontest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include <ql/instruments/swaption.hpp>
#include <ql/pricingengines/swaption/treeswaptionengine.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>
#include <ql/pricingengines/swaption/fdhullwhiteswaptionengine.hpp>
#include <ql/models/shortrate/onefactormodels/hullwhite.hpp>
#include <ql/cashflows/coupon.hpp>
#include <ql/time/daycounters/thirty360.hpp>
#include <ql/indexes/ibor/euribor.hpp>
#include <ql/time/schedule.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;

namespace {

    enum {
#if defined CL_GRAPH_GEN
        // Number of points for dependency graphics.
        pointNo = 100,
        // Number of points for performance graphic.
        iterNo = 50,
        // Step for portfolio size for performance testing .
        step = 1,
#else
        // Number of points for dependency graphics.
        pointNo = 1,
        // Number of points for performance graphic.
        iterNo = 1,
        // Step for portfolio size for performance testing .
        step = 3,
#endif
        // Defines performance accuracy. Its value is a minimum number
        // of calling of O(1) complexity methods per one performance test.
        iterNumFactor = 1,
    };

    enum EngineType { tree, fdHullWhite };

    struct CommonVars
    {
        CommonVars()
            : backup_()
            , startYears_(1)
            , length_(5)
            , type_(VanillaSwap::Payer)
            , nominal_(1000.0)
            , fixedConvention_(Unadjusted)
            , floatingConvention_(ModifiedFollowing)
            , fixedFrequency_(Annual)
            , floatingFrequency_(Semiannual)
            , fixedDayCount_(Thirty360())
            , termStructure_()
            , index_(boost::shared_ptr<IborIndex>(new Euribor6M(termStructure_)))
            , settlementDays_(2)
            , calendar_(index_->fixingCalendar())
            , today_(calendar_.adjust(Date::todaysDate()))
            , settlement_(calendar_.advance(today_, settlementDays_, Days))
        { }

        // Provides swaps that are different only by their fixed rate strike.
        boost::shared_ptr<VanillaSwap> makeSwap(Rate fixedRate)
        {
            Date start = calendar_.advance(settlement_, startYears_, Years);
            Date maturity = calendar_.advance(start, length_, Years);
            Schedule fixedSchedule(start, maturity,
                Period(fixedFrequency_),
                calendar_,
                fixedConvention_,
                fixedConvention_,
                DateGeneration::Forward, false);
            Schedule floatSchedule(start, maturity,
                Period(floatingFrequency_),
                calendar_,
                floatingConvention_,
                floatingConvention_,
                DateGeneration::Forward, false);
            boost::shared_ptr<VanillaSwap> swap(
                new VanillaSwap(type_, nominal_,
                fixedSchedule, fixedRate, fixedDayCount_,
                floatSchedule, index_, 0.0,
                index_->dayCounter()));
            swap->setPricingEngine(boost::shared_ptr<PricingEngine>(
                new DiscountingSwapEngine(termStructure_)));
            return swap;
        }

        // cleanup
        SavedSettings backup_;

        // underlying swap parameters
        Integer startYears_;
        Integer length_;
        VanillaSwap::Type type_;
        Real nominal_;
        BusinessDayConvention fixedConvention_;
        BusinessDayConvention floatingConvention_;
        Frequency fixedFrequency_;
        Frequency floatingFrequency_;
        DayCounter fixedDayCount_;
        RelinkableHandle<YieldTermStructure> termStructure_;
        boost::shared_ptr<IborIndex> index_;
        Natural settlementDays_;

        // global data
        Calendar calendar_;
        Date today_;
        Date settlement_;
    };

    std::string toString(EngineType engine)
    {
        switch (engine)
        {
        case tree:
            return "TreeSwaptionEngine";
        case fdHullWhite:
            return "FdHullWhiteSwaptionEngine";
        default:
            return "";
        }
    }

    struct SwaptionsData : public CommonVars
    {
        SwaptionsData(EngineType engine)
        {
            today_ = Date(15, February, 2002);
            Settings::instance().evaluationDate() = today_;
            settlement_ = Date(19, February, 2002);

            // flat yield term structure impling 1x5 swap at 5%
            termStructure_.linkTo(flatRate(settlement_,
                0.04875825,
                Actual365Fixed()));
            atmRate_ = makeSwap(0.0)->fairRate();
            Real a = 0.048696;
            Real sigma = 0.0058904;
            boost::shared_ptr<VanillaSwap> atmSwap = makeSwap(atmRate_);
            const Leg& leg = atmSwap->fixedLeg();
            model_ = boost::shared_ptr<HullWhite>(new HullWhite(termStructure_, a, sigma));
            std::vector<Date> exerciseDates(leg.size());
            for (Size i = 0; i < leg.size(); i++)
            {
                boost::shared_ptr<Coupon> coupon =
                    boost::dynamic_pointer_cast<Coupon>(leg[i]);
                exerciseDates[i] = coupon->accrualStartDate();
            }
            exercise_ = boost::shared_ptr<Exercise>(new BermudanExercise(exerciseDates));
            switch (engine)
            {
            case tree:
                engine_ = boost::shared_ptr<PricingEngine>(new TreeSwaptionEngine(model_, 200));
                break;
            case fdHullWhite:
                engine_ = boost::shared_ptr<PricingEngine>(new FdHullWhiteSwaptionEngine(model_, 100, 100, 0, 1e-7));
                break;
            }
        }

        // Returns net present value (NPV) of swaption with given fixed rate strike and other
        // parameters determined from class instance.
        Real swaptionNVP(Rate fixedRate)
        {
            boost::shared_ptr<VanillaSwap> swap(makeSwap(fixedRate));
            Swaption swaption(swap, exercise_);
            swaption.setPricingEngine(engine_);
            return swaption.NPV();
        }

        Rate atmRate_;
        boost::shared_ptr<HullWhite> model_;
        boost::shared_ptr<Exercise> exercise_;
        boost::shared_ptr<PricingEngine> engine_;
    };

    // Class that provides derivative calculation and testing methods. 
    class BSTestHelper
    {
    public:
        BSTestHelper(Size size, EngineType engine)
            : size_(size)
            , data_(engine)
            , strikes_(size)
            , f_(0)
            , prices_()
            , totalPrice_()
            , reverseResults_()
            , modelResults_()
        {
            if (size_ != 1)
            {
                for (Size i = 0; i < size_; i++)
                {
                    strikes_[i] = (0.4 + (1.2 * i) / (size_ - 1)) * data_.atmRate_;
                }
            }
            else
            {
                strikes_[0] = data_.atmRate_;
            }
        }

        // Tests adjoint derivative calculation in reverse mode only.
        bool testReverse()
        {
            recordTape();
            calcReverse();
            calcModel();
            return checkReverse();
        }

        // Records performance for reverse mode and for analytical model (to calculate derivatives).
        // Return true if result check passed.
        bool recordPerformance(PerformanceTime& perfTime)
        {
            perfTime.indepVarNumber_ = indepVarNumber();
            perfTime.timeTapeRecording_ = testTapePerformance();
            perfTime.timeAdjoint_ = testReversePerformance();
            perfTime.timeAnalytical_ = testModelPerformance();
            return checkReverse();
        }

        // Returns performance of tape recording.
        inline double testTapePerformance(Size testNo = 0)
        {
            if (!testNo)
            {
                testNo = minPerfIteration() / indepVarNumber() + 1;
            }
            boost::timer timer;
            for (Size i = 0; i < testNo; i++)
            {
                recordTape();
            }
            return timer.elapsed() / testNo;
        }

        // Returns performance of reverse mode derivatives calculations.
        inline double testReversePerformance(Size testNo = 0)
        {
            if (!testNo)
            {
                testNo = minPerfIteration() / indepVarNumber() + 1;
            }
            if (f_.get() == 0)
            {
                recordTape();
            }
            boost::timer timer;
            for (Size i = 0; i < testNo; i++)
            {
                calcReverse();
            }
            return timer.elapsed() / testNo;
        }

        // Returns performance of derivatives calculations with finite difference.
        inline double testModelPerformance(Size testNo = 0)
        {
            if (!testNo)
            {
                testNo = minPerfIteration() / indepVarNumber() + 1;
            }
            boost::timer timer;
            for (Size i = 0; i < testNo; i++)
            {
                calcModel();
            }
            return timer.elapsed() / testNo;
        }

        // Calculate derivatives with reverse mode only.
        inline void calcWithReverse()
        {
            recordTape();
            calcReverse();
        }

        // Returns derivatives calculated in reverse mode.
        inline std::vector<double> const& reverse()
        {
            return reverseResults_;
        }

        // Returns prices of every swaption in portfolio.
        inline std::vector<Real> const& prices()
        {
            return prices_;
        }

        // Returns fixed rate strikes.
        inline std::vector<Real> const& strikes()
        {
            return strikes_;
        }

        // Returns tape memory size.
        Size memory()
        {
            if (f_.get() == 0)
            {
                return 0;
            }
            return f_->Memory();
        }

        Size indepVarNumber() { return size_; }

        Size minPerfIteration() { return iterNumFactor; }

    protected:
        // Calculates derivatives in reverse mode.
        void calcReverse()
        {
            reverseResults_ = f_->Reverse(1, std::vector<double>(1, 1));
        }

        // Calculates price and records tape.
        void recordTape()
        {
            Independent(strikes_);
            calcPrice();
            f_.reset(new cl::TapeFunction<double>(strikes_, totalPrice_));
        }

        // Calculates price of portfolio and each swaption.
        inline void calcPrice()
        {
            totalPrice_.resize(1, 0);
            prices_.resize(size_);
            for (Size i = 0; i < size_; i++)
            {
                prices_[i] = data_.swaptionNVP(strikes_[i]);
                totalPrice_.front() += prices_[i];
            }
        }

    private:
        // Checks reverse mode and finite difference result consistency.
        bool checkReverse(double relativeTol = 1e-2, double absTol = 1e-10)
        {
            bool result = true;
            for (Size i = 0; i < size_; i++)
            {
                Real maxabs = std::max(std::abs(modelResults_[i]), std::abs(reverseResults_[i]));
                Real tol = std::max(relativeTol * maxabs, absTol);
                if (std::abs(reverseResults_[i] - modelResults_[i]) > tol)
                {
                    result = false;
                    BOOST_ERROR("\nFailed to reproduce expected derivative[" << i << "] at reverse mode."
                        << "\n    calculated: " << reverseResults_[i]
                        << "\n    expected:   " << modelResults_[i]
                        << "\n    tolerance:  " << tol);
                }
            }
            return result;
        }

        // Calculates derivatives using finite difference.
        inline void calcModel()
        {
            Real h = 1e-6 * data_.atmRate_;  // shift for finite diff. method
            modelResults_.resize(size_);
            for (Size i = 0; i < size_; i++)
            {
                modelResults_[i] = (data_.swaptionNVP(strikes_[i] + h) - data_.swaptionNVP(strikes_[i] - h)) / (2 * h);;
            }
        }

    protected:
        Size size_;
        SwaptionsData data_;
        std::vector<cl::TapeDouble> strikes_;
        std::auto_ptr<cl::TapeFunction<double> > f_;
        std::vector<Real> prices_;
        std::vector<cl::TapeDouble> totalPrice_;
        std::vector<double> reverseResults_;
        std::vector<Real> modelResults_;
    };

    // Makes graphics for performance time and tape size.
    // Returns true if calculated results checking passed.
    inline bool recordPerformance(EngineType engine, Size iterNo = 1)
    {
        bool result = true;
        std::vector<PerformanceTime> perfTimes(iterNo);
        std::vector<TapeSize> tapeMemory(iterNo);
        auto p = perfTimes.begin();
        auto m = tapeMemory.begin();
        for (Size i = 1; i <= iterNo; i++)
        {
            BSTestHelper helper(i * step, engine);
            result &= helper.recordPerformance(*p++);
            *m++ = { helper.indepVarNumber(), helper.memory() };
        }
        cl::AdjointTestOutput outPerform("AdjointBermudanSwaption//" + toString(engine), {
            { "title", "Bermudan swaption NPV differentiation performance with respect to strike rate" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of strike rates" }
            , { "line_box_width", "-5" }
        });
        outPerform << perfTimes;
        cl::AdjointTestOutput outSize("AdjointBermudanSwaption//" + toString(engine), {
            { "title", "Tape size dependence on number of strike rates" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
        });
        outSize << tapeMemory;
        return result;
    }

    // Struct for graphics recording.
    struct StrikeSens
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "Strike rate", ""
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, StrikeSens& v)
        {
                stm << v.strike_
                    << ";" << v.sensitivity_
                << std::endl;
            return stm;
        }

        Real strike_;
        Real sensitivity_;
    };

    // Makes graphics for strike sensitivity dependence.
    void recordDependence(EngineType engine)
    {
        std::vector<StrikeSens> outData(pointNo);
        BSTestHelper helper(pointNo, engine);
        helper.calcWithReverse();
        for (Size i = 0; i < pointNo; i++)
        {
            outData[i] = { helper.strikes()[i], helper.reverse()[i] };
        }

        cl::AdjointTestOutput out("AdjointBermudanSwaption//" + toString(engine) + "//output", {
            { "filename", "StrikeDependence" }
            , { "not_clear", "Not" }
            , { "title", "Strike rate sensitivity dependence on strike rate" }
            , { "ylabel", "Strike rate sensitivity" }
            , { "xlabel", "Strike rate" }
        });
        out << outData;
    }
}


bool AdjointBermudanSwaptionTest::testTreeSwaptionEngine()
{

    BOOST_TEST_MESSAGE("Testing fixed rate strike sensitivity of Bermudan swaption (with TreeSwaptionEngine) ...");

#ifdef CL_TAPE_CPPAD

    recordDependence(tree);
    return recordPerformance(tree);

#endif
    return true;
}


bool AdjointBermudanSwaptionTest::testFdHullWhiteSwaptionEngine()
{

    BOOST_TEST_MESSAGE("Testing fixed rate strike sensitivity of Bermudan swaption (with FdHullWhiteSwaptionEngine) ...");

#ifdef CL_TAPE_CPPAD

    recordDependence(fdHullWhite);
    return recordPerformance(fdHullWhite);

#endif
    return true;
}


test_suite* AdjointBermudanSwaptionTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Bermudan swaption tests");
    suite->add(QUANTLIB_TEST_CASE(&AdjointBermudanSwaptionTest::testTreeSwaptionEngine));
    suite->add(QUANTLIB_TEST_CASE(&AdjointBermudanSwaptionTest::testFdHullWhiteSwaptionEngine));
    return suite;
}


#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_bermudan_swaption)

BOOST_AUTO_TEST_CASE(testBermudanSwaptionTreeSwaptionEngine)
{
    BOOST_CHECK(AdjointBermudanSwaptionTest::testTreeSwaptionEngine());
}
BOOST_AUTO_TEST_CASE(testBermudanSwaptionFdHullWhiteSwaptionEngine)
{
    BOOST_CHECK(AdjointBermudanSwaptionTest::testFdHullWhiteSwaptionEngine());
}

BOOST_AUTO_TEST_SUITE_END()

#endif