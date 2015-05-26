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

#include "adjointeuropeanoptionportfoliotest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include <ql/quantlib.hpp>
#include <boost/timer.hpp>
#include <iostream>

using namespace QuantLib;
using namespace boost::unit_test_framework;


#ifdef CL_TAPE_CPPAD
namespace {
    enum {
#if defined CL_GRAPH_GEN
        // Number of points for performance graphic.
        iterNo = 100,
        // Step for portfolio size for performance testing .
        step = 1,
        // Defines performance accuracy. Its value is a minimum number
        // of calling of O(1) complexity methods per one performance test.
        // With iterNumFactor = 100000 one GreekTestHelper::recordPerformance()
        // call runs about 1 second in release mode.
        iterNumFactor = 1000000,
#else
        iterNo = 1,
        step = 1000,
        iterNumFactor = 1,
#endif
        // Number of points for dependency graphics.
        pointNo = 1000
    };

    // Parameters of single option.
    enum DataType {
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

    enum Greek {
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

    std::string greekName(Greek g)
    {
        static const std::string greekNames[] = {
            "Delta"
            , "Vega"
            , "Theta"
            , "Rho"
            , "Gamma"
            , "Vomma"
            , "Vanna"
            , "Charm"
            , "Veta"
            , "Vera"
            , "Speed"
            , "Ultima"
            , "Zomma"
            , "Color"
        };
        return greekNames[g];
    }

    // Returns BlackCalculator for option with parameters determinated by
    // optionData vector in format { stock, strike, sigma, time, rate}.
    inline BlackCalculator getBC(Option::Type type, const std::vector<Real>& optionData)
    {
        return BlackCalculator(type,
            optionData[strike],
            optionData[stock] * std::exp(optionData[rate] * optionData[time]),
            optionData[sigma] * std::sqrt(optionData[time]),
            std::exp(-optionData[rate] * optionData[time]));
    }


    // Abstract class that determinate method of adjoint result checking.
    class CheckingModel
    {
    public:
        // Relative tolerance and absolute tolerance settings.
        CheckingModel(double relativeTol, double absTol, std::string const& description = "BlackCalculator")
            : relativeTol_(relativeTol)
            , absTol_(absTol)
            , description_(description)
        { }

        virtual ~CheckingModel() {}

        // Should calculate Greek without adjoint.
        virtual Real operator()(Option::Type type, const std::vector<Real>& optionData) const = 0;

        std::string description_;
        double relativeTol_;
        double absTol_;
    };

    // Provides Greek value from given function.
    class SimpleCheckingModel
        : public CheckingModel
    {
    public:
        typedef Real(*func_type)(Option::Type, const std::vector<Real>&);
        // Constructor parameters:
        // function - function that calculates Greek.
        SimpleCheckingModel(func_type function,
            double relativeTol = 1e-8, double absTol = 1e-10)
            : CheckingModel(relativeTol, absTol)
            , function_(function)
        { }

        // Returns Greek value.
        inline Real operator()(Option::Type type, const std::vector<Real>& optionData) const
        {
            return (*function_)(type, optionData);
        }

    protected:
        func_type function_;
    };

    // Provides Greek value with finite difference derivative calculation.
    class FdCheckingModel
            : public CheckingModel
    {
    public:
        typedef Real (*func_type)(Option::Type, const std::vector<Real>&);
        // Constructor parameters:
        // function - function that is differentiated,
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
        inline Real operator()(Option::Type type, const std::vector<Real>& optionData) const
        {
            std::vector<Real> optionData_left(optionData);
            optionData_left[param_] -= shift_;
            std::vector<Real> optionData_right(optionData);
            optionData_right[param_] += shift_;
            return ((*function_)(type, optionData_right) - (*function_)(type, optionData_left)) / (2 * shift_);
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
        typedef Real(*func_type)(Option::Type, const std::vector<Real>&);
        // Constructor parameters:
        // function - function that is differentiated,
        // param - function parameter to differentiate,
        // shift - shift for finite difference calculation.
        Fd2CheckingModel(func_type function, DataType param,
            double shift = 1e-5,
            double relativeTol = 1e-3, double absTol = 1e-10,
            std::string const& description = "Finite diff.")
            : CheckingModel(relativeTol, absTol, description)
                , function_(function)
                , param_(param)
                , shift_(shift)
        { }

        // Returns second order finite difference derivative of function_ with respect to param_.
        inline Real operator()(Option::Type type, const std::vector<Real>& optionData) const
        {
            std::vector<Real> optionData_left(optionData);
            optionData_left[param_] -= shift_;
            std::vector<Real> optionData_right(optionData);
            optionData_right[param_] += shift_;
            return ((*function_)(type, optionData_right) + (*function_)(type, optionData_left) - 2 * (*function_)(type, optionData)) / (shift_ * shift_);
        }

    protected:
        func_type function_;
        DataType param_;
        double shift_;
    };


    // Abstract class that provides methods for Greek calculation and testing.
    class GreekTestHelper
    {
    public:
        // Abstract class for option portfolio data generating.
        class DataGenerator
        {
        public:
            virtual void generate(std::vector<std::vector<cl::TapeDouble>>&, Size) = 0;
            virtual ~DataGenerator() {}
        };

        // Generate option portfolio data with all parameters equal, exept
        // one that is uniform distributed.
        class UniformDistr
            : public DataGenerator
        {
        public:
            // Constructor parameters:
            // var - uniform distributed option parameter,
            // left, rigth - bounds of distribution.
            UniformDistr(DataType var, Real left, Real right)
                : var_(var)
                , left_(left)
                , right_(right)
            { }

            void generate(std::vector<std::vector<cl::TapeDouble>>& data, Size size)
            {
                data[rate  ].resize(size, 0.05);
                data[stock ].resize(size, 100.0);
                data[strike].resize(size, 120.0);
                data[sigma ].resize(size, 0.3);
                data[time  ].resize(size, 1.0);
                if (size != 1)
                {
                    for (Size i = 0; i < size; i++)
                    {
                        data[var_][i] = left_ + ((right_ - left_) * i) / (size - 1);
                    }
                }
                else
                {
                    data[var_][0] = left_;
                }
            }

        private:
            DataType var_;
            Real left_;
            Real right_;
        };

        // Generate option portfolio data with pseudorandom distributed parameters.
        class PseudoRandom
            : public DataGenerator
        {
        public:
            void generate(std::vector<std::vector<cl::TapeDouble>>& data, Size size)
            {
                data[rate].resize(size, 0.05);
                data[stock].resize(size);
                data[strike].resize(size);
                data[sigma].resize(size);
                data[time].resize(size);
                // Pseudorandom values.
                for (Size i = 0; i < size; i++)
                {
                    data[stock][i] = 100.0 + std::log(i + 1);
                    data[strike][i] = data[stock][i] * (1 + 0.1 * (0.87 * std::cos(3.26 * i) + 1.23 * std::sin(7.48 * i)));
                    data[sigma][i] = 0.3 + 0.2 * std::sin(i);
                    data[time][i] = 2.0 + std::cos(2 * i);
                }
            }
        };
        
        // Constructor parameters:
        // size - size of portfolio,
        // model - adjoint result checking method,
        // signFactor - the sign with which derivatives should be taken,
        // dataGen - data generator, default - PseudoRandom(),
        // type - type of options.
        GreekTestHelper(Size size, boost::shared_ptr<CheckingModel> model, int signFactor = 1, boost::shared_ptr<DataGenerator> dataGen = 0, Option::Type type = Option::Call)
            : size_(size)
            , model_(model)
            , signFactor_(signFactor)
            , type_(type)
            , data_(5)
            , f_(0)
            , prices_()
            , totalPrice_()
            , forwardResults_()
            , reverseResults_()
            , modelResults_()
            , log_(0, *this)
        {
            if (dataGen == 0)
            {
                PseudoRandom().generate(data_, size_);
            }
            else
            {
                dataGen->generate(data_, size_);
            }
        }

        virtual ~GreekTestHelper() {}

        // Tests adjoint derivative calculation in forward and reverse mode.
        bool test()
        {
            recordTape();
            calcForward();
            calcReverse();
            calcModel();
            return check();
        }

        // Tests adjoint derivative calculation in reverse mode only.
        bool testReverse()
        {
            recordTape();
            calcReverse();
            calcModel();
            return checkReverse();
        }

        // Prints performance for all calculations.
        bool printPerformance(std::ostream& out = std::cout)
        {
            out << "\nNumber of European option in portfolio: " << size_ << std::endl;
            out << "Time for tape recording:\t" << testTapePerformance() << std::endl;
            out << "Time for Forward mode:  \t" << testForwardPerformance() << std::endl;
            out << "Time for Reverse mode:  \t" << testReversePerformance() << std::endl;
            out << "Time for " << model_->description_ << ":  \t" << testModelPerformance() << std::endl;
            out << std::endl;
            return check();
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
            if (testNo == 0)
            {
                // Setting default number of tests.
                // recordTape() has O(size_) complexity, so divide by size_.
                testNo = iterNumFactor / size_ + 1;
            }
            log_.beforeTape();
            boost::timer timer;
            for (Size i = 0; i < testNo; i++)
            {
                recordTape();
            }
            log_.afterTape();
            return timer.elapsed() / testNo;
        }

        // Returns performance of forward mode derivatives calculations.
        inline double testForwardPerformance(Size testNo = 0)
        {
            if (!testNo)
            {
                // Setting default number of tests.
                // calcForward() has O(size_^2) complexity, so divide by size_^2.
                testNo = iterNumFactor / size_ / size_ + 1;
            }
            if (f_.get() == 0)
            {
                recordTape();
            }
            log_.beforeForward();
            boost::timer timer;
            for (Size i = 0; i < testNo; i++)
            {
                calcForward();
            }
            log_.afterForward();
            return timer.elapsed() / testNo;
        }

        // Returns performance of reverse mode derivatives calculations.
        inline double testReversePerformance(Size testNo = 0)
        {
            if (!testNo)
            {
                // Setting default number of tests.
                // calcReverse() has O(size_) complexity, so divide by size_.
                testNo = iterNumFactor / size_ + 1;
            }
            if (f_.get() == 0)
            {
                recordTape();
            }
            log_.beforeReverse();
            boost::timer timer;
            for (Size i = 0; i < testNo; i++)
            {
                calcReverse();
            }
            log_.afterReverse();
            return timer.elapsed() / testNo;
        }

        // Returns performance of Greek calculations with checking model.
        inline double testModelPerformance(Size testNo = 0)
        {
            if (!testNo)
            {
                // Setting default number of tests.
                // calcModel() has O(size_) complexity, so divide by size_.
                testNo = iterNumFactor / size_ + 1;
            }
            log_.beforeModel();
            boost::timer timer;
            for (Size i = 0; i < testNo; i++)
            {
                calcModel();
            }
            log_.afterModel();
            return timer.elapsed() / testNo;
        }

        // Calculate Greeks with reverse mode only.
        inline void calcWithReverse()
        {
            log_.beforeTape();
            recordTape();
            log_.afterTape();
            log_.beforeReverse();
            calcReverse();
            log_.afterReverse();
        }

        // Returns Greeks calculated in reverse mode.
        inline std::vector<double> const& reverse()
        {
            return reverseResults_;
        }

        // Returns prices of every option in portfolio.
        inline std::vector<Real> const& prices()
        {
            return prices_;
        }

        // Returns vector of one of the option parameters for portfolio.
        inline std::vector<cl::TapeDouble> const& data(DataType t)
        {
            return data_[t];
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

        void setLogger(cl::AdjointTestOutput* logger)
        {
            log_.logger_ = logger;
        }

        // Returns number of independent variables.
        virtual Size indepVarNumber() = 0;

    protected:
        // Calculates Greeks in forward mode.
        virtual void calcForward() = 0;
        // Calculates Greeks in reverse mode.
        virtual void calcReverse() = 0;
        // Calculates price and records tape.
        virtual void recordTape() = 0;

        // Calculates price of portfolio and each option.
        inline void calcPrice()
        {
            totalPrice_.resize(1, 0);
            prices_.resize(size_);
            for (Size i = 0; i < size_; i++)
            {
                prices_[i] = getBC(type_, optionData(i)).value();
                totalPrice_.front() += prices_[i];
            }
        }

    private:
        // Checks forward mode, reverse mode and checking model result consistency.
        bool check(double relativeTol = 0, double absTol = 0)
        {
            bool result = true;
            if (!relativeTol)
            {
                relativeTol = model_->relativeTol_;
            }
            if (!absTol)
            {
                absTol = model_->absTol_;
            }
            log_ << "Derivatives" << std::endl;
            log_ << "Forward mode || Reverse mode || " << model_->description_ << std::endl;
            for (Size i = 0; i < size_; i++)
            {
                log_ << std::setw(12) << forwardResults_[i] << " ";
                log_ << std::setw(15) << reverseResults_[i] << " ";
                log_ << std::setw(15) << modelResults_[i] << std::endl;

                // Maximum of absolut values of results.
                Real maxabs = std::max(std::abs(modelResults_[i]), std::max(std::abs(forwardResults_[i]), std::abs(reverseResults_[i])));
                Real tol = std::max(relativeTol * maxabs, absTol);
                if (std::abs(forwardResults_[i] - modelResults_[i]) > tol)
                {
                    result = false;
                    BOOST_ERROR("\nFailed to reproduce expected derivative[" << i << "] at forward mode."
                        << "\n    calculated: " << forwardResults_[i]
                        << "\n    expected:   " << modelResults_[i]
                        << "\n    tolerance:  " << tol);
                }
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

        // Checks reverse mode and checking model result consistency.
        bool checkReverse(double relativeTol = 0, double absTol = 0)
        {
            bool result = true;
            if (!relativeTol)
            {
                relativeTol = model_->relativeTol_;
            }
            if (!absTol)
            {
                absTol = model_->absTol_;
            }
            log_ << "Derivatives" << std::endl;
            log_ << "Reverse mode || " << model_->description_ << std::endl;
            for (Size i = 0; i < size_; i++)
            {
                log_ << std::setw(12) << reverseResults_[i] << " ";
                log_ << std::setw(15) << modelResults_[i] << std::endl;

                // Maximum of absolut values of results.
                Real maxabs = std::max(std::abs(modelResults_[i]),std::abs(reverseResults_[i]));
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

        // Calculates Greeks using checking model.
        inline void calcModel()
        {
            modelResults_.resize(size_);
            for (Size i = 0; i < size_; i++)
            {
                modelResults_[i] = (*model_)(type_, optionData(i));
            }
        }

        // Returns parameters for one option in format { stock, strike, sigma, time, rate}.
        inline std::vector<Real> optionData(Size index)
        {
            std::vector<Real> optionData(5);
            for (Size dt = stock; dt <= rate; dt++)
            {
                optionData[dt] = data_[dt][index];
            }
            return optionData;
        }

        struct Log
        {
            Log(cl::AdjointTestOutput* p, const GreekTestHelper& helper)
            : logger_(p)
            , this_(&helper)
            { }

            void beforeTape()
            {
                *this << "Start of tape recording: " << currentTime() << std::endl;
                timer_.restart();
            }
            void afterTape()
            {
                *this << "End of tape recording." << std::endl;
                *this << "Time for tape recording: " << timer_.elapsed() << std::endl;
                timer_.restart();
            }

            void beforeForward()
            {
                *this << "Start of differentiation in Forward mode: " << currentTime() << std::endl;
                timer_.restart();
            }
            void afterForward()
            {
                *this << "End of differentiation in Forward mode." << std::endl;
                *this << "Time for Forward mode: " << timer_.elapsed() << std::endl;
                timer_.restart();
            }

            void beforeReverse()
            {
                *this << "Start of differentiation in Reverse mode: " << currentTime() << std::endl;
                timer_.restart();
            }
            void afterReverse()
            {
                *this << "End of differentiation in Reverse mode." << std::endl;
                *this << "Time for Reverse mode: " << timer_.elapsed() << std::endl;
                timer_.restart();
            }

            void beforeModel()
            {
                *this << "Start of differentiation using " + this_->model_->description_ + " method: " << currentTime() << std::endl;
                timer_.restart();
            }
            void afterModel()
            {
                *this << "End of differentiation using " + this_->model_->description_ + " method." << std::endl;
                *this << "Time for " + this_->model_->description_ + " method: " << timer_.elapsed() << std::endl;
                timer_.restart();
            }

            template<class T>
            Log& operator<<(T&& t)
            {
                if (logger_)
                {
                    logger_->log() << std::forward<T>(t);
                }
                return *this;
            }

            // this is the type of std::cout
            typedef std::basic_ostream<char, std::char_traits<char> > cout_manip_type;

            // this is the function signature of std::endl
            typedef cout_manip_type& (*standard_manipulator)(cout_manip_type&);

            // define an operator<< to take in std::endl
            Log& operator<<(standard_manipulator mnp)
            {
                logger_->log() << mnp;
                return *this;
            }

            cl::AdjointTestOutput* logger_;
            boost::timer timer_;
            const GreekTestHelper* this_;
        };

    protected:
        Size size_;
        boost::shared_ptr<CheckingModel> model_;
        int signFactor_;
        Option::Type type_;
        std::vector<std::vector<cl::TapeDouble>> data_;
        std::auto_ptr<cl::TapeFunction<double> > f_;
        std::vector<Real> prices_;
        std::vector<cl::TapeDouble> totalPrice_;
        std::vector<double> forwardResults_;
        std::vector<double> reverseResults_;
        std::vector<Real> modelResults_;
        Log log_;
    };

    // Helper for first order Greeks.
    class Greek1Helper
        : public GreekTestHelper
    {
    public:
        // Constructor parameters:
        // size - size of portfolio,
        // independent - parameter of option to differentiate,
        // model - adjoint result checking method,
        // signFactor - the sign with which derivatives should be taken,
        // dataGen - data generator, default - PseudoRandom(),
        // type - type of options.
        Greek1Helper(Size size, DataType independent,
            boost::shared_ptr<CheckingModel> model, int signFactor = 1,
            boost::shared_ptr<DataGenerator> dataGen = 0,
            Option::Type type = Option::Call)
            : GreekTestHelper(size, model, signFactor, dataGen, type)
                , independent_(independent)
        { }

    protected:
        void calcForward()
        {
            forwardResults_.resize(size_);
            // Direction for directional derivative.
            std::vector<double> dX(size_, 0);
            for (Size i = 0; i < size_; i++)
            {
                dX[i] = 1;
                forwardResults_[i] = signFactor_ * f_->Forward(1, dX)[0];
                dX[i] = 0;
            }
        }

        void calcReverse()
        {
            reverseResults_ = f_->Reverse(1, std::vector<double>(1, 1));
            for (Size i = 0; i < size_; i++)
            {
                reverseResults_[i] *= signFactor_;
            }
        }

        void recordTape()
        {
            Independent(data_[independent_]);
            calcPrice();
            f_.reset(new cl::TapeFunction<double>(data_[independent_], totalPrice_));
        }

        inline Size indepVarNumber()
        {
            return size_;
        }

    private:
        DataType independent_;
    };

    // Helper for second order Greeks.
    class Greek2Helper
        : public GreekTestHelper
    {
    public:
        // Constructor parameters:
        // size - size of portfolio,
        // independent - parameter of option to differentiate,
        // model - adjoint result checking method,
        // signFactor - the sign with which derivatives should be taken,
        // dataGen - data generator, default - PseudoRandom(),
        // type - type of options.
        Greek2Helper(Size size, DataType independent,
            boost::shared_ptr<CheckingModel> model, int signFactor = 1,
            boost::shared_ptr<DataGenerator> dataGen = 0,
            Option::Type type = Option::Call)
            : GreekTestHelper(size, model, signFactor, dataGen, type)
                , independent_(independent)
        { }

    protected:
        void calcForward()
        {
            forwardResults_.resize(size_);
            // Direction for directional derivative.
            std::vector<double> dX(size_, 0);
            for (Size i = 0; i < size_; i++)
            {
                dX[i] = 1;
                f_->Forward(1, dX);
                dX[i] = 0;
                forwardResults_[i] = signFactor_ * 2 * f_->Forward(2, dX)[0];
            }
        }

        void calcReverse()
        {
            f_->Forward(1, std::vector<double>(size_, 1));
            std::vector<double> dw = f_->Reverse(2, std::vector<double>(1, 1));
            reverseResults_.resize(size_);
            for (Size i = 0; i < size_; i++)
            {
                reverseResults_[i] = signFactor_ * dw[2 * i + 1];
            }
        }

        void recordTape()
        {
            Independent(data_[independent_]);
            calcPrice();
            f_.reset(new cl::TapeFunction<double>(data_[independent_], totalPrice_));
        }

        inline Size indepVarNumber()
        {
            return size_;
        }

    private:
        DataType independent_;
    };

    // Helper for mixed second order Greeks.
    class Greek2MixedHelper
        : public GreekTestHelper
    {
    public:
        // Constructor parameters:
        // size - size of portfolio,
        // independent1, independent2 - parameters of option to differentiate,
        // model - adjoint result checking method,
        // signFactor - the sign with which derivatives should be taken,
        // dataGen - data generator, default - PseudoRandom(),
        // type - type of options.
        Greek2MixedHelper(Size size, DataType independent1, DataType independent2,
            boost::shared_ptr<CheckingModel> model, int signFactor = 1,
            boost::shared_ptr<DataGenerator> dataGen = 0,
            Option::Type type = Option::Call)
            : GreekTestHelper(size, model, signFactor, dataGen, type)
                , independent1_(independent1)
                , independent2_(independent2)
        { }

    protected:
        void calcForward()
        {
            forwardResults_.resize(size_);
            // Direction for directional derivative.
            std::vector<double> dX(2 * size_, 0);
            for (Size i = 0; i < size_; i++)
            {
                // Mixed derivative calculated from second order derivatives like this:
                // z = x + y
                // w = x - y
                // ddF/dzdz = ddF/dxdx + 2 * ddF/dxdy + ddF/dydy
                // ddF/dwdw = ddF/dxdx - 2 * ddF/dxdy + ddF/dydy
                // ddF/dxdy = (ddF/dzdz - ddF/dwdw) / 4
                dX[size_ + i] = 1;
                dX[i] = 1;
                f_->Forward(1, dX);
                dX[size_ + i] = 0;
                dX[i] = 0;
                double dzdz = 2 * f_->Forward(2, dX)[0];
                dX[size_ + i] = -1;
                dX[i] = 1;
                f_->Forward(1, dX);
                dX[size_ + i] = 0;
                dX[i] = 0;
                double dwdw = 2 * f_->Forward(2, dX)[0];
                forwardResults_[i] = signFactor_ * (dzdz - dwdw) / 4;
            }
        }

        void calcReverse()
        {
            reverseResults_.resize(size_);
            std::vector<double> dX(2 * size_, 0);
            std::fill(dX.begin(), dX.begin() + size_, 1);
            f_->Forward(1, dX);
            std::vector<double> dw = f_->Reverse(2, std::vector<double>(1, 1));
            for (Size i = 0; i < size_; i++)
            {
                reverseResults_[i] = signFactor_ * dw[2 * (size_ + i) + 1];
            }
        }

        void recordTape()
        {
            std::vector<cl::TapeDouble> ind(data_[independent1_]);
            ind.reserve(2 * size_);
            ind.insert(ind.end(), data_[independent2_].begin(), data_[independent2_].end());
            Independent(ind);
            std::copy(ind.begin(), ind.begin() + size_, data_[independent1_].begin());
            std::copy(ind.begin() + size_, ind.end(), data_[independent2_].begin());
            calcPrice();
            f_.reset(new cl::TapeFunction<double>(ind, totalPrice_));
        }

        inline Size indepVarNumber()
        {
            return 2 * size_;
        }

    private:
        DataType independent1_;
        DataType independent2_;
    };

    // Helper for third order Greeks.
    class Greek3Helper
        : public GreekTestHelper
    {
    public:
        // Constructor parameters:
        // size - size of portfolio,
        // independent - parameter of option to differentiate,
        // model - adjoint result checking method,
        // signFactor - the sign with which derivatives should be taken,
        // dataGen - data generator, default - PseudoRandom(),
        // type - type of options.
        Greek3Helper(Size size, DataType independent,
            boost::shared_ptr<CheckingModel> model, int signFactor = 1,
            boost::shared_ptr<DataGenerator> dataGen = 0,
            Option::Type type = Option::Call)
            : GreekTestHelper(size, model, signFactor, dataGen, type)
                , independent_(independent)
        { }

    protected:
        void calcForward()
        {
            forwardResults_.resize(size_);
            // Direction for directional derivative.
            std::vector<double> dX(size_, 0);
            for (Size i = 0; i < size_; i++)
            {
                dX[i] = 1;
                f_->Forward(1, dX);
                dX[i] = 0;
                f_->Forward(2, dX);
                forwardResults_[i] = signFactor_ * 6 * f_->Forward(3, dX)[0];
            }
        }

        void calcReverse()
        {
            reverseResults_.resize(size_);
            f_->Forward(1, std::vector<double>(size_, 1));
            f_->Forward(2, std::vector<double>(size_, 0));
            std::vector<double> dw = f_->Reverse(3, std::vector<double>(1, 1));
            for (Size i = 0; i < size_; i++)
            {
                reverseResults_[i] = signFactor_ * 2 * dw[3 * i + 2];
            }
        }

        void recordTape()
        {
            Independent(data_[independent_]);
            calcPrice();
            f_.reset(new cl::TapeFunction<double>(data_[independent_], totalPrice_));
        }

        inline Size indepVarNumber()
        {
            return size_;
        }

    private:
        DataType independent_;
    };

    // Helper for mixed third order Greeks.
    class Greek3MixedHelper
        : public GreekTestHelper
    {
    public:
        // Constructor parameters:
        // size - size of portfolio,
        // independent1 - parameter of option for first order derivative,
        // independent2 - parameter of option for secont order derivative,
        // model - adjoint result checking method,
        // signFactor - the sign with which derivatives should be taken,
        // dataGen - data generator, default - PseudoRandom(),
        // type - type of options.
        Greek3MixedHelper(Size size, DataType independent1, DataType independent2,
            boost::shared_ptr<CheckingModel> model, int signFactor = 1,
            boost::shared_ptr<DataGenerator> dataGen = 0,
            Option::Type type = Option::Call)
            : GreekTestHelper(size, model, signFactor, dataGen, type)
                , independent1_(independent1)
                , independent2_(independent2)
        { }

    protected:
        void calcForward()
        {
            forwardResults_.resize(size_);
            std::vector<double> dX(2 * size_, 0);
            for (Size i = 0; i < size_; i++)
            {
                // Mixed derivative calculated from third order derivatives like this:
                // z = x + y
                // w = x - y
                // d3F/dz3 = d3F/dx3 + 3 * d3F/dx2dy + 3 * d3F/dxdy2 + d3F/dy3
                // d3F/dw3 = d3F/dx3 - 3 * d3F/dx2dy + 3 * d3F/dxdy2 - d3F/dy3
                // d3F/dx2dy = (d3F/dz3 - d3F/dw3 - 2 * d3F/dy3) / 6
                dX[size_ + i] = 1;
                dX[i] = 1;
                f_->Forward(1, dX);
                dX[size_ + i] = 0;
                dX[i] = 0;
                f_->Forward(2, dX);
                double dz3 = 6 * f_->Forward(3, dX)[0];
                dX[size_ + i] = 1;
                dX[i] = -1;
                f_->Forward(1, dX);
                dX[size_ + i] = 0;
                dX[i] = 0;
                f_->Forward(2, dX);
                double dw3 = 6 * f_->Forward(3, dX)[0];
                dX[i] = 1;
                f_->Forward(1, dX);
                dX[i] = 0;
                f_->Forward(2, dX);
                double dy3 = 6 * f_->Forward(3, dX)[0];
                forwardResults_[i] = signFactor_ * (dz3 - dw3 - 2 * dy3) / 6;
            }
        }

        void calcReverse()
        {
            reverseResults_.resize(size_);
            std::vector<double> dX(2 * size_, 0);
            std::fill(dX.begin() + size_, dX.end(), 1);
            f_->Forward(1, dX);
            std::fill(dX.begin() + size_, dX.end(), 0);
            f_->Forward(2, std::vector<double>(2 * size_, 0));
            std::vector<double> dw = f_->Reverse(3, std::vector<double>(1, 1));
            for (Size i = 0; i < size_; i++)
            {
                reverseResults_[i] = signFactor_ * 2 * dw[3 * i + 2];
            }
        }

        void recordTape()
        {
            std::vector<cl::TapeDouble> ind(data_[independent1_]);
            ind.reserve(2 * size_);
            ind.insert(ind.end(), data_[independent2_].begin(), data_[independent2_].end());
            Independent(ind);
            std::copy(ind.begin(), ind.begin() + size_, data_[independent1_].begin());
            std::copy(ind.begin() + size_, ind.end(), data_[independent2_].begin());
            calcPrice();
            f_.reset(new cl::TapeFunction<double>(ind, totalPrice_));
        }

        inline Size indepVarNumber()
        {
            return 2 * size_;
        }

    private:
        DataType independent1_;
        DataType independent2_;
    };

    // Provides GreekTestHelper classes for setted Greek.
    class TestHelperFactory
    {
    public:
        TestHelperFactory(Greek greek) : greek_(greek) {}

        void setGreek(Greek greek)
        {
            greek_ = greek;
        }

        // Provides GreekTestHelper class for setted Greek with default portfolio data.
        inline boost::shared_ptr<GreekTestHelper> get(Size portfolioSize)
        {
            return get(portfolioSize,
                boost::shared_ptr<GreekTestHelper::DataGenerator>(
                    new GreekTestHelper::PseudoRandom()));
        }

        // Provides GreekTestHelper class for setted Greek with
        // uniform distributed on var parameter portfolio data.
        inline boost::shared_ptr<GreekTestHelper> getUniform(Size portfolioSize, DataType var, Real left, Real right)
        {
            return get(portfolioSize,
                boost::shared_ptr<GreekTestHelper::DataGenerator>(
                new GreekTestHelper::UniformDistr(var, left, right)));
        }

        // Provides GreekTestHelper class for setted Greek with given portfolio data.
        inline boost::shared_ptr<GreekTestHelper> get(Size portfolioSize, boost::shared_ptr<GreekTestHelper::DataGenerator> dataGen)
        {
            switch (greek_)
            {
            case Delta:
                // Differentiating with respect to stock price,
                // Checking with BlackCalculator method.
                return boost::shared_ptr<GreekTestHelper>(
                    new Greek1Helper(
                        portfolioSize, stock,
                        boost::shared_ptr<CheckingModel>(new SimpleCheckingModel(&delta)),
                        1, dataGen));
                break;
            case Vega:
                // Differentiating with respect volatility,
                // Checking with BlackCalculator method.
                return boost::shared_ptr<GreekTestHelper>(
                    new Greek1Helper(
                        portfolioSize, sigma,
                        boost::shared_ptr<CheckingModel>(new SimpleCheckingModel(&vega)),
                        1, dataGen));
                break;
            case Theta:
                // Differentiating with respect to maturity time,
                // Checking with BlackCalculator method.
                return boost::shared_ptr<GreekTestHelper>(
                    new Greek1Helper(
                        portfolioSize, DataType::time,
                        boost::shared_ptr<CheckingModel>(new SimpleCheckingModel(&theta)),
                        -1, dataGen));
                break;
            case Rho:
                // Differentiating with respect to interest rate,
                // Checking with BlackCalculator method.
                return boost::shared_ptr<GreekTestHelper>(
                    new Greek1Helper(
                        portfolioSize, rate,
                        boost::shared_ptr<CheckingModel>(new SimpleCheckingModel(&rho)),
                        1, dataGen));
                break;
            case Gamma:
                // Differentiating with respect to stock price,
                // Checking with BlackCalculator method.
                return boost::shared_ptr<GreekTestHelper>(
                    new Greek2Helper(
                        portfolioSize, stock,
                        boost::shared_ptr<CheckingModel>(new SimpleCheckingModel(&gamma)),
                        1, dataGen));
                break;
            case Vomma:
                return boost::shared_ptr<GreekTestHelper>(
                    new Greek2Helper(
                        portfolioSize, sigma,
                        boost::shared_ptr<CheckingModel>(new FdCheckingModel(&vega, sigma, 3e-7, 1e-4)),
                        1, dataGen));
                break;
            case Vanna:
                return boost::shared_ptr<GreekTestHelper>(
                    new Greek2MixedHelper(
                        portfolioSize, sigma, stock,
                        boost::shared_ptr<CheckingModel>(new FdCheckingModel(&delta, sigma, 3e-7, 1e-4, 1e-10)),
                        1, dataGen));
                break;
            case Charm:
                return boost::shared_ptr<GreekTestHelper>(
                    new Greek2MixedHelper(
                        portfolioSize, stock, DataType::time,
                        boost::shared_ptr<CheckingModel>(new FdCheckingModel(&theta, stock, 1e-6, 1e-4)),
                        -1, dataGen));
                break;
            case Veta:
                return boost::shared_ptr<GreekTestHelper>(
                    new Greek2MixedHelper(
                        portfolioSize, sigma, DataType::time,
                        boost::shared_ptr<CheckingModel>(new FdCheckingModel(&vega, DataType::time)),
                        1, dataGen));
                break;
            case Vera:
                return boost::shared_ptr<GreekTestHelper>(
                    new Greek2MixedHelper(
                        portfolioSize, sigma, rate,
                        boost::shared_ptr<CheckingModel>(new FdCheckingModel(&rho, sigma, 3e-7, 1e-4)),
                        1, dataGen));
                break;
            case Speed:
                return boost::shared_ptr<GreekTestHelper>(
                    new Greek3Helper(
                        portfolioSize, stock,
                        boost::shared_ptr<CheckingModel>(new FdCheckingModel(&gamma, stock, 1e-5, 1e-7)),
                        1, dataGen));
                break;
            case Ultima:
                return boost::shared_ptr<GreekTestHelper>(
                    new Greek3Helper(
                        portfolioSize, sigma,
                        boost::shared_ptr<CheckingModel>(new Fd2CheckingModel(&vega, sigma, 1e-4)),
                        1, dataGen));
                break;
            case Zomma:
                return boost::shared_ptr<GreekTestHelper>(
                    new Greek3MixedHelper(
                        portfolioSize, sigma, stock,
                        boost::shared_ptr<CheckingModel>(new FdCheckingModel(&gamma, sigma, 1e-6, 1e-4)),
                        1, dataGen));
                break;
            case Color:
                return boost::shared_ptr<GreekTestHelper>(
                    new Greek3MixedHelper(
                        portfolioSize, DataType::time, stock,
                        boost::shared_ptr<CheckingModel>(new FdCheckingModel(&gamma, DataType::time)),
                        1, dataGen));
                break;
            default:
                return 0;
            }
        }

    private:
        static inline Real delta(Option::Type type, const std::vector<Real>& optionData)
        {
            return getBC(type, optionData).delta(optionData[stock]);
        }

        static inline Real gamma(Option::Type type, const std::vector<Real>& optionData)
        {
            return getBC(type, optionData).gamma(optionData[stock]);
        }

        static inline Real vega(Option::Type type, const std::vector<Real>& optionData)
        {
            return getBC(type, optionData).vega(optionData[DataType::time]);
        }

        static inline Real theta(Option::Type type, const std::vector<Real>& optionData)
        {
            return getBC(type, optionData).theta(optionData[stock], optionData[DataType::time]);
        }

        static inline Real rho(Option::Type type, const std::vector<Real>& optionData)
        {
            return getBC(type, optionData).rho(optionData[DataType::time]);
        }

        Greek greek_;
    };
    
    // Makes graphics for performance time and tape size.
    // Returns true if calculated results checking passed.
    inline bool recordPerformance(Greek greek)
    {
        cl::AdjointTestOutput outPerform("AdjointEuropeanOptionPortfolio//" + greekName(greek) + "CallPortfolio", {
            { "title", greekName(greek) + " calculation performance" }
            , { "not_clear", "Not"}
            , { "ylabel", "Time (s)" }
            , { "line_box_width", "-5" }
            , { "cleanlog", "false" }
        });
        cl::AdjointTestOutput outSize("AdjointEuropeanOptionPortfolio//" + greekName(greek) + "CallPortfolio", {
            { "title", "Tape size dependence on number of options in portfolio" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
            , { "cleanlog", "false" }
        });
        TestHelperFactory factory(greek);
        std::vector<PerformanceTime> perfTimes(iterNo);
        std::vector<TapeSize> tapeMemory(iterNo);
        auto p = perfTimes.begin();
        auto m = tapeMemory.begin();
        bool result = true;
        for (Size i = 1; i <= iterNo; i++, ++p)
        {
            auto helper = factory.get(step * i);
            helper->setLogger(&outPerform);
            result &= helper->recordPerformance(*p);
            *m++ = { helper->indepVarNumber(), helper->memory() };
        }
        outPerform << perfTimes;
        outSize << tapeMemory;
        return result;
    }

    // Struct for graphics recording.
    template <DataType var, Greek greek>
    struct OutputData
    {
        static std::deque<std::string > get_columns()
        {
            return{ toString(var, false, true), greekName(greek) };
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

    
    // Makes graphics for greek dependence.
    // Function parameters:
    // var - option parameter on which dependence we record,
    // left, right - bounds of variation of var parameter,
    // greek - Greek which dependency is recording,
    // filename - name of file for graphic, "" - auto,
    // title - title of graphic, "" - auto.
    template <DataType var, Greek greek>
    void recordDependence(Real left, Real right, std::string filename, std::string title)
    {
        cl::AdjointTestOutput out("AdjointEuropeanOptionPortfolio//" + greekName(greek) + "CallPortfolio//output", {
            { "filename", filename.empty() ? toString(var) + "Dependence" : filename}
            , { "not_clear", "Not" }
            , { "title", title.empty() ? greekName(greek) + " dependence on " + toString(var, true) : title }
            , { "ylabel", "Greek value" }
            , { "line_box_width", "2" }
            , { "cleanlog", "false" }
        });
        auto helper = TestHelperFactory(greek).getUniform(pointNo, var, left, right);
        helper->setLogger(&out);
        helper->calcWithReverse();
        std::vector<OutputData<var, greek>> outData(pointNo);
        for (Size i = 0; i < pointNo; i++)
        {
            outData[i].data_.push_back(helper->data(var)[i]);
            outData[i].data_.push_back(helper->reverse()[i]);
        }
        out << outData;
    }
    void cleanFolder(Greek greek)
    {
        cl::AdjointTestOutput("AdjointEuropeanOptionPortfolio//" + greekName(greek) + "CallPortfolio");
    }
}
#endif


bool AdjointEuropeanOptionPortfolioTest::testDeltaCallPortfolio()
{
    /*BOOST_TEST_MESSAGE(
        "\nDelta test"
        "\nPortfolio of European Call options with different parameters is examined to"
        "\nmeasure sensitivity of the value of portfolio to changes in the underlying"
        "\nasset price of European option in portfolio while holding the other"
        "\nparameters fixed"
        "\n");*/
    BOOST_TEST_MESSAGE("Testing Delta calculations for Europeam Option Portfolio...");

#ifdef CL_TAPE_CPPAD

    cleanFolder(Delta);
    recordDependence<stock, Delta>(30, 250, "", "");
    return recordPerformance(Delta);

#endif
    return true;
}


bool AdjointEuropeanOptionPortfolioTest::testVegaCallPortfolio()
{
    /*BOOST_TEST_MESSAGE(
        "\nVega test"
        "\nPortfolio of European Call options with different parameters is examined to"
        "\nmeasure sensitivity of the value of portfolio to changes in the volatility"
        "\nof European option in portfolio while holding the other parameters fixed"
        "\n");*/
    BOOST_TEST_MESSAGE("Testing Vega calculations for Europeam Option Portfolio...");

#ifdef CL_TAPE_CPPAD

    cleanFolder(Vega);
    recordDependence<sigma, Vega>(0.015, 0.5, "", "");
    return recordPerformance(Vega);

#endif
    return true;
}


bool AdjointEuropeanOptionPortfolioTest::testThetaCallPortfolio()
{
    /*BOOST_TEST_MESSAGE(
        "\nTheta test"
        "\nPortfolio of European Call options with different parameters is examined to"
        "\nmeasure sensitivity of the value of portfolio to changes in the maturity time"
        "\nof European option in portfolio while holding the other parameters fixed"
        "\n");*/
    BOOST_TEST_MESSAGE("Testing Theta calculations for Europeam Option Portfolio...");

#ifdef CL_TAPE_CPPAD

    cleanFolder(Theta);
    recordDependence<DataType::time, Theta>(0.001, 2, "", "");
    return recordPerformance(Theta);

#endif
    return true;
}


bool AdjointEuropeanOptionPortfolioTest::testRhoCallPortfolio()
{
    /*BOOST_TEST_MESSAGE(
        "\nRho test"
        "\nPortfolio of European Call options with different parameters is examined to"
        "\nmeasure sensitivity of the value of portfolio to changes in the interest rate"
        "\nwhile holding the other parameters fixed"
        "\n");*/
    BOOST_TEST_MESSAGE("Testing Rho calculations for Europeam Option Portfolio...");

#ifdef CL_TAPE_CPPAD

    cleanFolder(Rho);
    recordDependence<rate, Rho>(0.0001, 0.5, "", "");
    return recordPerformance(Rho);

#endif
    return true;
}


bool AdjointEuropeanOptionPortfolioTest::testGammaCallPortfolio()
{
    /*BOOST_TEST_MESSAGE(
        "\nGamma test"
        "\nPortfolio of European Put options with different parameters is examined to"
        "\nmeasure second order sensitivity of the value of portfolio to changes in the"
        "\nunderlying asset price of European option in portfolio while holding the"
        "\nother parameters fixed"
        "\n");*/
    BOOST_TEST_MESSAGE("Testing Gamma calculations for Europeam Option Portfolio...");

#ifdef CL_TAPE_CPPAD

    cleanFolder(Gamma);
    recordDependence<stock, Gamma>(10, 200, "", "");
    return recordPerformance(Gamma);

#endif
    return true;
}


bool AdjointEuropeanOptionPortfolioTest::testVommaCallPortfolio()
{
    /*BOOST_TEST_MESSAGE(
        "\nVomma test"
        "\nPortfolio of European Call options with different parameters is examined to"
        "\nmeasure second order sensitivity of the value of portfolio to changes in the"
        "\nvolatility of European option in portfolio while holding the other parameters"
        "\n");*/
    BOOST_TEST_MESSAGE("Testing Vomma calculations for Europeam Option Portfolio...");

#ifdef CL_TAPE_CPPAD

    cleanFolder(Vomma);
    recordDependence<sigma, Vomma>(0.015, 0.25, "", "");
    return recordPerformance(Vomma);

#endif
    return true;
}


bool AdjointEuropeanOptionPortfolioTest::testVannaCallPortfolio()
{
    /*BOOST_TEST_MESSAGE(
        "\nVanna test"
        "\nPortfolio of European Call options with different parameters is examined to"
        "\nmeasure second order sensitivity of the value of portfolio to changes in the"
        "\nvolatility and the underlying asset price of European option in portfolio"
        "\nwhile holding the other parameters fixed"
        "\n");*/
    BOOST_TEST_MESSAGE("Testing Vanna calculations for Europeam Option Portfolio...");

#ifdef CL_TAPE_CPPAD

    cleanFolder(Vanna);
    recordDependence<stock, Vanna>(10   , 300, "", "");
    recordDependence<sigma, Vanna>(0.015, 0.9, "", "");
    return recordPerformance(Vanna);

#endif
    return true;
}


bool AdjointEuropeanOptionPortfolioTest::testCharmCallPortfolio()
{
    /*BOOST_TEST_MESSAGE(
        "\nCharm test"
        "\nPortfolio of European Call options with different parameters is examined to"
        "\nmeasure the instantaneous rate of change of Delta (sensitivity of the value of"
        "\nportfolio to changes in the underlying asset price) over the passage of time"
        "\n");*/
    BOOST_TEST_MESSAGE("Testing Charm calculations for Europeam Option Portfolio...");

#ifdef CL_TAPE_CPPAD

    cleanFolder(Charm);
    recordDependence<stock          , Charm>(50   , 200, "", "");
    recordDependence<DataType::time, Charm>(0.001, 0.5, "", "");
    return recordPerformance(Charm);

#endif
    return true;
}


bool AdjointEuropeanOptionPortfolioTest::testVetaCallPortfolio()
{
    /*BOOST_TEST_MESSAGE(
        "\nVeta test"
        "\nPortfolio of European Call options with different parameters is examined to"
        "\nmeasure the rate of change in the Vega (sensitivity of the value of portfolio"
        "\nto changes in the volatility) with respect to the passage of time"
        "\n");*/
    BOOST_TEST_MESSAGE("Testing Veta calculations for Europeam Option Portfolio...");

#ifdef CL_TAPE_CPPAD

    cleanFolder(Veta);
    recordDependence<sigma          , Veta>(0.015, 0.9, "", "");
    recordDependence<DataType::time, Veta>(0.001, 3  , "", "");
    return recordPerformance(Veta);

#endif
    return true;
}


bool AdjointEuropeanOptionPortfolioTest::testVeraCallPortfolio()
{
    /*BOOST_TEST_MESSAGE(
        "\nVera test"
        "\nPortfolio of European Call options with different parameters is examined to"
        "\nmeasure the rate of change in Rho (sensitivity of the value of portfolio to"
        "\nchanges in the interest rate) with respect to volatility"
        "\n");*/
    BOOST_TEST_MESSAGE("Testing Vera calculations for Europeam Option Portfolio...");

#ifdef CL_TAPE_CPPAD

    cleanFolder(Vera);
    recordDependence<rate , Vera>(0.0001, 1  , "", "");
    recordDependence<sigma, Vera>(0.0003, 0.5, "", "");
    return recordPerformance(Vera);

#endif
    return true;
}


bool AdjointEuropeanOptionPortfolioTest::testSpeedCallPortfolio()
{
    /*BOOST_TEST_MESSAGE(
        "\nSpeed test"
        "\nPortfolio of European Call options with different parameters is examined to"
        "\nmeasure third order sensitivity of the value of portfolio to changes in the"
        "\nunderlying asset price of European option in portfolio while holding the"
        "\nother parameters fixed"
        "\n");*/
    BOOST_TEST_MESSAGE("Testing Speed calculations for Europeam Option Portfolio...");

#ifdef CL_TAPE_CPPAD

    cleanFolder(Speed);
    recordDependence<stock, Speed>(40, 180, "", "");
    return recordPerformance(Speed);

#endif
    return true;
}


bool AdjointEuropeanOptionPortfolioTest::testUltimaCallPortfolio()
{
    /*BOOST_TEST_MESSAGE(
        "\nUltima test"
        "\nPortfolio of European Call options with different parameters is examined to"
        "\nmeasure third order sensitivity of the value of portfolio to changes in the"
        "\nvolatility of European option in portfolio while holding the other parameters"
        "\nfixed"
        "\n");*/
    BOOST_TEST_MESSAGE("Testing Ultima calculations for Europeam Option Portfolio...");

#ifdef CL_TAPE_CPPAD

    cleanFolder(Ultima);
    recordDependence<sigma, Ultima>(0.015, 0.2, "", "");
    return recordPerformance(Ultima);

#endif
    return true;
}


bool AdjointEuropeanOptionPortfolioTest::testZommaCallPortfolio()
{
    /*BOOST_TEST_MESSAGE(
        "\nZomma test"
        "\nPortfolio of European Call options with different parameters is examined to"
        "\nmeasure sensitivity of the Gamma greek (second order sensitivity of the value"
        "\nof portfolio to changes in the underlying asset price) to changes in the"
        "\nvolatility of European option in portfolio while holding the other parameters"
        "\nfixed"
        "\n");*/
    BOOST_TEST_MESSAGE("Testing Zomma calculations for Europeam Option Portfolio...");

#ifdef CL_TAPE_CPPAD

    cleanFolder(Zomma);
    recordDependence<stock, Zomma>(10  , 200 , "", "");
    recordDependence<sigma, Zomma>(0.03, 0.15, "", "");
    return recordPerformance(Zomma);

#endif
    return true;
}


bool AdjointEuropeanOptionPortfolioTest::testColorCallPortfolio()
{
    /*BOOST_TEST_MESSAGE(
        "\nColor test"
        "\nPortfolio of European Call options with different parameters is examined to"
        "\nmeasure the rate of change of Gamma (second order sensitivity of the value of"
        "\nportfolio to changes in the underlying asset price) over the passage of time"
        "\n");*/
    BOOST_TEST_MESSAGE("Testing Color calculations for Europeam Option Portfolio...");

#ifdef CL_TAPE_CPPAD

    cleanFolder(Color);
    recordDependence<stock          , Color>(40   , 250, "", "");
    recordDependence<DataType::time, Color>(0.001, 0.5, "", "");
    return recordPerformance(Color);

#endif
    return true;
}


test_suite* AdjointEuropeanOptionPortfolioTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Adjoint with European option portfolio tests");
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testDeltaCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testVegaCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testThetaCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testRhoCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testGammaCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testVommaCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testVannaCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testCharmCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testVetaCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testVeraCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testSpeedCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testUltimaCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testZommaCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testColorCallPortfolio));
    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_european_option_portfolio)

BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioDelta)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testDeltaCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioVega)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testVegaCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioTheta)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testThetaCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioRho)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testRhoCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioGamma)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testGammaCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioVomma)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testVommaCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioVanna)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testVannaCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioCharm)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testCharmCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioVeta)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testVetaCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioVera)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testVeraCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioSpeed)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testSpeedCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioUltima)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testUltimaCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioZomma)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testZommaCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioColor)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testColorCallPortfolio());
}

BOOST_AUTO_TEST_SUITE_END()

#endif