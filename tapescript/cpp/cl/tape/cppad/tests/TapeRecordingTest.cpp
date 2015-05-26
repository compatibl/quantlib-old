﻿/*
Copyright (C) 2003-2015 CompatibL

This file is part of TapeScript, an open source library and tape encoding
standard for adjoint algorithmic differentiation (AAD), available from

    http://github.com/compatibl/tapescript (source)
    http://tapescript.org (documentation)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include <boost/test/unit_test.hpp>
#include <boost/timer.hpp>
#include <cl/tape/tape.hpp>
#include <cl/tape/util/testutil.hpp>

static int dimension = (int)1e+7;

static std::vector<double> doubleVector = getRandomVector(dimension);

// Return time of tape recording for AD<double>.
double ADTapePerformance(std::ofstream & log)
{
    std::vector<CppAD::AD<double>> ADvec = vectorCast<CppAD::AD<double>>(getRandomVector(dimension));

    CppAD::AD<double> ADResult = 0;

    CppAD::Independent(ADvec);
    log << "Start of tape recording : " << currentTime() << std::endl;
    boost::timer timer;

    for (CppAD::AD<double> item : ADvec)
    {
        ADResult += CppAD::sin(item * CppAD::pow(ADResult, CppAD::sqrt(CppAD::pow(item, 2.0) * CppAD::exp(CppAD::cos(ADResult * CppAD::tan(item))))));
    }
    CppAD::ADFun<double> f(ADvec, std::vector<CppAD::AD<double>>({ ADResult }));
    double time = timer.elapsed();
    log << "End of tape recording" << std::endl;
    return time;
}

// Return time of tape recording for TapeDouble.
double TapeDoubleTapePerformance(std::ofstream & log)
{
    std::vector<cl::TapeDouble> CppVec = vectorCast<cl::TapeDouble>(getRandomVector(dimension));

    cl::TapeDouble CppResult = 0;

    cl::Independent(CppVec);
    log << "Start of tape recording : " << currentTime() << std::endl;
    boost::timer timer;

    for (cl::TapeDouble item : CppVec)
    {
        CppResult += std::sin(item * std::pow(CppResult, std::sqrt(std::pow(item, 2.0) * std::exp(std::cos(CppResult * std::tan(item))))));
    }
    cl::TapeFunction<double> f(CppVec, std::vector<cl::TapeDouble>({ CppResult }));
    double time = timer.elapsed();
    log << "End of tape recording" << std::endl;
    return time;
}

BOOST_AUTO_TEST_SUITE(TapeRecordingTest)

    BOOST_AUTO_TEST_CASE(TapeRecordingPerformanceTest)
{
        std::ofstream out("..\\..\\..\\..\\..\\..\\..\\output\\en-us\\Tape\\CppAD\\Tests\\TapeRecordingTest\\Log.csv",
            std::ofstream::out | std::ofstream::app);
        out << std::endl << "Independent variable vectors size " << dimension << std::endl;
        double ADtime = ADTapePerformance(out);
        double CppTime = TapeDoubleTapePerformance(out);
        out << "\tTime for AD<double> " << ADtime << " s" << std::endl;
        out << "\tTime for TapeDouble " << CppTime << " s" << std::endl;
        out << "\tThe relative difference  " << 1.0 * std::abs(CppTime - ADtime) / std::max(CppTime, ADtime)
            << std::endl ;
}

BOOST_AUTO_TEST_SUITE_END()
