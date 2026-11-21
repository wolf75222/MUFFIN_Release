//
//  Limiter.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef Limiter_hpp
#define Limiter_hpp

#include <stdio.h>
#include <bit>
#include <cstdint>
#include "../../Component.hpp"

class Limiter : public Component
{
public:
    Limiter(string name) : Component(name) {}
    virtual ~Limiter() {}
    virtual double operator()(const double& theta) = 0;
};

// c ? x : y by a mask on the bits, without branch: GCC does not turn a floating-point ?: into a select in a loop
// (the comparison may trap), so a loop that has one is not vectorised.
inline double selectValue(bool c, double x, double y)
{
    const std::uint64_t m = -static_cast<std::uint64_t>(c);
    return std::bit_cast<double>((m & std::bit_cast<std::uint64_t>(x)) | (~m & std::bit_cast<std::uint64_t>(y)));
}
// std::min and std::max: the same comparison and the same result, NaN and signed zeros included.
inline double minValue(double a, double b) {return selectValue(b < a, b, a);}
inline double maxValue(double a, double b) {return selectValue(a < b, b, a);}

#endif /* Limiter_hpp */
