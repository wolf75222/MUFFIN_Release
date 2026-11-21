//
//  ThirdOrderLimiter.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 10/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef ThirdOrderLimiter_hpp
#define ThirdOrderLimiter_hpp

#include <stdio.h>
#include <algorithm>
#include "Limiter.hpp"

class ThirdOrderLimiter : public Limiter
{
public:
    ThirdOrderLimiter(string name) : Limiter(name){}
    ~ThirdOrderLimiter() {}
    virtual double operator()(const double& theta);
    static double value(double theta) {return std::max(0., std::min(2*theta,std::min((2 + theta)/3., 2.)));} // Cada and Torrilhon 2008 limiter
};

#endif /* ThirdOrderLimiter_hpp */
