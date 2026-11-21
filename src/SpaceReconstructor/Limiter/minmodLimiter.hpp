//
//  minmodLimiter.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 12/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef minmodLimiter_hpp
#define minmodLimiter_hpp

#include <stdio.h>
#include <algorithm>
#include "Limiter.hpp"

class minmodLimiter : public Limiter
{
public:
    minmodLimiter(string name) : Limiter(name){}
    ~minmodLimiter() {}
    virtual double operator()(const double& theta);
    static double value(double theta) {return std::max(0., std::min(1.,theta));} // minmod limiter
};


#endif /* minmodLimiter_hpp */
