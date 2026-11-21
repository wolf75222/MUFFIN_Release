//
//  ThirdOrderLimiter.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 10/10/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "ThirdOrderLimiter.hpp"

// Self-register with the factory
#include "LimiterRegistrar.hpp"
REGISTER_LIMITER("thirdOrder", ThirdOrderLimiter);
#include <cmath>

double ThirdOrderLimiter::operator()(const double& theta) {
    
    return value(theta);
}
