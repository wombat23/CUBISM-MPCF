#pragma once

#include <cstdio>

#include "common.h"
#include "SOA2D.h"
#include "Types.h"

class SourceTerm
{
public:
    void compute(
            const Real *const pos , Real h_gridpoint,
            Real * src, const int gptfloats, Real t,
            InputStructVals inputStructVals
    );
};