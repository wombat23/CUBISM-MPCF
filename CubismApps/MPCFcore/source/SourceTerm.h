#pragma once

#include <cstdio>

#include "common.h"
#include "SOA2D.h"

class SourceTerm
{
public:
    void compute(
                 const Real *const pos ,
                 const Real * const src, const int gptfloats, Real t,
                 Real tSE, Real tSM,
                 Real iI, Real iU,
                 Real arcX, Real arcY, Real arcZ,
                 Real arcW, Real arcH
                 );
};