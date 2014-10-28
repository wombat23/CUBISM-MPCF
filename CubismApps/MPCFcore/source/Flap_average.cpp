/*
 *  Convection_CPP.cpp
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 5/6/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <cmath>
#include <cassert>

#include <vector>
#include <algorithm>
#include <cstdio>
#include <iostream>

#include "Flap_average.h"

using namespace std;

Flap_average::Flap_average(const Real a, const Real dtinvh): a(a), dtinvh(dtinvh) { }

void Flap_average::compute(const Real * const src, const int gptfloats,  Real & pAvg_local, Real & rAvg_local, Real & uAvg_local, Real & tAvg_local, int & N_local, const int nblocksx, int bx)
{
        pAvg_local=0;
        rAvg_local=0;
        uAvg_local=0;
        tAvg_local=0;
        N_local=0;

        Real gamma = 1.4;
        Real R = 8.3;

        const int N = _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*gptfloats;

        if (bx == nblocksx-1)
            for(int i=0; i<N; i+=gptfloats)
                {
                    int icell = i/gptfloats;
                    int ix = icell%_BLOCKSIZE_;

                    if (ix==_BLOCKSIZE_-1)
                    {
                        const Real rho = src[i];
                        const Real ru  = src[i+1];
                        const Real e   = src[i+4];
                    rAvg_local += rho;
                    uAvg_local += ru/rho;
                    Real p = (e - .5*ru*ru/rho)/gamma;
                    pAvg_local += p;
                    tAvg_local += p/(R*rho);
                    }
                } 
}
