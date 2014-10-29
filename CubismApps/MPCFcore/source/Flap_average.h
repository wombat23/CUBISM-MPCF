/*
 *  Flap_average.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 5/6/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include "common.h"
#include "SOA2D.h"

class Flap_average
{
public:
	Flap_average();
	
    void compute(const Real * const src, const int gptfloats, Real & pAvg_local, Real & rAvg_local, Real & uAvg_local, Real & tAvg_local, int & N_local, const int nblocksx, int bx);	
	
    
};
