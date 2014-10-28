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
	
	const Real a, dtinvh; //LSRK3-related "a" factor, and "lambda"
	
	//the only constructor for this class
	Flap_average(const Real a, const Real dtinvh);
	//main method of the class, it evaluates the convection term of the RHS
	

    void compute(const Real * const src, const int gptfloats, Real & pAvg_local, Real & rAvg_local, Real & uAvg_local, Real & tAvg_local, int & N_local, const int nblocksx, int bx);	
	
    
};
