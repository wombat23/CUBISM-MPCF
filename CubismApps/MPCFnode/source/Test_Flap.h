/*
 *  Test_ShockBubble.h
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 6/16/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "Test_SteadyState.h"



class Test_Flap: public Test_SteadyState
{
    void _ic(FluidGrid& grid);

protected:
	Real pCrit,    pInit;
	Real pAmbient, tAmbient;
    Real TInit;
    Real tSM, tSE, iI, iU;
    Real arcX,arcY,arcZ, arcWidth, arcHeight;
    Real zetaGrid, threshP, flRho, flS, flL;
		Real gamma, R_star;


    void _setup_constants();
    void _dumpStatistics(FluidGrid& grid, const int counter, const Real t, const Real dt);
    void _analysis(FluidGrid& grid, const int stepid);

public:
    InputStructVals inputStructVals;

    Test_Flap(const int argc, const char ** argv): Test_SteadyState(argc, argv) { }
    
	void run();
	void setup();
};

static Real phi=0, dphidt=0;
static Real t_current=0;

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabFlap: public BlockLab<BlockType,allocator>
{		
	typedef typename BlockType::ElementType ElementTypeBlock;
    
protected:
	bool is_xperiodic() {return false;}
	bool is_yperiodic() {return false;}
	bool is_zperiodic() {return false;}

public:
	BlockLabFlap(): BlockLab<BlockType,allocator>(){}
	
	void _apply_bc(const BlockInfo& info, const Real t=0)
	{	
        const Real dt = t-t_current;
        t_current = t;

        const Real pAvg = 1;
        const Real rAvg = 1;
        const Real uAvg = 0;
        const Real gamma = 1.4;

        BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);
        
        if (info.index[0]==0)           bc.template applyBC_reflecting<0,0>();
	if (info.index[0]==this->NX-1)  bc.template applyBC_reflecting<0,1>();
//if (info.index[0]==this->NX-1)  bc.template applyBC_flap_closed<0,1>();//bc.template applyBC_reflecting<0,1>();	
//        if (info.index[0]==this->NX-1)  bc.template applyBC_flap<0,1>(pAvg,rAvg,uAvg,dt,gamma,phi,dphidt);
        if (info.index[1]==0)			bc.template applyBC_reflecting<1,0>();
        if (info.index[1]==this->NY-1)	bc.template applyBC_reflecting<1,1>();
        if (info.index[2]==0)			bc.template applyBC_reflecting<2,0>();
        if (info.index[2]==this->NZ-1)	bc.template applyBC_reflecting<2,1>();
        
        
        
        }
};
