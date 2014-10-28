/*
 *  Test_ShockBubble.cpp
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 6/16/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#include <sstream>
#include <algorithm>

#ifdef _USE_NUMA_
#include <numa.h>
#include <omp.h>
#endif

#include <Profiler.h>

#include "Test_Flap.h"
#include "Tests.h"

using namespace std;

void Test_Flap::_ic(FluidGrid& grid)
{
	cout << "Flap Initial condition..." ;
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    const double G1 = Simulation_Environment::GAMMA1-1;
    const double G2 = Simulation_Environment::GAMMA2-1;
    const double F1 = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
    const double F2 = Simulation_Environment::GAMMA2*Simulation_Environment::PC2;
    printf("----> %f",pInit / TInit / R_star);
#pragma omp parallel
	{	
#ifdef _USE_NUMA_
		const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
		const int mynode = omp_get_thread_num() / cores_per_node;
		numa_run_on_node(mynode);
#endif

#pragma omp for
		for(int i=0; i<(int)vInfo.size(); i++)
		{
            BlockInfo info = vInfo[i];
            FluidBlock& b = *(FluidBlock*)info.ptrBlock;
            
            for(int iz=0; iz<FluidBlock::sizeZ; iz++)
                for(int iy=0; iy<FluidBlock::sizeY; iy++)
                    for(int ix=0; ix<FluidBlock::sizeX; ix++)
                    {
                        Real p[3];
                        info.pos(p, ix, iy, iz);


												b(ix, iy, iz).rho      = pInit / TInit / R_star;
												b(ix, iy, iz).u        = 0.0;                        
												b(ix, iy, iz).v        = 0.0;
												b(ix, iy, iz).w        = 0.0;
												b(ix, iy, iz).energy   = 0.5*(b(ix, iy, iz).u*b(ix, iy, iz).u
															+ b(ix, iy, iz).v*b(ix, iy, iz).v
															+ b(ix, iy, iz).w*b(ix, iy, iz).w)/b(ix, iy, iz).rho + (gamma-1)*pInit;
												b(ix, iy, iz).G        = 1.0/(gamma-1.0);
												b(ix, iy, iz).P        = 0.0;
                    }
        }		
	}	
	
	cout << "done." << endl;
}

void Test_Flap::_dumpStatistics(FluidGrid& grid, const int step_id, const Real t, const Real dt)
{
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    double rInt=0., uInt=0., vInt=0., wInt=0., eInt=0., vol=0., ke=0., r2Int=0., mach_max=-HUGE_VAL, p_max=-HUGE_VAL;
    const double h = vInfo[0].h_gridpoint;
    const double h3 = h*h*h;
    double wall_p_max=-HUGE_VAL;
    
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
				{
                    rInt += b(ix, iy, iz).rho;
                    uInt += b(ix, iy, iz).u;
                    vInt += b(ix, iy, iz).v;
                    wInt += b(ix, iy, iz).w;
                    eInt += b(ix, iy, iz).energy;
                    vol  += b(ix,iy,iz).G>0.5*(1/(Simulation_Environment::GAMMA1-1)+1/(Simulation_Environment::GAMMA2-1))? 1:0;
                    r2Int += b(ix, iy, iz).rho*(1-min(max((b(ix,iy,iz).G-1/(Simulation_Environment::GAMMA1-2))/(1/(Simulation_Environment::GAMMA1-1)-1/(Simulation_Environment::GAMMA2-1)),(Real)0),(Real)1));
                    ke   += 0.5/b(ix, iy, iz).rho * (b(ix, iy, iz).u*b(ix, iy, iz).u+b(ix, iy, iz).v*b(ix, iy, iz).v+b(ix, iy, iz).w*b(ix, iy, iz).w);
     
                    const double pressure = (b(ix, iy, iz).energy - 0.5/b(ix, iy, iz).rho * (b(ix, iy, iz).u*b(ix, iy, iz).u+b(ix, iy, iz).v*b(ix, iy, iz).v+b(ix, iy, iz).w*b(ix, iy, iz).w) - b(ix, iy, iz).P)/b(ix,iy,iz).G;
                    const double c = sqrt((1/b(ix,iy,iz).G+1)*(pressure+b(ix,iy,iz).P/b(ix,iy,iz).G/(1/b(ix,iy,iz).G+1))/b(ix, iy, iz).rho);                    
                    const double velmag = sqrt(b(ix, iy, iz).u*b(ix, iy, iz).u+b(ix, iy, iz).v*b(ix, iy, iz).v+b(ix, iy, iz).w*b(ix, iy, iz).w)/b(ix, iy, iz).rho;
                    
                    mach_max = max(mach_max, velmag/c);
                    p_max = max(p_max, pressure);
                    
                    if (info.index[2]==0 && iz==0)
                        wall_p_max = max(wall_p_max, pressure);
                }
    }
    
    FILE * f = fopen("integrals.dat", "a");
    fprintf(f, "%d %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", step_id, t, dt, rInt*h3, uInt*h3, 
            vInt*h3, wInt*h3, eInt*h3, vol*h3, ke*h3, r2Int*h3, mach_max, p_max, pow(0.75*vol*h3/M_PI,1./3.), wall_p_max);
    fclose(f);
    
    FILE * f2 = fopen("centerline_velocities.dat", "a");
    FILE * f3 = fopen("p_center.dat", "a");
    FILE * f4 = fopen("interface.dat","a");
    
    Lab lab;
    const int ss[3]={0,0,0};
    const int se[3]={2,2,2};
    lab.prepare(grid, ss, se, false);
    
    vector<pair<Real,Real> > velocities;
    vector<Real> iso_gamma;
    Real p_wall=0;
    
    Real x[3];
    
    for(int i=0; i<(int)vInfo.size(); i++)
    {        
        BlockInfo info = vInfo[i];
        
        if (info.index[1]!=grid.getBlocksPerDimension(1)/2) continue;
        
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
        lab.load(info);
        
        for(int iz=0; iz<FluidBlock::sizeZ; iz++)
            for(int iy=0; iy<FluidBlock::sizeY; iy++)
                for(int ix=0; ix<FluidBlock::sizeX; ix++)
                {
                    if (iy!=0 || iz!=0) continue;
                    
                    info.pos(x,ix,iy,iz);
                    
                    if((ix+1)==FluidBlock::sizeX && (info.index[0]+1)==grid.getBlocksPerDimension(0))
                    {    
                        const double ke = 0.5*(pow(b(ix, iy, iz).u,2)+pow(b(ix, iy, iz).v,2)+pow(b(ix, iy, iz).w,2))/b(ix, iy, iz).rho;

                        p_wall = (b(ix, iy, iz).energy - ke -  b(ix, iy, iz).P)/b(ix, iy, iz).G;
                    }
                    
                    if (abs(b(ix,iy,iz).G-2.35)<0.1) {
                        iso_gamma.push_back(x[0]);
                    }
                    
                    if( Simulation_Environment::heaviside(lab(ix, iy, iz).G-0.5*(1/(Simulation_Environment::GAMMA2-1)+1/(Simulation_Environment::GAMMA1-1)))*                        Simulation_Environment::heaviside(lab(ix+1, iy, iz).G-0.5*(1/(Simulation_Environment::GAMMA2-1)+1/(Simulation_Environment::GAMMA1-1))) == 0 &&
                       !(Simulation_Environment::heaviside(lab(ix, iy, iz).G-0.5*(1/(Simulation_Environment::GAMMA2-1)+1/(Simulation_Environment::GAMMA1-1)))==0 && 
                         Simulation_Environment::heaviside(lab(ix+1, iy, iz).G-0.5*(1/(Simulation_Environment::GAMMA2-1)+1/(Simulation_Environment::GAMMA1-1)))==0 ) ) 
                    {                   
                        const double pressure = (b(ix, iy, iz).energy - 0.5/b(ix, iy, iz).rho * (b(ix, iy, iz).u*b(ix, iy, iz).u+b(ix, iy, iz).v*b(ix, iy, iz).v+b(ix, iy, iz).w*b(ix, iy, iz).w) - b(ix, iy, iz).P)/b(ix,iy,iz).G;

                        velocities.push_back(pair<Real,Real>(sqrt(6.59*b(ix,iy,iz).rho*pressure),b(ix, iy, iz).u/b(ix, iy, iz).rho));
                    }
                }
    }
    
    std::sort(velocities.begin(), velocities.end(), sort_pred());
    std::sort(iso_gamma.begin(),iso_gamma.end());

    if (velocities.size()>0) fprintf(f2, "%d %f %e %e %f\n", step_id, t, velocities.front().second, velocities.back().second, velocities.front().first);
    if (iso_gamma.size()>0) fprintf(f4, "%d %e %f %f\n", step_id, t, iso_gamma.front(), iso_gamma.back());

    fprintf(f3, "%d %e %e\n", step_id, t, p_wall);
    
    fclose(f4);
    fclose(f3);
    fclose(f2);
}

void Test_Flap::_analysis(FluidGrid& grid, const int step_id)
{
    std::stringstream streamer;
    streamer<<"centerline_pressure-"<<step_id<<".dat";
    FILE * f = fopen(streamer.str().c_str(), "w");
    streamer.str("");
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    Real x[3];
    
    for(int i=0; i<(int)vInfo.size(); i++)
    {        
        BlockInfo info = vInfo[i];
        
        if (info.index[1]!=grid.getBlocksPerDimension(1)/2) continue;
        
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
        
        for(int iz=0; iz<FluidBlock::sizeZ; iz++)
            for(int iy=0; iy<FluidBlock::sizeY; iy++)
                for(int ix=0; ix<FluidBlock::sizeX; ix++)
                {
                    if (iy!=0 || iz!=0) continue;
                    
                    info.pos(x,ix,iy,iz);
                    
                    const double ke = 0.5*(pow(b(ix, iy, iz).u,2)+pow(b(ix, iy, iz).v,2)+pow(b(ix, iy, iz).w,2))/b(ix, iy, iz).rho;
                    
                    const double pressure = (b(ix, iy, iz).energy - ke -  b(ix, iy, iz).P)/b(ix, iy, iz).G;

                    fprintf(f, "%e %e\n", x[0], pressure);
                }
    }    
    
    fclose(f); 
    
    streamer<<"wall_pressure-"<<step_id<<".dat";
    f = fopen(streamer.str().c_str(), "w");
    streamer.str("");
    
    for(int i=0; i<(int)vInfo.size(); i++)
    {        
        BlockInfo info = vInfo[i];
        
        if (info.index[0]!=grid.getBlocksPerDimension(0)-1 || info.index[2]!=grid.getBlocksPerDimension(2)/2) continue;
        
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
        
        for(int iz=0; iz<FluidBlock::sizeZ; iz++)
            for(int iy=0; iy<FluidBlock::sizeY; iy++)
                for(int ix=0; ix<FluidBlock::sizeX; ix++)
                {
                    if (ix!=FluidBlock::sizeX-1 || iz!=0) continue;
                    
                    info.pos(x,ix,iy,iz);
                    
                    const double ke = 0.5*(pow(b(ix, iy, iz).u,2)+pow(b(ix, iy, iz).v,2)+pow(b(ix, iy, iz).w,2))/b(ix, iy, iz).rho;
                    const double pressure = (b(ix, iy, iz).energy - ke -  b(ix, iy, iz).P)/b(ix, iy, iz).G;
                    
                    fprintf(f, "%e %e\n", x[1], pressure);
                }
    }
    
    fclose(f);
}

void Test_Flap::run()
{	
    Real dt=0;
	bool bLoop = (NSTEPS>0) ? (step_id<NSTEPS) : (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1);
    
	while (bLoop)
	{
		cout << "time is " << t << endl;
		cout << "step_id is " << step_id << endl;
        
		if(step_id%DUMPPERIOD == 0)
		{
			profiler.push_start("DUMP");            
			std::stringstream streamer;
			streamer<<"data-"<<step_id;
			_dump(streamer.str());
			_vp(*grid);			
			profiler.pop_stop();
        }
        if(step_id%ANALYSISPERIOD == 0)
        {
            profiler.push_start("DUMP ANALYSIS");
            _analysis(*grid, step_id);
            profiler.pop_stop();
		}
        
		if (step_id%SAVEPERIOD == 0) _save();
        
        if(step_id%10 == 0)
        {
            profiler.push_start("DUMP STATISTICS");
            _dumpStatistics(*grid, step_id, t, dt);
            profiler.pop_stop();
        }
        
		profiler.push_start("EVOLVE");

        stepper->set_current_time(t);
		dt = (*stepper)(TEND-t);

		profiler.pop_stop();
		
		if(step_id%10 == 0)
			profiler.printSummary();			
                
		t+=dt;
		step_id++;
		bLoop = (NSTEPS>0) ? (step_id<NSTEPS) : (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1);
        if (dt==0) break;
	}
    
    cout << "Finishing run ...";
        
	std::stringstream streamer;
	streamer<<"data-"<<step_id;
	if (DUMPPERIOD < 1e5) _dump(streamer.str());
    
    cout << "done" << endl;
}

void Test_Flap::_setup_constants()
{
    parser.mute();
    
    bRESTART = parser("-restart").asBool();
    
    parser.set_strict_mode();
    TEND = parser("-tend").asDouble();
    DUMPPERIOD = parser("-dumpperiod").asInt();
    SAVEPERIOD = parser("-saveperiod").asInt();
    CFL = parser("-cfl").asDouble();
    BPDX = parser("-bpdx").asInt();
    
    parser.unset_strict_mode();
    
    pInit = parser("-pInit").asDouble(100000);
    pCrit = parser("-pCrit").asDouble(6000);
    TInit = parser("-tInit").asDouble(300);
    
    tSM  = parser("-tSM").asDouble(0.015);
    tSE  = parser("-tSE").asDouble(0.015);
    
    iI   = parser("-iI").asDouble(1250000);
    iU   = parser("-iU").asDouble(240);
    
    arcX = parser("-arcX").asDouble(0.1);
    arcY = parser("-arcY").asDouble(arcX);
    arcZ = parser("-arcZ").asDouble(arcX);
    
    arcWidth = parser("-arcWidth").asDouble(0.1);
    arcHeight = parser("-arcHeight").asDouble(arcWidth);
    
    zetaGrid = parser("-zetaGrid").asDouble(0.219);
    threshP  = parser("-threshP").asDouble(6000);
    
    flRho = parser("-flRho").asDouble(3500);
    flS   = parser("-flS").asDouble(0.022);
    flL   = parser("-flL").asDouble(0.850);

		gamma   = parser("-gamma").asDouble(1.4);
		R_star  = parser("-Rstar").asDouble(287); // k_Boltzman/molecular weight of air

    bASCIIFILES = parser("-ascii").asBool(false);

    BPDY = parser("-bpdy").asInt(BPDX);
    BPDZ = parser("-bpdz").asInt(BPDX);
    
    Simulation_Environment::GAMMA1 = parser("-g1").asDouble(1.4);
    Simulation_Environment::GAMMA2 = parser("-g2").asDouble(1.4);
    
    bVP = parser("-vp").asBool(0);
    VERBOSITY = parser("-verb").asInt(0);
    NSTEPS = parser("-nsteps").asInt(0);
    bAWK = parser("-awk").asBool(false);
    ANALYSISPERIOD = parser("-analysisperiod").asInt(std::numeric_limits<int>::max());
    
    assert(TEND >= 0.0);
    assert(BPDX >= 1);
    assert(BPDY >= 1);
    assert(BPDZ >= 1);
    assert(DUMPPERIOD > 0);
    assert(CFL > 0 && CFL<1);
    
    Simulation_Environment::PC1 = parser("-pc1").asDouble(0);
    Simulation_Environment::PC2 = parser("-pc2").asDouble(0);
    REPORT_FREQ = parser("-report").asInt(10);
}

void Test_Flap::setup()
{
	printf("////////////////////////////////////////////////////////////\n");
	printf("////////////         TEST FLAP       ///////////////\n");
	printf("////////////////////////////////////////////////////////////\n");
	
	_setup_constants();
	parser.mute();
	
	if (parser("-morton").asBool(0))
		grid = new GridMorton<FluidGrid>(BPDX, BPDY, BPDZ);
	else
		grid = new FluidGrid(BPDX, BPDY, BPDZ);
	
	assert(grid != NULL);
	
	stepper = new FlowStep_LSRK3(*grid, CFL, Simulation_Environment::GAMMA1, Simulation_Environment::GAMMA2, parser, VERBOSITY, &profiler, Simulation_Environment::PC1, Simulation_Environment::PC2, bAWK);
	
	if(bRESTART)
	{
		_restart();
		_dump("restartedcondition.vti");
	}
	else
		_ic(*grid);
}
