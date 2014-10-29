/*
 *  FlowStep_LSRK3.cpp
 *  MPCFnode
 *
 *  Created by Diego Rossinelli on 6/15/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#include <vector>

#include <utility>
#include <iostream>
#include <limits>

#include <omp.h>

#include <Timer.h>
#include <Profiler.h>
#include <Convection_CPP.h>

#if defined(_QPX_) || defined(_QPXEMU_)
#include <Convection_QPX.h>
#include <Update_QPX.h>
#include <MaxSpeedOfSound_QPX.h>
#ifdef _QPX_
//here are for cycle counting
#include <ucontext.h>
#include <signal.h>
#include <sys/time.h>
#include <errno.h>
#include "spi/include/upci/upci.h"
#endif
#endif

#include <Update.h>
#include <MaxSpeedOfSound.h>
#include <Flap_average.h>

#ifdef _USE_HPM_
#include <mpi.h>
extern "C" void HPM_Start(char *);
extern "C" void HPM_Stop(char *);
#else
#define HPM_Start(x)
#define HPM_Stop(x)
#endif

using namespace std;

#include "FlowStep_LSRK3.h"
#include "Tests.h"

namespace LSRK3data
{
	Real gamma1 = -1;
	Real gamma2 = -1;
	Real smoothlength = -1;
	Real pc1 = 0;
    Real pc2 = 0;
   
	int verbosity;
	int  NCORES;
	float PEAKPERF_CORE, PEAKBAND;
	
	string dispatcher;
	
	int step_id = 0;
    int ReportFreq = 1;
}

template<typename Lab, typename Kernel>
void _process(const Real a, const Real dtinvh, vector<BlockInfo>& myInfo, FluidGrid& grid, const Real t=0, bool tensorial=false)
{
	const int stencil_start[3] = {-3,-3,-3};
	const int stencil_end[3] = {4,4,4};
	BlockInfo * ary = &myInfo.front();
	const int N = myInfo.size();
	
	const int NTH = omp_get_max_threads();
	double total_time[NTH];

	static Lab * labs = NULL;

	if (labs == NULL)
	{
		printf("allocating %d labs\n", NTH); fflush(0);

		labs = new Lab[NTH];
		for(int i = 0; i < NTH; ++i)
			labs[i].prepare(grid, stencil_start, stencil_end, tensorial);
	}
	
#pragma omp parallel
	{
#ifdef _USE_NUMA_
        const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
        const int mynode = omp_get_thread_num() / cores_per_node;
	numa_run_on_node(mynode);
#endif
		
		const int tid = omp_get_thread_num();
		total_time[tid] = 0;
		
		Timer timer;
		Kernel kernel(a, dtinvh);
		
		Lab& mylab = labs[tid];
//		Lab mylab;
//		mylab.prepare(grid, stencil_start, stencil_end, tensorial);

#pragma omp for schedule(runtime)
		for(int i=0; i<N; i++)
		{
			//we want to measure the time spent in ghost reconstruction
			timer.start();
			mylab.load(ary[i], t);
			total_time[tid] += timer.stop();
			
            const Real * const srcfirst = &mylab(-3,-3,-3).rho;
            const int labSizeRow = mylab.template getActualSize<0>();
            const int labSizeSlice = labSizeRow*mylab.template getActualSize<1>();
			
			Real * const destfirst = &((FluidBlock*)ary[i].ptrBlock)->tmp[0][0][0][0];
			
			kernel.compute(srcfirst, FluidBlock::gptfloats, labSizeRow, labSizeSlice, 
						   destfirst, FluidBlock::gptfloats, FluidBlock::sizeX, FluidBlock::sizeX*FluidBlock::sizeY);
		}
		
#pragma omp single
		{
			double min_val = total_time[0], max_val = total_time[0], sum = total_time[0];
			
			for(int i=1; i<NTH; i++)
			{
				min_val = min(min_val, total_time[i]);
				max_val = max(max_val, total_time[i]);
				sum += total_time[i];
			}
			
			if (LSRK3data::verbosity >= 1)
				printf("(min,max,avg) of lab.load() is (%5.10e, %5.10e, %5.10e)\n", min_val, max_val, sum/NTH);
		}
	}
}

template < typename TSOS>
Real _computeSOS_OMP(FluidGrid& grid,  bool bAwk)
{
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    const int N = vInfo.size();
    const BlockInfo * const ary = &vInfo.front();

#if (_OPENMP < 201107)   
    Real * tmp = NULL;
    int error = posix_memalign((void**)&tmp, std::max(8, _ALIGNBYTES_), sizeof(Real) * N);
    assert(error == 0);
    
    Real * const local_sos = tmp;
    
#ifdef _USE_HPM_
    HPM_Start("dt_for");
#endif
	
#pragma omp parallel
    {
        
#ifdef _USE_NUMA_
		const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
        const int mynode = omp_get_thread_num() / cores_per_node;
        numa_run_on_node(mynode);
#endif
        
        TSOS kernel;
#pragma omp for schedule(runtime)
        for (int i=0; i<N; ++i)
        {
            FluidBlock & block = *(FluidBlock *)ary[i].ptrBlock;
            local_sos[i] =  kernel.compute(&block.data[0][0][0].rho, FluidBlock::gptfloats);
        }
    }
	
#ifdef _USE_HPM_
    HPM_Stop("dt_for");
#endif
	
#ifdef _USE_HPM_
	HPM_Start("dt_reduce");
#endif
    
#if defined(_QPX_) 
//|| defined(_QPXEMU_)    
	const int N8 = 8 * (N / 8);
    vector4double sos4A = vec_splats(0);
    vector4double sos4B = vec_splats(0);
	
    for(int i = 0; i < N8; i += 8)   
	{
		const vector4double candidate0 = vec_lda(0L, local_sos+i);
		const vector4double candidate1 = vec_lda(sizeof(Real) * 4, local_sos+i);
		
		const vector4double flag0 = vec_cmplt(sos4A, candidate0);
		const vector4double flag1 = vec_cmplt(sos4B, candidate1);
		
		sos4A = vec_sel(sos4A, candidate0, flag0);
		sos4B = vec_sel(sos4B, candidate1, flag1);
	}
    sos4A = mymax(sos4A, sos4B);
    sos4A = mymax(sos4A, vec_perm(sos4A, sos4A, vec_gpci(2323)));
	
    Real global_sos = std::max(vec_extract(sos4A, 0), vec_extract(sos4A, 1));
	
    for(int i=N8; i<N; ++i)
		global_sos = max(global_sos, local_sos[i]);
#else
	
    Real global_sos = local_sos[0];
    
#pragma omp parallel
    {
	    Real mymax = local_sos[0];
		
#pragma omp for schedule(runtime)
	    for(int i=0; i<N; ++i)
		    mymax = max(local_sos[i], mymax);

#pragma omp critical
	    {
		    global_sos = max(global_sos, mymax);
	    }
    }
#endif
	
#ifdef _USE_HPM_
	HPM_Stop("dt_reduce");
#endif
	
    free(tmp);
    
#else   /* FUSED_MAXSOS */

    Real global_sos = 0;
#ifdef _USE_HPM_
    HPM_Start("dt_for_reduce");
#endif

#pragma omp parallel
    {
      TSOS kernel;
        Real mymax = 0;

#pragma omp for schedule(runtime) reduction(max:global_sos)
        for (size_t i=0; i<N; ++i)
        {
            FluidBlock & block = *(FluidBlock *)ary[i].ptrBlock;
            global_sos =  kernel.compute(&block.data[0][0][0].rho, FluidBlock::gptfloats);
        }
    }

#ifdef _USE_HPM_
        HPM_Stop("dt_for_reduce");
#endif

#endif
	return global_sos;
}

Real FlowStep_LSRK3::_computeSOS(bool bAwk)
{
    Real sos = -1;
	
    const string kernels = parser("-kernels").asString("cpp");
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
	Timer timer;
    
	timer.start();
#if defined(_QPX_) || defined(_QPXEMU_)	
	if (kernels == "qpx")
		sos = _computeSOS_OMP<MaxSpeedOfSound_QPX>(grid,  bAwk);
	else
#endif
		sos = _computeSOS_OMP<MaxSpeedOfSound_CPP>(grid,  bAwk);
	
    const Real time = timer.stop();
    
    if (LSRK3data::verbosity >= 1 && LSRK3data::step_id % LSRK3data::ReportFreq == 0)
    {
		MaxSpeedOfSound_CPP::printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, 1, vInfo.size(), time);
        
        cout << "MAXSOS: " << time << "s (per substep), " << time/vInfo.size()*1e3 << " ms (per block)" << endl;
    }
    
    return sos;
}

template < typename TFLAP>
Real _computeFLAP_omp(FluidGrid& grid)
{
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    const int N = vInfo.size();
    const BlockInfo * const ary = &vInfo.front();

    Real global_rAvg = 0;
    Real global_pAvg = 0;
    Real global_uAvg = 0;
    Real global_tAvg = 0;
 	int global_n = 0;


#pragma omp parallel
    {
      TFLAP kernel;
      
#pragma omp for schedule(runtime) reduction(+:global_rAvg) reduction(+:global_pAvg) reduction(+:global_uAvg) reduction(+:global_tAvg) reduction(+:global_n) 

        for (size_t i=0; i<N; ++i)
        {
            FluidBlock & block = *(FluidBlock *)ary[i].ptrBlock;
            kernel.compute(&block.data[0][0][0].rho, FluidBlock::gptfloats, global_pAvg, global_rAvg, global_uAvg, global_tAvg, global_n, grid.getBlocksPerDimension(0),vInfo[i].index[0]);
        }
    }
    
    global_rAvg /= global_n;
    global_pAvg /= global_n;
    global_uAvg /= global_n;
    global_tAvg /= global_n;

	return global_pAvg;
}

Real FlowStep_LSRK3::_computeFLAP()
{
   Real pAvg = -1;

    vector<BlockInfo> vInfo = grid.getBlocksInfo();

	pAvg  = _computeFLAP_omp<Flap_average>(grid);

    return pAvg;
}

template<typename Kflow, typename Kupdate>
struct LSRKstep
{
    template<typename T>
    void check_error(const double tol, T ref[], T val[], const int N)
    {
        for(int i=0; i<N; ++i)
        {
            assert(!std::isnan(ref[i]));
            
            assert(!std::isnan(val[i]));
            
            const double err = ref[i] - val[i];
            
            const double relerr = err/std::max(1e-6, (double)std::max(fabs(val[i]), fabs(ref[i])));
            
            if (LSRK3data::verbosity>= 1) printf("+%1.1e,", relerr);
            
            if (fabs(relerr) >= tol && fabs(err) >= tol)
                printf("\n%d: %e %e -> %e %e\n", i, ref[i], val[i], err, relerr);
            
            //assert(fabs(relerr) < tol || fabs(err) < tol);
        }
        if (LSRK3data::verbosity >=1) printf("\t");
    }

    template<typename T>
    void check_error(const double tol, T ref, T val)
    {
        check_error(tol, &ref, &val, 1);
    }

    template<typename T>
    void check_error(const double tol, T ref[], T val[], const int N, const int bidx[], const int idx[], const int comp)
    {
        for(int i=0; i<N; ++i)
        {
            assert(!std::isnan(ref[i]));
            
            assert(!std::isnan(val[i]));
            
            const double err = fabs(ref[i]) - fabs(val[i]);
            
            const double relerr = err/std::max(1e-6, (double)std::max(fabs(val[i]), fabs(ref[i])));
            
            if (LSRK3data::verbosity>= 1) printf("+%1.1e,", relerr);
            
            if (fabs(relerr) >= tol && fabs(err) >= tol)
                printf("\n%d: %e %e -> %e %e, located at mirror block %d,%d,%d and point %d,%d,%d, comp %d \n", i, ref[i], val[i], err, relerr, bidx[0], bidx[1], bidx[2], idx[0], idx[1], idx[2], comp);
            
            
            //assert(fabs(relerr) < tol || fabs(err) < tol);
        }
        if (LSRK3data::verbosity >=1) printf("\t");
    }
    
    void _check_symmetry(FluidGrid& grid)
    {
        vector<BlockInfo> vInfo = grid.getBlocksInfo();
        
        for(int i=0; i<(int)vInfo.size(); i++)
        {
            BlockInfo info = vInfo[i];
            
            if (info.index[1]>grid.getBlocksPerDimension(1)/2-1) continue;
            
            FluidBlock& b = *(FluidBlock*)info.ptrBlock;
            
            const int bidx_mirror[3] = {info.index[0], grid.getBlocksPerDimension(1)-info.index[1]-1, info.index[2]};
            const int i_mirror = bidx_mirror[0] + (bidx_mirror[1] + bidx_mirror[2]*grid.getBlocksPerDimension(1))*grid.getBlocksPerDimension(0);
            assert(i_mirror<grid.getBlocksPerDimension(0)*grid.getBlocksPerDimension(1)*grid.getBlocksPerDimension(2) && i_mirror>=0);
            
            BlockInfo info_mirror = vInfo[i_mirror];
            FluidBlock& b_mirror = *(FluidBlock*)info_mirror.ptrBlock;
            
            for(int iz=0; iz<FluidBlock::sizeZ; iz++)
                for(int iy=0; iy<FluidBlock::sizeY; iy++)
                    for(int ix=0; ix<FluidBlock::sizeX; ix++)
                    {
                        const int idx_mirror[3] = {ix, FluidBlock::sizeY-iy-1, iz};
                        
                        check_error(std::numeric_limits<Real>::epsilon()*50, &b(ix,iy,iz).rho, &b_mirror(idx_mirror[0],idx_mirror[1],idx_mirror[2]).rho, 7, bidx_mirror, idx_mirror, 0);
//                        check_error(std::numeric_limits<Real>::epsilon()*50, &b.tmp[iz][iy][ix][0], &b_mirror.tmp[idx_mirror[2]][idx_mirror[1]][idx_mirror[0]][0], 7, bidx_mirror, idx_mirror, 1);
                    }
        }
        
        cout << "Symmetry check done" << endl;
    }
    
    LSRKstep(FluidGrid& grid, Real dtinvh, const Real current_time, bool bAwk)
    {
        vector<BlockInfo> vInfo = grid.getBlocksInfo();
        
        vector< vector<double> > timings;

        //_check_symmetry(grid);
        timings.push_back(step(grid, vInfo, 0      , 1./4, dtinvh, current_time));
        timings.push_back(step(grid, vInfo, -17./32, 8./9, dtinvh, current_time));
        timings.push_back(step(grid, vInfo, -32./27, 3./4, dtinvh, current_time));
        
        const double avg1 = ( timings[0][0] + timings[1][0] + timings[2][0] )/3;
        const double avg2 = ( timings[0][1] + timings[1][1] + timings[2][1] )/3;
        
        if (LSRK3data::verbosity >= 1 && LSRK3data::step_id % LSRK3data::ReportFreq == 0)
        {
            cout << "FLOWSTEP: " << avg1 << "s (per substep), " << avg1/vInfo.size()*1e3 << " ms (per block)" << endl;
            cout << "UPDATE: " << avg2 << "s (per substep), " << avg2/vInfo.size()*1e3 << " ms (per block)" << endl;
			//const float PEAKPERF_CORE, const float PEAKBAND, const int NCORES, const int NTIMES, const int NBLOCKS, const float MEASUREDTIME            
            Kflow::printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, 1, vInfo.size(), avg1);
            Kupdate::printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, 1, vInfo.size(), avg2);
        }
    }
    
    vector<double> step(FluidGrid& grid, vector<BlockInfo>& vInfo, Real a, Real b, Real dtinvh, const Real current_time)
    {
        Timer timer;
        vector<double> res;
        
        LSRK3data::FlowStep<Kflow, Lab> rhs(a, dtinvh);
	
#ifdef _USE_HPM_
	HPM_Start("RHS");	
#endif
        timer.start();     
        _process<Lab, Kflow>(a, dtinvh, vInfo, grid, current_time);
        const double t1 = timer.stop();
#ifdef _USE_HPM_
        HPM_Stop("RHS");
#endif
        LSRK3data::Update<Kupdate> update(b, &vInfo.front());
        
#ifdef _USE_HPM_
	HPM_Start("Update");
#endif
        timer.start();
        update.omp(vInfo.size());
        const double t2 = timer.stop();
#ifdef _USE_HPM_
        HPM_Stop("Update");
#endif
		res.push_back(t1);
		res.push_back(t2);
		
		return res;
    }
};

void FlowStep_LSRK3::set_constants()
{
    LSRK3data::gamma1 = gamma1;
    LSRK3data::gamma2 = gamma2;
    LSRK3data::smoothlength = smoothlength;
    LSRK3data::verbosity = verbosity;
    LSRK3data::PEAKBAND = PEAKBAND;
    LSRK3data::PEAKPERF_CORE = PEAKPERF_CORE;
    LSRK3data::NCORES = parser("-ncores").asInt(1);
    LSRK3data::pc1 = pc1;
    LSRK3data::pc2 = pc2;
    LSRK3data::dispatcher = blockdispatcher;
    LSRK3data::ReportFreq = parser("-report").asInt(20);
}

Real FlowStep_LSRK3::operator()(const Real max_dt)
{
    set_constants();
	
	Timer timer;
	timer.start();
#ifdef _USE_HPM_
	HPM_Start("dt");
#endif
    
    const Real maxSOS = _computeSOS(bAwk);
	
#ifdef _USE_HPM_
	HPM_Stop("dt");
#endif
	
	const double t_sos = timer.stop();
	
	cout << "sos take " << t_sos << " sec" << endl;
	
    double dt = min(max_dt, CFL*h/maxSOS);
    cout << "sos max is " << setprecision(8) << maxSOS << ", " << "dt is "<< dt << "\n";
    
    if (maxSOS > 1e6)
    {
        cout << "Speed of sound is too high. Is it realistic?" << endl;
        abort();
    }
    
    if (dt<std::numeric_limits<double>::epsilon()*1e1)
    {
        cout << "Last time step encountered." << endl;
        return 0;
    }

    if (LSRK3data::verbosity >= 1)
        cout << "Dispatcher is " << LSRK3data::dispatcher << endl;
    
    if (parser("-kernels").asString("cpp")=="cpp")
        LSRKstep<Convection_CPP, Update_CPP>(grid, dt/h, current_time, bAwk);
#if defined(_QPX_) || defined(_QPXEMU_)    
	else if (parser("-kernels").asString("cpp")=="qpx")
		LSRKstep<Convection_QPX, Update_QPX>(grid, dt/h, current_time, bAwk);
#endif
    else
    {
        cout << "combination not supported yet" << endl;
        abort();
    }
    
    	const Real pressure_avg=_computeFLAP();    
	cout << "Average pressure at flap " << pressure_avg << endl;    

    LSRK3data::step_id++;
    
    return dt;
}
