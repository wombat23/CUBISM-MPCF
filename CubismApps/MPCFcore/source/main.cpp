/*
 *  main.cpp
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 2/27/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#if defined(_SSE_) && _BLOCKSIZE_%4!=0
#error BLOCKSIZE NOT GOOD FOR SSE
#elif defined(_AVX_) && (_BLOCKSIZE_%8!=0 || !defined(_SP_COMP_))
#error BLOCKSIZE NOT GOOD FOR AVX
#endif

#include <iostream>

#include <ArgumentParser.h>

#include "TestTypes.h"
#include "FlowStep_Test.h"
#include "SurfaceTension_Test.h"
#include "LocalKernel_Test.h"
#include "Diffusion_Test.h"

#include "Convection_CPP.h"
#include "Update.h"
#include "MaxSpeedOfSound.h"
#include "SurfaceTension_CPP.h"
#include "Diffusion_CPP.h"
#ifdef _SSE_
#include "Convection_SSE.h"
#include "Diffusion_SSE.h"
#include "SurfaceTension_SSE.h"
#endif
#ifdef _AVX_
#include "Convection_AVX.h"
#include "SurfaceTension_AVX.h"
#include "Diffusion_AVX.h"
#endif

using namespace std;

int main (int argc, const char ** argv) 
{	
	ArgumentParser parser(argc, argv);
	
	int N = parser("-n").asInt(1000);
	const double accuracy = parser("-accuracy").asDouble(1e-4);
	const int NB = parser("-nblocks").asInt(512);
	const string kernel = parser("-kernel").asString("all");
	const bool bAwk = parser("-awk").asBool(false);
	const bool bFlush2Zero = parser("-f2z").asBool(false);
	const bool bProfile = parser("-profile").asBool(false);
	
	if (bFlush2Zero)
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	
#ifdef _SP_COMP_
	const double PP = parser("-pp").asDouble(21);
#else
	const double PP = parser("-pp").asDouble(21./2);
#endif //_SP_COMP_
	const double PB = parser("-pb").asDouble(4.5);
	
	LocalKernel_Test local_test;
	
	if (kernel == "SurfaceTension_CPP" || kernel == "all")
	{		
		SurfaceTension_Test surf_test;
		SurfaceTension_CPP st;
		
		printKernelName("SurfaceTension_CPP:");
		if(!bProfile)
		{
			surf_test.accuracy(st, accuracy, bAwk);
			surf_test.performance(st, PP*1e9, PB*1e9, NB, N, bAwk);	
		}
		else
			surf_test.profile(st, PP*1e9, PB*1e9, NB, N, bAwk);
	}
	
	if (kernel == "SurfaceTension_SSE" || kernel == "all")
	{		
		SurfaceTension_Test surf_test;
		SurfaceTension_SSE st;
		
		printKernelName("SurfaceTension_SSE:");
		if(!bProfile)
		{
			surf_test.accuracy(st, accuracy, bAwk);
			surf_test.performance(st, PP*1e9, PB*1e9, NB, N, bAwk);	
		}
		else
			surf_test.profile(st, PP*1e9, PB*1e9, NB, N, bAwk);
	}

#if defined(_AVX_) 
	if (kernel == "SurfaceTension_AVX" || kernel == "all")
	{		
		SurfaceTension_Test surf_test;
		SurfaceTension_AVX st;
		
		printKernelName("SurfaceTension_AVX:");
		if(!bProfile)
		{
			surf_test.accuracy(st, accuracy, bAwk);
			surf_test.performance(st, PP*1e9, PB*1e9, NB, N, bAwk);	
		}
		else
			surf_test.profile(st, PP*1e9, PB*1e9, NB, N, bAwk);
	}
#endif
	
	if (kernel == "Diffusion_CPP" || kernel == "all")
	{	
		Diffusion_CPP diffusion;
		Diffusion_Test diffusion_test;
		printKernelName("Diffusion_CPP:");
		
		if(!bProfile)
		{
			diffusion_test.accuracy(diffusion, accuracy, bAwk);
			diffusion_test.performance(diffusion, PP*1e9, PB*1e9, NB, N, bAwk);
		}
		else
			diffusion_test.profile(diffusion, PP*1e9, PB*1e9, NB, N, bAwk);
	}
	
	if (kernel == "Diffusion_SSE" || kernel == "all")
	{	
		Diffusion_SSE diffusion;
		Diffusion_Test diffusion_test;
		printKernelName("Diffusion_SSE:");
		
		if(!bProfile)
		{
			diffusion_test.accuracy(diffusion, accuracy, bAwk);
			diffusion_test.performance(diffusion, PP*1e9, PB*1e9, NB, N, bAwk);
		}
		else
			diffusion_test.profile(diffusion, PP*1e9, PB*1e9, NB, N, bAwk);	
	}
	
#if defined(_AVX_) 
	if (kernel == "Diffusion_AVX" || kernel == "all")
	{	
		Diffusion_AVX diffusion;
		Diffusion_Test diffusion_test;
		printKernelName("Diffusion_AVX:");
		
		if(!bProfile)
		{
			diffusion_test.accuracy(diffusion, accuracy, bAwk);
			diffusion_test.performance(diffusion, PP*1e9, PB*1e9, NB, N, bAwk);
		}
		else
			diffusion_test.profile(diffusion, PP*1e9, PB*1e9, NB, N, bAwk);	
	}
#endif
	
#ifdef _SSE_
	if (kernel == "MaxSOS_SSE" || kernel == "all")
	{
		MaxSpeedOfSound_CPP maxsos_cpp(2.5, 2.1, 1, 0, 0);
		MaxSpeedOfSound_SSE maxsos_sse(2.5, 2.1, 1, 0, 0);
		printKernelName("MaxSpeedOfSound_SSE:");
		local_test.accuracy(maxsos_sse, maxsos_cpp, accuracy, bAwk);
		local_test.performance(maxsos_sse, maxsos_cpp, PP*1e9, PB*1e9, NB, N*10, bAwk);
		printEndKernelTest();
	}
	
	if (kernel == "Update_SSE" | kernel == "all")
	{
		Update_CPP update_cpp;
		Update_SSE update_sse;
		printKernelName("Update_SSE:");
		local_test.accuracy(update_sse, update_cpp, accuracy, bAwk);
		local_test.performance(update_sse, update_cpp, PP*1e9, PB*1e9, NB, N*10, bAwk);
		printEndKernelTest();
	}
#endif //_SSE_
	
	FlowStep_Test test;
	
#ifdef _AVX_
#if _ALIGNBYTES_ % 32 == 0
	if (kernel == "Convection_AVX" || kernel == "all")
	{
		Convection_AVX convection(0, 1, 2.5, 2.1, 1, 0, 0);
		//test.profile(flowstep_avxdiego, PP*1e9, PB*1e9, N); 
		printKernelName("Convection_AVX:");
		test.accuracy(convection, accuracy, bAwk);
		test.performance(convection, PP*1e9, PB*1e9, NB, N, bAwk);
		printEndKernelTest();
	}
#endif //_ALIGNBYTES_
#endif // _AVX_
	
#ifdef _SSE_
#if _ALIGNBYTES_ % 16 == 0
	if (kernel == "Convection_SSE" || kernel == "all")
	{
		Convection_SSE convection(0, 1, 2.5, 2.1, 1, 0, 0);
		//test.profile(flowstep_diego, PP*1e9, PB*1e9, N); 
		printKernelName("Convection_SSE:");
		//test.profile(flowstep_diego, PP*1e9, PB*1e9, N); 
		test.accuracy(convection, accuracy, bAwk);
		test.performance(convection, PP*1e9, PB*1e9, NB, N, bAwk);
		printEndKernelTest();
	}
#endif //_ALIGNBYTES_
	
#endif //_SSE_
	
	if (kernel == "Convection_CPP" || kernel == "all")
	{
		Convection_CPP convection(0, 1, 2.5, 2.1, 1, 0, 0);
		printKernelName("FlowStep_CPP:");
		test.accuracy(convection, accuracy, bAwk);
		test.performance(convection, PP*1e9, PB*1e9, NB, N, bAwk);
		printEndKernelTest();
	}
	
	return 0;
}
