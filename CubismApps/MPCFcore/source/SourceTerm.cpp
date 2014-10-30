#include "SourceTerm.h"


/*
      src - entrance for block points
      gptfloats - number of "properties" - need for offset
 */

void SourceTerm::compute(
             const Real *const  pos ,
             const Real * const src, const int gptfloats,
             Real t,
             Real tSE, Real tSM,
             Real iI, Real iU,
             Real arcX, Real arcY, Real arcZ,
             Real arcW, Real arcH
                         )
{
    
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