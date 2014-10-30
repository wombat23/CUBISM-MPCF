#include "SourceTerm.h"



bool valueInRange(Real value, Real min, Real max){
    return (value >= min) && (value <= max);
}

bool rectOverlap(Real x1, Real y1 , Real width_1, Real height_1, Real x2, Real y2, Real width_2, Real height_2)
{
    bool xOverlap = valueInRange(x1, x2, x2 + width_2) ||
            valueInRange(x2, x1, x1 + width_1);

    bool yOverlap = valueInRange(y1, y2, y2 + height_2) ||
            valueInRange(y2, y1, y1 + height_1);

    return xOverlap && yOverlap;
}

bool pointInside(Real x, Real y, Real x2, Real y2, Real width, Real height){
    bool xInside = valueInRange(x, x2, x2 + width);
    bool yInside = valueInRange(y, y2, y2 + height);

    return xInside && yInside;
}



/*
      src - entrance for block points
      gptfloats - number of "properties" - need for offset
 */

void SourceTerm::compute(
        const Real *const pos , Real h_gridpoint,
        Real * src, const int gptfloats, Real t,
        InputStructVals inputStructVals )
{


    const int N = _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*gptfloats; // number of properties! 7 per point and we have _blockSize_^3 points


//    if (rectOverlap(
//            pos[0], pos[1], _BLOCKSIZE_*h_gridpoint, _BLOCKSIZE_*h_gridpoint,
//            inputStructVals.arcX, inputStructVals.arcY, inputStructVals.arcWidth, inputStructVals.arcHeight)
//            ) {
        // So this block has intersection with source. So we can iterate through points
        // and add source term to them.

        for (int i = 0; i < N; i += gptfloats) { // go through all points
            int icell = i / gptfloats;

            int ix = icell % _BLOCKSIZE_;
            int iy = (icell /_BLOCKSIZE_ ) % _BLOCKSIZE_ ;

            Real posX = pos[0] + ix * h_gridpoint;
            Real posY = pos[1] + iy * h_gridpoint;

            if (pointInside(
                    posX, posY,
                    inputStructVals.arcX, inputStructVals.arcY,
                    inputStructVals.arcWidth, inputStructVals.arcHeight)
                    ) {

		//src[i]   += 1000;
                src[i+4] += 1000;
            }

        }

//    }
}
