#ifndef __AJALSOQJAJSDFASD_H__
#define __AJALSOQJAJSDFASD_H__

#include <stdio.h>

// This file is a stub, to make TRWS compile.


struct TypeTruncatedQuadratic2D
{
    typedef double REAL;

    struct Edge
    {
	void DistanceTransformL2(int /*K*/, int /*stride*/, REAL /*alpha*/, REAL* /*source*/, REAL* /*dest*/,
				 int* /*parabolas*/, int* /*intersections*/)
	{
	    printf("\n\
+-------------------------------------------------------------+\n\
|   In order to run TRW-S with truncted L2 terms,             |\n\
|   you need to download the implementation from              |\n\
|      http://pub.ist.ac.at/~vnk/papers/TRW-S.html            |\n\
|   and copy file  typeTruncatedQuadratic2D.h                 |\n\
|   to the main directory (i.e. replace the existing file)    |\n\
+-------------------------------------------------------------+\n\
			");
	    exit(1);
	}
    };
};

#endif
