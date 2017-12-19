#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <new>
#include "BP-S.h"

#define private public
#include "typeTruncatedQuadratic2D.h"
#undef private


#define m_D(pix,l)  m_D[(pix)*m_nLabels+(l)]
#define m_V(l1,l2)  m_V[(l1)*m_nLabels+(l2)]


#define MIN(a,b)  (((a) < (b)) ? (a) : (b))
#define MAX(a,b)  (((a) > (b)) ? (a) : (b))
#define TRUNCATE_MIN(a,b) { if ((a) > (b)) (a) = (b); }
#define TRUNCATE_MAX(a,b) { if ((a) < (b)) (a) = (b); }
#define TRUNCATE TRUNCATE_MIN

/////////////////////////////////////////////////////////////////////////////
//                  Operations on vectors (arrays of size K)               //
/////////////////////////////////////////////////////////////////////////////

inline void CopyVector(BPS::REAL* to, MRF::CostVal* from, int K)
{
    BPS::REAL* to_finish = to + K;
    do
	{
	    *to ++ = *from ++;
	} while (to < to_finish);
}

inline void AddVector(BPS::REAL* to, BPS::REAL* from, int K)
{
    BPS::REAL* to_finish = to + K;
    do
	{
	    *to ++ += *from ++;
	} while (to < to_finish);
}

inline BPS::REAL SubtractMin(BPS::REAL *D, int K)
{
    int k;
    BPS::REAL delta;
	
    delta = D[0];
    for (k=1; k<K; k++) TRUNCATE(delta, D[k]);
    for (k=0; k<K; k++) D[k] -= delta;

    return delta;
}

// Functions UpdateMessageTYPE (see the paper for details):
//
// - Set Di[ki] := gamma*Di_hat[ki] - M[ki]
// - Set M[kj] := min_{ki} (Di[ki] + V[ki,kj])
// - Normalize message: 
//        delta := min_{kj} M[kj]
//        M[kj] := M[kj] - delta
//        return delta
//
// If dir = 1, then the meaning of i and j is swapped.

///////////////////////////////////////////
//                  L1                   //
///////////////////////////////////////////

inline BPS::REAL UpdateMessageL1(BPS::REAL* M, BPS::REAL* Di_hat, int K, BPS::REAL gamma, MRF::CostVal lambda, MRF::CostVal smoothMax)
{
    int k;
    BPS::REAL delta;

    delta = M[0] = gamma*Di_hat[0] - M[0];
    for (k=1; k<K; k++)
	{
	    M[k] = gamma*Di_hat[k] - M[k];
	    TRUNCATE(delta, M[k]);
	    TRUNCATE(M[k], M[k-1] + lambda);
	}

    M[--k] -= delta;
    TRUNCATE(M[k], lambda*smoothMax);
    for (k--; k>=0; k--)
	{
	    M[k] -= delta;
	    TRUNCATE(M[k], M[k+1] + lambda);
	    TRUNCATE(M[k], lambda*smoothMax);
	}

    return delta;
}

////////////////////////////////////////
//               L2                   //
////////////////////////////////////////

inline BPS::REAL UpdateMessageL2(BPS::REAL* M, BPS::REAL* Di_hat, int K, BPS::REAL gamma, MRF::CostVal lambda, MRF::CostVal smoothMax, void *buf)
{
    BPS::REAL* Di = (BPS::REAL*) buf;
    int* parabolas = (int*) ((char*)buf + K*sizeof(BPS::REAL));
    int* intersections = parabolas + K;
    TypeTruncatedQuadratic2D::REAL* Di_tmp = (TypeTruncatedQuadratic2D::REAL*) (intersections + K + 1);
    TypeTruncatedQuadratic2D::REAL* M_tmp = Di_tmp + K;
    TypeTruncatedQuadratic2D::Edge* tmp = NULL;

    int k;
    BPS::REAL delta;

    assert(lambda >= 0);

    Di[0] = gamma*Di_hat[0] - M[0];
    delta = Di[0];
    for (k=1; k<K; k++)
	{
	    Di[k] = gamma*Di_hat[k] - M[k];
	    TRUNCATE(delta, Di[k]);
	}

    if (lambda == 0)
	{
	    for (k=0; k<K; k++) M[k] = 0;
	    return delta;
	}

    for (k=0; k<K; k++) Di_tmp[k] = Di[k];
    tmp->DistanceTransformL2(K, 1, lambda, Di_tmp, M_tmp, parabolas, intersections);
    for (k=0; k<K; k++) M[k] = (BPS::REAL) M_tmp[k];

    for (k=0; k<K; k++)
	{
	    M[k] -= delta;
	    TRUNCATE(M[k], lambda*smoothMax);
	}

    return delta;
}


//////////////////////////////////////////////////
//                FIXED_MATRIX                  //
//////////////////////////////////////////////////

inline BPS::REAL UpdateMessageFIXED_MATRIX(BPS::REAL* M, BPS::REAL* Di_hat, int K, BPS::REAL gamma, MRF::CostVal lambda, MRF::CostVal* V, void* buf)
{
    BPS::REAL* Di = (BPS::REAL*) buf;
    int ki, kj;
    BPS::REAL delta;

    for (ki=0; ki<K; ki++)
	{
	    Di[ki] = gamma*Di_hat[ki] - M[ki];
	}

    for (kj=0; kj<K; kj++)
	{
	    M[kj] = Di[0] + lambda*V[0]; 
	    V ++;
	    for (ki=1; ki<K; ki++)
		{
		    TRUNCATE(M[kj], Di[ki] + lambda*V[0]);
		    V ++;
		}
	}

    delta = M[0];
    for (kj=1; kj<K; kj++) TRUNCATE(delta, M[kj]);
    for (kj=0; kj<K; kj++) M[kj] -= delta;

    return delta;
}

/////////////////////////////////////////////
//                GENERAL                  //
/////////////////////////////////////////////

inline BPS::REAL UpdateMessageGENERAL(BPS::REAL* M, BPS::REAL* Di_hat, int K, BPS::REAL gamma, int dir, MRF::CostVal* V, void* buf)
{
    BPS::REAL* Di = (BPS::REAL*) buf;
    int ki, kj;
    BPS::REAL delta;

    for (ki=0; ki<K; ki++)
	{
	    Di[ki] = (gamma*Di_hat[ki] - M[ki]);
	}

    if (dir == 0)
	{
	    for (kj=0; kj<K; kj++)
		{
		    M[kj] = Di[0] + V[0]; 
		    V ++;
		    for (ki=1; ki<K; ki++)
			{
			    TRUNCATE(M[kj], Di[ki] + V[0]);
			    V ++;
			}
		}
	}
    else
	{
	    for (kj=0; kj<K; kj++)
		{
		    M[kj] = Di[0] + V[0];
		    V += K;
		    for (ki=1; ki<K; ki++)
			{
			    TRUNCATE(M[kj], Di[ki] + V[0]);
			    V += K;
			}
		    V -= K*K - 1;
		}
	}

    delta = M[0];
    for (kj=1; kj<K; kj++) TRUNCATE(delta, M[kj]);
    for (kj=0; kj<K; kj++) M[kj] -= delta;

    return delta;
}

inline BPS::REAL UpdateMessageGENERAL(BPS::REAL* M, BPS::REAL* Di_hat, int K, BPS::REAL gamma, BPS::SmoothCostGeneralFn fn, int i, int j, void* buf)
{
    BPS::REAL* Di = (BPS::REAL*) buf;
    int ki, kj;
    BPS::REAL delta;

    for (ki=0; ki<K; ki++)
	{
	    Di[ki] = (gamma*Di_hat[ki] - M[ki]);
	}

    for (kj=0; kj<K; kj++)
	{
	    M[kj] = Di[0] + fn(i, j, 0, kj); 
	    for (ki=1; ki<K; ki++)
		{
		    delta = Di[ki] + fn(i, j, ki, kj);
		    TRUNCATE(M[kj], delta);
		}
	}

    delta = M[0];
    for (kj=1; kj<K; kj++) TRUNCATE(delta, M[kj]);
    for (kj=0; kj<K; kj++) M[kj] -= delta;

    return delta;
}






BPS::BPS(int width, int height, int nLabels,EnergyFunction *eng):MRF(width,height,nLabels,eng)
{
    Allocate();
}
BPS::BPS(int nPixels, int nLabels,EnergyFunction *eng):MRF(nPixels,nLabels,eng)
{
    Allocate();
}

BPS::~BPS()
{ 
    delete[] m_answer;
    if ( m_needToFreeD ) delete [] m_D;
    if ( m_needToFreeV ) delete [] m_V;
    if ( m_messages ) delete [] m_messages;
    if ( m_DBinary ) delete [] m_DBinary;
    if ( m_horzWeightsBinary ) delete [] m_horzWeightsBinary;
    if ( m_vertWeightsBinary ) delete [] m_vertWeightsBinary;
}


void BPS::Allocate()
{
    m_type = NONE;
    m_needToFreeV = false;
    m_needToFreeD = false;

    m_D = NULL;
    m_V = NULL;
    m_horzWeights = NULL;
    m_vertWeights = NULL;
    m_horzWeightsBinary = NULL;
    m_vertWeightsBinary = NULL;

    m_DBinary = NULL;
    m_messages = NULL;
    m_messageArraySizeInBytes = 0;

    m_answer = new Label[m_nPixels];
}

void BPS::clearAnswer()
{
    memset(m_answer, 0, m_nPixels*sizeof(Label));
    if (m_messages)
	{
	    memset(m_messages, 0, m_messageArraySizeInBytes);
	}
}


MRF::EnergyVal BPS::smoothnessEnergy()
{
    EnergyVal eng = (EnergyVal) 0;
    EnergyVal weight;
    int x,y,pix;

    if ( m_grid_graph )
	{
	    if ( m_smoothType != FUNCTION  )
		{
		    for ( y = 0; y < m_height; y++ )
			for ( x = 1; x < m_width; x++ )
			    {
				pix    = x+y*m_width;
				weight = m_varWeights ? m_horzWeights[pix-1] :  1;
				eng = eng + m_V(m_answer[pix],m_answer[pix-1])*weight;
			    }

		    for ( y = 1; y < m_height; y++ )
			for ( x = 0; x < m_width; x++ )
			    {
				pix = x+y*m_width;
				weight = m_varWeights ? m_vertWeights[pix-m_width] :  1;
				eng = eng + m_V(m_answer[pix],m_answer[pix-m_width])*weight;
			    }
		}
	    else
		{
		    for ( y = 0; y < m_height; y++ )
			for ( x = 1; x < m_width; x++ )
			    {
				pix = x+y*m_width;
				eng = eng + m_smoothFn(pix,pix-1,m_answer[pix],m_answer[pix-1]);
			    }

		    for ( y = 1; y < m_height; y++ )
			for ( x = 0; x < m_width; x++ )
			    {
				pix = x+y*m_width;
				eng = eng + m_smoothFn(pix,pix-m_width,m_answer[pix],m_answer[pix-m_width]);
			    }
		}
	}
    else
	{
	    // not implemented
	}

    return(eng);
}



MRF::EnergyVal BPS::dataEnergy()
{
    EnergyVal eng = (EnergyVal) 0;

    
    if ( m_dataType == ARRAY) 
	{
	    for ( int i = 0; i < m_nPixels; i++ )
		eng = eng + m_D(i,m_answer[i]);
	}
    else
	{
	    for ( int i = 0; i < m_nPixels; i++ )
		eng = eng + m_dataFn(i,m_answer[i]);
	}
    return(eng);

}


void BPS::setData(DataCostFn dcost)
{
    int i, k;

    m_dataFn = dcost;
    CostVal* ptr;
    m_D = new CostVal[m_nPixels*m_nLabels];

    for (ptr=m_D, i=0; i<m_nPixels; i++)
	for (k=0; k<m_nLabels; k++, ptr++)
	    {
		*ptr = m_dataFn(i,k);
	    }
    m_needToFreeD = true;
}

void BPS::setData(CostVal* data)
{
    m_D = data;
    m_needToFreeD = false;
}


void BPS::setSmoothness(SmoothCostGeneralFn cost)
{
    assert(m_horzWeights == NULL && m_vertWeights == NULL && m_V == NULL);

    int x, y, i, ki, kj;
    CostVal* ptr;

    m_smoothFn = cost;
    m_type = GENERAL;

    if (!m_allocateArrayForSmoothnessCostFn) return;

    // try to cache all the function values in an array for efficiency
    m_V = new(std::nothrow) CostVal[2*m_nPixels*m_nLabels*m_nLabels];
    if (!m_V) {
	fprintf(stderr, "not caching smoothness cost values (not enough memory)\n");
	return; // if not enough space, just call the function directly
    }

    m_needToFreeV = true;

    for (ptr=m_V,i=0,y=0; y<m_height; y++)
	for (x=0; x<m_width; x++, i++)
	    {
		if (x < m_width-1)
		    {
			for (kj=0; kj<m_nLabels; kj++)
			    for (ki=0; ki<m_nLabels; ki++)
				{
				    *ptr++ = cost(i,i+1,ki,kj);
				}
		    }
		else ptr += m_nLabels*m_nLabels;

		if (y < m_height-1)
		    {
			for (kj=0; kj<m_nLabels; kj++)
			    for (ki=0; ki<m_nLabels; ki++)
				{
				    *ptr++ = cost(i,i+m_width,ki,kj);
				}
		    }
		else ptr += m_nLabels*m_nLabels;
	    }
}
void BPS::setSmoothness(CostVal* V)
{
    m_type = FIXED_MATRIX;
    m_V = V;
}


void BPS::setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda)
{
    assert(smoothExp == 1 || smoothExp == 2);
    assert(lambda >= 0);

    m_type = (smoothExp == 1) ? L1 : L2;

    int ki, kj;
    CostVal cost;

    m_needToFreeV = true;

    m_V = new CostVal[m_nLabels*m_nLabels];

    for (ki=0; ki<m_nLabels; ki++)
	for (kj=ki; kj<m_nLabels; kj++)
	    {
		cost = (CostVal) ((smoothExp == 1) ? kj - ki : (kj - ki)*(kj - ki));
		if (cost > smoothMax) cost = smoothMax;
		m_V[ki*m_nLabels + kj] = m_V[kj*m_nLabels + ki] = cost*lambda;
	    }

    m_smoothMax = smoothMax;
    m_lambda = lambda;
}


void BPS::setCues(CostVal* hCue, CostVal* vCue)
{
    m_horzWeights = hCue;
    m_vertWeights  = vCue;
}


void BPS::initializeAlg()
{
    assert(m_type != NONE);

    int i;

    // determine type
    if (m_type == L1 && m_nLabels == 2)
	{
	    m_type = BINARY;
	}

    // allocate messages
    int messageNum = (m_type == BINARY) ? 4*m_nPixels : 4*m_nPixels*m_nLabels;
    m_messageArraySizeInBytes = messageNum*sizeof(REAL);
    m_messages = new REAL[messageNum];
    memset(m_messages, 0, messageNum*sizeof(REAL));

    if (m_type == BINARY)
	{
	    assert(m_DBinary == NULL && m_horzWeightsBinary == NULL && m_horzWeightsBinary == NULL);
	    m_DBinary = new CostVal[m_nPixels];
	    m_horzWeightsBinary = new CostVal[m_nPixels];
	    m_vertWeightsBinary = new CostVal[m_nPixels];

	    if ( m_dataType == ARRAY)
		{
		    for (i=0; i<m_nPixels; i++)
			{
			    m_DBinary[i] = m_D[2*i+1] - m_D[2*i];
			}
		}
	    else
		{
		    for (i=0; i<m_nPixels; i++)
			{
			    m_DBinary[i] = m_dataFn(i,1) - m_dataFn(i,0);
			}
		}

	    assert(m_V[0] == 0 && m_V[1] == m_V[2] && m_V[3] == 0);
	    for (i=0; i<m_nPixels; i++)
		{
		    m_horzWeightsBinary[i] = (m_varWeights) ? m_V[1]*m_horzWeights[i] : m_V[1];
		    m_vertWeightsBinary[i] = (m_varWeights) ? m_V[1]*m_vertWeights[i] : m_V[1];
		}
	}
}

void BPS::optimizeAlg(int nIterations)
{
    assert(m_type != NONE);

    if (m_grid_graph)
	{
	    switch (m_type)
		{
		case L1:            optimize_GRID_L1(nIterations);           break;
		case L2:            optimize_GRID_L2(nIterations);           break;
		case FIXED_MATRIX:  optimize_GRID_FIXED_MATRIX(nIterations); break;
		case GENERAL:       optimize_GRID_GENERAL(nIterations);      break;
		case BINARY:        optimize_GRID_BINARY(nIterations);       break;
		default: assert(0); exit(1);
		}
	}
    else {printf("\nNot implemented for general graphs yet, exiting!");exit(1);}

    //	printf("lower bound = %f\n", m_lowerBound);

    ////////////////////////////////////////////////
    //          computing solution                //
    ////////////////////////////////////////////////

    if (m_type != BINARY)
	{
	    int x, y, n, K = m_nLabels;
	    CostVal* D_ptr;
	    REAL* M_ptr;
	    REAL* Di;
	    REAL delta;
	    int ki, kj;

	    Di = new REAL[K];

	    n = 0;
	    D_ptr = m_D;
	    M_ptr = m_messages;

	    for (y=0; y<m_height; y++)
		for (x=0; x<m_width; x++, D_ptr+=K, M_ptr+=2*K, n++)
		    {
			CopyVector(Di, D_ptr, K);

			if (m_type == GENERAL)
			    {
				if (m_V)
				    {
					CostVal* ptr = m_V + 2*(x+y*m_width-1)*K*K;
					if (x > 0)
					    {
						kj = m_answer[n-1];
						for (ki=0; ki<K; ki++)
						    {
							Di[ki] += ptr[kj + ki*K];
						    }
					    }
					ptr -= (2*m_width-3)*K*K;
					if (y > 0)
					    {
						kj = m_answer[n-m_width];
						for (ki=0; ki<K; ki++)
						    {
							Di[ki] += ptr[kj + ki*K];
						    }
					    }
				    }
				else
				    {
					if (x > 0)
					    {
						kj = m_answer[n-1];
						for (ki=0; ki<K; ki++)
						    {
							Di[ki] += m_smoothFn(n, n-1, ki, kj);
						    }
					    }
					if (y > 0)
					    {
						kj = m_answer[n-m_width];
						for (ki=0; ki<K; ki++)
						    {
							Di[ki] += m_smoothFn(n, n-m_width, ki, kj);
						    }
					    }
				    }
			    }
			else // m_type == L1, L2 or FIXED_MATRIX
			    {
				if (x > 0)
				    {
					kj = m_answer[n-1];
					CostVal lambda = (m_varWeights) ? m_horzWeights[n-1] : 1;
					for (ki=0; ki<K; ki++)
					    {
						Di[ki] += lambda*m_V[kj*K + ki];
					    }
				    }
				if (y > 0)
				    {
					kj = m_answer[n-m_width];
					CostVal lambda = (m_varWeights) ? m_vertWeights[n-m_width] : 1;
					for (ki=0; ki<K; ki++)
					    {
						Di[ki] += lambda*m_V[kj*K + ki];
					    }
				    }
			    }

			if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
			if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

			// compute min
			delta = Di[0];
			m_answer[n] = 0;
			for (ki=1; ki<K; ki++)
			    {
				if (delta > Di[ki])
				    {
					delta = Di[ki];
					m_answer[n] = ki;
				    }
			    }
		    }

	    delete [] Di;
	}
    else // m_type == BINARY
	{
	    int x, y, n;
	    REAL* M_ptr;
	    REAL Di;

	    n = 0;
	    M_ptr = m_messages;

	    for (y=0; y<m_height; y++)
		for (x=0; x<m_width; x++, M_ptr+=2, n++)
		    {
			Di = m_DBinary[n];
			if (x > 0) Di += (m_answer[n-1] == 0)       ? m_horzWeightsBinary[n-1]       : -m_horzWeightsBinary[n-1];
			if (y > 0) Di += (m_answer[n-m_width] == 0) ? m_vertWeightsBinary[n-m_width] : -m_vertWeightsBinary[n-m_width];

			if (x < m_width-1)  Di += M_ptr[0]; // message (x+1,y)->(x,y)
			if (y < m_height-1) Di += M_ptr[1]; // message (x,y+1)->(x,y)

			// compute min
			m_answer[n] = (Di >= 0) ? 0 : 1;
		    }
	}
}

void BPS::optimize_GRID_L1(int nIterations)
{
    int x, y, n, K = m_nLabels;
    CostVal* D_ptr;
    REAL* M_ptr;
    REAL* Di;

    Di = new REAL[K];

    for ( ; nIterations > 0; nIterations --)
	{
	    // forward pass
	    n = 0;
	    D_ptr = m_D;
	    M_ptr = m_messages;

	    for (y=0; y<m_height; y++)
		for (x=0; x<m_width; x++, D_ptr+=K, M_ptr+=2*K, n++)
		    {
			CopyVector(Di, D_ptr, K);
			if (x > 0) AddVector(Di, M_ptr-2*K, K); // message (x-1,y)->(x,y)
			if (y > 0) AddVector(Di, M_ptr-(2*m_width-1)*K, K); // message (x,y-1)->(x,y)
			if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
			if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

			if (x < m_width-1) 
			    {
				CostVal lambda = (m_varWeights) ? m_lambda*m_horzWeights[n] : m_lambda;
				UpdateMessageL1(M_ptr, Di, K, 1, lambda, m_smoothMax);
			    }
			if (y < m_height-1) 
			    {
				CostVal lambda = (m_varWeights) ? m_lambda*m_vertWeights[n] : m_lambda;
				UpdateMessageL1(M_ptr+K, Di, K, 1, lambda, m_smoothMax);
			    }
		    }

	    // backward pass
	    n --;
	    D_ptr -= K;
	    M_ptr -= 2*K;

	    for (y=m_height-1; y>=0; y--)
		for (x=m_width-1; x>=0; x--, D_ptr-=K, M_ptr-=2*K, n--)
		    {
			CopyVector(Di, D_ptr, K);
			if (x > 0) AddVector(Di, M_ptr-2*K, K); // message (x-1,y)->(x,y)
			if (y > 0) AddVector(Di, M_ptr-(2*m_width-1)*K, K); // message (x,y-1)->(x,y)
			if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
			if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

			SubtractMin(Di, K);

			if (x > 0) 
			    {
				CostVal lambda = (m_varWeights) ? m_lambda*m_horzWeights[n-1] : m_lambda;
				UpdateMessageL1(M_ptr-2*K, Di, K, 1, lambda, m_smoothMax);
			    }
			if (y > 0) 
			    {
				CostVal lambda = (m_varWeights) ? m_lambda*m_vertWeights[n-m_width] : m_lambda;
				UpdateMessageL1(M_ptr-(2*m_width-1)*K, Di, K, 1, lambda, m_smoothMax);
			    }
		    }
	}

    delete [] Di;
}

void BPS::optimize_GRID_L2(int nIterations)
{
    int x, y, n, K = m_nLabels;
    CostVal* D_ptr;
    REAL* M_ptr;
    REAL* Di;
    void* buf;

    Di = new REAL[K];
    buf = new char[2*K*sizeof(TypeTruncatedQuadratic2D::REAL) + (2*K+1)*sizeof(int) + K*sizeof(REAL)];

    for ( ; nIterations > 0; nIterations --)
	{
	    // forward pass
	    n = 0;
	    D_ptr = m_D;
	    M_ptr = m_messages;

	    for (y=0; y<m_height; y++)
		for (x=0; x<m_width; x++, D_ptr+=K, M_ptr+=2*K, n++)
		    {
			CopyVector(Di, D_ptr, K);
			if (x > 0) AddVector(Di, M_ptr-2*K, K); // message (x-1,y)->(x,y)
			if (y > 0) AddVector(Di, M_ptr-(2*m_width-1)*K, K); // message (x,y-1)->(x,y)
			if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
			if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

			if (x < m_width-1) 
			    {
				CostVal lambda = (m_varWeights) ? m_lambda*m_horzWeights[n] : m_lambda;
				UpdateMessageL2(M_ptr, Di, K, 1, lambda, m_smoothMax, buf);
			    }
			if (y < m_height-1) 
			    {
				CostVal lambda = (m_varWeights) ? m_lambda*m_vertWeights[n] : m_lambda;
				UpdateMessageL2(M_ptr+K, Di, K, 1, lambda, m_smoothMax, buf);
			    }
		    }

	    // backward pass
	    n --;
	    D_ptr -= K;
	    M_ptr -= 2*K;

	    for (y=m_height-1; y>=0; y--)
		for (x=m_width-1; x>=0; x--, D_ptr-=K, M_ptr-=2*K, n--)
		    {
			CopyVector(Di, D_ptr, K);
			if (x > 0) AddVector(Di, M_ptr-2*K, K); // message (x-1,y)->(x,y)
			if (y > 0) AddVector(Di, M_ptr-(2*m_width-1)*K, K); // message (x,y-1)->(x,y)
			if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
			if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

			SubtractMin(Di, K);

			if (x > 0) 
			    {
				CostVal lambda = (m_varWeights) ? m_lambda*m_horzWeights[n-1] : m_lambda;
				UpdateMessageL2(M_ptr-2*K, Di, K, 1, lambda, m_smoothMax, buf);
			    }
			if (y > 0) 
			    {
				CostVal lambda = (m_varWeights) ? m_lambda*m_vertWeights[n-m_width] : m_lambda;
				UpdateMessageL2(M_ptr-(2*m_width-1)*K, Di, K, 1, lambda, m_smoothMax, buf);
			    }
		    }
	}

    delete [] Di;
    delete [] (REAL *)buf;
}


void BPS::optimize_GRID_BINARY(int nIterations)
{
    int x, y, n;
    REAL* M_ptr;
    REAL Di;

    for ( ; nIterations > 0; nIterations --)
	{
	    // forward pass
	    n = 0;
	    M_ptr = m_messages;

	    for (y=0; y<m_height; y++)
		for (x=0; x<m_width; x++, M_ptr+=2, n++)
		    {
			Di = m_DBinary[n];
			if (x > 0) Di += M_ptr[-2]; // message (x-1,y)->(x,y)
			if (y > 0) Di += M_ptr[-2*m_width+1]; // message (x,y-1)->(x,y)
			if (x < m_width-1) Di += M_ptr[0]; // message (x+1,y)->(x,y)
			if (y < m_height-1) Di += M_ptr[1]; // message (x,y+1)->(x,y)

			REAL DiScaled = Di * 1;
			if (x < m_width-1) 
			    {
				Di = DiScaled - M_ptr[0];
				CostVal lambda = m_horzWeightsBinary[n];
				if (lambda < 0) { Di = -Di; lambda = -lambda; }
				if (Di > lambda) M_ptr[0] = lambda;
				else             M_ptr[0] = (Di < -lambda) ? -lambda : Di;
			    }
			if (y < m_height-1) 
			    {
				Di = DiScaled - M_ptr[1];
				CostVal lambda = m_vertWeightsBinary[n];
				if (lambda < 0) { Di = -Di; lambda = -lambda; }
				if (Di > lambda) M_ptr[1] = lambda;
				else             M_ptr[1] = (Di < -lambda) ? -lambda : Di;
			    }
		    }

	    // backward pass
	    n --;
	    M_ptr -= 2;

	    for (y=m_height-1; y>=0; y--)
		for (x=m_width-1; x>=0; x--, M_ptr-=2, n--)
		    {
			Di = m_DBinary[n];
			if (x > 0) Di += M_ptr[-2]; // message (x-1,y)->(x,y)
			if (y > 0) Di += M_ptr[-2*m_width+1]; // message (x,y-1)->(x,y)
			if (x < m_width-1) Di += M_ptr[0]; // message (x+1,y)->(x,y)
			if (y < m_height-1) Di += M_ptr[1]; // message (x,y+1)->(x,y)

			REAL DiScaled = Di * 1;
			if (x > 0) 
			    {
				Di = DiScaled - M_ptr[-2];
				CostVal lambda = m_horzWeightsBinary[n-1];
				if (lambda < 0) { Di = -Di; lambda = -lambda; }
				if (Di > lambda) M_ptr[-2] = lambda;
				else             M_ptr[-2] = (Di < -lambda) ? -lambda : Di;
			    }
			if (y > 0) 
			    {
				Di = DiScaled - M_ptr[-2*m_width+1];
				CostVal lambda = m_vertWeightsBinary[n-m_width];
				if (lambda < 0) { Di = -Di; lambda = -lambda; }
				if (Di > lambda) M_ptr[-2*m_width+1] = lambda;
				else             M_ptr[-2*m_width+1] = (Di < -lambda) ? -lambda : Di;
			    }
		    }
	}
}

void BPS::optimize_GRID_FIXED_MATRIX(int nIterations)
{
    int x, y, n, K = m_nLabels;
    CostVal* D_ptr;
    REAL* M_ptr;
    REAL* Di;
    void* buf;

    Di = new REAL[K];
    buf = new REAL[K];

    for ( ; nIterations > 0; nIterations --)
	{
	    // forward pass
	    n = 0;
	    D_ptr = m_D;
	    M_ptr = m_messages;

	    for (y=0; y<m_height; y++)
		for (x=0; x<m_width; x++, D_ptr+=K, M_ptr+=2*K, n++)
		    {
			CopyVector(Di, D_ptr, K);
			if (x > 0) AddVector(Di, M_ptr-2*K, K); // message (x-1,y)->(x,y)
			if (y > 0) AddVector(Di, M_ptr-(2*m_width-1)*K, K); // message (x,y-1)->(x,y)
			if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
			if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

			if (x < m_width-1) 
			    {
				CostVal lambda = (m_varWeights) ? m_horzWeights[n] : 1;
				UpdateMessageFIXED_MATRIX(M_ptr, Di, K, 1, lambda, m_V, buf);
			    }
			if (y < m_height-1) 
			    {
				CostVal lambda = (m_varWeights) ? m_vertWeights[n] : 1;
				UpdateMessageFIXED_MATRIX(M_ptr+K, Di, K, 1, lambda, m_V, buf);
			    }
		    }

	    // backward pass
	    n --;
	    D_ptr -= K;
	    M_ptr -= 2*K;

	    for (y=m_height-1; y>=0; y--)
		for (x=m_width-1; x>=0; x--, D_ptr-=K, M_ptr-=2*K, n--)
		    {
			CopyVector(Di, D_ptr, K);
			if (x > 0) AddVector(Di, M_ptr-2*K, K); // message (x-1,y)->(x,y)
			if (y > 0) AddVector(Di, M_ptr-(2*m_width-1)*K, K); // message (x,y-1)->(x,y)
			if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
			if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

			SubtractMin(Di, K);

			if (x > 0) 
			    {
				CostVal lambda = (m_varWeights) ? m_horzWeights[n-1] : 1;
				UpdateMessageFIXED_MATRIX(M_ptr-2*K, Di, K, 1, lambda, m_V, buf);
			    }
			if (y > 0) 
			    {
				CostVal lambda = (m_varWeights) ? m_vertWeights[n-m_width] : 1;
				UpdateMessageFIXED_MATRIX(M_ptr-(2*m_width-1)*K, Di, K, 1, lambda, m_V, buf);
			    }
		    }
	}

    delete [] Di;
    delete [] (REAL *)buf;
}

void BPS::optimize_GRID_GENERAL(int nIterations)
{
    int x, y, n, K = m_nLabels;
    CostVal* D_ptr;
    REAL* M_ptr;
    REAL* Di;
    void* buf;

    Di = new REAL[K];
    buf = new REAL[K];

    for ( ; nIterations > 0; nIterations --)
	{
	    // forward pass
	    n = 0;
	    D_ptr = m_D;
	    M_ptr = m_messages;
	    CostVal* V_ptr = m_V;

	    for (y=0; y<m_height; y++)
		for (x=0; x<m_width; x++, D_ptr+=K, M_ptr+=2*K, V_ptr+=2*K*K, n++)
		    {
			CopyVector(Di, D_ptr, K);
			if (x > 0) AddVector(Di, M_ptr-2*K, K); // message (x-1,y)->(x,y)
			if (y > 0) AddVector(Di, M_ptr-(2*m_width-1)*K, K); // message (x,y-1)->(x,y)
			if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
			if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

			if (x < m_width-1) 
			    {
				if (m_V) UpdateMessageGENERAL(M_ptr, Di, K, 1, /* forward dir*/ 0, V_ptr, buf);
				else     UpdateMessageGENERAL(M_ptr, Di, K, 1,   m_smoothFn, n, n+1,      buf);
			    }
			if (y < m_height-1) 
			    {
				if (m_V) UpdateMessageGENERAL(M_ptr+K, Di, K, 1, /* forward dir*/ 0, V_ptr+K*K, buf);
				else     UpdateMessageGENERAL(M_ptr+K, Di, K, 1,   m_smoothFn, n, n+m_width,    buf);
			    }
		    }

	    // backward pass
	    n --;
	    D_ptr -= K;
	    M_ptr -= 2*K;
	    V_ptr -= 2*K*K;

	    for (y=m_height-1; y>=0; y--)
		for (x=m_width-1; x>=0; x--, D_ptr-=K, M_ptr-=2*K, V_ptr-=2*K*K, n--)
		    {
			CopyVector(Di, D_ptr, K);
			if (x > 0) AddVector(Di, M_ptr-2*K, K); // message (x-1,y)->(x,y)
			if (y > 0) AddVector(Di, M_ptr-(2*m_width-1)*K, K); // message (x,y-1)->(x,y)
			if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
			if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

			// normalize Di, update lower bound
			SubtractMin(Di, K);

			if (x > 0) 
			    {
				if (m_V) UpdateMessageGENERAL(M_ptr-2*K, Di, K, 1, /* backward dir */ 1, V_ptr-2*K*K, buf);
				else     UpdateMessageGENERAL(M_ptr-2*K, Di, K, 1,   m_smoothFn, n, n-1,              buf);
			    }
			if (y > 0)
			    {
				if (m_V) UpdateMessageGENERAL(M_ptr-(2*m_width-1)*K, Di, K, 1, /* backward dir */ 1, V_ptr-(2*m_width-1)*K*K, buf);
				else     UpdateMessageGENERAL(M_ptr-(2*m_width-1)*K, Di, K, 1,   m_smoothFn, n, n-m_width,                    buf);

			    }
		    }
	}

    delete [] Di;
    delete [] (REAL *)buf;
}
