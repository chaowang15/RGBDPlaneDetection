#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "ICM.h"

#define m_D(pix,l)  m_D[(pix)*m_nLabels+(l)]
#define m_V(l1,l2)  m_V[(l1)*m_nLabels+(l2)]


        


ICM::ICM(int width, int height, int nLabels,EnergyFunction *eng):MRF(width,height,nLabels,eng)
{
    m_needToFreeV = 0;
    initializeAlg();
}
ICM::ICM(int nPixels, int nLabels,EnergyFunction *eng):MRF(nPixels,nLabels,eng)
{
    m_needToFreeV = 0;
    initializeAlg();
}

ICM::~ICM()
{ 
    delete[] m_answer;
    if (!m_grid_graph) delete[] m_neighbors;
    if ( m_needToFreeV ) delete[] m_V;
}


void ICM::initializeAlg()
{
    m_answer = (Label *) new Label[m_nPixels];
    if ( !m_answer ){printf("\nNot enough memory, exiting");exit(0);}

    if (!m_grid_graph)
    {
        m_neighbors = (LinkedBlockList *) new LinkedBlockList[m_nPixels];
        if (!m_neighbors) {printf("Not enough memory,exiting");exit(0);};
    }
}

void ICM::clearAnswer()
{
    memset(m_answer, 0, m_nPixels*sizeof(Label));
}

void ICM::setNeighbors(int pixel1, int pixel2, CostVal weight)
{
    assert(!m_grid_graph);
    assert(pixel1 < m_nPixels && pixel1 >= 0 && pixel2 < m_nPixels && pixel2 >= 0);


    Neighbor *temp1 = (Neighbor *) new Neighbor;
    Neighbor *temp2 = (Neighbor *) new Neighbor;

    if ( !temp1 || ! temp2 ) {printf("\nNot enough memory, exiting");exit(0);}

    temp1->weight  = weight;
    temp1->to_node = pixel2;

    temp2->weight  = weight;
    temp2->to_node = pixel1;

    m_neighbors[pixel1].addFront(temp1);
    m_neighbors[pixel2].addFront(temp2);
}


MRF::EnergyVal ICM::smoothnessEnergy()
{
    EnergyVal eng = (EnergyVal) 0;
    EnergyVal weight;
    int x,y,pix,i;

    if ( m_grid_graph )
    {
        if ( m_smoothType != FUNCTION  )
        {
            for ( y = 0; y < m_height; y++ )
                for ( x = 1; x < m_width; x++ )
                {
                    pix    = x+y*m_width;
                    weight = m_varWeights ? m_horizWeights[pix-1] :  1;
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
        
        Neighbor *temp; 

        if ( m_smoothType != FUNCTION  )
        {
            for ( i = 0; i < m_nPixels; i++ )
                if ( !m_neighbors[i].isEmpty() )
                {
                    m_neighbors[i].setCursorFront();
                    while ( m_neighbors[i].hasNext() )
                    {
                        temp = (Neighbor *) m_neighbors[i].next();

                        if ( i < temp->to_node )
                            eng = eng + m_V(m_answer[i],m_answer[temp->to_node])*(temp->weight);
                    }
                }
        }
        else
        {
            for ( i = 0; i < m_nPixels; i++ )
                if ( !m_neighbors[i].isEmpty() )
                {
                    m_neighbors[i].setCursorFront();
                    while ( m_neighbors[i].hasNext() )
                    {
                        temp = (Neighbor *) m_neighbors[i].next();
                        if ( i < temp->to_node )
                            eng = eng + m_smoothFn(i,temp->to_node, m_answer[i],m_answer[temp->to_node]);
                }
            }
        }
        
    }

    return(eng);
}



MRF::EnergyVal ICM::dataEnergy()
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


void ICM::setData(DataCostFn dcost)
{
    m_dataFn = dcost;
}

void ICM::setData(CostVal* data)
{
    m_D = data;
}


void ICM::setSmoothness(SmoothCostGeneralFn cost)
{
    m_smoothFn = cost;
}
void ICM::setSmoothness(CostVal* V)
{
    m_V = V;
}


void ICM::setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda)
{
    int i, j;
    CostVal cost;


    m_needToFreeV = 1;

    m_V = (CostVal *) new CostVal[m_nLabels*m_nLabels*sizeof(CostVal)];
    if (!m_V) { fprintf(stderr, "Not enough memory!\n"); exit(1); }


    for (i=0; i<m_nLabels; i++)
        for (j=i; j<m_nLabels; j++)
        {
            cost = (CostVal) ((smoothExp == 1) ? j - i : (j - i)*(j - i));
            if (cost > smoothMax) cost = smoothMax;
            m_V[i*m_nLabels + j] = m_V[j*m_nLabels + i] = cost*lambda;
        }

}


void ICM::setCues(CostVal* hCue, CostVal* vCue)
{
    m_horizWeights = hCue;
    m_vertWeights  = vCue;
}


void ICM::optimizeAlg(int nIterations)
{
    int x, y, i, j, n;
    Label* l;
    CostVal* dataPtr;
    CostVal *D = (CostVal *) new CostVal[m_nLabels];
    if ( !D ) {printf("\nNot enough memory, exiting");exit(0);}

    if ( !m_grid_graph) {printf("\nICM is not implemented for nongrids yet!");exit(1);}

    for ( ; nIterations > 0; nIterations --)
    {
        n = 0;
        l = m_answer;
        dataPtr = m_D;

        for (y=0; y<m_height; y++)
        for (x=0; x<m_width; x++, l++, dataPtr+=m_nLabels, n++)
        {
            // set array D
            if (m_dataType == FUNCTION)
            {
                for (i=0; i<m_nLabels; i++)
                {
                    D[i] = m_dataFn(x+y*m_width, i);
                }
            }
            else memcpy(D, dataPtr, m_nLabels*sizeof(CostVal));
            

            // add smoothness costs
            if (m_smoothType == FUNCTION)
            {
                if (x > 0)
                {
                    j = *(l-1);
                    for (i=0; i<m_nLabels; i++) D[i] += m_smoothFn(x+y*m_width-1, x+y*m_width, j, i);
                }
                if (y > 0)
                {
                    j = *(l-m_width);
                    for (i=0; i<m_nLabels; i++) D[i] += m_smoothFn(x+y*m_width-m_width,x+y*m_width , j, i);
                }
                if (x < m_width-1)
                {
                    j = *(l+1);
                    for (i=0; i<m_nLabels; i++) D[i] += m_smoothFn(x+y*m_width+1, x+y*m_width, i, j);
                }
                if (y < m_height-1)
                {
                    j = *(l+m_width);
                    for (i=0; i<m_nLabels; i++) D[i] += m_smoothFn(x+y*m_width+m_width, x+y*m_width, i, j);
                }
            }
            else
            {
                if (x > 0)
                {
                    j = *(l-1);
                    CostVal lambda = (m_varWeights) ? m_horizWeights[n-1] : 1;
                    for (i=0; i<m_nLabels; i++) D[i] += lambda * m_V[j*m_nLabels + i];
                }
                if (y > 0)
                {
                    j = *(l-m_width);
                    CostVal lambda = (m_varWeights) ? m_vertWeights[n-m_width] : 1;
                    for (i=0; i<m_nLabels; i++) D[i] += lambda * m_V[j*m_nLabels + i];
                }
                if (x < m_width-1)
                {
                    j = *(l+1);
                    CostVal lambda = (m_varWeights) ? m_horizWeights[n] : 1;
                    for (i=0; i<m_nLabels; i++) D[i] += lambda * m_V[j*m_nLabels + i];
                }
                if (y < m_height-1)
                {
                    j = *(l+m_width);
                    CostVal lambda = (m_varWeights) ? m_vertWeights[n] : 1;
                    for (i=0; i<m_nLabels; i++) D[i] += lambda * m_V[j*m_nLabels + i];
                }
            }

            // compute minimum of D, set new label for (x,y)
            CostVal D_min = D[0];
            *l = 0;
            for (i=1; i<m_nLabels; i++)
            {
                if (D_min > D[i])
                {
                    D_min = D[i];
                    *l = i;
                }
            }
        }
    }

    delete[] D;
}

