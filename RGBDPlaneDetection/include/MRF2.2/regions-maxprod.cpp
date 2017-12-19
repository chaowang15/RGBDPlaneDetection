// (C) 2002 Marshall Tappen, MIT AI Lab mtappen@mit.edu

#include <limits>
#include <stdio.h>
#include "MaxProdBP.h"
#include "assert.h"
int numIterRun;
// Some of the GBP code has been disabled here

#define mexPrintf printf
#define mexErrMsgTxt printf
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3


OneNodeCluster::OneNodeCluster()
{
}

int OneNodeCluster::numStates;

FLOATTYPE vec_min(FLOATTYPE *vec, int length)
{

  FLOATTYPE min = vec[0];
  for(int i = 0; i < length; i++)
    if(vec[i] < min)
      min = vec[i];

  return min;
  }

FLOATTYPE vec_max(FLOATTYPE *vec, int length)
{

  FLOATTYPE max = vec[0];
  for(int i = 0; i < length; i++)
    if(vec[i] > max)
      max = vec[i];

  return max;
}

void getPsiMat(OneNodeCluster &/*cluster*/, FLOATTYPE *&destMatrix, 
	       int r, int c, MaxProdBP *mrf, int direction, FLOATTYPE &var_weight)
{
  int mrfHeight = mrf->getHeight();
  int mrfWidth = mrf->getWidth();
  int numLabels = mrf->getNLabels();
  int x=c;
  int y=r;
  int i;
  
  FLOATTYPE *currMatrix = mrf->getScratchMatrix();
  if(mrf->getSmoothType() != MRF::FUNCTION)
  {
    if(((direction==UP) &&(r==0)) ||
       ((direction==DOWN) &&(r==(mrfHeight-1)))||
       ((direction==LEFT) &&(c==0))||
       ((direction==RIGHT) &&(c==(mrfWidth-1))))
    {
      for( i=0; i < numLabels * numLabels; i++)
      {
	currMatrix[i] = 0;
      }
      
    }
    else
    {
	MRF::CostVal weight_mod = 1;
      if(mrf->varWeights())
      {
	if(direction==LEFT)
	  weight_mod = mrf->getHorizWeight(r,c-1);
	else if (direction ==RIGHT)
	  weight_mod = mrf->getHorizWeight(r,c);
	else if (direction ==UP)
	  weight_mod = mrf->getVertWeight(r-1,c);
	else if (direction == DOWN)
	  weight_mod = mrf->getVertWeight(r,c);
	
      }
      for( i = 0; i < numLabels*numLabels; i++)
      {
	if(weight_mod!=1)
	{
	  currMatrix[i] = FLOATTYPE(mrf->m_V[i]*weight_mod);
	}
	else
	  currMatrix[i] = FLOATTYPE(mrf->m_V[i]);
      }
      destMatrix = currMatrix;
      var_weight = (float)weight_mod;
    }
    
  }
  else
  {
    if(((direction==UP) &&(r==0)) ||
       ((direction==DOWN) &&(r==(mrfHeight-1)))||
       ((direction==LEFT) &&(c==0))||
       ((direction==RIGHT) &&(c==(mrfWidth-1))))
    {
      for( i=0; i < numLabels * numLabels; i++)
      {
	currMatrix[i] = 0;
      }
    }
    else
    {
      for( i = 0; i < numLabels; i++)
      {
	for(int j = 0; j < numLabels; j++)
	{
	  MRF::CostVal cCost;
	  if(direction==LEFT)
	    cCost = mrf->m_smoothFn(x+y*mrf->m_width,
				    x+y*mrf->m_width-1 , j, i);
	  
	  else if (direction ==RIGHT)
	    cCost = mrf->m_smoothFn(x+y*mrf->m_width,
				    x+y*mrf->m_width+1 , i, j);
	  else if (direction ==UP)
	    cCost = mrf->m_smoothFn(x+y*mrf->m_width,
				    x+(y-1)*mrf->m_width , j, i);
	  else if (direction == DOWN)
	    cCost = mrf->m_smoothFn(x+y*mrf->m_width,
				    x+(y+1)*mrf->m_width , i, j);
	  else 
	  {
	    cCost = mrf->m_smoothFn(x+y*mrf->m_width,
				    x+(y+1)*mrf->m_width-1 , j, i);
	    assert(0);
	  }
	  
	  currMatrix[i*numLabels+j] = (float)cCost;
	}
      }
    }
  }
  destMatrix = currMatrix;
}


void getVarWeight(OneNodeCluster &/*cluster*/,  int r, int c, MaxProdBP *mrf, int direction, FLOATTYPE &var_weight)
{
  MRF::CostVal weight_mod = 1;
  if(mrf->varWeights())
  {
    if(direction==LEFT)
      weight_mod = mrf->getHorizWeight(r,c-1);
    else if (direction ==RIGHT)
      weight_mod = mrf->getHorizWeight(r,c);
    else if (direction ==UP)
      weight_mod = mrf->getVertWeight(r-1,c);
    else if (direction == DOWN)
      weight_mod = mrf->getVertWeight(r,c);
    
  }
  var_weight = (FLOATTYPE) weight_mod;
  //  printf("%d\n",weight_mod);
}

void initOneNodeMsgMem(OneNodeCluster *nodeArray, FLOATTYPE *memChunk, 
		       const int numNodes, const int msgChunkSize)
{
  FLOATTYPE *currPtr = memChunk;
  OneNodeCluster *currNode = nodeArray;
  FLOATTYPE *nextRoundChunk = new FLOATTYPE[nodeArray[1].numStates];
  // MEMORY LEAK? where does this ever get deleted??
  for(int i = 0; i < numNodes; i++)
  {

    currNode->receivedMsgs[0] = currPtr; currPtr+=msgChunkSize;
    currNode->receivedMsgs[1] = currPtr; currPtr+=msgChunkSize;
    currNode->receivedMsgs[2] = currPtr; currPtr+=msgChunkSize;
    currNode->receivedMsgs[3] = currPtr; currPtr+=msgChunkSize;

    currNode->nextRoundReceivedMsgs[0] = nextRoundChunk;
    currNode->nextRoundReceivedMsgs[1] = nextRoundChunk;
    currNode->nextRoundReceivedMsgs[2] = nextRoundChunk;
    currNode->nextRoundReceivedMsgs[3] = nextRoundChunk;

    currNode++;
  }
}

inline void l1_dist_trans_comp(FLOATTYPE smoothMax, FLOATTYPE c, FLOATTYPE* tmpMsgDest, FLOATTYPE * msgProd, int numStates)
{
  int q;
  for(int i=0; i < numStates; i++)
    tmpMsgDest[i] = msgProd[i];

  for (q = 1; q <= numStates-1; q++)
  {
    if (tmpMsgDest[q]  > tmpMsgDest[q-1]+c)
      tmpMsgDest[q] = tmpMsgDest[q-1]+c;
  }

  for (q = numStates-2; q >= 0; q--)
  {
    if (tmpMsgDest[q]  > tmpMsgDest[q+1]+c)
      tmpMsgDest[q] = tmpMsgDest[q+1]+c;
  }

  FLOATTYPE minPotts = msgProd[0] + smoothMax;
  for(q = 0; q <= numStates -1; q++)
  {
    if((msgProd[q]+smoothMax) < minPotts)
      minPotts = msgProd[q]+smoothMax;
  }
   for(q = 0; q <= numStates -1; q++)
   {
     if((tmpMsgDest[q]) > minPotts)
       tmpMsgDest[q] = minPotts;
     tmpMsgDest[q] = -tmpMsgDest[q];
   }
   //   printf("%f %f %f\n",smoothMax,c,minPotts);
}



inline void l2_dist_trans_comp(FLOATTYPE smoothMax, FLOATTYPE c, FLOATTYPE* tmpMsgDest, FLOATTYPE * msgProd, int numStates)
{


  FLOATTYPE *z = new FLOATTYPE[numStates];
  int       *v = new int[numStates];
  int j=0;
  FLOATTYPE INFINITY_ =std::numeric_limits<float>::infinity();

  z[0]=-1*INFINITY_;
  z[1]=INFINITY_;
  v[0]=0;
  int q;
  if(c==0)
  {
    FLOATTYPE minVal = msgProd[0];

    for (q = 0; q < numStates; q++)
    {
      if(msgProd[q] < minVal)
	minVal = msgProd[q];
    }
    for (q = 0; q < numStates; q++)
    {
      tmpMsgDest[q]=-minVal;
    }
    delete [] z;
    delete [] v;
    return;
  }

  for(q = 1; q <= numStates -1; q++)
  {
    FLOATTYPE s;
    while(   (s = ((msgProd[q] + c*q*q) - (msgProd[v[j]] + c*v[j]*v[j]))/
	    (2*c*q-2*c*v[j])) <= z[j])
    { 
      j -=1;
    }
    
    j += 1;
    v[j] = q;
    z[j] = s;
    z[j+1] = INFINITY_;
    
    
  }
  j = 0;
  FLOATTYPE minPotts = msgProd[0] + smoothMax;
  for(q = 0; q <= numStates -1; q++)
  {
    while(z[j+1] < q)
    {
      j +=1;
    }
    tmpMsgDest[q] = c*(q-v[j])*(q-v[j]) + msgProd[v[j]];
    if((msgProd[q]+smoothMax) < minPotts)
      minPotts = msgProd[q]+smoothMax;
  }
  for(q = 0; q <= numStates -1; q++)
  {
    if((tmpMsgDest[q]) > minPotts)
      tmpMsgDest[q] = minPotts;
    tmpMsgDest[q] = -tmpMsgDest[q];
  }
  delete [] z;
  delete [] v;

}

void OneNodeCluster::ComputeMsgRight(FLOATTYPE *msgDest, int r, int c, MaxProdBP *mrf)
{


  FLOATTYPE *nodeLeftMsg = receivedMsgs[LEFT],
    *nodeDownMsg = receivedMsgs[DOWN],
    *nodeUpMsg =    receivedMsgs[UP];

  FLOATTYPE weight_mod;
  getVarWeight(*this,r,c,mrf,RIGHT,weight_mod);
  
  FLOATTYPE     *tmpMsgDest =msgDest;

  
  if(mrf->m_type==MaxProdBP::L1 || mrf->m_type==MaxProdBP::L2)
  {
    FLOATTYPE *msgProd = new FLOATTYPE[numStates];
    const FLOATTYPE lambda = (FLOATTYPE)mrf->m_lambda;
    const FLOATTYPE smoothMax = (FLOATTYPE)mrf->m_smoothMax;
    for(int leftNodeInd = 0; leftNodeInd < numStates; leftNodeInd++)
    {
      msgProd[leftNodeInd]  = -nodeLeftMsg[leftNodeInd] + 
	-nodeUpMsg[leftNodeInd] + 
	-nodeDownMsg[leftNodeInd] + localEv[leftNodeInd];
    }
    if(mrf->m_type==MaxProdBP::L1)
    {
      l1_dist_trans_comp( weight_mod*smoothMax*lambda, lambda*weight_mod, tmpMsgDest, msgProd, numStates);
    }
    else
    {
      l2_dist_trans_comp( weight_mod*smoothMax*lambda, lambda*weight_mod, tmpMsgDest, msgProd, numStates);
    }
    delete [] msgProd;
  }
  else if ((mrf->getSmoothType()==MRF::FUNCTION)||(mrf->getSmoothType()==MRF::ARRAY))
  {
    FLOATTYPE *psiMat, var_weight;
  
    getPsiMat(*this,psiMat,r,c,mrf,RIGHT, var_weight);
    FLOATTYPE *cmessage = msgDest;    
    for(int rightNodeInd = 0; rightNodeInd < numStates; rightNodeInd++)
    {
      
      *cmessage = 0;
      for(int leftNodeInd = 0; leftNodeInd < numStates; leftNodeInd++)
      {
	FLOATTYPE tmp = 	nodeLeftMsg[leftNodeInd] + 
	  nodeUpMsg[leftNodeInd] + 
	  nodeDownMsg[leftNodeInd]
	  - localEv[leftNodeInd] 
	  - psiMat[leftNodeInd * numStates + rightNodeInd];

	if((tmp > *cmessage)||(leftNodeInd==0))
	  *cmessage = tmp;
	
      }
      cmessage++;
    }
  }
  else {
      fprintf(stderr, "not implemented!\n");
      exit(1);
  }

    
  FLOATTYPE max = msgDest[0];
  for(int i=0; i < numStates; i++)
    msgDest[i] -= max;
}


// This means, "Compute the message to send left."

void OneNodeCluster::ComputeMsgLeft(FLOATTYPE *msgDest, int r, int c, MaxProdBP *mrf)
{

  FLOATTYPE *nodeRightMsg = receivedMsgs[RIGHT],
    *nodeDownMsg = receivedMsgs[DOWN],
    *nodeUpMsg =    receivedMsgs[UP];


  int do_dist=(int)(mrf->getSmoothType()==MRF::THREE_PARAM);
  FLOATTYPE *tmpMsgDest=msgDest;
  if(do_dist)
  {
    FLOATTYPE weight_mod;
    getVarWeight(*this,r,c,mrf,LEFT, weight_mod);
    
    FLOATTYPE *msgProd = new FLOATTYPE[numStates];
      
    const FLOATTYPE lambda = (FLOATTYPE)mrf->m_lambda;

    const FLOATTYPE smoothMax = (FLOATTYPE)mrf->m_smoothMax;

    for(int rightNodeInd = 0; rightNodeInd < numStates; rightNodeInd++)
    {
      msgProd[rightNodeInd]  = -nodeRightMsg[rightNodeInd] + 
	-nodeUpMsg[rightNodeInd] + 
	-nodeDownMsg[rightNodeInd] 
	+localEv[rightNodeInd] ;
    }
    if(mrf->m_smoothExp==1)
    {      l1_dist_trans_comp( weight_mod*smoothMax*lambda, lambda*weight_mod, tmpMsgDest, msgProd, numStates);
    }
    else
      l2_dist_trans_comp( weight_mod*smoothMax*lambda, lambda*weight_mod, tmpMsgDest, msgProd, numStates);
    
    delete [] msgProd;
  }
  else if ((mrf->getSmoothType()==MRF::FUNCTION)||(mrf->getSmoothType()==MRF::ARRAY))
  {
    FLOATTYPE *psiMat, var_weight;
    
    getPsiMat(*this,psiMat,r,c,mrf,LEFT, var_weight);
  
    FLOATTYPE *cmessage = msgDest;
    
    for(int leftNodeInd = 0; leftNodeInd < numStates; leftNodeInd++)
    {
      
      *cmessage = 0;
      for(int rightNodeInd = 0; rightNodeInd < numStates; rightNodeInd++)
      {
	FLOATTYPE tmp = 	nodeRightMsg[rightNodeInd] + 
	  nodeUpMsg[rightNodeInd] + 
	  nodeDownMsg[rightNodeInd] 
	  -localEv[rightNodeInd] 
	  - psiMat[leftNodeInd * numStates + rightNodeInd] ; 
       
	if((tmp > *cmessage)||(rightNodeInd==0))
	  *cmessage = tmp;

      }
      cmessage++;
    }
    
  }
  else
    assert(0);

  

//   FLOATTYPE max = vec_max(msgDest,numStates);
  FLOATTYPE max = msgDest[0];

  for(int i=0; i < numStates; i++)
    msgDest[i] -= max;

}

void OneNodeCluster::ComputeMsgUp(FLOATTYPE *msgDest, int r, int c, MaxProdBP *mrf)
{
  FLOATTYPE *nodeRightMsg = receivedMsgs[RIGHT],
    *nodeDownMsg = receivedMsgs[DOWN],
    *nodeLeftMsg =    receivedMsgs[LEFT];
 



  int do_dist=(int)(mrf->getSmoothType()==MRF::THREE_PARAM);
  if(do_dist)
  {
    FLOATTYPE weight_mod;
    getVarWeight(*this,r,c,mrf,UP,weight_mod);
    
    FLOATTYPE *tmpMsgDest = msgDest;
    FLOATTYPE *msgProd = new FLOATTYPE[numStates];
    
    const FLOATTYPE lambda = (FLOATTYPE)mrf->m_lambda;
    
    const FLOATTYPE smoothMax = (FLOATTYPE)mrf->m_smoothMax;
    
    for(int downNodeInd = 0; downNodeInd < numStates; downNodeInd++)
    {
      msgProd[downNodeInd]  = -nodeRightMsg[downNodeInd] +
	-nodeLeftMsg[downNodeInd] + 
	-nodeDownMsg[downNodeInd] + 
	+localEv[downNodeInd] ;
    }
    //       printf("%f %f %f %f\n",nodeLeftMsg[leftNodeInd] ,
    // 	     nodeUpMsg[leftNodeInd] , 
    // 	     nodeDownMsg[leftNodeInd] ,localEv[leftNodeInd]);
    if(mrf->m_smoothExp==1)
    {      l1_dist_trans_comp( weight_mod*smoothMax*lambda, lambda*weight_mod, tmpMsgDest, msgProd, numStates);
    }
    else
      l2_dist_trans_comp( weight_mod*smoothMax*lambda, lambda*weight_mod, tmpMsgDest, msgProd, numStates);
    delete [] msgProd;
  }
  else if ((mrf->getSmoothType()==MRF::FUNCTION)||(mrf->getSmoothType()==MRF::ARRAY))
  {
    FLOATTYPE *psiMat, var_weight;
    
    getPsiMat(*this,psiMat,r,c,mrf,UP, var_weight);
    
    FLOATTYPE *cmessage = msgDest;
    
    for(int upNodeInd = 0; upNodeInd < numStates; upNodeInd++)
    {
      
      *cmessage = 0;
      for(int downNodeInd = 0; downNodeInd < numStates; downNodeInd++)
      {
	FLOATTYPE tmp = nodeRightMsg[downNodeInd] +
	  nodeLeftMsg[downNodeInd] + 
	  nodeDownMsg[downNodeInd] + 
	  -localEv[downNodeInd] 
	  -psiMat[upNodeInd * numStates + downNodeInd] ; 
	
	if((tmp > *cmessage)||(downNodeInd==0))
	  *cmessage = tmp;
	

      }
      cmessage++;
    }
  }
  else
    assert(0);
  FLOATTYPE max = msgDest[0];
  //  FLOATTYPE max = vec_max(msgDest,numStates);
  for(int i=0; i < numStates; i++)
    msgDest[i]  -=max;
}

void OneNodeCluster::ComputeMsgDown(FLOATTYPE *msgDest, int r, int c, MaxProdBP *mrf)
{

  FLOATTYPE *nodeRightMsg = receivedMsgs[RIGHT],
    *nodeUpMsg = receivedMsgs[UP],
    *nodeLeftMsg =    receivedMsgs[LEFT];
 
  int do_dist=(int)(mrf->getSmoothType()==MRF::THREE_PARAM);
  if(do_dist)
  {

    FLOATTYPE weight_mod;
    getVarWeight(*this,r,c,mrf,DOWN,weight_mod);
    
    FLOATTYPE *tmpMsgDest = msgDest;
    FLOATTYPE *msgProd = new FLOATTYPE[numStates];
    
    const FLOATTYPE lambda = (FLOATTYPE)mrf->m_lambda;
    
    const FLOATTYPE smoothMax = (FLOATTYPE)mrf->m_smoothMax;
    
    for(int upNodeInd = 0; upNodeInd < numStates; upNodeInd++)
    {
      msgProd[upNodeInd]  = -nodeRightMsg[upNodeInd] +
	-nodeLeftMsg[upNodeInd] + 
      -nodeUpMsg[upNodeInd] +
      +localEv[upNodeInd] ;
    }

    if(mrf->m_smoothExp==1)
    {      l1_dist_trans_comp( weight_mod*smoothMax*lambda, lambda*weight_mod, tmpMsgDest, msgProd, numStates);
    }
    else
      l2_dist_trans_comp( weight_mod*smoothMax*lambda, lambda*weight_mod, tmpMsgDest, msgProd, numStates);
    delete [] msgProd;

  }
  else if((mrf->getSmoothType()==MRF::FUNCTION)||(mrf->getSmoothType()==MRF::ARRAY))
  {
    FLOATTYPE *psiMat, var_weight;
    
    getPsiMat(*this,psiMat,r,c,mrf,DOWN, var_weight);
    
    FLOATTYPE *cmessage = msgDest;
    
    for(int downNodeInd = 0; downNodeInd < numStates; downNodeInd++)
    {
      
      *cmessage = 0;
      for(int upNodeInd = 0; upNodeInd < numStates; upNodeInd++)
      {
	FLOATTYPE tmp = 	nodeRightMsg[upNodeInd] +
	  nodeLeftMsg[upNodeInd] + 
	  nodeUpMsg[upNodeInd] +
	  -localEv[upNodeInd] 
	-psiMat[upNodeInd * numStates + downNodeInd] ; 
	
	if((tmp > *cmessage)||(upNodeInd==0))
	  *cmessage = tmp;
	
      }
      cmessage++;
    }
    
  }
  else
    assert(0);

  FLOATTYPE max = msgDest[0];
  //  FLOATTYPE max = vec_max(msgDest,numStates);
  for(int i=0; i < numStates; i++)
    msgDest[i] -=max;
  

}





void OneNodeCluster::getBelief(FLOATTYPE *beliefVec)
{
	for(int i = 0; i < numStates; i++)
	{
		beliefVec[i] = receivedMsgs[UP][i] + receivedMsgs[DOWN][i] +
		  receivedMsgs[LEFT][i] + receivedMsgs[RIGHT][i] - localEv[i];
	}

}

int OneNodeCluster::getBeliefMaxInd()
{
  FLOATTYPE currBelief,bestBelief;
  int bestInd = 0;
  {
    int i = 0;
    bestBelief = receivedMsgs[UP][i] + receivedMsgs[DOWN][i] +
      receivedMsgs[LEFT][i] + receivedMsgs[RIGHT][i] - localEv[i];
  }
  for(int i = 1; i < numStates; i++)
  {
    currBelief = receivedMsgs[UP][i] + receivedMsgs[DOWN][i] +
      receivedMsgs[LEFT][i] + receivedMsgs[RIGHT][i] - localEv[i];
    if(currBelief > bestBelief)
    {
      bestInd=i;
      bestBelief = currBelief;
    }

  }
  return bestInd;

}


void computeMessagesLeftRight(OneNodeCluster *nodeArray, const int numCols, const int /*numRows*/, const int currRow, const FLOATTYPE alpha, MaxProdBP *mrf)
{
  const int numStates = OneNodeCluster::numStates;
  const FLOATTYPE omalpha = 1.0f - alpha;
  int i;
  int col;
  for( col = 0; col < numCols-1; col++)
  {
    nodeArray[currRow * numCols + col].ComputeMsgRight(nodeArray[currRow * numCols + col+1].nextRoundReceivedMsgs[LEFT],currRow, col, mrf);
    for(i = 0; i < numStates; i++)
    {
      nodeArray[currRow * numCols + col+1].receivedMsgs[LEFT][i] = 
	omalpha * nodeArray[currRow * numCols + col+1].receivedMsgs[LEFT][i] + 
	alpha * nodeArray[currRow * numCols + col+1].nextRoundReceivedMsgs[LEFT][i];
    }
  } 
  for( col = numCols-1; col > 0; col--)
  {
    nodeArray[currRow * numCols + col].ComputeMsgLeft(nodeArray[currRow * numCols + col-1].nextRoundReceivedMsgs[RIGHT],currRow, col, mrf);
    for(i = 0; i < numStates; i++)
    {
      nodeArray[currRow * numCols + col-1].receivedMsgs[RIGHT][i] = 
	omalpha * nodeArray[currRow * numCols + col-1].receivedMsgs[RIGHT][i] + 
	alpha * nodeArray[currRow * numCols + col-1].nextRoundReceivedMsgs[RIGHT][i];
    }
  } 

}

void computeMessagesUpDown(OneNodeCluster *nodeArray, const int numCols, const int numRows, const int currCol, const FLOATTYPE alpha, MaxProdBP *mrf)
{
  const int numStates = OneNodeCluster::numStates;
  const FLOATTYPE omalpha = 1.0f - alpha;
  int i;
  int row;
  for(row = 0; row < numRows-1; row++)
  {
    nodeArray[row * numCols + currCol].ComputeMsgDown(nodeArray[(row+1) * numCols + currCol].nextRoundReceivedMsgs[UP],row, currCol, mrf);
    for(i = 0; i < numStates; i++)
    {
      nodeArray[(row+1) * numCols + currCol].receivedMsgs[UP][i] = 
	omalpha * nodeArray[(row+1) * numCols + currCol].receivedMsgs[UP][i] + 
	alpha * nodeArray[(row+1) * numCols + currCol].nextRoundReceivedMsgs[UP][i];
    }
  } 
  for( row = numRows-1; row > 0; row--)
  {
    nodeArray[row * numCols + currCol].ComputeMsgUp(nodeArray[(row-1) * numCols + currCol].nextRoundReceivedMsgs[DOWN], row, currCol, mrf);
    for(i = 0; i < numStates; i++)
    {
      nodeArray[(row-1) * numCols + currCol].receivedMsgs[DOWN][i] = 
	omalpha * nodeArray[(row-1) * numCols + currCol].receivedMsgs[DOWN][i] + 
	alpha * nodeArray[(row-1) * numCols + currCol].nextRoundReceivedMsgs[DOWN][i];
    }
  } 

}


