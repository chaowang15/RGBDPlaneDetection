// (C) 2002 Marshall Tappen, MIT AI Lab

#ifndef _reg_h
#define _reg_h

#define FLOATTYPE float
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3
#include "MaxProdBP.h"

class MaxProdBP;
class OneNodeCluster
{
public:
  OneNodeCluster();
  
  //static const int UP = 0, DOWN = 1, LEFT = 2, RIGHT = 3;
  static int numStates;
  
  FLOATTYPE   *receivedMsgs[4],
              *nextRoundReceivedMsgs[4],
              *localEv;


//   FLOATTYPE *psi[4];
  void ComputeMsgRight(FLOATTYPE *msgDest, int r, int c, MaxProdBP *mrf);
  void ComputeMsgUp(FLOATTYPE *msgDest, int r, int c, MaxProdBP *mrf);

  void ComputeMsgLeft(FLOATTYPE *msgDest, int r, int c, MaxProdBP *mrf);

  void ComputeMsgDown(FLOATTYPE *msgDest, int r, int c, MaxProdBP *mrf);

  void getBelief(FLOATTYPE *beliefVec);
  int getBeliefMaxInd();
  
  int msgChange(FLOATTYPE thresh);

  void deliverMsgs();
    
};

void initOneNodeMsgMem(OneNodeCluster *nodeArray, FLOATTYPE *memChunk, const int numNodes, 
                       const int msgChunkSize);

void computeMessagesUpDown(OneNodeCluster *nodeArray, const int numCols, const int numRows,
                           const int currCol, const FLOATTYPE alpha, MaxProdBP *mrf);
void computeMessagesLeftRight(OneNodeCluster *nodeArray, const int numCols, const int numRows,
                              const int currRow, const FLOATTYPE alpha, MaxProdBP *mrf);

void computeOneNodeMessagesPeriodic(OneNodeCluster *nodeTopArray, OneNodeCluster *nodeBotArray,
                                    const int numCols, const FLOATTYPE alpha);


#endif
