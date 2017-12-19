#ifndef __ICM_H__
#define __ICM_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "mrf.h"
#include "LinkedBlockList.h"


class ICM : public MRF{
public:
    ICM(int width, int height, int nLabels, EnergyFunction *eng);
    ICM(int nPixels, int nLabels,EnergyFunction *eng);
    ~ICM();
    void setNeighbors(int pix1, int pix2, CostVal weight);
    Label getLabel(int pixel){return(m_answer[pixel]);};
    void setLabel(int pixel,Label label){m_answer[pixel] = label;};
    Label* getAnswerPtr(){return(m_answer);};
    void clearAnswer();
    void setParameters(int /*numParam*/, void * /*param*/){printf("No optional parameters to set"); exit(1);}
    EnergyVal smoothnessEnergy();
    EnergyVal dataEnergy();

protected:
    void setData(DataCostFn dcost); 
    void setData(CostVal* data);    
    void setSmoothness(SmoothCostGeneralFn cost);
    void setSmoothness(CostVal* V);
    void setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda);
    void setCues(CostVal* hCue, CostVal* vCue); 
    void initializeAlg();
    void optimizeAlg(int nIterations);

private:
    Label *m_answer;
    CostVal *m_V;
    CostVal *m_D;
    CostVal *m_horizWeights;
    CostVal *m_vertWeights;
    DataCostFn m_dataFn;
    SmoothCostGeneralFn m_smoothFn;
    bool m_needToFreeV;

    typedef struct NeighborStruct{
        int     to_node;
        CostVal weight;
    } Neighbor;

    LinkedBlockList *m_neighbors;
};


#endif /*  __ICM_H__ */


