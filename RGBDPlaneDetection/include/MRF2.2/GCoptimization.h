/* GCoptimization.h */
/* Olga Veksler, 2005, olga@csd.uwo.ca */
/* Copyright 2005 Olga Veksler (olga@csd.uwo.ca)*/
/* email olga@csd.uwo.ca for any questions, suggestions and bug reports */

/* See file README.txt for instructions on using this library */

/****************************************************************************************************/
/*

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/



#ifndef __GCOPTIMIZATION_H__
#define __GCOPTIMIZATION_H__

#include "LinkedBlockList.h"
#include <assert.h>
#include "graph.h"
#include "energy.h"
#define m_datacost(pix,lab)     (m_datacost[(pix)*m_nLabels+(lab)] )
#define m_smoothcost(lab1,lab2) (m_smoothcost[(lab1)+(lab2)*m_nLabels] )
#define USE_MEMBER_FUNCTION 0
#define PASS_AS_PARAMETER   1


class GCoptimization:public MRF
{
public:

    ///////////////////////////////////////////////////////////////////////////////////////
    //    First define needed data types                                                 //
    ///////////////////////////////////////////////////////////////////////////////////////
    /* Type for the total energy calculation. Change it in Graph.h   */
    typedef Graph::flowtype EnergyType;

    /* Type for the individual terms in the energy function.Change it in Graph.h */
    typedef Graph::captype EnergyTermType;

    /* Type of label. Can be set to char, short, int, long */
    typedef int LabelType;

    /* Type for pixel. Can be set to  short, int, long */ 
    typedef int PixelType;

    
    ///////////////////////////////////////////////////////////////////////////////////////
    //    Next define constructors                                                       //
    ///////////////////////////////////////////////////////////////////////////////////////


    /* Use this constructor only for grid of size width by height.  If you use this constructor,  */
    /* 4 connected grid neigbhorhood structure is assumed.  num_labels is the number of labels,   */
    /* Labels are assumed to be in  {0, 1, 2,....num_labels-1}                                    */
    /* dataSetup and smoothSetup can take only 2 values, namely USE_MEMBER_FUNCTION (0) and       */
    /* PASS_AS_PARAMETER(1)                                                                       */ 
    /* In case dataSetup = USE_MEMBER_FUNCTION, to set the data cost you must use the             */
    /* member function setDataCost(pixel,label,cost).  If dataSetup = PASS_AS_PARAMETER,          */
    /* to set the data costs you must use any of the setData() functions.                         */
    /* In case smoothSetup = USE_MEMBER_FUNCTION, to set the smooth cost you must use the         */
    /* member function setSmoothCost(label1,label2,cost).  If smoothSetup = PASS_AS_PARAMETER,    */
    /* to set the smoothness costs you must use any of the setSmoothness() functions.             */
    GCoptimization(PixelType width,PixelType height,int num_labels,EnergyFunction *eng);


    /* This is the constructor for non-grid graphs. Everything is the same as in the constructor  */
    /* above, but neighborhood structure must  be specified by any of the two setNeighbors()      */
    /* functions */
    GCoptimization(PixelType num_pixels,int nLabels, EnergyFunction *eng);

    ~GCoptimization();

    /* Returns current label assigned to input pixel */
    inline LabelType getLabel(PixelType pixel){assert(pixel >= 0 && pixel < m_nPixels);return(m_labeling[pixel]);};

    /* Sets data cost of pixel to label */

   /* returns pointer to the answer */
    Label* getAnswerPtr(){return(m_labeling);}

    // Sets answer to zeros
    void clearAnswer();

    /* Makes pixel1 and pixel2 neighbors of each other. Can be called only 1 time for each         */
    /* unordered pair of pixels. Parameter weight can be used to set spacially varying terms       */
    /* If the desired penalty for neighboring pixels pixel1 and pixel2 is                          */
    /*  V(label1,label2) = weight*SmoothnessPenalty(label1,label2), then                           */
    /* member function setLabel should be called as: setLabel(pixel1,pixel2,weight)                */
    void setNeighbors(PixelType pixel1, PixelType pixel2, EnergyTermType weight);


    /* This function can be used to change the label of any pixel at any time      */
    inline void setLabel(PixelType pixel, LabelType label){
        assert(label >= 0 && label < m_nLabels && pixel >= 0 && pixel < m_nPixels);m_labeling[pixel] = label;};
    
    /* By default, the labels are visited in random order for both the swap and alpha-expansion moves */
    /* Use this function with boolean argument 0 to fix the order to be not random                    */
    /* Use this function with argumnet 1 to fix the order back to random                              */
    void setLabelOrder(bool RANDOM_LABEL_ORDER);


     void setParameters(int numParam, void *param);

    /* Returns Data Energy of current labeling */
    EnergyType dataEnergy();

    /* Returns Smooth Energy of current labeling */
    EnergyType smoothnessEnergy();

protected:
	void initializeAlg() {};

	/* This function is used to set the data term, and it can be used only if dataSetup = PASS_AS_PARAMETER */
    /* DataCost is an array s.t. the data cost for pixel p and  label l is stored at                        */
    /* DataCost[pixel*num_labels+l].  If the current neighborhood system is a grid, then                    */
    /* the data term for label l and pixel with coordinates (x,y) is stored at                              */ 
    /* DataCost[(x+y*width)*num_labels + l]. Thus the size of array DataCost is num_pixels*num_labels       */
     void setData(EnergyTermType *DataCost);

    /* This function is used to set the data term, and it can be used only if dataSetup = PASS_AS_PARAMETER */
    /* dataFn is a pointer to a function  f(Pixel p, Label l), s.t. the data cost for pixel p to have       */
    /* label l  is given by f(p,l) */
     void setData(DataCostFn dataFn);

    /* This function is used to set the smoothness term, and it can be used only if                         */
     /* smoothSetup = PASS_AS_PARAMETER                                                                      */
     /*  V is an array of costs, such that V(label1,label2)  is stored at V[label1+num_labels*label2]        */
     /* If graph is a grid, then using this  function only if the smooth costs are not spacially varying     */
     /* that is the smoothness penalty V depends only on labels, but not on pixels                           */
     void setSmoothness(EnergyTermType* V);

    void setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda);

     /* This function is used to set the smoothness term, and it can be used only if                         */
     /* smoothSetup = PASS_AS_PARAMETER AND the graph is a grid                                              */
     /* hCue(x,y,label1,label2) is the smoothness penalty for pixels (x,y) and (x+1,y) to have labels        */
     /* label1 and label2, correspondingly.                                                                  */
     /* vCue(x,y,label1,label2) is the smoothness penalty for pixels (x,y) and (x,y+1) to have labels        */
     /* label1 and label2, correspondingly.                                                                  */
     void setCues(CostVal* hCue, CostVal* vCue); // hCue and vCue must be arrays of size width*height in row major order. 
 

     /* This function is used to set the smoothness term, and it can be used only if                         */
     /* smoothSetup = PASS_AS_PARAMETER AND the graph is a grid                                              */
     /* horz_cost is a function f(x,y,label1,label2) such that smoothness penalty for neigboring pixels      */
     /* (x,y) and (x+1,y) to have labels, respectively, label1 and label2 is f(x,y,label1,label2)            */
     /* vert_cost is a function f(x,y,label1,label2) such that smoothness penalty for neigboring pixels      */
     /* (x,y) and (x,y+1) to have labels, respectively, label1 and label2 is f(x,y,label1,label2)            */
     void setSmoothness(SmoothCostGeneralFn cost);


    virtual void optimizeAlg(int nIterations) = 0;

    typedef struct NeighborStruct{
        PixelType  to_node;
        EnergyTermType weight;
    } Neighbor;


    LabelType *m_labeling;
    bool m_random_label_order;
    bool m_needToFreeV;
    EnergyTermType *m_datacost;
    EnergyTermType *m_smoothcost;
    EnergyTermType *m_vertWeights;
    EnergyTermType *m_horizWeights;
    LinkedBlockList *m_neighbors;

    LabelType *m_labelTable;
    PixelType *m_lookupPixVar;
    

    EnergyTermType m_weight;

    /* Pointers to function for energy terms */
    DataCostFn m_dataFnPix;
    SmoothCostGeneralFn m_smoothFnPix;

    void commonGridInitialization( PixelType width, PixelType height, int nLabels);
    void commonNonGridInitialization(PixelType num_pixels, int num_labels);
    void commonInitialization();    

    void scramble_label_table();

    EnergyType giveDataEnergyArray();
    EnergyType giveDataEnergyFnPix();

    EnergyType giveSmoothEnergy_G_ARRAY_VW();
    EnergyType giveSmoothEnergy_G_ARRAY();
    EnergyType giveSmoothEnergy_G_FnPix();
    EnergyType giveSmoothEnergy_NG_ARRAY();
    EnergyType giveSmoothEnergy_NG_FnPix();

    void add_t_links_ARRAY(Energy *e,Energy::Var *variables,int size,LabelType alpha_label);
    void add_t_links_FnPix(Energy *e,Energy::Var *variables,int size,LabelType alpha_label);
            

    void initialize_memory();
    void terminateOnError(bool error_condition,const char *message);

};


class Swap: public GCoptimization
{
public:
    Swap(PixelType width,PixelType height,int num_labels, EnergyFunction *eng);
    Swap(PixelType nPixels, int num_labels, EnergyFunction *eng);
    ~Swap();
    
    /* Peforms swap algorithm. Runs it the specified number of iterations */
    EnergyType swap(int max_num_iterations);

    /* Peforms swap algorithm. Runs it until convergence */
    EnergyType swap();

    /* Peforms one iteration (one pass over all labels)  of swap algorithm.*/
    EnergyType oneSwapIteration();

    /* Peforms  swap on a pair of labels, specified by the input parameters alpha_label, beta_label */
    EnergyType alpha_beta_swap(LabelType alpha_label, LabelType beta_label);

protected:
    void optimizeAlg(int nIterations);

private:
    PixelType *m_pixels;
    void set_up_swap_energy_G_ARRAY_VW(int size, LabelType alpha_label,LabelType beta_label,PixelType *pixels,Energy* e, Energy::Var *variables);
    void set_up_swap_energy_G_ARRAY(int size, LabelType alpha_label,LabelType beta_label,PixelType *pixels,Energy* e, Energy::Var *variables);
    void set_up_swap_energy_G_FnPix(int size, LabelType alpha_label,LabelType beta_label,PixelType *pixels,Energy* e, Energy::Var *variables);
    void set_up_swap_energy_NG_ARRAY(int size, LabelType alpha_label,LabelType beta_label,PixelType *pixels,Energy* e, Energy::Var *variables);     
    void set_up_swap_energy_NG_FnPix(int size, LabelType alpha_label,LabelType beta_label,PixelType *pixels,Energy* e, Energy::Var *variables);     
    EnergyType start_swap(int max_iterations);
    void add_t_links_ARRAY_swap(Energy *e,Energy::Var *variables,int size,LabelType alpha_label,LabelType beta_label,PixelType *pixels);
    void add_t_links_FnPix_swap(Energy *e,Energy::Var *variables,int size,LabelType alpha_label,LabelType beta_label,PixelType *pixels);
    void perform_alpha_beta_swap(LabelType alpha_label, LabelType beta_label);

};

class Expansion: public GCoptimization
{
public:
    Expansion(PixelType width,PixelType height,int num_labels,EnergyFunction *eng):GCoptimization(width,height,num_labels,eng){};
    Expansion(PixelType nPixels, int num_labels,EnergyFunction *eng):GCoptimization(nPixels,num_labels,eng){};


    /* Peforms expansion algorithm. Runs the number of iterations specified by max_num_iterations */
    /* Returns the total energy of labeling   */
    EnergyType expansion(int max_num_iterations);

    /* Peforms expansion algorithm. Runs it until convergence */
    /* Returns the total energy of labeling   */
    EnergyType expansion();

    /* Peforms one iteration (one pass over all labels)  of expansion algorithm.*/
    EnergyType oneExpansionIteration();

    /* Peforms  expansion on one label, specified by the input parameter alpha_label */
    EnergyType alpha_expansion(LabelType alpha_label);
    
protected:
    void optimizeAlg(int nIterations);


private:
    void set_up_expansion_energy_G_ARRAY_VW(int size, LabelType alpha_label,Energy* e, Energy::Var *variables);
    void set_up_expansion_energy_G_ARRAY(int size, LabelType alpha_label,Energy* e, Energy::Var *variables);
    void set_up_expansion_energy_G_FnPix(int size, LabelType alpha_label,Energy* e, Energy::Var *variables);
    void set_up_expansion_energy_NG_ARRAY(int size, LabelType alpha_label,Energy* e, Energy::Var *variables);       
    void set_up_expansion_energy_NG_FnPix(int size, LabelType alpha_label,Energy* e, Energy::Var *variables);       
    void perform_alpha_expansion(LabelType label);  
    EnergyType start_expansion(int max_iterations);

};




#endif

