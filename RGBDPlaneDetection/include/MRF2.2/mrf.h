/* Copyright Olga Veksler, Ramin Zabih, Vladimir Kolmogorov, and Daniel Scharstein
 * Send any questions to schar@middlebury.edu
 */

#ifndef __MRF_H__
#define __MRF_H__
#include <stdio.h>

class EnergyFunction;

class MRF
{
public:

    // *********** CONSTRUCTORS/DESTRUCTOR
    // Constructor. After you call this, you must call setData and setSmoothness            
    // Use this constructor for 2D grid graphs of size width by height Standard 4-connected 
    // neighborhood system is assumed. Labels are in the range 0,1,...nLabels - 1           
    // Width is in  the range 0,1,...width-1 and height is in the range 0,1,...height-1     
    // Input parameter eng specifies the data and smoothness parts of the energy
    // For 2D grids, since 4 connected neighborhood structure is assumed, this 
    // fully specifies the energy
    MRF(int width, int height, int nLabels, EnergyFunction *eng);

    // Use this constructor for a general neighborhood system. Pixels are in the range      
    // 0,1,..nPixels-1, and labels are in the range 0,1,...,nLabels-1 
    // Input parameter eng specifies the data and smoothness parts of the energy
    // after this constructor you need to call setNeighbors() to specify the neighborhood system
    MRF(int nPixels, int nLabels, EnergyFunction *eng);

    virtual ~MRF() { }
    
    // Returns true if energy function has been specified, returns false otherwise
    // By default, it always returns true. Can be modified by the supplier of
    // optimization algorithm                                                      
    virtual int isValid(){return true;};  
    

    // *********** EVALUATING THE ENERGY
    typedef int Label;
	
	//  [12/14/2017 chao]
	// Original code uses int as type of energy and cost value
    //typedef int EnergyVal;        /* The total energy of a labeling */
    //typedef int CostVal;          /* costs of individual terms of the energy */
	
	// New code
	typedef double EnergyVal;
	typedef double CostVal;
 
    EnergyVal totalEnergy();      /* returns energy of current labeling */
    virtual EnergyVal dataEnergy() = 0;        /* returns the data part of the energy */
    virtual EnergyVal smoothnessEnergy() = 0;  /* returns the smoothness part of the energy */

    //Functional representation for data costs
    typedef CostVal (*DataCostFn)(int pix, Label l); 

    // Functional representation for the general cost function type 
    typedef CostVal (*SmoothCostGeneralFn)(int pix1, int pix2,  Label l1, Label l2); 

    // For general smoothness functions, some implementations try to cache all function values in an array
    // for efficiency.  To prevent this, call the following function before calling initialize():
    void dontCacheSmoothnessCosts() {m_allocateArrayForSmoothnessCostFn = false;}

    // Use this function only for non-grid graphs. Sets pix1 and pix2 to be neighbors 
    // with the specified weight.  Can be called ONLY once for each pair of pixels    
    // That is if pixel1 and pixel2 are neihbors, call  either setNeighbors(pixel1,pixel2,weight) 
    // or setNeighbors(pixel2,pixel1,weight), but NOT BOTH                                     
    virtual void setNeighbors(int pix1, int pix2, CostVal weight)= 0;


    void initialize();

    // Runs optimization for nIterations. Input parameter time returns the time it took to
    // perform nIterations of optimization
    void optimize(int nIterations, float& time);
    virtual void optimizeAlg(int nIterations)=0;

    // *********** ACCESS TO SOLUTION
    // Returns pointer to array of size nPixels. Client may then read/write solution (but not deallocate array).
    virtual Label* getAnswerPtr()= 0;
    // returns the label of the input pixel
    virtual Label getLabel(int pixel)= 0;
    // sets label of a pixel
    virtual void setLabel(int pixel,Label label)= 0;
    // sets all the labels to zero
    virtual void clearAnswer()  = 0;

    // use this function to pass any parameters to optimization algorithm. 
    // The first argument is the number of passed, parameters  and
    // the second argument is the pointer to the array of parameters
    virtual void setParameters(int numParam, void *param) = 0;

    // This function returns lower bound computed by the algorithm (if any)
    // By default, it returns 0.
    virtual double lowerBound(){return((double) 0);};
    // Returns 0 if the energy is not suitable for current optimization algorithm 
    // Returns 1 if the energy is suitable for current optimization algorithm
    // Returns 2 if current optimizaiton algorithm does not check the energy 
    virtual char checkEnergy();

    typedef enum
        {
            FUNCTION,
            ARRAY,
            THREE_PARAM,
            NONE
        } InputType;


protected:
    int  m_width, m_height;  // width and height of a grid,if graph is a grid
    int  m_nPixels;          // number of pixels, for both grid and non-grid graphs
    int  m_nLabels;          // number of labels, for both grid and non-grid graphs
    bool m_grid_graph;   // true if the graph is a 2D grid
    bool m_varWeights;   // true if weights are spatially varying. To be used only with 2D grids
    bool m_initialized;  // true if array m_V is allocated memory.  
    EnergyFunction *m_e;
    
    InputType m_dataType;     
    InputType m_smoothType;

    // *********** SET THE DATA COSTS 
    // Following 2 functions set the data costs
    virtual void setData(DataCostFn dcost)=0; 
    virtual void setData(CostVal* data)=0;   

    // *********** SET THE SMOOTHNESS COSTS 
    // following 3 functions set the smoothness costs 
    // there are 2 ways to represent the smoothness costs, one with array, one with function
    // In addition, for 2D grid graphs spacially varying weights can be specified by 2 arrays

    // Smoothness cost depends on labels  V(l1,l2) for all edges (except for a multiplier - see setCues() ).
    // V must be symmetric: V(l1,l2) = V(l2,l1)
    // V must be an array of size nLabels*nLabels. It is NOT copied into internal memory
    virtual void setSmoothness(CostVal* V)=0;
    
    // General smoothness cost can be specified by passing pointer to a function 
    virtual  void setSmoothness(SmoothCostGeneralFn cost)=0;

    // To prevent implementations from caching all general smoothness costs values, the flag below
    // can be set to false by calling dontCacheSmoothnessCosts() before calling initialize():
    bool m_allocateArrayForSmoothnessCostFn;

    // Use if the smoothness is V(l1,l2) = lambda * min ( |l1-l2|^m_smoothExp, m_smoothMax )
    // Can also add optional spatially varying weights for 2D grid graphs using setCues()
    virtual void setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda)=0;

    // You are not required to call setCues, in which case there is no multiplier.
    // Function below cannot be called for general cost function.
    // This function can be only used for a 2D grid graph
    // hCue and vCue must be arrays of size width*height in row major order. 
    // They are NOT copied into internal memory.
    // hCue(x,y) holds the variable weight for edge between pixels (x+1,y) and (x,y)
    // vCue(x,y) holds the variable weight for edge between pixels (x,y+1) and (x,y)
    virtual void setCues(CostVal* hCue, CostVal* vCue)=0; 
    
    virtual void initializeAlg()=0; // called by initialize()

    void commonInitialization(EnergyFunction *e);
    void checkArray(CostVal *V);
};

// *********** This class is for data costs
// Data costs can be specified eithe by an array or by a pointer to a function 
// If specified by an array, use constructor DataCost(cost) where 
// cost is the array of type CostVal. The cost of pixel p and label l is
// stored at cost[p*nLabels+l] where nLabels is the number of labels
// If data costs are to be specified by a function, pass
// a pointer to a function 
// CostVal costFn(int pix, Label lab)
// which returns the
// data cost of pixel pix to be assigned label lab


class DataCost
{
friend class MRF;
public:
    typedef MRF::CostVal CostVal;
    typedef MRF::DataCostFn DataCostFn;
    DataCost(CostVal *cost){m_costArray = cost;m_type = MRF::ARRAY; };
    DataCost(DataCostFn costFn){m_costFn = costFn;m_type = MRF::FUNCTION;};
private:
    MRF::CostVal *m_costArray;
    MRF::DataCostFn m_costFn;
    MRF::InputType m_type;     
};

// ***************** This class represents smoothness costs 
// If the smoothness is V(l1,l2) = lambda * min ( |l1-l2|^m_smoothExp, m_smoothMax )
// use constructor SmoothnessCost(smoothExp,smoothMax,lambda)
// If, in addition,  there are spacially varying weights use constructor
// SmoothnessCost(smoothExp,smoothMax,lambda,hWeights,vWeights)
// hWeights and vWeights can be only used for a 2D grid graph
// hWeights and vWeights must be arrays of size width*height in row major order. 
// They are NOT copied into internal memory.
// hWeights(x,y) holds the variable weight for edge between pixels (x+1,y) and (x,y)
// vWeights(x,y) holds the variable weight for edge between pixels (x,y+1) and (x,y)
//  If the smoothness costs are specified by input array V of type CostVal and
// size nLabels*nLabels, use consructor SmoothnessCost(V). 
// If in addition, there are 
// are spacially varying weights use constructor SmoothnessCost(V,hWeights,vWeights)
// Note that array V must be of size nLabels*nLabels, and be symmetric. 
// That is V[i*nLabels+j] = V[j*nLabels+i]
// Finally, if the smoothness term is specified by a general function, use
// constructor SmoothnessCost(costFn)


class SmoothnessCost
{
    friend class MRF;
public:
    typedef MRF::CostVal CostVal;
    // Can be used for  2D grids and for general graphs
    // In case if used for 2D grids, the smoothness term WILL NOT be spacially varying
    SmoothnessCost(int smoothExp,CostVal smoothMax,CostVal lambda)
        {m_type=MRF::THREE_PARAM;m_smoothMax = smoothMax;m_smoothExp = smoothExp;m_lambda=lambda;m_varWeights=false;};

    // Can be used only for 2D grids 
    // the smoothness term WILL BE be spacially varying
    SmoothnessCost(int smoothExp,CostVal smoothMax,CostVal lambda,CostVal *hWeights, CostVal *vWeights)
        {m_type=MRF::THREE_PARAM;m_smoothMax = smoothMax;m_smoothExp = smoothExp;m_lambda=lambda;
        m_varWeights = true;m_hWeights = hWeights; m_vWeights = vWeights;};
    
    // Can be used 2D grids and for general graphs
    // In case if used for 2D grids, the smoothness term WILL NOT be spacially varying
    SmoothnessCost(CostVal *V){m_V = V;m_type = MRF::ARRAY;m_varWeights=false;};

    // Can be used only for 2D grids 
    // the smoothness term WILL BE be spacially varying
    SmoothnessCost(CostVal *V,CostVal *hWeights, CostVal *vWeights )
        {m_V = V;m_hWeights = hWeights; m_vWeights = vWeights; m_varWeights = true; m_type=MRF::ARRAY;};

    // Can be used 2D grids and for general graphs
    SmoothnessCost(MRF::SmoothCostGeneralFn costFn){m_costFn = costFn;m_type = MRF::FUNCTION;m_varWeights=false;};

private:
    CostVal *m_V,*m_hWeights, *m_vWeights;
    MRF::SmoothCostGeneralFn m_costFn;
    MRF::InputType m_type;
    int m_smoothExp;
    CostVal m_smoothMax,m_lambda;
    bool m_varWeights;
    EnergyFunction *m_eng;
};


class EnergyFunction
{
public:

    EnergyFunction(DataCost *dataCost,SmoothnessCost *smoothCost)
        {m_dataCost = dataCost;m_smoothCost = smoothCost;};
    DataCost *m_dataCost;
    SmoothnessCost *m_smoothCost;
};


#endif /*  __MRF_H__ */

/*
    virtual EnergyVal dataEnergy() = 0;        
    virtual EnergyVal smoothnessEnergy() = 0;  
    virtual void setNeighbors(int pix1, int pix2, CostVal weight)= 0;
    virtual void optimizeAlg(int nIterations)=0;
    virtual Label* getAnswerPtr()= 0;
    virtual Label getLabel(int pixel)= 0;
    virtual void setLabel(int pixel,Label label)= 0;
    virtual void clearAnswer()  = 0;
    virtual void setParameters(int numParam, void *param) = 0;
    virtual void setData(DataCostFn dcost)=0; 
    virtual void setData(CostVal* data)=0;   
    virtual void setSmoothness(CostVal* V)=0;
    virtual void setSmoothness(SmoothCostGeneralFn cost)=0;
    virtual void setCues(CostVal* hCue, CostVal* vCue)=0; 
    virtual void setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda);
    virtual EnergyVal lowerBound(){return((EnergyVal) 0);};
*/
