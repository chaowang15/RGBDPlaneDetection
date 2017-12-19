**************************************************
* MRF energy minimization software               *
* Version 2.2                                    *
* August 30, 2012                                *
**************************************************

* Changes since version 2.0: 
* added code to count GC truncations (see energy.h)

* Changes since version 2.1: 
* minor changes to reduce warnings and allow compiling on 64-bit machines with gcc 4.6
* In the Makefile, 64-bit pointers are now selected by default


This directory contains the MRF energy minimization software accompanying
the paper

  [1] A Comparative Study of Energy Minimization Methods for Markov
      Random Fields.
      R. Szeliski, R. Zabih, D. Scharstein, O. Veksler, V. Kolmogorov,
      A. Agarwala, M. Tappen, and C. Rother.
      In Ninth European Conference on Computer Vision (ECCV 2006),
      volume 2, pages 16-29, Graz, Austria, May 2006.

Optimization methods included:
1) ICM
2) Graph Cuts
3) Max-Product Belief Propagation
4) Sequential tree-reweighted message passing (TRW-S)

Instructions for compiling and using the software are included below.

CREDITS:
* MRF code, graph cut interface, and example code by Olga Veksler
* Graph cut library by Yuri Boykov and Vladimir Kolmogorov
* Belief propagation code by Marshall Tappen
* TRW-S by Vladimir Kolmogorov


=========================== INSTRUCTIONS FOR CITATIONS ========================

If you use this software, you should cite the paper [1].  In addition,
since this software builds on the algorithms and libraries developed by
many different people, you should also follow the instructions about proper
citing below.

(a) If you use the GraphCuts optimization (written by Olga Veksler, using
    the libraries provided by Yuri Boykov and Vladimir Kolmogorov), you
    should cite the following papers:
      
    [2] Fast Approximate Energy Minimization via Graph Cuts.
        Y. Boykov, O. Veksler, and R. Zabih.
        In IEEE Transactions on Pattern Analysis and Machine Intelligence
        (PAMI), vol. 23, no. 11, pages 1222-1239, November 2001.  

    [3] What Energy Functions can be Minimized via Graph Cuts?
        V. Kolmogorov and R. Zabih. 
        In IEEE Transactions on Pattern Analysis and Machine Intelligence
        (PAMI), vol. 26, no. 2, pages 147-159, February 2004. 
        An earlier version appeared in European Conference on Computer
        Vision (ECCV), May 2002.

    [4] An Experimental Comparison of Min-Cut/Max-Flow Algorithms for
        Energy Minimization in Vision. 
        Y. Boykov and Vladimir Kolmogorov.
        In IEEE Transactions on Pattern Analysis and Machine Intelligence
        (PAMI), vol. 26, no. 9, pages 1124-1137, September 2004. 

    The graph-cuts software can be used only for research purposes. If you
    wish to use the graph-cuts software for commercial purposes, you should
    be aware that there is a US patent:

        R. Zabih, Y. Boykov, O. Veksler.
        "System and method for fast approximate energy minimization via
        graph cuts",  
        United Stated Patent 6,744,923, June 1, 2004.

    Citation [2] develops the graph-cuts algorithms (alpha-expansion and
    the swap algorithms), citation [3] gives a simpler graph construction
    for the methods described in [2], and citation [4] gives an efficient
    min-cut/max-flow algorithm for computing the minimum graph cut.
  
    The file energy.h included here is a slight modification of the file
    energy.h in energy-v1.1.src by Vladimir Kolmogorov, available at
    http://www.adastral.ucl.ac.uk/~vladkolm/software.html

    The files block.h, graph.cpp, graph.h, and maxflow.cpp included here
    are slight modifications of the files in maxflow-v2.2.src/forward_star,
    by Vladimir Kolmogorov and Yuri Boykov, available at
    http://www.adastral.ucl.ac.uk/~vladkolm/software.html

(b) If you are using the Belief Propagation software (provided by Marshall
    Tappen), you should cite 

    [5] M. F. Tappen and W. T. Freeman.
        Comparison of Graph Cuts with Belief Propagation for Stereo, using
        Identical MRF Parameters. 
        In Proceedings of the Ninth IEEE International Conference on
        Computer Vision (ICCV), pages 900-907, 2003

(c) If you are using the TRW-S software (provided by Vladimir Kolmogorov),
    you should cite 

    [6] M. J. Wainwright, T. S. Jaakkola and A. S. Willsky.
        MAP estimation via agreement on trees: Message-passing and
        linear-programming approaches. 
        In IEEE Transactions on Information Theory, vol. 51, no. 11, pages
        3697--3717, November 2005.

    [7] V. Kolmogorov.
        Convergent Tree-reweighted Message Passing for Energy Minimization.
        In IEEE Transactions on Pattern Analysis and Machine Intelligence
        (PAMI), vol. 28, no. 10, October 2006.  


============================= INSTRUCTIONS FOR USAGE ==========================

A simple Makefile is included - typing "make" will compile the library and
a sample driver file "example.cpp" into an executable "example".

Here is a brief description of how to use the code; see also example.cpp
and mrf.h.

The default type for energy is int.  If you want to change this type to
other numerical type, such as float, double, etc., in mrf.h, change the
following 2 lines to the desired type:

typedef int EnergyVal;        /* The total energy of a labeling */
typedef int CostVal;          /* costs of individual terms of the energy */

Note that CostVal is the type of individual energy terms (such as the data
cost for pixel p and label l), and EnergyVal is the type for the total
energy value. Therefore, the sum of all individual energy terms should not
"overflow" the EnergyVal type.

Step 1: Set up an energy function
   An energy function is specified by setting data costs and smoothness
   costs. Data costs and smoothness costs can be specified by an array or
   function. In addition, smoothness cost can be specified by parameters
   lambda, k, V_max, as explained in section 2 of paper [1], in which case
   the smoothness cost V(|l1-l2|) = lambda*min(V_max,|l1-l2|^k). Spatially
   varying per-pairing weights w_pq can be also specified.
   
   (a) Setting up Data Terms
       There are 2 ways to set up data terms: 
       (i)  you can allocate an array dCost of size (numberofPixels)*(numberofLabels) of the CostVal type,
            where dCost[pix*numberofLabels + l] is the cost of assigning label l to pixel pix
       (ii) or you can set up a function 
            MRF::CostVal dCost(int pix, int i), which takes pixel pix and label i and returns the data cost

       After you initialized data cost array, or setup a function data cost, to initialize data cost
       for optimization, call: 
              
       DataCost *data = new DataCost(dCost);

   (b) Setting up Smoothness Terms
       There are 3 ways to set up smoothness terms:
       (i)   Set up a symmetric array smooth of size numberofLabels*numberofLabels, 
             where V(l1,l2) = smooth[l1+numberofLabels*l2] = smooth[l2+numberofLabels*l1]
       (ii)  Specify V_max and k 
       (iii) The most general way is to set up afunction 
             MRF::CostVal fnCost(int pix1, int pix2, int i, int j) which takes pixels pix1 and pix2, 
             labels i and j and returns the smoothness penalty V for assigning labels i to pix1 and j to pix2. 
 
      In addition, if you used (i) or (ii) to specify the smoothness penalty V, and if the neighborhood system
      is the standard 4 connected grid, you can specify the spacial varying weights w_pq in arrays vCue, hCue,
      each of which has size width*height, where width is the width of the grid and height is the height of
      the set of the grid.  hCue[x+y*width] holds the variable weight for edge between pixels (x+1,y) and (x,y)
      and vCue[x+y*width] holds the variable weight for edge between pixels (x,y+1) and (x,y)
    
      After you've done steps (i), or (ii), or (iii), to initialize smoothness cost
      for optimization, call: 
         In case of (i) and NO spacially varying weights w_pq:  SmoothnessCost *smooth = new SmoothnessCost(smooth);
         In case of (i) and spacially varying weights w_pq:     SmoothnessCost *smooth = new SmoothnessCost(smooth,hCue,vCue);
             Note that in this case, V(pix1,pix2,l1,l2) = smooth[l1+l2*numofLabels]*hCue[min(pix1,pix2)] if pix1 and pix2
             are horizontal neighbors, and  V(pix1,pix2,l1,l2) = smooth[l1+l2*numofLabels]*vCue[min(pix1,pix2)] 
             if pix1 and pix2 are vertical neighbors
 
         In case of (ii) and NO spacially varying weights w_pq: SmoothnessCost *smooth = new SmoothnessCost(k,V_max,lambda);
         In case of  (ii) and spacially varying weights w_pq:   SmoothnessCost *smooth = new SmoothnessCost(k,V_max,lambda,hCue,vCue);
             Note that in this case, V(pix1,pix2,l1,l2) = lambda*min(V_max,|l1-l2]^k)*hCue[min(pix1,pix2)] if pix1 and pix2
             are horizontal neighbors, and  V(pix1,pix2,l1,l2) = lambda*min(V_max,|l1-l2]^k)*vCue[min(pix1,pix2)] 
             if pix1 and pix2 are vertical neighbors

         In case of (iii): SmoothnessCost *smooth = new SmoothnessCost(fnCost); 

     After the data and smoothness terms are set-up, the energy can be specified by calling,
     EnergyFunction *eng    = new EnergyFunction(data,smooth);

      Note that we directly use all the arrays that you make, that is we do not copy them. So don't free any memory
      in the arrays that you allocate until you are done with optimization.      

Step 2: Invoke an optimization algorithm, assuming the grid has dimensions width by height, and number of labels is numberOfLabels, and
eng is set up as above:

For ICM:
    float t;
    MRF* mrf = new ICM(width,height,numberOfLabels,eng);
    mrf->initialize();  
    mrf->clearAnswer();
    mrf->optimize(5,t);  // run for 5 iterations, store time t it took 
    MRF::EnergyVal E_smooth = mrf->smoothnessEnergy();
    MRF::EnergyVal E_data   = mrf->dataEnergy();
    printf("Total Energy = %d (Smoothness energy %d, Data Energy %d)\n", E_smooth+E_data,E_smooth,E_data);
    for (int pix =0; pix < width*height; pix++ ) printf("Label of pixel %d is %d",pix, mrf->getLabel(pix));
    delete mrf;


For Graph-cuts with expansion algorithm:
    float t;
    MRF* mrf = new Expansion(width,height,numberOfLabels,eng);
    mrf->initialize();  
    mrf->clearAnswer();
    mrf->optimize(5,t);  // run for 5 iterations, store time t it took 
    MRF::EnergyVal E_smooth = mrf->smoothnessEnergy();
    MRF::EnergyVal E_data   = mrf->dataEnergy();
    printf("Total Energy = %d (Smoothness energy %d, Data Energy %d)\n", E_smooth+E_data,E_smooth,E_data);
    for (int pix =0; pix < width*height; pix++ ) printf("Label of pixel %d is %d",pix, mrf->getLabel(pix));
    delete mrf;

For Graph-cuts with swap algorithm:
    float t;
    MRF* mrf = new Swap(width,height,numberOfLabels,eng);
    mrf->initialize();  
    mrf->clearAnswer();
    mrf->optimize(5,t);  // run for 5 iterations, store time t it took 
    MRF::EnergyVal E_smooth = mrf->smoothnessEnergy();
    MRF::EnergyVal E_data   = mrf->dataEnergy();
    printf("Total Energy = %d (Smoothness energy %d, Data Energy %d)\n", E_smooth+E_data,E_smooth,E_data);
    for (int pix =0; pix < width*height; pix++ ) printf("Label of pixel %d is %d",pix, mrf->getLabel(pix));
    delete mrf;

For Max-product Belief Propagation:
    float t;
    MRF* mrf = new MaxProdBP(width,height,numberOfLabels,eng);
    mrf->initialize();  
    mrf->clearAnswer();
    mrf->optimize(5,t);  // run for 5 iterations, store time t it took 
    MRF::EnergyVal E_smooth = mrf->smoothnessEnergy();
    MRF::EnergyVal E_data   = mrf->dataEnergy();
    printf("Total Energy = %d (Smoothness energy %d, Data Energy %d)\n", E_smooth+E_data,E_smooth,E_data);
    for (int pix =0; pix < width*height; pix++ ) printf("Label of pixel %d is %d",pix, mrf->getLabel(pix));
    delete mrf;


The above code works for grid graphs. Currently, the graph cuts algorithms (swap and expansion)  are also 
implemented for general neighborhood systems, in which case 
after specifying energy you also have to specify the neighborhood system using function: 
    mrf->setNeighbors(int pix1, int pix2, CostVal weight); 
which makes pix1 and pix2 neighbors with spacially varying weight. Note that in this case, you can't use hCue and vCue, 
when setting up the smoothness terms, since they only work for 4-connected grid. Also, if you use general neighborhood 
system, you have to use a different constructor:

For Graph-cuts with expansion algorithm:
    float t;
    MRF* mrf = new Expansion(numberofPixels,numberOfLabels,eng);
    
    for (int i = 0; i <numberofPixels; i++ )     // this nested loop sets up a full neighborhood system
        for ( j = 0; j < numberofPixels; j++ )   // with weights w_pq= (p-q)^2        
            mrf->setNeighbors(i,j, (i-j)*(i-j));

    mrf->initialize();  
    mrf->clearAnswer();
    mrf->optimize(5,t);  // run for 5 iterations, store time t it took 
    MRF::EnergyVal E_smooth = mrf->smoothnessEnergy();
    MRF::EnergyVal E_data   = mrf->dataEnergy();
    printf("Total Energy = %d (Smoothness energy %d, Data Energy %d)\n", E_smooth+E_data,E_smooth,E_data);
    for (int pix =0; pix < width*height; pix++ ) printf("Label of pixel %d is %d",pix, mrf->getLabel(pix));
    delete mrf;
     

For Graph-cuts with swap algorithm:
    float t;
    MRF* mrf = new Swap(numberofPixels,numberOfLabels,eng);
    
    for (int i = 0; i <numberofPixels; i++ )     // this nested loop sets up a full neighborhood system
        for ( j = 0; j < numberof Pixels; j++ )  // with weights w_pq= (p-q)^2        
            mrf->setNeighbors(i,j, (i-j)*(i-j));

    mrf->initialize();  
    mrf->clearAnswer();
    mrf->optimize(5,t);  // run for 5 iterations, store time t it took 
    MRF::EnergyVal E_smooth = mrf->smoothnessEnergy();
    MRF::EnergyVal E_data   = mrf->dataEnergy();
    printf("Total Energy = %d (Smoothness energy %d, Data Energy %d)\n", E_smooth+E_data,E_smooth,E_data);
    for (int pix =0; pix < width*height; pix++ ) printf("Label of pixel %d is %d",pix, mrf->getLabel(pix));
    delete mrf;
