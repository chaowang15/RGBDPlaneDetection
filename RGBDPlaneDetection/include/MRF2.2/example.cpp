// example.cpp -- illustrates calling the MRF code

static const char *usage = "usage: %s [energyType] (a number between 0 and 3)\n";

// uncomment "#define COUNT_TRUNCATIONS" in energy.h to enable counting of truncations

#include "mrf.h"
#include "ICM.h"
#include "GCoptimization.h"
#include "MaxProdBP.h"
#include "TRW-S.h"
#include "BP-S.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <new>

const int sizeX = 50;
const int sizeY = 50;
const int numLabels = 20;

MRF::CostVal D[sizeX * sizeY * numLabels];
MRF::CostVal V[numLabels * numLabels];
MRF::CostVal hCue[sizeX * sizeY];
MRF::CostVal vCue[sizeX * sizeY];

#ifdef COUNT_TRUNCATIONS
int truncCnt, totalCnt;
#endif

EnergyFunction* generate_DataARRAY_SmoothFIXED_FUNCTION()
{
	int i, j;

	// generate function
	for (i = 0; i < numLabels; i++) {
		for (j = i; j < numLabels; j++) {
			V[i * numLabels + j] = V[j * numLabels + i] = (i == j) ? 0 : (MRF::CostVal)2.3;
		}
	}
	MRF::CostVal* ptr;
	for (ptr = &D[0]; ptr < &D[sizeX * sizeY * numLabels]; ptr++) *ptr = ((MRF::CostVal)(rand() % 100)) / 10 + 1;
	for (ptr = &hCue[0]; ptr < &hCue[sizeX * sizeY]; ptr++) *ptr = rand() % 3 + 1;
	for (ptr = &vCue[0]; ptr < &vCue[sizeX * sizeY]; ptr++) *ptr = rand() % 3 + 1;

	// allocate energy
	DataCost *data         = new DataCost(D);
	SmoothnessCost *smooth = new SmoothnessCost(V, hCue, vCue);
	EnergyFunction *energy    = new EnergyFunction(data, smooth);

	return energy;
}

EnergyFunction* generate_DataARRAY_SmoothTRUNCATED_LINEAR()
{
	// generate function
	MRF::CostVal* ptr;
	for (ptr = &D[0]; ptr < &D[sizeX * sizeY * numLabels]; ptr++) *ptr = ((MRF::CostVal)(rand() % 100)) / 10 + 1;
	for (ptr = &hCue[0]; ptr < &hCue[sizeX * sizeY]; ptr++) *ptr = rand() % 3;
	for (ptr = &vCue[0]; ptr < &vCue[sizeX * sizeY]; ptr++) *ptr = rand() % 3;
	MRF::CostVal smoothMax = (MRF::CostVal)25.5, lambda = (MRF::CostVal)2.7;

	// allocate energy
	DataCost *data         = new DataCost(D);
	SmoothnessCost *smooth = new SmoothnessCost(1, smoothMax, lambda, hCue, vCue);
	EnergyFunction *energy    = new EnergyFunction(data, smooth);

	return energy;
}

EnergyFunction* generate_DataARRAY_SmoothTRUNCATED_QUADRATIC()
{

	// generate function
	MRF::CostVal* ptr;
	for (ptr = &D[0]; ptr < &D[sizeX * sizeY * numLabels]; ptr++) *ptr = ((MRF::CostVal)(rand() % 100)) / 10 + 1;
	for (ptr = &hCue[0]; ptr < &hCue[sizeX * sizeY]; ptr++) *ptr = rand() % 3;
	for (ptr = &vCue[0]; ptr < &vCue[sizeX * sizeY]; ptr++) *ptr = rand() % 3;
	MRF::CostVal smoothMax = (MRF::CostVal)5.5, lambda = (MRF::CostVal)2.7;

	// allocate energy
	DataCost *data         = new DataCost(D);
	SmoothnessCost *smooth = new SmoothnessCost(2, smoothMax, lambda, hCue, vCue);
	EnergyFunction *energy    = new EnergyFunction(data, smooth);

	return energy;
}


MRF::CostVal dCost(int pix, int i)
{
	return ((pix * i + i + pix) % 30) / ((MRF::CostVal) 3);
}

MRF::CostVal fnCost(int pix1, int pix2, int i, int j)
{
	if (pix2 < pix1) { // ensure that fnCost(pix1, pix2, i, j) == fnCost(pix2, pix1, j, i)
		int tmp;
		tmp = pix1; pix1 = pix2; pix2 = tmp;
		tmp = i; i = j; j = tmp;
	}
	MRF::CostVal answer = (pix1 * (i + 1) * (j + 2) + pix2 * i * j * pix1 - 2 * i * j * pix1) % 100;
	return answer / 10;
}


EnergyFunction* generate_DataFUNCTION_SmoothGENERAL_FUNCTION()
{
	DataCost *data         = new DataCost(dCost);
	SmoothnessCost *smooth = new SmoothnessCost(fnCost);
	EnergyFunction *energy = new EnergyFunction(data, smooth);

	return energy;
}

int main(int argc, char **argv)
{
	MRF* mrf;
	EnergyFunction *energy;
	MRF::EnergyVal E;
	double lowerBound;
	float t, tot_t;
	int iter;

	int seed = 1124285485;
	srand(seed);

	int Etype = 0;

	if (argc != 2) {
		fprintf(stderr, usage, argv[0]);
		exit(1);
	}

	if (argc > 1)
		Etype = atoi(argv[1]);

	try {
		switch (Etype) {
		// Here are 4 sample energies to play with.
		case 0:
			energy = generate_DataARRAY_SmoothFIXED_FUNCTION();
			fprintf(stderr, "using fixed (array) smoothness cost\n");
			break;
		case 1:
			energy = generate_DataARRAY_SmoothTRUNCATED_LINEAR();
			fprintf(stderr, "using truncated linear smoothness cost\n");
			break;
		case 2:
			energy = generate_DataARRAY_SmoothTRUNCATED_QUADRATIC();
			fprintf(stderr, "using truncated quadratic smoothness cost\n");
			break;
		case 3:
			energy = generate_DataFUNCTION_SmoothGENERAL_FUNCTION();
			fprintf(stderr, "using general smoothness functions\n");
			break;
		default:
			fprintf(stderr, usage, argv[0]);
			exit(1);
		}

		bool runICM       = true;
		bool runExpansion = true;
		bool runSwap      = true;
		bool runMaxProdBP = true;
		bool runTRWS      = true;
		bool runBPS       = true;

		////////////////////////////////////////////////
		//                     ICM                    //
		////////////////////////////////////////////////
		if (runICM) {
			printf("\n*******Started ICM *****\n");

			mrf = new ICM(sizeX, sizeY, numLabels, energy);
			mrf->initialize();
			mrf->clearAnswer();

			E = mrf->totalEnergy();
			printf("Energy at the Start= %g (%g,%g)\n", (float)E,
			       (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

			tot_t = 0;
			for (iter = 0; iter < 6; iter++) {
				mrf->optimize(10, t);

				E = mrf->totalEnergy();
				tot_t = tot_t + t ;
				printf("energy = %g (%f secs)\n", (float)E, tot_t);
			}

			delete mrf;
		}

		////////////////////////////////////////////////
		//          Graph-cuts expansion              //
		////////////////////////////////////////////////
		if (runExpansion) {
			printf("\n*******Started graph-cuts expansion *****\n");
			mrf = new Expansion(sizeX, sizeY, numLabels, energy);
			mrf->initialize();
			mrf->clearAnswer();

			E = mrf->totalEnergy();
			printf("Energy at the Start= %g (%g,%g)\n", (float)E,
			       (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

#ifdef COUNT_TRUNCATIONS
			truncCnt = totalCnt = 0;
#endif
			tot_t = 0;
			for (iter = 0; iter < 6; iter++) {
				mrf->optimize(1, t);

				E = mrf->totalEnergy();
				tot_t = tot_t + t ;
				printf("energy = %g (%f secs)\n", (float)E, tot_t);
			}
#ifdef COUNT_TRUNCATIONS
			if (truncCnt > 0)
				printf("***WARNING: %d terms (%.2f%%) were truncated to ensure regularity\n",
				       truncCnt, (float)(100.0 * truncCnt / totalCnt));
#endif

			delete mrf;
		}

		////////////////////////////////////////////////
		//          Graph-cuts swap                   //
		////////////////////////////////////////////////
		if (runSwap) {
			printf("\n*******Started graph-cuts swap *****\n");
			mrf = new Swap(sizeX, sizeY, numLabels, energy);
			mrf->initialize();
			mrf->clearAnswer();

			E = mrf->totalEnergy();
			printf("Energy at the Start= %g (%g,%g)\n", (float)E,
			       (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

#ifdef COUNT_TRUNCATIONS
			truncCnt = totalCnt = 0;
#endif
			tot_t = 0;
			for (iter = 0; iter < 8; iter++) {
				mrf->optimize(1, t);

				E = mrf->totalEnergy();
				tot_t = tot_t + t ;
				printf("energy = %g (%f secs)\n", (float)E, tot_t);
			}
#ifdef COUNT_TRUNCATIONS
			if (truncCnt > 0)
				printf("***WARNING: %d terms (%.2f%%) were truncated to ensure regularity\n",
				       truncCnt, (float)(100.0 * truncCnt / totalCnt));
#endif


			delete mrf;
		}

		////////////////////////////////////////////////
		//          Belief Propagation                //
		////////////////////////////////////////////////
		if (runMaxProdBP) {
			printf("\n*******  Started MaxProd Belief Propagation *****\n");
			mrf = new MaxProdBP(sizeX, sizeY, numLabels, energy);
			mrf->initialize();
			mrf->clearAnswer();

			E = mrf->totalEnergy();
			printf("Energy at the Start= %g (%g,%g)\n", (float)E,
			       (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

			tot_t = 0;
			for (iter = 0; iter < 10; iter++) {
				mrf->optimize(1, t);

				E = mrf->totalEnergy();
				tot_t = tot_t + t ;
				printf("energy = %g (%f secs)\n", (float)E, tot_t);
			}


			delete mrf;
		}

		////////////////////////////////////////////////
		//                  TRW-S                     //
		////////////////////////////////////////////////
		if (runTRWS) {
			printf("\n*******Started TRW-S *****\n");
			mrf = new TRWS(sizeX, sizeY, numLabels, energy);

			// can disable caching of values of general smoothness function:
			//mrf->dontCacheSmoothnessCosts();

			mrf->initialize();
			mrf->clearAnswer();


			E = mrf->totalEnergy();
			printf("Energy at the Start= %g (%g,%g)\n", (float)E,
			       (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

			tot_t = 0;
			for (iter = 0; iter < 10; iter++) {
				mrf->optimize(10, t);

				E = mrf->totalEnergy();
				lowerBound = mrf->lowerBound();
				tot_t = tot_t + t ;
				printf("energy = %g, lower bound = %f (%f secs)\n", (float)E, lowerBound, tot_t);
			}

			delete mrf;
		}

		////////////////////////////////////////////////
		//                  BP-S                     //
		////////////////////////////////////////////////
		if (runBPS) {
			printf("\n*******Started BP-S *****\n");
			mrf = new BPS(sizeX, sizeY, numLabels, energy);

			// can disable caching of values of general smoothness function:
			//mrf->dontCacheSmoothnessCosts();

			mrf->initialize();
			mrf->clearAnswer();

			E = mrf->totalEnergy();
			printf("Energy at the Start= %g (%g,%g)\n", (float)E,
			       (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

			tot_t = 0;
			for (iter = 0; iter < 10; iter++) {
				mrf->optimize(10, t);

				E = mrf->totalEnergy();
				tot_t = tot_t + t ;
				printf("energy = %g (%f secs)\n", (float)E, tot_t);
			}

			delete mrf;
		}
	}
	catch (std::bad_alloc) {
		fprintf(stderr, "*** Error: not enough memory\n");
		exit(1);
	}

	return 0;
}
