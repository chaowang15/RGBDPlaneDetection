#include "plane_detection.h"

PlaneDetection plane_detection;

//-----------------------------------------------------------------
// MRF energy functions
MRF::CostVal dCost(int pix, int label)
{
	return plane_detection.dCost(pix, label);
}

MRF::CostVal fnCost(int pix1, int pix2, int i, int j)
{
	return plane_detection.fnCost(pix1, pix2, i, j);
}

void runMRFOptimization()
{
	DataCost *data = new DataCost(dCost);
	SmoothnessCost *smooth = new SmoothnessCost(fnCost);
	EnergyFunction *energy = new EnergyFunction(data, smooth);
	int width = kDepthWidth, height = kDepthHeight;
	MRF* mrf = new Expansion(width * height, plane_detection.plane_num_ + 1, energy);
	// Set neighbors for the graph
	for (int row = 0; row < height; row++)
	{
		for (int col = 0; col < width; col++)
		{
			int pix = row * width + col;
			if (col < width - 1) // horizontal neighbor
				mrf->setNeighbors(pix, pix + 1, 1);
			if (row < height - 1) // vertical
				mrf->setNeighbors(pix, pix + width, 1);
			if (row < height - 1 && col < width - 1) // diagonal
				mrf->setNeighbors(pix, pix + width + 1, 1);
		}
	}
	mrf->initialize();
	mrf->clearAnswer();
	float t;
	mrf->optimize(5, t);  // run for 5 iterations, store time t it took 
	MRF::EnergyVal E_smooth = mrf->smoothnessEnergy();
	MRF::EnergyVal E_data = mrf->dataEnergy();
	cout << "Optimized Energy: smooth = " << E_smooth << ", data = " << E_data << endl;
	cout << "Time consumed in MRF: " << t << endl;

	// Get MRF result
	for (int row = 0; row < height; row++)
	{
		for (int col = 0; col < width; col++)
		{
			int pix = row * width + col;
			plane_detection.opt_seg_img_.at<cv::Vec3b>(row, col) = plane_detection.plane_colors_[mrf->getLabel(pix)];
			plane_detection.opt_membership_img_.at<int>(row, col) = mrf->getLabel(pix);
		}
	}
	delete mrf;
	delete energy;
	delete smooth;
	delete data;
}
//-----------------------------------------------------------------


void printUsage()
{
	cout << "Usage: RGBDPlaneDetection <-o> color_image depth_image" << endl;
	cout << "-o: run MRF-optimization based plane refinement" << endl;
}

int main(int argc, char** argv)
{
	if (argc != 3 && argc != 4)
	{
		printUsage();
		return -1;
	}
	bool run_mrf = string(argv[1]) == "-o" ? true : false;
	string color_filename = run_mrf ? string(argv[2]) : string(argv[1]);
	string depth_filename = run_mrf ? string(argv[3]) : string(argv[2]);
	
	plane_detection.readDepthImage(depth_filename);
	plane_detection.readColorImage(color_filename);
	plane_detection.runPlaneDetection();

	if (run_mrf)
	{
		plane_detection.prepareForMRF();
		runMRFOptimization();
	}

	// The filename of a depth image is like '/depth-image-path/frame-XXXXXX-depth.png'
	string output_prefix = depth_filename.substr(0, depth_filename.length() - 10);
	plane_detection.writeOutputFiles(output_prefix, run_mrf);
	return 0;
}