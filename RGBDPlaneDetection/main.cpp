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
	cout << "Usage: RGBDPlaneDetection <-o> color_image depth_image output_folder" << endl;
	cout << "-o: run MRF-optimization based plane refinement" << endl;
}

int main(int argc, char** argv)
{
	if (argc != 4 && argc != 5)
	{
		printUsage();
		return -1;
	}
	bool run_mrf = string(argv[1]) == "-o" ? true : false;
	string color_filename = run_mrf ? string(argv[2]) : string(argv[1]);
	string depth_filename = run_mrf ? string(argv[3]) : string(argv[2]);
	string output_folder = run_mrf ? string(argv[4]) : string(argv[3]);
	
	plane_detection.readDepthImage(depth_filename);
	plane_detection.readColorImage(color_filename);
	plane_detection.runPlaneDetection();

	if (run_mrf)
	{
		plane_detection.prepareForMRF();
		runMRFOptimization();
	}
	int pos = color_filename.find_last_of("/\\");
	string frame_name = color_filename.substr(pos + 1);
	frame_name = frame_name.substr(0, frame_name.length() - 10);
	plane_detection.writeOutputFiles(output_folder, frame_name, run_mrf);
	return 0;
}