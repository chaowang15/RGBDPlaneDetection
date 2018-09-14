# RGBDPlaneDetection
This code is an implementation of RGB-D plane detection and color-based plane refinement with MRF optimization. 

## Reference
**General pipeline**: 

[1] Huang, Jingwei, et al. *3DLite: Towards Commodity 3D Scanning for Content Creation*. ACM Transactions on Graphics, 2017.

[2] Wang, Chao, et al. "Plane-Based Optimization of Geometry and Texture for RGB-D Reconstruction of Indoor Scenes". International conference on 3D Vision (3DV) 2018, pp. 533-541.

**Plane detection**: 

Feng, Chen, Yuichi Taguchi, and Vineet R. Kamat. *Fast plane extraction in organized point clouds using agglomerative hierarchical clustering*. Robotics and Automation (ICRA), 2014 IEEE International Conference on. IEEE, 2014.

**MRF**:

See the readme file in [MRF 2.2](http://vision.middlebury.edu/MRF/code/) link for details about reference papers.

## Dependencies
- OpenCV
- Eigen 3
- [MRF 2.2](http://vision.middlebury.edu/MRF/code/) (already included)
- Fast single-image plane detection code: [PEAC](http://www-personal.umich.edu/~cforrest/research.html) (already included)

## Usage
```
RGBDPlaneDetection <-o> color_image depth_image
```
- `-o` is running MRF optimization to refine planes.

Two scripts `demo_win.sh` and `demo_linux.sh` are provided to run the code on a RGB-D sequence. Note to modify the corresponding paths.

## Build
In Windows, use Visual Studio to open sln file and compile and build the code. It tests successfully in Visual Studio 2010 on Windows 7 and 8.1. 

In Linux, run `build_linux.sh` to build the code.


## Running time
Without MRF optimization, the execution code by Visual Studio 2010 runs at about 20 FPS (including data I/O) on RGBD images with resolution 640x480 in the computer with 12GB RAM and intel i7 processor. With MRF optimization, the same code runs much slower at about 7 seconds per frame on the same data.  

## Output
- Plane segmentation image in PNG
- Plane label file in TXT: label of the plane which each pixel belongs to (starting from 0 to #planes - 1). If a pixel is not on any plane, then its label value is #planes.
- Plane data file in TXT. Each line represents one plane with format like this:
```
#plane_index(starting from 0 to #planes - 1) number_of_points_on_the_plane plane_color_in_png_image(r,g,b between [0,255]) plane_normal(1x3) plane_center(1x3) sx sy sz sxx syy szz sxy syz sxz
```
Here `(sx sy sz)` are average of sum of all 3D points `(x, y, z)` on the plane, `(sxx syy szz sxy syz sxz)` are the average of sum of `x*x, y*y, z*z, x*y, y*z, z*z` of all 3D points on the plane, respectively.

## Note
- Currently the code only works on [BundleFusion](http://graphics.stanford.edu/projects/bundlefusion/) or [3DLite](http://graphics.stanford.edu/projects/3dlite/) RGBD data. If you want to use other kinds of RGBD data, you need to rewrite the part of reading color and depth images, and reset the camera intrinsic parameters in `plane_detection.h`.
- Note for the scale factor for depth images in `plane_detection.h`.
- Sometimes the MRF 2.2 source code crashes in Visual Studio in Windows due to some kind of memory leak problem, but the code works well in Linux. If you meet the problem, just have a try in Linux, or try to implement the graph-cut/max-flow code by yourself. One suggestion is to use *boykov_kolmogorov_max_flow* in Boost library. 
