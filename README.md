# RGBDPlaneDetection
RGBD plane detection and color-based plane refinement with MRF optimization.

## Reference
- Feng, Chen, Yuichi Taguchi, and Vineet R. Kamat. **Fast plane extraction in organized point clouds using agglomerative hierarchical clustering**. Robotics and Automation (ICRA), 2014 IEEE International Conference on. IEEE, 2014.
- Huang, Jingwei, et al. **3DLite: Towards Commodity 3D Scanning for Content Creation**. ACM Transactions on Graphics 2017 (2017).

## Dependencies
- OpenCV
- Eigen 3
- [MRF 2.2](http://vision.middlebury.edu/MRF/code/) (already included)
- [PEAC](http://www-personal.umich.edu/~cforrest/research.html) (already included)

## Usage
```
RGBDPlaneDetection <-o> color_image depth_image
```
- `-o` is running MRF optimization to refine planes.

## Output
- Plane segmentation image in PNG
- Plane label image in PNG: plane label each pixel belongs to (starting from 0 to `number_of_points_on_the_plane` - 1). If a pixel is not on any plane, then its label value is `number_of_points_on_the_plane`.
- Plane data file in TXT. Each line represents one plane with format like this:
```
#plane_index(starting from 0) number_of_points_on_the_plane plane_color_in_png_image(r,g,b between [0,255]) plane_normal(1x3) plane_center(1x3) sx sy sz sxx syy szz sxy syz sxz
```
Here `(sx sy sz)` are average of sum of all 3D points `(x, y, z)` on the plane, `(sxx syy szz sxy syz sxz)` are the average of sum of `x*x, y*y, z*z, x*y, y*z, z*z` of all 3D points on the plane, respectively.

## Note
- Currently the code only works on [BundleFusion](http://graphics.stanford.edu/projects/bundlefusion/) or [3DLite](http://graphics.stanford.edu/projects/3dlite/) RGBD data. If you want to use other kinds of RGBD data, you need to rewrite the part of reading color and depth images, and reset the camera intrinsic parameters in `plane_detection.h`.
- Note for the scale factor for depth images in `plane_detection.h`.
