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

## Note
- Camera intrinsic parameters are set in `plane_detection.h`
- 
