# CoastCams
CoastCams is an open-source collection of existing MATLAB scripts to quantify key wave parameters (e.g., wave height, wave period), mean water levels, and morphology (e.g., shoreline positions) in the nearshore environment. The repository performs the analysis on oblique orthorectified timestack images from land-based coastal monitoring systems. The proposed approach is a combination of several key parameters that aims to get a better understanding of nearshore processes by leveraging the strength of existing codes. CoastCams provides a unified and simplified method that is accessible for coastal managers, engineers, and scientists with a user-friendly and practical method to monitor and identify key drivers in coastal zone. In this paper, we present the standalone remote video-based method and validate the estimated hydro parameters with sensors deployed in the nearshore on a rocky platform in Socoa, France. 

## Usage
### Overview
CoastCams builds upon the foundation laid by [CIRN](https://github.com/Coastal-Imaging-Research-Network), while expanding on the capabilities by making the codes accessible to estimate nearshore processes, mean water levels, and morpholigcal changes in a unified and simplified manner that is accessible to a wide range of uers, i.e., from experts to novices. More information can be found [in this paper](https://www.sciencedirect.com/science/article/pii/S136481522300186X). 

The input of CoastCams are georectified timestack images from coastal video cameras. The creation of timestack images and georectification can be achieved with the [Quantitative Coastal Imaging Toolbox](https://github.com/Coastal-Imaging-Research-Network/CIRN-Quantitative-Coastal-Imaging-Toolbox). 

### User Inputs
The whole repository can be downloaded and added to your MATLAB path.
The only user inputs are required in [S01_AnalysisTimestackImages](https://github.com/NuytsSiegmund/CoastCams/edit/main/UserScripts/S01_AnalysisTimestackImages) and are dependent are the dimensions and acquisition of your timestack images.

Follow these steps:
1. Open S01_AnalysisTimestackImages
2. Add all paths to your MATLAB workspace in section B;
3. Select the timestack images in section C;
4. Add your specific parameters for image processing in section D:
   * ⋅⋅⋅D1⋅⋅
   * dt = Frequency acquisition of the camera e.g., freq = 2 (2 images per second);
   * H_camera = Camera height above MSL im metre;
   * res = Size of each pixel on timestack image in metre;
   * rotation = Waves in the timestack image should come from top-left corner - rotate the timestack image accordingly




## Contributing and Issues
Having a problem? Post an issue in the [Issues Page](https://github.com/NuytsSiegmund/CoastCams/issues)

If you're willing to contribute: 

1. Fork the repository. A fork is a copy on which you can make your changes.
2. Create a new branch on your fork
3. Commit your changes and push them to your branch
4. When the branch is ready to be merged, create a Pull Request (how to make a clean pull request explained [here](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request))

## Publication

More information can be found in the following publication: 

Nuyts, S., Almar, R., Morichon, D., Dealbera, S., Abalia, A., Muñoz, J. M., Abessolo, G. O., & Regard, V. (2023). CoastCams: A MATLAB toolbox making accessible estimations of nearshore processes, mean water levels, and morphology from timestack images. Environmental Modelling & Software, 168, 105800. https://doi.org/https://doi.org/10.1016/j.envsoft.2023.105800 

