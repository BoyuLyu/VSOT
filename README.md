# VSOT

### Table of contents
- [Introduction](#Introduction)
- [System requirements](#System-requirements)
- [How to use](#How-to-use)
- [Reference](#Reference)
- [License](#License)
- [Acknowledgments](#Acknowledgments)

## Introduction
Morphological analysis of dendrites is critical to understanding the function and plasticity of neural circuits. The nanoscale details of recent public Electron Microscopy (EM) datasets with segmentations provide unprecedented opportunities to investigate subcellular structures and generate hypotheses regarding the underlying biological mechanisms. However, the methods for segmenting and quantifying the EM subcellular structures lag far behind the impressive data generation capability. We introduce VSOT, a novel method based on Volume-Surface Optimization, designed for the accurate structural analysis of dendrites in volumetric EM data. VSOT accurately segments dendritic compartments and quantifies subcellular structures by leveraging advanced optimization techniques that combine local surface and global volumetric information. Our computational analysis pipeline detects dendritic spines and further subdivides them into spine heads and necks. Through applications testing on public datasets of dendritic spines, as well as on a first-of-its-kind dataset of dendritic spine heads that we constructed, we demonstrated that our method generates more accurate segmentations quantitively compared to state-of-the-art approaches. By applying our methods to a large-scale EM dataset containing multiple layers of the brain, we effectively reveal how the nanostructure of dendrites varies in different brain areas. Furthermore, we explored the structural relationships between neurons and astrocytes in tripartite synaptic structures. With our newly developed computation methods, neuroscientists can exploit the newly generated large-scale volumetric EM data to address various scientific questions and advance the understanding of neural circuits.


## System requirements
- The code has been tested on Ubuntu 20.04 & 22.04
- MATLAB 2023b

## Run on code ocean
[code ocean compute capsule](https://codeocean.com/capsule/3574450/tree)
Free to run VSOT, with our example data as well as your own data.
## How to use locally
### data
- Test data can be downloaded from [here](https://virginiatech-my.sharepoint.com/:u:/g/personal/boyu93_vt_edu/EWu28K-ZjStNognr4h5H26wBErK5nY1EtrKRrnStCjhWfA?e=Nr27RH)
- Multi-layer cortex data constructed from MICrONS can be requested from the author Boyu Lyu [boyu93@vt.edu](mailto:boyu93@vt.edu)
### dependencies
- gcc/g++
- build-essential
- [tinymesh](https://github.com/tatsy/tinymesh.git) please download and extract to ./resources
- [iso2mesh](https://github.com/fangq/iso2mesh.git) please download and extract to ./resources
### run
After finishing downloading all the dependencies, run main.m for
- (1) Run the whole pipeline on an example dataset
- (2) Reproduce the performance testing for level-1 & level-2 segmentation
- (3) Plot the results of structural quantification.


## License
This project is licensed under the GNU license - see the [LICENSE](LICENSE) file for details.
