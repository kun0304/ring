# Computed Tomography Ring Artifact Correction Method with Super-Pixel and Adaptive Relative Total Variation

This repository contains the MATLAB implementation of the research paper titled "Computed Tomography Ring Artifact Correction Method with Super-Pixel and Adaptive Relative Total Variation" by Na Li and Xiaokun Liang (E-mail: xk.liang@siat.ac.cn).

## Overview

The method presented in this repository addresses the significant challenge of removing ring artifacts in computed tomography (CT) images. It specifically targets the complexities introduced by photon counting CT (PCCT) technology, using a combination of super-pixel segmentation and an adaptive form of Relative Total Variation (RTV). This dual-stage approach not only ensures the efficient correction of artifacts but also maintains the integrity and balance of the image details, which is crucial in medical imaging applications.

## Getting Started

### Prerequisites

- MATLAB (Tested on MATLAB R2020b and later)
- Image Processing Toolbox

### Installation

Clone the repository:

```
git clone https://github.com/kun0304/ring.git
```

### Usage

1. Load your CT image in `.mat` format.
2. Run `main_run.m` in MATLAB, which utilizes functions in `ring_remove.m`.
3. Processed images will be displayed, showing the corrected artifacts.

## Files in the Repository

- `main_run.m`: Main MATLAB script to run the correction process.
- `ring_remove.m`: Contains the functions for the artifact correction.
- `ring_img_01.mat`: Sample image for testing the algorithm.

## Contributing

Contributions to improve the algorithm are welcome. Please fork the repository, make your changes, and submit a pull request.

## License

This project is released under the MIT License - see [LICENSE.md](LICENSE.md) for details.

## Acknowledgments

This work has been significantly aided by the code and methods developed in the following research articles:

- Xu L, Yan Q, Xia Y, et al. "Structure extraction from texture via relative total variation." ACM Transactions on Graphics (TOG), 2012, 31(6): 1-10.
- Garcia D. "Robust smoothing of gridded data in one and higher dimensions with missing values." Computational Statistics & Data Analysis, 2010, 54(4): 1167-1178.

We extend our gratitude to the authors of these papers for their invaluable contributions to the field and for the insights their work has provided, which were instrumental in the development of our method.
