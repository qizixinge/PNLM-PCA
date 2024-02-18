# [The source code for the new paper is expected to be updated before February 21st](https://arxiv.org/abs/2308.14145)
These are the code files for this article. This study improves the pre-filtered rotationally invariant non-local PCA method, ensuring its efficiency while improving algorithm performance. This denoising algorithm has two obvious advantages: first, it has the noise estimator built-in, and can even handle spatially varying Rician noise situations well; the second reason is that it has a weak dependence on parameters, and this algorithm has a certain degree of robustness in denoising performance for different types of images under the same set of parameters. We can extract an auxiliary tool from the new synthesis algorithm obtained in this study, which further reduces its dependence on parameters. It is highly promising to further improve their performance by combining it with state-of-the-art methods.

DOI:&#x20;

## 1.The description of files in this repository

The contents contained in the MATLAB files are the following:&#x20;

- 1.`t2_ai_msles2_1mm_pn0_rf0.rawb`: An example data file of T2w MS lesion image from the BrainWeb dataset

- 2.`t1_icbm_normal_1mm_pn0_rf0.rawb`: An example data file of T1w normal image from the BrainWeb dataset

- 3.`precomputation.mat`: A MATLAB data file which provides pre-computed relevant data for images' Rician correction

- 4.`optimal_parameter.m`:A MATLAB script file which demonstrates the process of finding the optimal parameters

- 5.`Experiments.m`: A MATLAB script file which completely presents the code used in the whole process of the experiments

- 6.`ssim_index3d.m`: A MATLAB function file which calculate the SSIM of a denoised image

- 7.`RicianSTD.m`: A MATLAB function file which is used to estimate the Rician noise in the MR image using the object-based method

- 8.`ricernd.m`: A MATLAB function file which add Rician noise on a noise-free image

- 9.`RI_NLM.m`: A MATLAB function file which display the principle of rotationally invariant non-local mean filter

- 10.`psnrallM`: A MATLAB function which is used to find the optimal value of \$\$ and T in the NL-PCA filter&#x20;

- 11.`psnrall`: A MATLAB function which is used to find the optimal value of (d, M, w, \$\$, T) in the NL-PCA filter&#x20;

- 12.`PD`: A MATLAB function which serve as the powerful auxiliary tool mentioned in the paper

- 13.`PCA_PRI_PCAr`: A MATLAB function which shows the new collaborative algorithm proposed by us

- 14.`NLPCApso`: A MATLAB function which is the adaptation of NL-PCA filter for the PSO algorithm

- 15.`NLPCAmy`: A MATLAB function which tries to improve the efficiency of the NL-PCA filter under the optimal parameter set obtained by the PSO algorithm

- 16.`NLPCA.m`: A MATLAB function file which demonstrates the code for realizing the NL-PCA filter

- 17.`kmeans`: for image segmentation

- 18.`farras.m`: nearly symmetric filters for orthogonal 2-channel perfect reconstruction filter bank

- 19.`epsi.m`: A MATLAB function file which calculate the correction factor \$\$

- 20.`dwt3D`: 3-D Discrete Wavelet Transform

- 21.`cshift3D`: 3D Circular Shift

- 22.`afb3D`: 3D Analysis Filter Bank

- 23.`OAS1_0001_MR1_mpr-1_anon.img`: An example data file of T1w image from the OASIS dataset

- 24.`OAS1_0001_MRI_mpr-1_anon.hdr`: The above file cannot be read without it.

- 25.`T1w_acpc_dc_restore.nii.gz`: An example data file of T1w image from the HCP dataset



