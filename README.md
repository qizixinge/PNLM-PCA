# [New non-local mean methods for MRI denoising based on global self-similarity between values(Modifying in progress)](https://arxiv.org/abs/2308.14145)
These are the code files for this article. This study improves the pre-filtered rotationally invariant non-local PCA method, ensuring its efficiency while improving algorithm performance. This work has two main contributions. First, the proposed new algorithm PNLM-PCA exhibits superior denoising performance in different types of MRI images over some state-of-the-art algorithms. Second, a novel structure of NLM is proposed, which simultaneously improves speed and accuracy. The proposed new structure is conducive to promoting the development of NLM based denoising algorithms for natural and MR images.

DOI:&#x20;

## The description of files in this repository

The contents contained in the MATLAB files are the following:&#x20;

- 1.`t2_ai_msles2_1mm_pn0_rf0.rawb`: An example data file of T2w MS lesion image from the BrainWeb dataset.

- 2.`pd_icbm_normal_1mm_pn0_rf0.rawb`: An example data file of PDw normal image from the BrainWeb dataset.

- 3.`precomputation.mat`: A MATLAB data file which provides pre-computed relevant data for images' Rician correction.

- 4.`optimal_parameter.m`:A MATLAB script file which demonstrates the process of finding the optimal parameters of NLPCA.

- 5.`Experiments.m`: A MATLAB script file which completely presents the code used in the whole process of the experiments.

- 6.`ssim_index3d.m`: A MATLAB function file that calculates the SSIM of a denoised image.

- 7.`RicianSTD.m`: A MATLAB function file which is used to estimate the Rician noise in the MR image using the object-based method.

- 8.`ricernd.m`: A MATLAB function file that adds Rician noise on a noise-free image.

- 9.`RI_NLM.m`: A MATLAB function file that displays a possible version of rotationally invariant non-local mean filter.

- 10.`psnrallM.m`: A MATLAB function which is used to find the optimal value of \$\tau\beta\$ and T in the NL-PCA filter.

- 11.`psnrall.m`: A MATLAB function which is used to find the optimal value of (d, M, w, \$\tau\beta\$, T) in the NL-PCA filter.

- 12.`bm4dw.m`: A MATLAB function which serves as the powerful auxiliary tool mentioned in the paper.

- 13.`bwp.m`: A MATLAB function that shows the new collaborative algorithm proposed by us.

- 14.`NLPCApso.m`: A MATLAB function which is the adaptation of NL-PCA filter for the PSO algorithm.

- 15.`NLPCAnp.m`: A MATLAB function that demonstrates the code for realizing the NL-PCA filter with exact noise level.

- 16.`NLPCA.m`: A MATLAB function file that demonstrates the code for realizing the NL-PCA filter.

- 17.`kmeans.m`: for image segmentation.

- 18.`farras.m`: nearly symmetric filters for orthogonal 2-channel perfect reconstruction filter bank.

- 19.`epsi.m`: A MATLAB function file which calculate the correction factor \$\xi\$

- 20.`dwt3D.m`: 3-D Discrete Wavelet Transform.

- 21.`cshift3D.m`: 3D Circular Shift.

- 22.`afb3D.m`: 3D Analysis Filter Bank.

- 23.`OAS1_0001_MR1_mpr-1_anon.img`: An example data file of T1w image from the OASIS dataset.

- 24.`OAS1_0001_MRI_mpr-1_anon.hdr`: The above file cannot be read without it.

- 25.`T1w_acpc_dc_restore.nii.gz`: An example data file of T1w image from the HCP dataset.

- 26.`RINLMmy.m`:

- 27.`RiC.m`:

- 28.`truncateslice.m`:

- 29.`cPRI_NL_PCA.mexw64`:

- 30.`dis255.m`:

- 31.`getTransfMatrix.m`:

- 32.`optimizationRINLMmy.m`:

- 33.`precomputation.m`:

- 34.`RiceOptVST`:

- 35.`BM4D`: 






