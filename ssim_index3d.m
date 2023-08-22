function [ssim_mean, ssim_map] = ssim_index3d(image1, image2, sw, index)


%----------------------------------------------------------------------
%Permission to use, copy, or modify this software and its documentation
%for educational and research purposes only and without fee is hereby
%granted, provided that this copyright notice and the original authors'
%names appear on all copies and supporting documentation. This program
%shall not be used, rewritten, or adapted as the basis of a commercial
%software or hardware product without first obtaining permission of the
%authors. The authors make no representations about the suitability of
%this software for any purpose. It is provided "as is" without express
%or implied warranty.
%----------------------------------------------------------------------
%
%This is an implementation of the algorithm for calculating the
%Structural SIMilarity (SSIM) index between two images. Please refer
%to the following paper:
%
%Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
%quality assessment: From error measurement to structural similarity"
%IEEE Transactios on Image Processing, vol. 13, no. 4, Apr. 2004.
%
%Kindly report any suggestions or corrections to zhouwang@ieee.org
%========================================================================


if (nargin < 2 | nargin > 5)
   ssim_mean = -Inf;
   ssim_map = -Inf;
   return;
end

if (size(image1) ~= size(image2))
   ssim_mean = -Inf;
   ssim_map = -Inf;
   return;
end

s = size(image1);
% sw = [2 2 2];

if (nargin == 2)
   if ((s(1) < sw(1)) | (s(2) < sw(2)) | (s(3) < sw(3)))
	   ssim_mean = -Inf;
	   ssim_map = -Inf;
      return
   end
   sw = [2 2 2];
    index = find(image1 ~=0);
%    window = fspecial('gaussian', 11, 1.5);	%
   window = gkernel(1.5,sw);	%
   K(1) = 0.01;								      % default settings
   K(2) = 0.03;								      %
   L = 255;                                  %
end

if (nargin == 3)
   if ((s(1) < sw(1)) | (s(2) < sw(2)) | (s(3) < sw(3)))
	   ssim_mean = -Inf;
	   ssim_map = -Inf;
      return
   end
%    window = fspecial('gaussian', 11, 1.5);
     window = gkernel(1.5,sw);	%
     index = find(image1 ~=0);
   K(1) = 0.01;								      % default settings
   K(2) = 0.03;	  
   L = 255;
   if (length(K) == 2)
      if (K(1) < 0 | K(2) < 0)
		   ssim_mean = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   ssim_mean = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end

if (nargin == 4)
   window = gkernel(1.5,sw);	%
   K(1) = 0.01;								      % default settings
   K(2) = 0.03;	
   L = 255;
   if (length(K) == 2)
      if (K(1) < 0 | K(2) < 0)
		   ssim_mean = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   ssim_mean = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end


C1 = (K(1)*L)^2;
C2 = (K(2)*L)^2;
window = window/sum(window(:));
image1 = double(image1);
image2 = double(image2);

mu1   = convn( image1,window, 'same');
mu2   = convn( image2,window, 'same');
mu1gs = mu1.*mu1;
mu2gs = mu2.*mu2;
mu1_mu2 = mu1.*mu2;
sigma1gs = convn( image1.*image1,window, 'same') - mu1gs;
sigma2gs = convn( image2.*image2,window, 'same') - mu2gs;
sigma12 = convn( image1.*image2,window, 'same') - mu1_mu2;

if (C1 > 0 & C2 > 0)
   ssim_map = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))./((mu1gs + mu2gs + C1).*(sigma1gs + sigma2gs + C2));
else
   numerator1 = 2*mu1_mu2 + C1;
   numerator2 = 2*sigma12 + C2;
	denominator1 = mu1gs + mu2gs + C1;
   denominator2 = sigma1gs + sigma2gs + C2;
   ssim_map = ones(size(mu1));
   index = (denominator1.*denominator2 > 0);
   ssim_map(index) = (numerator1(index).*numerator2(index))./(denominator1(index).*denominator2(index));
   index = (denominator1 ~= 0) & (denominator2 == 0);
   ssim_map(index) = numerator1(index)./denominator1(index);
end

temp = zeros(size(image1));
temp(index) = 1;
iind = find(temp ==0);
ssim_map(iind)=1;
ssim_mean = mean(ssim_map(index));

return

function [gaussianKer]=gkernel(sigma,sw)

% Pierrick Coupe - pierrick.coupe@gmail.com                                                                         
% Brain Imaging Center, Montreal Neurological Institute.                     
% Mc Gill University                                                         
%                                                                            
% Copyright (C) 2008 Pierrick Coupe             

for x = 1:(2*sw(1)+1)
    for y=1:(2*sw(2)+1)
        for z=1:(2*sw(3)+1)
            centerr2 = (x-(sw(1)+1))^2 + (y-(sw(2)+1))^2 + (z-(sw(3)+1))^2;
            gaussianKer(x, y, z) = exp(-centerr2/(2*sigma^2));
        end
    end
end