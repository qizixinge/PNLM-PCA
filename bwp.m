function [y_est] = bwp(y, z, sigma, is_rician, alpha, do_wiener)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Specify the std. dev. of the corrupting noise
if (exist('sigma','var') ~= 1),
    sigma = 0; %% default std of the AWGN (enables noise estimation)
end
if (exist('is_rician','var') ~= 1),
    is_rician = 0; %% default noise distribution
end
if (exist('alpha','var') ~= 1),
    alpha = 1.0; %% default sharpening parameter
end
if (exist('do_wiener','var') ~= 1),
    do_wiener = 1; %% 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Following are the parameters for the Normal Profile.
%%%%

%%%% Select transforms ('dct', 'dst', 'hadamard', or anything that is listed by 'help wfilters'):
transform_3D_HT_name     = 'bior1.5'; %% transform used for the HT filtering
transform_3D_Wiener_name = 'dct';     %% transform used for the Wiener filtering
transform_4D_name        = 'haar';    %% transform used in the grouping dim, 
                                      %  the same for HT and Wiener filt.

%%%% Hard-thresholding (HT) parameters:
N1                  = 4;    %% cube has size (N1 x N1 x N3)
N3                  = N1;   %% cube has size (N1 x N1 x N3)
Nstep               = 3;    %% sliding step to process every next reference cube
N2                  = 16;   %% maximum number of similar cubes
Ns                  = 11;   %% length of the side of the search neighborhood for full-search cube-matching
tau_match           = 3000;%2.9; %% threshold for the cube-distance (d-distance)
lambda_thr4D        = 2.7;  %% threshold parameter for the hard-thresholding in 4D transform domain

%%%% Wiener filtering parameters:
N1_wiener           = N1;
N3_wiener           = N1_wiener;
Nstep_wiener        = Nstep;
N2_wiener           = 32;
Ns_wiener           = Ns;
tau_match_wiener    = 400;%0.4;

%%%% Cube-matching parameters:
synchronous = 0;  %% if 1, the grouped cubes have coordinates lying in the same slice
decLevel    = 0;  %% decimation levels of the dyadic wavelet 2D transform 
                  %  0 means full decimation, higher values decrease the dec. number



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Create transform matrices, etc.
%%%%
[Tfor, Tinv]     = getTransfMatrix(N1, transform_3D_HT_name, decLevel);
[TforW, TinvW]   = getTransfMatrix(N1_wiener, transform_3D_Wiener_name, 0);
[Tfor3, Tinv3]   = getTransfMatrix(N3, transform_3D_HT_name, decLevel);
[Tfor3W, Tinv3W] = getTransfMatrix(N3_wiener, transform_3D_Wiener_name, 0);

Tfor4 = cell(1,max(N2,N2_wiener));
Tinv4 = cell(1,max(N2,N2_wiener));
for hpow = 0:ceil(log2(max(N2,N2_wiener))),
    h = 2^hpow;
    [Tfor4rd, Tinv4rd] = getTransfMatrix(h, transform_4D_name, 0);
    Tfor4{h} = double(Tfor4rd);
    Tinv4{h} = double(Tinv4rd');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% If needed, read images, generate noise, etc.
%%%%
if (exist('y','var') && ndims(y)==3 && ~exist('z','var'))
    if sigma~=0
        randn('seed', 0);
        if is_rician
            z = sqrt( (y+sigma.*randn(size(y))).^2 + (sigma.*randn(size(y))).^2 );
        else
            z = y + sigma*randn(size(y));
        end
    else
        error('Standard deviation sigma must be greater than 0.')
    end
else
    % convert z and y to double precision if needed
    z = double(z);
    y = double(y);
end
if (size(z,4) ~= 1) || (size(z,3) == 1),
    error('BM4D accepts only grayscale volumetric images.');
end
% Check if the true image y is a valid one; if not, then we cannot compute PSNR, etc.
y_is_invalid_image = (length(size(z)) ~= length(size(y))) | ...
    (size(z,1) ~= size(y,1)) | (size(z,2) ~= size(y,2)) | (size(z,3) ~= size(y,3));
if (y_is_invalid_image),
    dump_output_information = 0;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Rician data forward stabilization
%%%%
if is_rician && ((sigma~=0 && ~exist('riceVST.m','file')) || ...
        (sigma==0 && ~exist('ricePairInversion.m','file')))
    error(['VST framework for Rician-distributed data not found. ',...
        'Please download it from http://www.cs.tut.fi/~foi/RiceOptVST/.'])
end
if is_rician && sigma~=0
    smoothVST = 'A';
    sigma_rice = sigma;
    sigma = 1;
    z = riceVST(z,sigma_rice,smoothVST);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Normalizing input
%%%%
maxtransformed = max(z(:));
mintransformed = min(z(:));
scale_range = 0.7;
scale_shift = (1-scale_range)/2;
z = (z-mintransformed) / (maxtransformed-mintransformed);
z = z*scale_range+scale_shift;
sigma = sigma/(maxtransformed-mintransformed);
sigma = sigma*scale_range;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 1. Produce the basic estimate by HT filtering
%%%%
tic;
if is_rician && sigma==0
    [y_hat, sigma_hat] = bm4d_thr_rice_mex(z, Nstep, N1, N2, N3,...
        lambda_thr4D, tau_match*N1*N1*N1/(255*255), (Ns-1)/2, synchronous, ...
        0, alpha, Tfor, Tinv', Tfor3, Tinv3', Tfor4, Tinv4 );
else
    % [y_hat, sigma_hat] = bm4d_thr_mex(z, Nstep, N1, N2, N3,...
    %     lambda_thr4D, tau_match*N1*N1*N1/(255*255), (Ns-1)/2, synchronous, ...
    %     sigma, alpha, Tfor, Tinv', Tfor3, Tinv3', Tfor4, Tinv4 );
    y_hat=NLPCAnp(z,sigma,2.2);

end
estimate_elapsed_time = toc;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 2. Produce the final estimate by Wiener filtering (using the
%%%%  hard-thresholding initial estimate)
%%%%
tic;
if alpha~=1.0 || ~do_wiener
    y_est     = y_hat;
    sigma_est = sigma_hat;
else
    if is_rician && sigma==0
        [y_est, sigma_est] = bm4d_wie_rice_mex(z, y_hat, Nstep_wiener, N1_wiener, N2_wiener, N3_wiener, ...
            tau_match_wiener*N1_wiener*N1_wiener*N1_wiener/(255*255), (Ns_wiener-1)/2, ...
            synchronous, 0, TforW, TinvW', Tfor3W, Tinv3W', Tfor4, Tinv4 );
    else
        [y_est, sigma_est] = bm4d_wie_mex(z, y_hat, Nstep_wiener, N1_wiener, N2_wiener, N3_wiener, ...
            tau_match_wiener*N1_wiener*N1_wiener*N1_wiener/(255*255), (Ns_wiener-1)/2, ...
            synchronous, sigma, TforW, TinvW', Tfor3W, Tinv3W', Tfor4, Tinv4 );
    end
end
wiener_elapsed_time = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Restoring original range
%%%%
y_est     = (y_est-scale_shift)/scale_range;
y_est     = y_est*(maxtransformed-mintransformed)+mintransformed;
sigma_est = sigma_est/scale_range;
sigma_est = sigma_est*(maxtransformed-mintransformed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Rician data inverse stabilization
%%%%
if is_rician && sigma~=0
    y_est = riceVST_EUI(y_est,sigma_rice,smoothVST);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Calculate the final estimate's PSNR, print it, and show the
%%%% denoised image next to the noisy one
%%%%
PSNR = 0; %% Remains 0 if the true image y is not available
if ~y_is_invalid_image
    ind  = y>background;
    PSNR = 10*log10(1/mean((y(ind)-y_est(ind)).^2)); % y is valid
end

if dump_output_information
    fprintf('Final estimate, PSNR: %.2f dB \nTotal execution time: %.1f sec \n', ...
        PSNR, wiener_elapsed_time + estimate_elapsed_time);
end

return;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [Transf_Matrix, invTransf_Matrix] = getTransfMatrix (N, transform_type, Nden)
%
% Create forward and inverse transform matrices, which allow for perfect
% reconstruction. The forward transform matrix is normalized so that the
% l2-norm of each basis element is 1.
%
% [Tforward, Tinverse] = getTransfMatrix (N, transform_type, dec_levels)
%
%  INPUTS:
%
%   N           --> Size of the transform (for wavelets, must be 2^K)
%
%   transform_type  --> 'dct', 'dst', 'hadamard', or anything that is
%                       listed by 'help wfilters' (bi-orthogonal wavelets)
%                       'DCrand' -- an orthonormal transform with a DC and all
%                       the other basis elements of random nature
%
%   dec_levels      --> If a wavelet transform is generated, this is the
%                       desired decomposition level. Must be in the
%                       range [0, log2(N)-1], where "0" implies
%                       full decomposition.
%
%  OUTPUTS:
%
%   Tforward        --> (N x N) Forward transform matrix
%
%   Tforward        --> (N x N) Inverse transform matrix
%

if ~exist('Nden','var')
    Nden = 0;
end

if N == 1,
    Transf_Matrix = 1;
elseif strcmp(transform_type, 'dct') == 1,
    Transf_Matrix    = dct(eye(N));
elseif strcmp(transform_type, 'dst') == 1,
    Transf_Matrix    = dst(eye(N));
elseif strcmp(transform_type, 'DCrand') == 1,
    x = randn(N); x(1:end,1) = 1; [Q,R] = qr(x);
    if (Q(1) < 0),
        Q = -Q;
    end;
    Transf_Matrix = Q';
elseif strcmp(transform_type, 'hadamard') == 1,
    Transf_Matrix    = hadamard(N);
else %% wavelet transform

    %%% Set periodic boundary conditions, to preserve bi-orthogonality
    dwtmode('per','nodisp');

    [LO_D, HI_D, LO_R, HI_R] = wfilters(transform_type);
    for i = 1:N
        Transf_Matrix(i,:)=waverec(circshift([1 zeros(1,N-1)],[Nden i-1]), ...
            2.^[Nden Nden:log2(N)], LO_D, -HI_D);  %% construct transform matrix
    end
end

%%% Normalize the basis elements
Transf_Matrix = (Transf_Matrix' * diag(sqrt(1./sum(Transf_Matrix.^2,2))))';

%%% Compute the inverse transform matrix
invTransf_Matrix = inv(Transf_Matrix);

return;
