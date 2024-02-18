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

end