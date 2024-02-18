function [y_est] = bm4dw(z,y_hat,sigma)
%z与y_hat不用提前归一
% if dtr==8
%     ma=1;%255
% elseif dtr == 12
%     ma=4096/255;
% end
z = double(z);
y_hat = double(y_hat);
%[Tfor, Tinv]     = getTransfMatrix(4, 'bior1.5', 0);
[TforW, TinvW]   = getTransfMatrix(4, 'dct', 0);
%[Tfor3, Tinv3]   = getTransfMatrix(4, 'bior1.5', 0);
[Tfor3W, Tinv3W] = getTransfMatrix(4, 'dct', 0);

Tfor4 = cell(1,max(16,32));
Tinv4 = cell(1,max(16,32));
for hpow = 0:ceil(log2(max(16,32)))
    h = 2^hpow;
    [Tfor4rd, Tinv4rd] = getTransfMatrix(h, 'haar', 0);
    Tfor4{h} = double(Tfor4rd);
    Tinv4{h} = double(Tinv4rd');
end


    maxtransformed = max(z(:));
    mintransformed = min(z(:));
    scale_range = 0.7;
    scale_shift = (1-scale_range)/2;
    z = (z-mintransformed) / (maxtransformed-mintransformed);
    z = z*scale_range+scale_shift;
    y_hat = (y_hat-mintransformed) / (maxtransformed-mintransformed);
    y_hat = y_hat*scale_range+scale_shift;
    sigma = sigma/(maxtransformed-mintransformed);
    sigma = sigma*scale_range;

% [y_est, sigma_est] = bm4d_wie_mex(z, y_hat, Nstep_wiener, N1_wiener, N2_wiener, N3_wiener, ...
%             tau_match_wiener*N1_wiener*N1_wiener*N1_wiener/(255*255), (Ns_wiener-1)/2, ...
%             synchronous, sigma, TforW, TinvW', Tfor3W, Tinv3W', Tfor4, Tinv4 );
    [y_est, ~] = bm4d_wie_mex(z, y_hat, 3, 4, 32, 4, ...
            400*4*4*4/(255*255), (11-1)/2, ...
            0, sigma, TforW, TinvW', Tfor3W, Tinv3W', Tfor4, Tinv4 );

y_est     = (y_est-scale_shift)/scale_range;
    y_est     = y_est*(maxtransformed-mintransformed)+mintransformed;


end

