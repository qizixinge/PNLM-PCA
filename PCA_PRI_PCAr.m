function [fim] = PCA_PRI_PCAr(nim, x1, noise_map, option)
%external noise map won't be use when option = 0
%noise map should have the same size as nim
global etta
global fai
[dnim,noise_mapr]=NLPCA(nim,x1(1),x1(2),x1(3));
yichang=find(isnan(dnim));
dnim(yichang)=nim(yichang);
A1=size(nim,1);B1=size(nim,2);C1=size(nim,3);
if option==0
    cfai=(0.9953*noise_mapr-1.716)./(noise_mapr-1.787);%proposed by the author
    idx=(noise_mapr<1.913); cfai(idx)=0;
    %For Manjon et al. (2015)
    %cfai=(0.9846*noise_mapr-1.6331)./(noise_mapr-1.86+0.1175);
    %idx=(noise_mapr<1.86);cfai(idx)=0;
    noise_map=cfai.*noise_mapr;
end
%The fixed point formula of SNR in Koay and Baser (2006) 
%establishes the analytical relationship between the magnitude 
%signal-to-noise ratio and the correction factor. We can sample 
%SNR at dense intervals to obtain the corresponding $\xi$ and
%$\gamma$ array, and then points of ($\gamma$, $\frac{1}{\sqrt{\xi}}$) 
%are fitted to obtain a model closer to the analytical solution. The 
%difference between the above two mapping functions is large at low SNR.
        eetta = dnim ./ noise_map;
        eetta=reshape(eetta,1,A1*B1*C1);
        idxe = (eetta<min(etta)); eetta(idxe) = min(etta);%avoid outliers
        idxd = (eetta>max(etta)); eetta(idxd) = max(etta);%avoid outliers
        temp = interp1(etta, fai, eetta, 'linear');%arrayfun(@ettainv, eetta);
        temp=reshape(temp,A1,B1,C1);
        gim = noise_map.* temp;
        NNim = padarray(nim,[5, 5, 5],'symmetric');
        GGim = padarray(gim, [5, 5, 5], 'symmetric');
        Noise_map=padarray(noise_map,[5, 5, 5],'symmetric');%sigmago*ones((A1+10),(B1+10),(C1+10));%
        miu = imboxfilt3( GGim  , 'padding' , 'symmetric' ) ;%3D mean matrix of ggim      
        PRINLMfim = RI_NLM(GGim, NNim, Noise_map, miu,0.4);%~5min

[fim,~]=NLPCA(PRINLMfim,x1(1),x1(2),x1(3));
yichang=find(isnan(fim));
fim(yichang)=PRINLMfim(yichang);


end

