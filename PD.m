function [fim] = PD(nim, gim, noise_map, T)
%noise map can be obtained by a noise estimator
%noise map should have the same size as nim
NNim1 = padarray(nim,[5, 5, 5],'symmetric');
GGim = padarray(gim,[5, 5, 5],'symmetric');
Noise_map=padarray(noise_map,[5, 5, 5],'symmetric');
miu = imboxfilt3( GGim  , 'padding' , 'symmetric' ) ;%3D mean matrix of ggim
PRINLMfim = RI_NLM(GGim, NNim1, Noise_map, miu,0.4);%~5min

[fim,~]=NLPCA(PRINLMfim, 1, T, T);
%when tau * beta is constant, the performance of PD remains unchanged
yichang=find(isnan(fim));
fim(yichang)=PRINLMfim(yichang);


end

