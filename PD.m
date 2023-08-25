function [fim] = PD(nim, gim, noise_map, x1)
%noise map can be obtained by a noise estimator
%noise map should have the same size as nim
NNim1 = padarray(nim,[5, 5, 5],'symmetric');
GGim = padarray(gim,[5, 5, 5],'symmetric');
Noise_map=padarray(noise_map,[5, 5, 5],'symmetric');
miu = imboxfilt3( GGim  , 'padding' , 'symmetric' ) ;%3D mean matrix of ggim
PRINLMfim = RI_NLM(GGim, NNim1, Noise_map, miu,0.4);%~5min

[fim,~]=NLPCA(PRINLMfim,x1(1), x1(2), x1(3));
%when x1(1) * x1(2) is constant, the performance of PD remains unchanged
%we could let x1(3) = x1(1) * x1(2) according to the discussion in the
%paper
yichang=find(isnan(fim));
fim(yichang)=PRINLMfim(yichang);


end

