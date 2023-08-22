function [fim] = PD(nim, gim, noise_map, x1)
NNim1 = padarray(nim,[5, 5, 5],'symmetric');
GGim = padarray(gim,[5, 5, 5],'symmetric');
Noise_map=padarray(noise_map,[5, 5, 5],'symmetric');
miu = imboxfilt3( GGim  , 'padding' , 'symmetric' ) ;%3D mean matrix of ggim
PRINLMfim = RI_NLM(GGim, NNim1, Noise_map, miu,0.4);%~5min

[fim,~]=NLPCA(PRINLMfim,x1(1), x1(2), x1(3));

yichang=find(isnan(fim));
fim(yichang)=PRINLMfim(yichang);


end

