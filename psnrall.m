function [psnr] = psnrall(x)

global im
global nnim
% global etta 
% global fai
% global Nim
x=x';
psnr=zeros(50,1);
x(1,:)=round(x(1,:));x(2,:)=round(x(2,:));x(3,:)=round(x(3,:));
for i=1:50
    if ( x(2,i)-(2*x(3,i)+2-x(1,i) )^3>0  )||( x(1,i)-x(3,i)-1>0  )||( -x(2,i)+x(1,i)^3>0  )%||(x(4,i)-x(5,i)>0) 
        psnr(i)=0;
    else
        dnim=NLPCApso(single(nnim),x(1,i),x(2,i),x(3,i),x(4,i),x(5,i));
%         sigma2 = imboxfilt3(Nim.^2 , 3 , 'padding' , 'symmetric')-imboxfilt3(Nim , 3 , 'padding' , 'symmetric').^2;
%         sigma = sqrt(sigma2);
%         logamma =  imboxfilt3(Nim , 3 , 'padding' , 'symmetric')./sigma;
%         cfai =  (0.9953*logamma-1.716)./(logamma-1.787) ; idx = (logamma<1.913); cfai(idx) = 0;
%         noise_mapg = cfai.*noise_mapr;
%         noise_mapg=imboxfilt3(noise_mapg , 15 , 'padding' , 'symmetric');
%         Cov = std( noise_mapg , 1 , 'all' )/mean( noise_mapg , 'all');
%         if Cov <0.15
%             noise_map = mean( noise_mapg , 'all')*ones(181,217,181);
%         else
%             noise_map = noise_mapg;
%         end
%        noise_map=2.55*ones(181,217,181);
%        eetta = dnim ./ noise_map;
%        eetta=reshape(eetta,1,181*217*181);
%        idxe = (eetta<min(etta)); eetta(idxe) = min(etta);%avoid outliers
%        temp = interp1(etta, fai, eetta, 'linear');%arrayfun(@ettainv, eetta);
%        temp=reshape(temp,181,217,181);
%        gim = noise_map.* temp;
%         NNim = padarray(Nim,[5, 5, 5],'symmetric');
%         GGim = padarray(gim,[5, 5, 5],'symmetric');
%         Noise_map=padarray(noise_map,[5, 5, 5],'symmetric');
%         miu = imboxfilt3( GGim  , 'padding' , 'symmetric' ) ;%3D mean matrix of ggim
%         alpha=0.4;
%         PRINLMfim = RI_NLM(GGim, NNim, Noise_map, miu,alpha);%~5min
        %PRINLMfim = RI_NLM(GGim, GGim, Noise_map, miu);%~2min
        index = find(im>0);
        psnr(i) = -20*log10(255/sqrt(mean((im(index)-dnim(index)).^2)));
%        psnr(i) = -20*log10(255/sqrt(mean((im(index)-gim(index)).^2)));
    end
end
end

