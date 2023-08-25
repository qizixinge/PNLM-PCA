%%
%a series of experiemnts
%We have obtained the set of parameters that make the algorithm have the best denoising performance
% d=3;M=38;w=3;tauopt=1.1418;betaopt=1.7557;
% taurtn=1.397 ;betartn=1.4350;
% x1=[d,M,w,tauopt,betaopt];x2=[d,M,w,taurtn,betartn];
noisesm=255*0.01*[1 3 5 7 9];
noiselr=255*0.01*[1 7 13 19 25];
x1=[2.12 1.16 2.46];
alpha=0.4;
%x1=[4 64 4 2.2 1.16];
%%
%precomputation
%ptheta=linspace(0,50,1000);
%effSNR=sqrt((2+ptheta.^2)./kexi(ptheta)-1);
global etta
global fai
load('precomputation.mat')
%%
%load in data
% read the file
imsize=[181,217,181];
file1 ='t1_icbm_normal_1mm_pn0_rf0.rawb';
%file1='brainwebmag.rawb';
fid1 = fopen(file1,'r');    
im1=zeros(imsize(1:3));
for z=1:imsize(3)    
    im1(:,:,z) = fread(fid1,imsize(1:2));
end
fclose(fid1);

file2 ='t2_ai_msles2_1mm_pn0_rf0.rawb';
%file2='brainwebimaginary.rawb';
fid2 = fopen(file2,'r');    
im2=zeros(imsize(1:3));
for z=1:imsize(3)    
    im2(:,:,z) = fread(fid2,imsize(1:2));
end
fclose(fid2);

% file3 ='T1w_acpc_dc_restore.nii';
% fid3 = fopen(file3,'r');    
% im3=zeros(imsize(1:3));
% for z=1:imsize(3)    
%     im3(:,:,z) = fread(fid3,imsize(1:2));
% end
% fclose(fid3);
ker=[1 1 1];
index=find(im1>0);
%%
%First check the correctness of the algorithm qualitatively
%opt vs rtn, stationary
%This section is estimated to take 10 hours to run
% PSNRopt1=zeros(1,5);PSNRrtn1=zeros(1,5);
% PSNRopt2=zeros(1,5);PSNRrtn2=zeros(1,5);
% ssimopt1=zeros(1,5);ssimrtn1=zeros(1,5);
% ssimopt2=zeros(1,5);ssimrtn2=zeros(1,5);
% x1=[d,M,w,tauopt,betaopt];x2=[d,M,w,taurtn,betartn];
% for i=1:5
%     Nim1=ricernd(im1, noisesm(i)*ones(imsize(1:3)));
%     nnim1 = Nim1;
%     Nim2=ricernd(im2, noisesm(i)*ones(imsize(1:3)));
%     nnim2 = Nim2;
%     [PSNRopt1(i),ssimopt1(i),PRINLMfim15]=psnrallnv(im1,nnim1,Nim1,x1,0.4);
%     [PSNRrtn1(i),ssimrtn1(i),PRINLMfim14]=psnrallnv(im1,nnim1,Nim1,x2,0.5);
%     [PSNRopt2(i),ssimopt2(i),PRINLMfim25]=psnrallnv(im2,nnim2,Nim2,x1,0.4);
%     [PSNRrtn2(i),ssimrtn2(i),PRINLMfim24]=psnrallnv(im2,nnim2,Nim2,x2,0.5);
% end
%%
%timing algorithm
%NL-PCA
tic
Nim1=ricernd(im1, 2.55*ones(imsize(1:3)));
[dnim,~]=NLPCA(single(Nim1),x1(1),x1(2),x1(3));
toc

% sigma2 = imboxfilt3(Nim1.^2 , 3 , 'padding' , 'symmetric')-imboxfilt3(Nim1 , 3 , 'padding' , 'symmetric').^2;
%         sigma = sqrt(sigma2);
%         logamma =  imboxfilt3(Nim1 , 3 , 'padding' , 'symmetric')./sigma;
%         cfai =  (0.9953*logamma-1.716)./(logamma-1.787) ; idx = (logamma<1.913); cfai(idx) = 0;
%         noise_mapg = cfai.*noise_mapr;
%         noise_mapg=imboxfilt3(noise_mapg , 15 , 'padding' , 'symmetric');
%         Cov = std( noise_mapg , 1 , 'all' )/mean( noise_mapg , 'all');
%         if Cov <0.15
%             noise_map = mean( noise_mapg , 'all')*ones(181,217,181);
%         else
%             noise_map = noise_mapg;
%         end
        noise_map=2.55*ones(181,217,181);
        eetta = dnim ./ noise_map;
        eetta=reshape(eetta,1,181*217*181);
        idxe = (eetta<min(etta)); eetta(idxe) = min(etta);%avoid outliers
        temp = interp1(etta, fai, eetta, 'linear');%arrayfun(@ettainv, eetta);
        temp=reshape(temp,181,217,181);
        gim = noise_map.* temp;
NNim1 = padarray(Nim1,[5, 5, 5],'symmetric');
        GGim = padarray(gim,[5, 5, 5],'symmetric');
        Noise_map=padarray(noise_map,[5, 5, 5],'symmetric');
        miu = imboxfilt3( GGim  , 'padding' , 'symmetric' ) ;%3D mean matrix of ggim
        
%PRI-NLM
tic
        PRINLMfim = RI_NLM(GGim, NNim1, Noise_map, miu,0.4);%~5min
toc
%

%psnr = 20*log10(255/sqrt(mean((im1(index)-PRINLMfim(index)).^2)));
%%
%noisy image
psnr1=zeros(1,5); psnr2=zeros(1,5);ssim1=zeros(1,5);ssim2=zeros(1,5);
index1=find(im1>0);index2=find(im2>0);
for i=1:5
    nnim1=ricernd(im1,noisesm(i)*ones(181,217,181));
    psnr1(i)=20*log10(255/sqrt(mean((im1(index1)-nnim1(index1)).^2)));
    ssim1(i)=ssim_index3d(im1,nnim1,[1 1 1],index1);
end    
for i=1:5
    nnim2=ricernd(im2,noisesm(i)*ones(181,217,181));
    psnr2(i)=20*log10(255/sqrt(mean((im2(index2)-nnim2(index2)).^2)));
    ssim2(i)=ssim_index3d(im2,nnim2,[1 1 1],index2);
end    

%%
%draw map 9%
% residual1=im1-PRINLMfim15;
% residual2=im2-PRINLMfim25;
% colormap(gray);
% n=round(size(im1,3)/2);
% subplot(4,4,1),imagesc(imrotate(im1(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,4,2),imagesc(imrotate(Nim1(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,4,3),imagesc(imrotate(PRINLMfim15(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,4,4),imagesc(imrotate(residual1(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,4,5),imagesc(imrotate(im2(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,4,6),imagesc(imrotate(Nim2(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,4,7),imagesc(imrotate(PRINLMfim25(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,4,8),imagesc(imrotate(residual2(:,:,n),90));set(gca,'xtick',[],'ytick',[]);


%%
%median filter
PSNRmed=zeros(1,5);ssimed=zeros(1,5);
PSNRnm=zeros(1,5);ssimnm=zeros(1,5);
for i=1:5
        Nim1=ricernd(im1, noiselr(i)*ones(imsize(1:3)));
        nnim1 = medfilt3(Nim1);
        [dnim1,~]=NLPCA(Nim1,x1(1),x1(2),x1(3));
        [dnim2,~]=NLPCA(nnim1,x1(1),x1(2),x1(3));
        PSNRmed(i)=20*log10(255/sqrt(mean((im1(index)-dnim2(index)).^2)));
        ssimed(i)=ssim_index3d( dnim2 , im1 , ker , index );
        PSNRnm(i)=20*log10(255/sqrt(mean((im1(index)-dnim1(index)).^2)));
        ssimnm(i)=ssim_index3d( dnim1 , im1 , ker , index );
%         if i==5
%             xuyao1=PRINLMfimm;
%             xuyao2=PRINLMfimnm;
%         end
end

% Nim1=ricernd(im1, 255*0.09*ones(imsize(1:3)));
% 
% residualm=im1-xuyao1;
% residualnm=im1-xuyao2;
% colormap(gray);
% n=round(size(im1,3)/2);
% subplot(4,4,9),imagesc(imrotate(im1(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,4,10),imagesc(imrotate(Nim1(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,4,11),imagesc(imrotate(xuyao1(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,4,12),imagesc(imrotate(residualm(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,4,13),imagesc(imrotate(xuyao2(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,4,14),imagesc(imrotate(residualnm(:,:,n),90));set(gca,'xtick',[],'ytick',[]);


%%
%Universality test ——PRINLM
% im0= ;
% psnrexam=zeros(1,5);ssimexam=zeros(1,5);
% for i=1:5
%     noise3D=noisesm(i)*ones(imsize(1:3)); 
%     Nim1=ricernd(im1, noise3D);
%     Noise_map= noisesm(i)*ones( 191,227,191);
%     cim=normrnd(im1, RMSE(i)*ones(imsize(1:3)));%construct a pre-filtered image
%     cgim=padarray(cim ,[5, 5, 5],'symmetric');
%     miu=imboxfilt3( cgim  , 'padding' , 'symmetric' ) ;%3D mean matrix of ggim
%     alpha=0.4;NNim1=padarray(Nim1,[5, 5, 5],'symmetric');
%     PRINLMfim2 = RI_NLM(cgim, NNim1, Noise_map, miu,alpha);%~5min
%     PRINLMfim = RI_NLM(GGim, NNim1, Noise_map, miu,alpha);%~5min
%     psnrexam(i)= 20*log10(255/sqrt(mean((im0(index)-PRINLMfim(index)).^2)));
%     ssimexam(i)= ssim_index3d( PRINLMfim , im0 , ker , index );  
% end    
%%
%Horizontal comparison im1
psnrd=zeros(1,5);ssimd=zeros(1,5);
psnrdgd=zeros(1,5);ssimdgd=zeros(1,5);
psnrdg=zeros(1,5);ssimdg=zeros(1,5);
psnrdd=zeros(1,5);ssimdd=zeros(1,5);
psnrdgp=zeros(1,5);ssimdgp=zeros(1,5);
psnrdgpp=zeros(1,5);ssimdgpp=zeros(1,5);
psnrdgpd=zeros(1,5);ssimdgpd=zeros(1,5);
for i=1:5
        nnim1=ricernd(im1, noisesm(i)*ones(imsize(1:3)));
        [dnim,~]=NLPCA(single(nnim1),x1(1),x1(2),x1(3));
        noise_map=noisesm(i)*ones(181,217,181);
        eetta = dnim ./ noise_map;
        eetta=reshape(eetta,1,181*217*181);
        idxe = (eetta<min(etta)); eetta(idxe) = min(etta);%avoid outliers
        temp = interp1(etta, fai, eetta, 'linear');%arrayfun(@ettainv, eetta);
        temp=reshape(temp,181,217,181);
        gim = noise_map.* temp;
        [dgdim,~]=NLPCA(single(gim),x1(1),x1(2),x1(3));
        yichang=find(isnan(dgdim));
        dgdim(yichang)=gim(yichang);
        [ddim,~]=NLPCA(single(dnim),x1(1),x1(2),x1(3));
        index = find(im1>0);
        psnrd(i) = 20*log10(255/sqrt(mean((im1(index)-dnim(index)).^2)));
        ssimd(i)=ssim_index3d( dnim , im1 , ker , index );
        psnrdg(i) = 20*log10(255/sqrt(mean((im1(index)-gim(index)).^2)));
        ssimdg(i)=ssim_index3d( gim , im1 , ker , index );
        psnrdgd(i) = 20*log10(255/sqrt(mean((im1(index)-dgdim(index)).^2)));
        ssimdgd(i)=ssim_index3d( dgdim , im1 , ker , index );
        psnrdd(i)=20*log10(255/sqrt(mean((im1(index)-ddim(index)).^2)));
        ssimdd(i)=ssim_index3d( ddim , im1 , ker , index );  
        
        NNim1 = padarray(nnim1,[5, 5, 5],'symmetric');
        GGim = padarray(gim,[5, 5, 5],'symmetric');
        Noise_map=noisesm(i)*ones(191,227,191);%padarray(noise_map,[5, 5, 5],'symmetric');
        miu = imboxfilt3( GGim  , 'padding' , 'symmetric' ) ;%3D mean matrix of ggim
        PRINLMfim = RI_NLM(GGim, NNim1, Noise_map, miu,alpha);%~5min
        PPim=padarray(PRINLMfim,[5 5 5],'symmetric'  );
        miup= imboxfilt3( PPim  , 'padding' , 'symmetric' ) ;
        P2fim=RI_NLM(PPim, NNim1, noisesm(i)*ones(191,227,191), miup,0.4);%~5min
        psnrdgpp(i)=20*log10(255/sqrt(mean((im1(index)-P2fim(index)).^2)));
        ssimdgpp(i)=ssim_index3d( P2fim , im1 , ker , index );
        psnrdgp(i)=20*log10(255/sqrt(mean((im1(index)-PRINLMfim(index)).^2)));
        ssimdgp(i)=ssim_index3d( PRINLMfim , im1 , ker , index );
        [dpim,~]=NLPCA(single(PRINLMfim),x1(1),x1(2),x1(3));
        yichang=find(isnan(dpim));
        dpim(yichang)=PRINLMfim(yichang);
        psnrdgpd(i)=20*log10(255/sqrt(mean((im1(index)-dpim(index)).^2)));
        ssimdgpd(i)=ssim_index3d( dpim , im1 , ker , index );
end

%%
%Horizontal comparison im2
psnrd2=zeros(1,5);ssimd2=zeros(1,5);
psnrdgd2=zeros(1,5);ssimdgd2=zeros(1,5);
psnrdg2=zeros(1,5);ssimdg2=zeros(1,5);
psnrdd2=zeros(1,5);ssimdd2=zeros(1,5);ker=[1 1 1];
psnrdgpp2=zeros(1,5);ssimdgpp2=zeros(1,5);
psnrdgp2=zeros(1,5);ssimdgp2=zeros(1,5);
psnrdgpd2=zeros(1,5);ssimdgpd2=zeros(1,5);
for i=1:5
        nnim2=ricernd(im2, noisesm(i)*ones(imsize(1:3)));

        [dnim,~] =NLPCA(single(nnim2),x1(1),x1(2),x1(3));
        noise_map=noisesm(i)*ones(181,217,181);
        eetta = dnim ./ noise_map;
        eetta=reshape(eetta,1,181*217*181);
        idxe = (eetta<min(etta)); eetta(idxe) = min(etta);%avoid outliers
        temp = interp1(etta, fai, eetta, 'linear');%arrayfun(@ettainv, eetta);
        temp=reshape(temp,181,217,181);
        gim = noise_map.* temp;
        [dgdim,~] =NLPCA(single(gim),x1(1),x1(2),x1(3));
        yichang=find(isnan(dgdim));
        dgdim(yichang)=gim(yichang);
        [ddim,~]=NLPCA(single(dnim),x1(1),x1(2),x1(3));
        index = find(im2>0);
        psnrd2(i) = 20*log10(255/sqrt(mean((im2(index)-dnim(index)).^2)));
        ssimd2(i)=ssim_index3d( dnim , im2 , ker , index );
        psnrdg2(i) = 20*log10(255/sqrt(mean((im2(index)-gim(index)).^2)));
        ssimdg2(i)=ssim_index3d( gim , im2 , ker , index );
        psnrdgd2(i) = 20*log10(255/sqrt(mean((im2(index)-dgdim(index)).^2)));
        ssimdgd2(i)=ssim_index3d( dgdim , im2 , ker , index );
        psnrdd2(i)=20*log10(255/sqrt(mean((im2(index)-ddim(index)).^2)));
        ssimdd2(i)=ssim_index3d( ddim , im2 , ker , index );
%     miut=miu(2:180,2:216,2:180);
%     T1=miut(1,1,:);T2=miut(1,215,:);
%     T3=miut(:,1,1);T4=miut(1,:,179);
%     T5=miut(1,:,1);T6=miut(179,1,:);
%     T7=miut(179,215,:);T8=miut(:,1,179);
%     T9=miut(:,215,179);T10=miut(179,:,1);
%     T11=miut(179,:,179);T12=miut(:,215,1);
%     T=[T1(:);T2(:);T3(:);T4(:);T5(:);T6(:);T7(:);T8(:);T9(:);T10(:);T11(:);T12(:)];
%     [T, ~]=sort(T); 
%     med=median(T);
%     temp = sum(T<1.5*med); medT = median(T(1:temp ));
%     noisestimate=sqrt(2/pi)*medT;
        NNim2 = padarray(nnim2,[5, 5, 5],'symmetric');
        GGim = padarray(gim,[5, 5, 5],'symmetric');
        Noise_map=noisesm(i)*ones(191,227,191);%Noise_map=padarray(noise_map,[5, 5, 5],'symmetric');
        miu = imboxfilt3( GGim  , 'padding' , 'symmetric' ) ;%3D mean matrix of ggim
        PRINLMfim = RI_NLM(GGim, NNim2, Noise_map, miu,alpha);%~5min
        PPim=padarray(PRINLMfim,[5 5 5],'symmetric'  );
        miup= imboxfilt3( PPim  , 'padding' , 'symmetric' ) ;
        P2fim=RI_NLM(PPim, NNim2, noisesm(i)*ones(191,227,191), miup,0.4);%~5min
        psnrdgpp2(i)=20*log10(255/sqrt(mean((im2(index)-P2fim(index)).^2)));
        ssimdgpp2(i)=ssim_index3d( P2fim , im2 , ker , index );
        psnrdgp2(i)=20*log10(255/sqrt(mean((im2(index)-PRINLMfim(index)).^2)));
        ssimdgp2(i)=ssim_index3d( PRINLMfim , im2 , ker , index );
        [dpim,~] =NLPCA(single(PRINLMfim),x1(1),x1(2),x1(3));
        yichang=find(isnan(dpim));
        dpim(yichang)=PRINLMfim(yichang);
        psnrdgpd2(i)=20*log10(255/sqrt(mean((im2(index)-dpim(index)).^2)));
        ssimdgpd2(i)=ssim_index3d( dpim , im2 , ker , index );
end

% colormap(gray);
% n=round(size(im2,3)/2);
% subplot(3,3,1),imagesc(imrotate(im2(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(3,3,2),imagesc(imrotate(Nim2(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(3,3,3),imagesc(imrotate(dnim(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(3,3,4),imagesc(imrotate(ddim(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(3,3,5),imagesc(imrotate(gim(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(3,3,6),imagesc(imrotate(dgdim(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(3,3,7),imagesc(imrotate(PRINLMfim(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(3,3,8),imagesc(imrotate(P2fim(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(3,3,9),imagesc(imrotate(dpim(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% psnrn=zeros(1,5);
% for i=1:5
%     Nim2=ricernd(im2, noisesm(i)*ones(imsize(1:3)));
%     psnrn(i)=20*log10(255/sqrt(mean((im2(index)-Nim2(index)).^2)));
% end
%%
%theoretical limit T1w
psnrth=zeros(1,5);ssimth=zeros(1,5);
%psnrthth=zeros(1,5);ssimthth=zeros(1,5);

for i=1:5
        noise3D=noisesm(i)*ones(imsize(1:3));
        Nim1=ricernd(im1, noise3D);
        Noise_map= noisesm(i)*ones( 191,227,191);
        tgim=padarray(im1,[5, 5, 5],'symmetric');
        miu=imboxfilt3( tgim  , 'padding' , 'symmetric' ) ;%3D mean matrix of ggim
        alpha=0.4;NNim1=padarray(Nim1,[5, 5, 5],'symmetric');
        PRINLMfim2 = RI_NLM(tgim, NNim1, Noise_map, miu,alpha);%~5min
%         PP=padarray(PRINLMfim2,[5, 5, 5],'symmetric');
%         miumiu=imboxfilt3( PP  , 'padding' , 'symmetric' ) ;%3D mean matrix of PRINLMfim2
%         PRINLMfim22=RI_NLM(PP, NNim1, Noise_map, miumiu,alpha);%~5min
        index = find(im1>0);
        psnrth(i)= 20*log10(255/sqrt(mean((im1(index)-PRINLMfim2(index)).^2)));
        ker=[1 1 1];
        ssimth(i)=ssim_index3d(PRINLMfim2,im1,ker,index);
%         psnrthth(i)= 20*log10(255/sqrt(mean((im1(index)-PRINLMfim22(index)).^2)));
%         ssimthth(i)=ssim_index3d(PRINLMfim22,im1,ker,index);
end
%%
% %simulation large psnr
% %p, pp, pd, pdp (scheme 2) im1
% psnrc=zeros(1,5);ssimc=zeros(1,5);
% 
% psnrcp=zeros(1,5);ssimcp=zeros(1,5);
% psnrcpp=zeros(1,5);ssimcpp=zeros(1,5);
% psnr_dncnn=[46.32, 40.47, 37.82, 36.20, 34.71];
% RMSE=255./10.^(psnr_dncnn/20);
% psnrcpd=zeros(1,5);ssimcpd=zeros(1,5);
% psnrcpdp=zeros(1,5);ssimcpdp=zeros(1,5);
% %ndncnn=RMSE/mean(sqrt( kexi(  ) ), 'all' );% difficult to solve, Non-monotonic variation
% %ndncnn=[1.24 2.40 3.23 3.87 4.56]; %derive from precomputation 
% for i=1:5
%         noise3D=noisesm(i)*ones(imsize(1:3)); 
%         Nim1=ricernd(im1, noise3D);
%         Noise_map= noisesm(i)*ones( 191,227,191);
%         %tgim=padarray(im1,[5, 5, 5],'symmetric');
%         %cim=ricernd(im1, RMSE(i)*ones(imsize(1:3)));%construct a pre-filtered image
%         cim=normrnd(im1, RMSE(i)*ones(imsize(1:3)));%construct a pre-filtered image
%         cgim=padarray(cim ,[5, 5, 5],'symmetric');
%         miu=imboxfilt3( cgim  , 'padding' , 'symmetric' ) ;%3D mean matrix of ggim
%         alpha=0.4;NNim1=padarray(Nim1,[5, 5, 5],'symmetric');
%         PRINLMfim2 = RI_NLM(cgim, NNim1, Noise_map, miu,alpha);%~5min
%         PP=padarray(PRINLMfim2,[5, 5, 5],'symmetric');
%         miumiu=imboxfilt3( PP  , 'padding' , 'symmetric' ) ;%3D mean matrix of PRINLMfim2
%         PRINLMfim22=RI_NLM(PP, NNim1, Noise_map, miumiu,alpha);%~5min
%         index = find(im1>0);
%         psnrcp(i)= 20*log10(255/sqrt(mean((im1(index)-PRINLMfim2(index)).^2)));
%         ker=[1 1 1];
%         ssimcp(i)=ssim_index3d(PRINLMfim2,im1,ker,index);
%         psnrcpp(i)= 20*log10(255/sqrt(mean((im1(index)-PRINLMfim22(index)).^2)));
%         ssimcpp(i)=ssim_index3d(PRINLMfim22,im1,ker,index);    
%         
%         [dpim,~]=NLPCA(single(PRINLMfim2),x1(1),x1(2),x1(3));
%         yichang=find(isnan(dpim));
%         dpim(yichang)=PRINLMfim2(yichang);
%         DP=padarray(dpim,[5, 5, 5],'symmetric');
%         dmiu=imboxfilt3( DP  , 'padding' , 'symmetric' ) ;%3D mean matrix of PRINLMfim2
%         pdpim=RI_NLM(DP, NNim1, Noise_map, dmiu,alpha);%~5min
%         index = find(im1>0);
%         psnrc(i)=20*log10(255/sqrt(mean((im1(index)-cim(index)).^2)));
%         psnrcpd(i)= 20*log10(255/sqrt(mean((im1(index)-dpim(index)).^2)));
%         ker=[1 1 1];
%         ssimc(i)=ssim_index3d(cim,im1,ker,index);
%         ssimcpd(i)=ssim_index3d(dpim,im1,ker,index);
%         psnrcpdp(i)= 20*log10(255/sqrt(mean((im1(index)-pdpim(index)).^2)));
%         ssimcpdp(i)=ssim_index3d(pdpim,im1,ker,index);    
% end
% colormap(gray);
% n=round(size(im1,3)/2);
% subplot(2,2,1),imagesc(imrotate(PRINLMfim2(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(2,2,2),imagesc(imrotate(PRINLMfim22(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(2,2,3),imagesc(imrotate(dpim(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(2,2,4),imagesc(imrotate(pdpim(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
%%
%simulation large psnr, but impact NLM this time
% psnrci=zeros(1,5);ssimci=zeros(1,5);
% psnrcpi=zeros(1,5);ssimcpi=zeros(1,5);
% psnrcppi=zeros(1,5);ssimcppi=zeros(1,5);
% psnr_dncnn=[46.32, 40.47, 37.82, 36.20, 34.71];
% RMSE=255./10.^(psnr_dncnn/20);
% psnrcpdi=zeros(1,5);ssimcpdi=zeros(1,5);
% psnrcpdpi=zeros(1,5);ssimcpdpi=zeros(1,5);
% %ndncnn=RMSE/mean(sqrt( kexi(  ) ), 'all' );% difficult to solve, Non-monotonic variation
% %ndncnn=[1.24 2.40 3.23 3.87 4.56]; %derive from precomputation 
% for i=1:5
%         noise3D=noisesm(i)*ones(imsize(1:3)); 
%         Nim1=ricernd(im1, noise3D);
%         Noise_map= noisesm(i)*ones( 191,227,191);
%         %tgim=padarray(im1,[5, 5, 5],'symmetric');
%         %cim=ricernd(im1, RMSE(i)*ones(imsize(1:3)));%construct a pre-filtered image
%         cim=normrnd(im1, RMSE(i)*ones(imsize(1:3)));%construct a pre-filtered image
%         cgim=padarray(cim ,[5, 5, 5],'symmetric');
%         miu=imboxfilt3( cgim  , 'padding' , 'symmetric' ) ;%3D mean matrix of ggim
%         alpha=0.4;NNim1=padarray(Nim1,[5, 5, 5],'symmetric');
%         PRINLMfim2 = RI_NLM(cgim, NNim1, Noise_map, miu,alpha);%~5min
%         PP=padarray(PRINLMfim2,[5, 5, 5],'symmetric');
%         miumiu=imboxfilt3( PP  , 'padding' , 'symmetric' ) ;%3D mean matrix of PRINLMfim2
%         PRINLMfim22=RI_NLM(PP, NNim1, Noise_map, miumiu,alpha);%~5min
%         index = find(im1>0);
%         psnrcpi(i)= 20*log10(255/sqrt(mean((im1(index)-PRINLMfim2(index)).^2)));
%         ker=[1 1 1];
%         ssimcpi(i)=ssim_index3d(PRINLMfim2,im1,ker,index);
%         psnrcppi(i)= 20*log10(255/sqrt(mean((im1(index)-PRINLMfim22(index)).^2)));
%         ssimcppi(i)=ssim_index3d(PRINLMfim22,im1,ker,index);    
%         
%         [dpim , noise_mapr]=NLPCA(single(PRINLMfim2),x1(1),x1(2),x1(3),x1(4),x1(5), 1);
%         yichang=find(isnan(dpim));
%         dpim(yichang)=PRINLMfim2(yichang);
%         DP=padarray(dpim,[5, 5, 5],'symmetric');
%         dmiu=imboxfilt3( DP  , 'padding' , 'symmetric' ) ;%3D mean matrix of PRINLMfim2
%         pdpim=RI_NLM(DP, NNim1, Noise_map, dmiu,alpha);%~5min
%         index = find(im1>0);
%         psnrci(i)=20*log10(255/sqrt(mean((im1(index)-cim(index)).^2)));
%         psnrcpdi(i)= 20*log10(255/sqrt(mean((im1(index)-dpim(index)).^2)));
%         ker=[1 1 1];
%         ssimci(i)=ssim_index3d(cim,im1,ker,index);
%         ssimcpdi(i)=ssim_index3d(dpim,im1,ker,index);
%         psnrcpdpi(i)= 20*log10(255/sqrt(mean((im1(index)-pdpim(index)).^2)));
%         ssimcpdpi(i)=ssim_index3d(pdpim,im1,ker,index);    
% end




%%
%simulation small psnr
%p, pp, pd, pdp (scheme 2) im1
psnrc=zeros(1,5);ssimc=zeros(1,5);

psnrcp=zeros(1,5);ssimcp=zeros(1,5);
psnrcpp=zeros(1,5);ssimcpp=zeros(1,5);
psnr_dncnn_brainweb=[45.00 39.61 36.75 34.93 33.43];
RMSE=255./10.^(psnr_dncnn_brainweb/20);
psnrcpd=zeros(1,5);ssimcpd=zeros(1,5);
psnrcpdp=zeros(1,5);ssimcpdp=zeros(1,5);
%ndncnn=RMSE/mean(sqrt( kexi(  ) ), 'all' );% difficult to solve, Non-monotonic variation
%ndncnn=[1.24 2.40 3.23 3.87 4.56]; %derive from precomputation 
for i=1:5
        noise3D=noisesm(i)*ones(imsize(1:3)); 
        Nim1=ricernd(im1, noise3D);
        Noise_map= noisesm(i)*ones( 191,227,191);
        %tgim=padarray(im1,[5, 5, 5],'symmetric');
        %cim=ricernd(im1, RMSE(i)*ones(imsize(1:3)));%construct a pre-filtered image
        cim=normrnd(im1, RMSE(i)*ones(imsize(1:3)));%construct a pre-filtered image
        cgim=padarray(cim ,[5, 5, 5],'symmetric');
        miu=imboxfilt3( cgim  , 'padding' , 'symmetric' ) ;%3D mean matrix of ggim
        alpha=0.4;NNim1=padarray(Nim1,[5, 5, 5],'symmetric');
        PRINLMfim1 = RI_NLM(cgim, NNim1, Noise_map, miu,alpha);%~5min
        PP=padarray(PRINLMfim1,[5, 5, 5],'symmetric');
        miumiu=imboxfilt3( PP  , 'padding' , 'symmetric' ) ;%3D mean matrix of PRINLMfim2
        PRINLMfim11=RI_NLM(PP, NNim1, Noise_map, miumiu,alpha);%~5min
        index = find(im1>0);
        psnrcp(i)= 20*log10(255/sqrt(mean((im1(index)-PRINLMfim1(index)).^2)));
        ker=[1 1 1];
        ssimcp(i)=ssim_index3d(PRINLMfim1,im1,ker,index);
        psnrcpp(i)= 20*log10(255/sqrt(mean((im1(index)-PRINLMfim11(index)).^2)));
        ssimcpp(i)=ssim_index3d(PRINLMfim11,im1,ker,index);    
        
        [dpim1,~] =NLPCA(single(PRINLMfim1),x1(1),x1(2),x1(3));
        yichang=find(isnan(dpim1));
        dpim1(yichang)=PRINLMfim1(yichang);
        DP=padarray(dpim1,[5, 5, 5],'symmetric');
        dmiu=imboxfilt3( DP  , 'padding' , 'symmetric' ) ;%3D mean matrix of PRINLMfim2
        pdpim1=RI_NLM(DP, NNim1, Noise_map, dmiu,alpha);%~5min
        index = find(im1>0);
        psnrc(i)=20*log10(255/sqrt(mean((im1(index)-cim(index)).^2)));
        psnrcpd(i)= 20*log10(255/sqrt(mean((im1(index)-dpim1(index)).^2)));
        ker=[1 1 1];
        ssimc(i)=ssim_index3d(cim,im1,ker,index);
        ssimcpd(i)=ssim_index3d(dpim1,im1,ker,index);
        psnrcpdp(i)= 20*log10(255/sqrt(mean((im1(index)-pdpim1(index)).^2)));
        ssimcpdp(i)=ssim_index3d(pdpim1,im1,ker,index);    
end
%%
colormap(gray);
n=round(size(im1,3)/2);
n1=round(size(im1,1)/3);nn1=2*n1;
n2=round(size(im1,2)/3);nn2=2*n2;
%line
%subplot(4,5,1),imagesc(imrotate(im1(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
%subplot(4,5,2),imagesc(imrotate(Nim1(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
%subplot(4,5,3),imagesc(imrotate(PRINLMfim1(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
%subplot(4,5,4),imagesc(imrotate(PRINLMfim11(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
%subplot(4,5,5),imagesc(imrotate(dpim1(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
%subplot(4,5,6),imagesc(imrotate(pdpim1(:,:,n),90));set(gca,'xtick',[],'ytick',[]);

%Partial Enlarged View, Center Right, bottom center
subplot(4,5,1),imagesc(imrotate(im1(nn1:end,n2:nn2,n),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,5,2),imagesc(imrotate(Nim1(nn1:end,n2:nn2,n),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,5,3),imagesc(imrotate(PRINLMfim1(nn1:end,n2:nn2,n),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,5,4),imagesc(imrotate(PRINLMfim11(nn1:end,n2:nn2,n),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,5,5),imagesc(imrotate(dpim1(nn1:end,n2:nn2,n),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,5,6),imagesc(imrotate(pdpim1(nn1:end,n2:nn2,n),90));set(gca,'xtick',[],'ytick',[]);

%%
%scheme 2 , im2
% psnrcpd=zeros(1,5);
% ssimcpd=zeros(1,5);
% psnrcpdp=zeros(1,5);
% ssimcpdp=zeros(1,5);
% for i=1:5
%     noise3D=noisesm(i)*ones(imsize(1:3)); 
%     Nim2=ricernd(im2, noise3D);
%     Noise_map= noisesm(i)*ones( 191,227,191);
%     %tgim=padarray(im1,[5, 5, 5],'symmetric');
%     cim=ricernd(im2, RMSE(i)*ones(imsize(1:3)));%construct a pre-filtered image
%     cgim=padarray(cim ,[5, 5, 5],'symmetric');
%     miu=imboxfilt3( cgim  , 'padding' , 'symmetric' ) ;%3D mean matrix of ggim
%     alpha=0.4;NNim2=padarray(Nim2,[5, 5, 5],'symmetric');
%     PRINLMfim2 = RI_NLM(cgim, NNim2, Noise_map, miu,alpha);%~5min
%     [dpim , noise_mapr]=NLPCA(single(PRINLMfim2),x1(1),x1(2),x1(3),x1(4),x1(5), 1);
%     yichang=find(isnan(dpim));
%     dpim(yichang)=PRINLMfim2(yichang);
%     DP=padarray(dpim,[5, 5, 5],'symmetric');
%     dmiu=imboxfilt3( DP  , 'padding' , 'symmetric' ) ;%3D mean matrix of PRINLMfim2
%     pdpim=RI_NLM(DP, NNim2, Noise_map, dmiu,alpha);%~5min
%     index = find(im1>0);
%     psnrcpd(i)= 20*log10(255/sqrt(mean((im2(index)-dpim(index)).^2)));
%     ker=[1 1 1];
%     ssimcpd(i)=ssim_index3d(PRINLMfim2,im2,ker,index);
%     psnrcpdp(i)= 20*log10(255/sqrt(mean((im2(index)-pdpim(index)).^2)));
%     ssimcpdp(i)=ssim_index3d(PRINLMfim22,im2,ker,index);    
% end
%%
%precomputation
% psnrtest=zeros(1,411);
% for i=100:510
%     sg=i/100;
%     noise3D=sg*ones(imsize(1:3)); 
%     Nim1=ricernd(im1, noise3D);
%     psnrtest(i-99)= 20*log10(255/sqrt(mean((im1(index)-Nim1(index)).^2)));
% end
% psnrtest=[psnrtest, zeros(1,100)];
% for i=511:610
%     sg=i/100;
%     noise3D=sg*ones(imsize(1:3)); 
%     Nim1=ricernd(im1, noise3D);
%     psnrtest(i-99)= 20*log10(255/sqrt(mean((im1(index)-Nim1(index)).^2)));
% end
%%
%scheme 1 based on more impact NLM
%First check the accuracy of the scheme
% an absolutely powerful method to obtain local noise estimation given im.
%when intensity of noise relates and only relates to signal intensity
% noise3D=noisesm(1)*ones(imsize(1:3)); 
% Nim1=ricernd(im1, noise3D);
% SIGr=zeros(imsize(1:3)); 
% for j=0:255
%     vtemp=find(im1==j);
%     %Nim1=ricernd(im1, noise3D);
%     ntemp=Nim1(vtemp);
%     %ntemp= PRINLMfim2(vtemp);
%     sigma2=mean(ntemp.^2 )-(mean(ntemp ))^2;
%     sigma=sqrt(sigma2);
%     SIGr(vtemp)= sigma;
% end
% SIGg=SIGr./sqrt(kexi( im1/2.55 )) ;
%fully correct !!

% implementation 
%f2 becomes weaker than f at high PSNR, we need to solve it
% psnr_dncnn=[46.32, 40.47, 37.82, 36.20, 34.71];
% RMSE=255./10.^(psnr_dncnn/20);
% ndncnn=[1.24 2.40 3.23 3.87 4.56]; %derive from precomputation - Rician
% psnrcthxian=zeros(1,5);ssimcthxian=zeros(1,5);
% for i=1:5
%     noise3D=noisesm(i)*ones(imsize(1:3)); 
%     Nim1=ricernd(im1, noise3D);
%     Noise_map= noisesm(i)*ones( 191,227,191);
%     %tgim=padarray(im1,[5, 5, 5],'symmetric');
%     cim=normrnd(im1, RMSE(i)*ones(imsize(1:3)));%construct a pre-filtered image
%     cgim=padarray(cim ,[5, 5, 5],'symmetric');
%     miu=imboxfilt3( cgim  , 'padding' , 'symmetric' ) ;%3D mean matrix of ggim
%     alpha=0.4;NNim1=padarray(Nim1,[5, 5, 5],'symmetric');
%     PRINLMfim = RI_NLM(cgim, NNim1, Noise_map, miu,alpha);%~5min
%     Ptemp = padarray(PRINLMfim,[5 5 5],'symmetric');
%     mu=imboxfilt3( Ptemp  , 'padding' , 'symmetric' );
%     nosemap=zeros(imsize(1:3)); 
%     %rounding off
%     Pfim = round( PRINLMfim);
%     maxPfim = max(Pfim,[],'all');
%     for j=0:maxPfim
%         vtemp=find(Pfim==j);
%         ntemp=cgim(vtemp);
% %         lower=mean( ntemp )-2*std(ntemp ,1) ;
% %         upper=mean( ntemp )+2*std(ntemp ,1) ;
% %         ntemp=ntemp( (ntemp<upper)&(ntemp>lower));
%         sigma2=mean(ntemp.^2 )-(mean(ntemp ))^2;
%         %More sophisticated algorithms are needed internally
%         sigma=sqrt(sigma2);
%         nosemap(vtemp)= sigma;
%     end
%     %max nosemap=67+
%     %cfai =  (0.9953*logamma-1.716)./(logamma-1.787) ; idx = (logamma<1.913); cfai(idx) = 0;
%     nosemap=padarray(nosemap,[5 5 5],'symmetric');
%     Fim = RI_NLMxian(Ptemp, cgim, nosemap, mu, alpha );
%     psnrcthxian(i)= 20*log10(255/sqrt(mean((im1(index)-Fim(index)).^2)));
%     ssimcthxian(i)=ssim_index3d(Fim,im1,ker,index);        
% end





%%
%The image itself was used to guide logamma
%spatial locality
% psnryc1=zeros(1,5);ssimyc1=zeros(1,5);
% pre=zeros(3,3,3);
% pre(:,:,1)=[3 2 3;2 1 2;3 2 3] ;pre(:,:,2)=[2 1 2;1 0 1;2 1 2] ;pre(:,:,3)= [3 2 3;2 1 2;3 2 3] ;
% 
% for i=1:5
%     Nim1=ricernd(im1, noiselr(i)*ones(imsize(1:3)));
%     %gaussian kernel
%     window=exp( -pre.^2/(2*noiselr(i)^2 ));
%     window=window/sum(window(:));
%     u1im  = convn( Nim1,window, 'same');
%     u2im = convn( Nim1.*Nim1,window, 'same') ;
%     
%     %u2im= imboxfilt3(Nim1.^2 , 3 , 'padding' , 'symmetric');
%     %u1im=imboxfilt3(Nim1, 3 , 'padding' , 'symmetric');
%     sigma2=u2im-u1im.^2;
%     sigma = sqrt(sigma2);
%     logamma=u1im./sigma;
%     logamma=reshape(logamma,1,181*217*181);
%     idxmin = (logamma<min(effSNR)); logamma(idxmin) = min(effSNR);%avoid outliers
%     idxmax=(logamma>max(effSNR));
%     theta = interp1(effSNR, ptheta, logamma, 'linear');theta(idxmax)=sqrt((logamma(idxmax)).^2-1);
%     theta=reshape(theta,181,217,181);
%     sigmag=sqrt(u2im./(2+theta.^2));
%     %         temp = interp1(fai, etta, theta, 'linear');%arrayfun(@ettainv, eetta);
%     %         temp=reshape(temp,181,217,181);
%     %         sigmag= imboxfilt3(Nim , 3 , 'padding' , 'symmetric')./temp;
%     %         theta=reshape(theta,181,217,181);
%     vim=theta.*sigmag;
%     sigmag=imboxfilt3(sigmag , 15 , 'padding' , 'symmetric');
%     Cov = std( sigmag , 1 , 'all' )/mean( sigmag , 'all');
%     if Cov <0.15
%         noise_map = mean( sigmag , 'all')*ones(181,217,181);
%     else
%         noise_map = sigmag;
%     end
%     NNim = padarray(Nim1,[5, 5, 5],'symmetric');
%     GGim = padarray(vim,[5, 5, 5],'symmetric');
%     Noise_map=padarray(noise_map,[5, 5, 5],'symmetric');
%     miu = imboxfilt3( GGim  , 'padding' , 'symmetric' ) ;%3D mean matrix of ggim
%     alpha=0.4;
%     PRINLMfim = RI_NLM(GGim, NNim, Noise_map, miu,alpha);%~5min
%     index = find(im1>0);
%     psnryc1(i) = 20*log10(255/sqrt(mean((im1(index)-PRINLMfim(index)).^2)));
%     ssimyc1(i) = ssim_index3d(im1,PRINLMfim,[1 1 1],index);
% end
%%
%get CIFTI files into MATLAB
%filename = 'D:\MATLAB\biyelunwen\T1w_acpc_dc_restore.nii ';
%wbcommand = 'D:\scientific_research\biyesheji\software\workbench-windows64-v1.5.0\workbench\bin_windows64\wb_command.exe ';
%cii = ciftiopen(filename,wbcommand);
%CIFITdata=cii.cdata;

%V=niftiread('T1w_acpc_dc.nii');
%%
%real clinical data
nii=load_nii('T1w_acpc_dc_restore.nii.gz');
imch=nii.img;
%reslice_nii('IXI002-Guys-0828-PD.nii.gz', '0828-PD.nii.gz');
%nii=load_nii('0828-PD.nii.gz');
%imci=nii.img;
nii=load_nii('OAS1_0001_MR1_mpr-1_anon.img');
imco=nii.img;
imco=single(imco);
%OASIS
tic
sigmago=RicianSTD(imco);
toc

noise_mapo=sigmago*ones(size(imco));
tic
[dnim,~]=NLPCA(imco,2.12,1.16,2.46);
toc
yichang=find(isnan(dnim));
dnim(yichang)=imco(yichang);
A1=size(imco,1);B1=size(imco,2);C1=size(imco,3);
        eetta = dnim ./ noise_mapo;
        eetta=reshape(eetta,1,A1*B1*C1);
        idxe = (eetta<min(etta)); eetta(idxe) = min(etta);%avoid outliers
        idxd = (eetta>max(etta)); eetta(idxd) = max(etta);%avoid outliers
        temp = interp1(etta, fai, eetta, 'linear');%arrayfun(@ettainv, eetta);
        temp=reshape(temp,A1,B1,C1);
        gimo = noise_mapo.* temp;
        NNim = padarray(imco,[5, 5, 5],'symmetric');
        GGim = padarray(gimo,[5, 5, 5],'symmetric');
        Noise_mapo=sigmago*ones((A1+10),(B1+10),(C1+10));%padarray(noise_map,[5, 5, 5],'symmetric');
        miu = imboxfilt3( GGim  , 'padding' , 'symmetric' ) ;%3D mean matrix of ggim
        tic
        PRINLMfimo = RI_NLM(GGim, NNim, Noise_mapo, miu,alpha);%~5min
        toc
        tic
        [dpimo,~]=NLPCA(PRINLMfimo,2.12,1.16,2.46);
        toc
        yichang=find(isnan(dpimo));
dpimo(yichang)=PRINLMfimo(yichang);
residual1=imco-dpimo;
residualg1=imco-gimo;

%HCP T1w
tic
sigmagh=RicianSTD(imch);
toc

noise_maph=sigmagh*ones(size(imch));
tic
[dnim,~]=NLPCA(imch,2.12,1.16,2.46);
toc
yichang=find(isnan(dnim));
dnim(yichang)=imch(yichang);
A2=size(imch,1);B2=size(imch,2);C2=size(imch,3);
        eetta = dnim ./ noise_maph;
        eetta=reshape(eetta,1,A2*B2*C2);
        idxe = (eetta<min(etta)); eetta(idxe) = min(etta);%avoid outliers
        temp = interp1(etta, fai, eetta, 'linear');%arrayfun(@ettainv, eetta);
        temp=reshape(temp,A2,B2,C2);
        gimh = noise_maph.* temp;
        NNim = padarray(imch,[5, 5, 5],'symmetric');
        GGim = padarray(gimh,[5, 5, 5],'symmetric');
        Noise_maph=sigmagh*ones((A2+10),(B2+10),(C2+10));%padarray(noise_map,[5, 5, 5],'symmetric');
        miu = imboxfilt3( GGim  , 'padding' , 'symmetric' ) ;%3D mean matrix of ggim
        tic
        PRINLMfimh = RI_NLM(GGim, NNim, Noise_maph, miu,alpha);%~5min
        toc
        tic
        [dpimh,~]=NLPCA(PRINLMfimh,2.12,1.16,2.46);
        toc
yichang=find(isnan(dpimh));
dpimh(yichang)=PRINLMfimh(yichang);
residual2=imch-dpimh;
residualg2=imch-gimh;
%%
%drawing
colormap(gray);
ni1=round(size(imco,1)/2); nh1=round(size(imch,1)/2);
ni2=round(size(imco,2)/2); nh2=round(size(imch,2)/2);
ni3=round(size(imco,3)/2); nh3=round(size(imch,3)/2);

%traverse OASIS dgpd
subplot(4,5,7),imagesc(imrotate(reshape(imco(:,ni2,:),256,128),180));set(gca,'xtick',[],'ytick',[]);
subplot(4,5,8),imagesc(imrotate(reshape(dpimo(:,ni2,:),256,128),180));set(gca,'xtick',[],'ytick',[]);
subplot(4,5,9),imagesc(imrotate(reshape(residual1(:,ni2,:),256,128),180));set(gca,'xtick',[],'ytick',[]);
%coronal OASIS dgpd
subplot(4,5,10),imagesc(imrotate(reshape(imco(ni1,:,:),256,128),180));set(gca,'xtick',[],'ytick',[]);
subplot(4,5,11),imagesc(imrotate(reshape(dpimo(ni1,:,:),256,128),180));set(gca,'xtick',[],'ytick',[]);
subplot(4,5,12),imagesc(imrotate(reshape(residual1(ni1,:,:),256,128),180));set(gca,'xtick',[],'ytick',[]);
%sagittal OASIS dgpd
subplot(4,5,13),imagesc(imrotate(imco(:,:,ni3),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,5,14),imagesc(imrotate(dpimo(:,:,ni3),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,5,15),imagesc(imrotate(residual1(:,:,ni3),90));set(gca,'xtick',[],'ytick',[]);
%sagittal HCP dgpd
subplot(4,5,16),imagesc(imrotate(reshape(imch(nh1,:,:),311,260),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,5,17),imagesc(imrotate(reshape(dpimh(nh1,:,:),311,260),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,5,18),imagesc(imrotate(reshape(residual2(nh1,:,:),311,260),90));set(gca,'xtick',[],'ytick',[]);

%sagittal HCP dg
%subplot(3,6,10),imagesc(imrotate(reshape(imch(nh1,:,:),311,260),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,5,19),imagesc(imrotate(reshape(gimh(nh1,:,:),311,260),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,5,20),imagesc(imrotate(reshape(residualg2(nh1,:,:),311,260),90));set(gca,'xtick',[],'ytick',[]);
%(:,:,nh3),260,311),90
%residual3=(gimo-dpimo);
%subplot(3,6,13),imagesc(imrotate(reshape(residual3(:,ni2,:),256,128),180));set(gca,'xtick',[],'ytick',[]);
% %traverse OASIS
% subplot(3,6,7),imagesc(imrotate(reshape(imco(:,ni2,:),256,128),180));set(gca,'xtick',[],'ytick',[]);
% subplot(3,6,8),imagesc(imrotate(reshape(dpimo(:,ni2,:),256,128),180));set(gca,'xtick',[],'ytick',[]);
% subplot(3,6,9),imagesc(imrotate(reshape(residual1(:,ni2,:),256,128),180));set(gca,'xtick',[],'ytick',[]);
% %coronal HCP
% subplot(3,6,10),imagesc(imrotate(reshape(imch(:,nh2,:),260,260),90));set(gca,'xtick',[],'ytick',[]);
% subplot(3,6,11),imagesc(imrotate(reshape(dpimh(:,nh2,:),260,260),90));set(gca,'xtick',[],'ytick',[]);
% subplot(3,6,12),imagesc(imrotate(reshape(residual2(:,nh2,:),260,260),90));set(gca,'xtick',[],'ytick',[]);
% %coronal OASIS
% subplot(3,6,13),imagesc(imrotate(reshape(imco(:,:,ni3),256,256),90));set(gca,'xtick',[],'ytick',[]);
% subplot(3,6,14),imagesc(imrotate(reshape(dpimo(:,:,ni3),256,256),90));set(gca,'xtick',[],'ytick',[]);
% subplot(3,6,15),imagesc(imrotate(reshape(residual1(:,:,ni3),256,256),90));set(gca,'xtick',[],'ytick',[]);
% %traverse HCP
% subplot(3,6,16),imagesc(imrotate(reshape(imch(:,:,nh3),260,311),90));set(gca,'xtick',[],'ytick',[]);
% subplot(3,6,17),imagesc(imrotate(reshape(dpimh(:,:,nh3),260,311),90));set(gca,'xtick',[],'ytick',[]);
% subplot(3,6,18),imagesc(imrotate(reshape(residual2(:,:,nh3),260,311),90));set(gca,'xtick',[],'ytick',[]);

%subplot(2,3,3),imagesc(imrotate(imch(87:173,104:207,n),90));set(gca,'xtick',[],'ytick',[]);
%subplot(2,3,4),imagesc(imrotate(imci(85:171,91:182,n),90));set(gca,'xtick',[],'ytick',[]);
%subplot(2,3,5),imagesc(imrotate(PRINLMfimh(87:173,104:207,n),90));set(gca,'xtick',[],'ytick',[]);
%subplot(2,3,6),imagesc(imrotate(PRINLMfimi(85:171,91:182,n),90));set(gca,'xtick',[],'ytick',[]);


