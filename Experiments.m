%%

noisesm=255*0.01*[1 3 5 7 9];

%x1=[2.12 1.16 2.46];
alpha=0.4;
imsizes=[181,217,10];
imsize=[181,217,181];

ker=[1 1 1];
%x1=[4 64 4 2.2 1.16];
%%
%precomputation

global etta
global fai
load('precomputation.mat')
%%
%load in data
% read the file

file1 ='pd_icbm_normal_1mm_pn0_rf0.rawb';
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

% data need to be convert from single to double
nii=load_nii('T1w_acpc_dc_restore.nii.gz');
imh=nii.img;
imh=double(imh);
ma3=max(imh,[],'all');
nii=load_nii('OAS1_0001_MR1_mpr-1_anon.img');
imo=nii.img;
imo= double(imo);

% reslice_nii('IXI132-HH-1415-T2.nii.gz', 'HH-1415-T2.nii.gz');
% nii=load_nii('HH-1415-T2.nii.gz');
% imix=nii.img;
% maix=max(imix,[],'all');

noimch=RicianSTD(imh);

noimo = RicianSTD(imo);

%noimix = RicianSTD(imix);


%%
%timing algorithm 9%
Nimbw=ricernd(im2, 2.55*9*ones(imsize(1:3)));
Nimh=imh;
% Nimbw=ricernd(im1, 2.55*9*ones(181,217,6));
% Nimh=ricernd(imh, 40.96*9*ones(260,311,6));
Nimo=imo;

%NL-PCA


tb(1,1)=datetime;
Dnimpbw=NLPCAnp( (Nimbw),2.55*9,2.2);
tb(2,1)=datetime;
th(1,1)=datetime;
Dnimph=NLPCAnp( (Nimh),noimch,2.2);
th(2,1)=datetime;
to(1,1)=datetime;
Dnimpo=NLPCAnp( (Nimo),noimo,2.2);
to(2,1)=datetime;



%BM4D
tb(1,2)=datetime;
[~,Dnimbbw,~]=bm4d(1, Nimbw/255, 0.09, 1,1,1,0,0);%0,1,1,0,0
tb(2,2)=datetime;
th(1,2)=datetime;
[~,Dnimbh,~]=bm4d(1, Nimh/1024, noimch/1024, 1,1,1,0,0);Dnimbh=Dnimbh*1024;
th(2,2)=datetime;
to(1,2)=datetime;
[~,Dnimbo,~]=bm4d(1, Nimo/1024, noimo/1024, 1,1,1,0,0);Dnimbo=Dnimbo*1024;
to(2,2)=datetime;
% %BM4D noise estimation only care for time, not accuracy
% [~,~,noim]=bm4d(1, Nimbw/255, 0, 1,1,1,0,0);%0,1,1,0,0
% toc;
% tic;
% [~,~,noim]=bm4d(1, Nimh/4096, 0, 1,1,1,0,0);
% toc;
% tic;
% [~,~,noim]=bm4d(1, Nimo/4096, 0, 1,1,1,0,0);
% toc;

%PRIsl-NL-PCA subscript
tb(1,3)=datetime;
%Ppim = RI_NLM(RiC(Dnimpbw,2.55*9), Nimbw, 2.55*9, alpha);
Ppimbw = cPRI_NL_PCA(Nimbw,9,1,2.55*9*ones(181,217,181),RiC(Dnimpbw,2.55*9),1); 
tb(2,3)=datetime;
th(1,3)=datetime;
%Ppim = RI_NLM(RiC(Dnimph,40.96*9), Nimh, 40.96*9, alpha);
Ppimh = cPRI_NL_PCA(Nimh,9,1,noimch*ones(260,311,260),RiC(Dnimph,noimch),1); 
th(2,3)=datetime;
to(1,3)=datetime;
Ppim = RI_NLM(RiC(Dnimpo,noimo), Nimo, noimo, alpha);
Ppimo = cPRI_NL_PCA(Nimo,9,1,noimo*ones(256,256,128),RiC(Dnimpo,noimo),1); 
to(2,3)=datetime;
%PRIvg-NL-PCA

tb(1,4)=datetime;
ppim = RINLMmy(RiC(Dnimpbw,2.55*9),Nimbw,2.55*9,alpha,8);

tb(2,4)=datetime;
th(1,4)=datetime;

rt=1024/255;
        nnimh=(Nimh/rt);dnnim=RiC(Dnimph,noimch)/rt;
        nnimh(nnimh>255)=255;dnnim(dnnim>255)=255;
        ppim = RINLMmy(dnnim, nnimh, noimch/rt,alpha,8);

% ppimh = RINLMmy(RiC(Dnimph,10.24*9)/rt,Nimh/rt,2.55*9,alpha,8);
% ppimh=ppimh*1024/255;
th(2,4)=datetime;

to(1,4)=datetime;
rt=1024/255;
        nnimo=(Nimo/rt);dnnim=RiC(Dnimpo,noimo)/rt;
        nnimo(nnimo>255)=255;dnnim(dnnim>255)=255;
ppim= RINLMmy(dnnim,nnimo,noimo/rt,alpha,8);
to(2,4)=datetime;


%BM4Dw-NL-PCA

tb(1,5)=datetime;
bwpim=bwp(1,Nimbw/255,0.09,1,1,1);

tb(2,5)=datetime;
th(1,5)=datetime;
bwpim=bwp(1,Nimh/1024,noimch/1024,1,1,1);
th(2,5)=datetime;
to(1,5)=datetime;
bwpim=bwp(1,Nimo/1024,noimo/1024,1,1,1);
to(2,5)=datetime;

%PRI-bh
tb(1,6)=datetime;
[~,bhimbw,~]=bm4d(1, Nimbw/255, 0.09, 1,1,0,0,0);bhimbw=bhimbw*255;
        Pbhim=cPRI_NL_PCA(Nimbw,9,1,2.55*9*ones(181,217,181),bhimbw,1);

tb(2,6)=datetime;
th(1,6)=datetime;
[~,bhimh,~]=bm4d(1, Nimh/1024, noimch/1024, 1,1,0,0,0);bhimh=bhimh*1024;
        Pbhim=cPRI_NL_PCA(Nimh,9,1,noimch*ones(260,311,260),bhimh,1);
th(2,6)=datetime;
to(1,6)=datetime;
[~,bhimo,~]=bm4d(1, Nimo/1024, noimo/1024, 1,1,0,0,0);bhimo=bhimo*1024;
        Pbhim=cPRI_NL_PCA(Nimo,9,1,noimo*ones(256,256,128),bhimo,1);
to(2,6)=datetime;




%pri-bh

tb(1,7)=datetime;




pbhim=RINLMmy(bhimbw,Nimbw,2.55*9,0.4,8);

tb(2,7)=datetime;
th(1,7)=datetime;

rt=1024/255;
        nnimh=(Nimh/rt);dnnim=RiC(bhimh,noimch)/rt;
        nnimh(nnimh>255)=255;dnnim(dnnim>255)=255;
pbhim=RINLMmy(dnnim,nnimh,noimch/rt,0.4,8);



th(2,7)=datetime;
to(1,7)=datetime;

rt=1024/255;
        nnimo=(Nimo/rt);dnnim=RiC(bhimo,noimo)/rt;
        nnimo(nnimo>255)=255;dnnim(dnnim>255)=255;
ppimo = RINLMmy(dnnim,nnimo,noimo/rt,alpha,8);


to(2,7)=datetime;







%%
%truncate for faster computation
% im1=truncateslice(im1,3);
% im2=truncateslice(im2,3);
% imh=truncateslice(imh,3);
im1=truncateslice(im1,5);
im2=truncateslice(im2,5);
imh=truncateslice(imh,5);
%imix=truncateslice(imix,5);
%imo=truncateslice(imo,3);
index1=find(im1>0);
index2=find(im2>0);
index3=find(imh>3*noimch);
%indexix=find(imix>3*noimix);









%%
%median filter NLPCA
PSNRmed=zeros(1,5);ssimed=zeros(1,5);
PSNRnm=zeros(1,5);ssimnm=zeros(1,5);
for i=1:5
        nimbw=ricernd(im1, noisesm(i)*ones(imsizes(1:3)));
        nnim = medfilt3(nimbw);
        dnim1=NLPCAnp(nimbw,noisesm(i),2.2);dnim1=RiC(dnim1,noisesm(i));
        dnim2=NLPCAnp(nnim,noisesm(i),2.2);dnim2=RiC(dnim2,noisesm(i));
        PSNRmed(i)=20*log10(255/sqrt(mean((im1(index1)-dnim2(index1)).^2)));
        ssimed(i)=ssim_index3d( dnim2 , im1 , ker , index1 );
        PSNRnm(i)=20*log10(255/sqrt(mean((im1(index1)-dnim1(index1)).^2)));
        ssimnm(i)=ssim_index3d( dnim1 , im1 , ker , index1 );
%         if i==5
%             xuyao1=PRINLMfimm;
%             xuyao2=PRINLMfimnm;
%         end
end









%%
%Horizontal comparison im1
psnrn=zeros(1,5);ssimn=zeros(1,5);
psnrp=zeros(1,5);ssimp=zeros(1,5);
psnrb=zeros(1,5);ssimb=zeros(1,5);
psnrPp=zeros(1,5);ssimPp=zeros(1,5);
psnrpp=zeros(1,5);ssimpp=zeros(1,5);
psnrbwp=zeros(1,5);ssimbwp=zeros(1,5);
psnrPbh=zeros(1,5);ssimPbh=zeros(1,5);
psnrpbh=zeros(1,5);ssimpbh=zeros(1,5);

for i=1:5
        nim1=round(ricernd(im1, noisesm(i)*ones(imsizes(1:3))));
        psnrn(i) = 20*log10(255/sqrt(mean((im1(index1)-nim1(index1)).^2)));
        ssimn(i)=ssim_index3d( nim1 , im1 , ker , index1 );

        %[dnim,~]=NLPCA( (nim1),2.2,1.29,2);
        dnim = NLPCAnp(nim1,noisesm(i),2.2);
        dnim = RiC(dnim,noisesm(i));
        
        psnrp(i) = 20*log10(255/sqrt(mean((im1(index1)-round(dnim(index1))).^2)));
        ssimp(i)=ssim_index3d( round(dnim) , im1 , ker , index1 );
        
        [~,bim,~] = bm4d(1, nim1/255, noisesm(i)/255, 1,1,1,0,0);bim=bim*255;
        psnrb(i) = 20*log10(255/sqrt(mean((im1(index1)-round(bim(index1))).^2)));
        ssimb(i)=ssim_index3d( round(bim) , im1 , ker , index1 );

        Ppim = cPRI_NL_PCA(nim1,9,1,noisesm(i)*ones(181,217,10),dnim, 1);
        ppim = RINLMmy(dnim, nim1, noisesm(i),0.4,8);
        psnrPp(i) = 20*log10(255/sqrt(mean((im1(index1)-round(Ppim(index1))).^2)));
        ssimPp(i)=ssim_index3d( round(Ppim) , im1 , ker , index1 );
        psnrpp(i) = 20*log10(255/sqrt(mean((im1(index1)-round(ppim(index1))).^2)));
        ssimpp(i)=ssim_index3d( round(ppim) , im1 , ker , index1 );

        % bwpim=bm4dw(nim1,dnim,noisesm(i));
        % bwpim=RiC(bwpim,noisesm(i));
        bwpim = bwp(1,nim1/255,noisesm(i)/255,1,1,1);bwpim=bwpim*255;
        psnrbwp(i) = 20*log10(255/sqrt(mean((im1(index1)-round(bwpim(index1))).^2)));
        ssimbwp(i)=ssim_index3d( round(bwpim) , im1 , ker , index1 );

        % bhim=bm4dh(nim1,noisesm(i));
        % bhim=RiC(bhim,noisesm(i));
        [~,bhim,~]=bm4d(1, nim1/255, noisesm(i)/255, 1,1,0,0,0);bhim=bhim*255;
        Pbhim=cPRI_NL_PCA(nim1,9,1,noisesm(i)*ones(181,217,10),bhim,1);
        pbhim=RINLMmy(bhim,nim1,noisesm(i),0.4,8);
        psnrPbh(i) = 20*log10(255/sqrt(mean((im1(index1)-round(Pbhim(index1))).^2)));
        ssimPbh(i)=ssim_index3d( round(Pbhim) , im1 , ker , index1 );
        psnrpbh(i) = 20*log10(255/sqrt(mean((im1(index1)-round(pbhim(index1))).^2)));
        ssimpbh(i)=ssim_index3d( round(pbhim) , im1 , ker , index1 );


end
% residualPp=im1-Ppim;
% residualpp=im1-ppim;
% residualb=im1-bim;
% colormap(gray);

% subplot(4,8,1),imagesc(imrotate(im1(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,8,2),imagesc(imrotate(nim1(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,8,3),imagesc(imrotate(Ppim(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,8,4),imagesc(imrotate(ppim(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,8,5),imagesc(imrotate(bim(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,8,6),imagesc(imrotate(residualPp(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,8,7),imagesc(imrotate(residualpp(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,8,8),imagesc(imrotate(residualb(:,:,3),90));set(gca,'xtick',[],'ytick',[]);


%%
%Horizontal comparison im2
psnrn2=zeros(1,5);ssimn2=zeros(1,5);
psnrp2=zeros(1,5);ssimp2=zeros(1,5);
psnrb2=zeros(1,5);ssimb2=zeros(1,5);
psnrPp2=zeros(1,5);ssimPp2=zeros(1,5);
psnrpp2=zeros(1,5);ssimpp2=zeros(1,5);
psnrbwp2=zeros(1,5);ssimbwp2=zeros(1,5);
psnrPbh2=zeros(1,5);ssimPbh2=zeros(1,5);
psnrpbh2=zeros(1,5);ssimpbh2=zeros(1,5);
for i=1:5
%for i=4
        nim2=round(ricernd(im2, noisesm(i)*ones(imsizes(1:3))));
        psnrn2(i) = 20*log10(255/sqrt(mean((im2(index2)-nim2(index2)).^2)));
        ssimn2(i)=ssim_index3d( nim2 , im2 , ker , index2 );

        dnim=NLPCAnp( (nim2),noisesm(i),2.2);
        dnim = RiC(dnim,noisesm(i));
        
        psnrp2(i) = 20*log10(255/sqrt(mean((im2(index2)-round(dnim(index2))).^2)));
        ssimp2(i)=ssim_index3d( round(dnim) , im2 , ker , index2 );
        
        [~,bim,~] = bm4d(1, nim2/255, noisesm(i)/255, 1,1,1,0,0);bim=bim*255;
        psnrb2(i) = 20*log10(255/sqrt(mean((im2(index2)-round(bim(index2))).^2)));
        ssimb2(i)=ssim_index3d( round(bim) , im2 , ker , index2 );

        Ppim = cPRI_NL_PCA( nim2,9,1, noisesm(i)*ones(181,217,10),dnim,1);
        ppim = RINLMmy(dnim, nim2, noisesm(i),alpha,8);
        psnrPp2(i) = 20*log10(255/sqrt(mean((im2(index2)-round(Ppim(index2))).^2)));
        ssimPp2(i)=ssim_index3d( round(Ppim) , im2 , ker , index2 );
        psnrpp2(i) = 20*log10(255/sqrt(mean((im2(index2)-round(ppim(index2))).^2)));
        ssimpp2(i)=ssim_index3d( round(ppim) , im2 , ker , index2 );

        % bwpim=bm4dw(nim2,dnim,noisesm(i));
        bwpim = bwp(1,nim2/255,noisesm(i)/255,1,1,1);bwpim=bwpim*255;
        psnrbwp2(i) = 20*log10(255/sqrt(mean((im2(index2)-round(bwpim(index2))).^2)));
        ssimbwp2(i)=ssim_index3d( round(bwpim) , im2 , ker , index2 );

        % bhim=bm4dh(nim2,noisesm(i));
        % bhim=RiC(bhim,noisesm(i));

        [~,bhim,~]=bm4d(1, nim2/255, noisesm(i)/255, 1,1,0,0,0);bhim=bhim*255;
        Pbhim=cPRI_NL_PCA(nim2,9,1,noisesm(i)*ones(181,217,10),bhim,1);
        pbhim=RINLMmy(bhim,nim2,noisesm(i),alpha,8);
        psnrPbh2(i) = 20*log10(255/sqrt(mean((im2(index2)-round(Pbhim(index2))).^2)));
        ssimPbh2(i)=ssim_index3d( round(Pbhim) , im2 , ker , index2 );
        psnrpbh2(i) = 20*log10(255/sqrt(mean((im2(index2)-round(pbhim(index2))).^2)));
        ssimpbh2(i)=ssim_index3d( round(pbhim) , im2 , ker , index2 );


end
residualPp=im2-Ppim;
residualpp=im2-ppim;
residualb=im2-bim;
residualbwp=im2-bwpim;
colormap(gray);

subplot(4,10,1),imagesc(imrotate(im2(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,2),imagesc(imrotate(nim2(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,3),imagesc(imrotate(Ppim(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,4),imagesc(imrotate(ppim(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,5),imagesc(imrotate(bim(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,6),imagesc(imrotate(bwpim(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,7),imagesc(imrotate(residualPp(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,8),imagesc(imrotate(residualpp(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,9),imagesc(imrotate(residualb(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,10),imagesc(imrotate(residualbwp(:,:,3),90));set(gca,'xtick',[],'ytick',[]);


%%
% %Partial Enlarged View, Center Right, bottom center
% colormap(gray);
% n=round(size(im1,3)/2);
% n1=round(size(im1,1)/3);nn1=2*n1;
% n2=round(size(im1,2)/3);nn2=2*n2;
% %line
% %subplot(4,5,1),imagesc(imrotate(im1(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% %subplot(4,5,2),imagesc(imrotate(Nim1(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% %subplot(4,5,3),imagesc(imrotate(PRINLMfim1(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% %subplot(4,5,4),imagesc(imrotate(PRINLMfim11(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% %subplot(4,5,5),imagesc(imrotate(dpim1(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% %subplot(4,5,6),imagesc(imrotate(pdpim1(:,:,n),90));set(gca,'xtick',[],'ytick',[]);
% 
% 
% subplot(4,5,1),imagesc(imrotate(im1(nn1:end,n2:nn2,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,5,2),imagesc(imrotate(Nimbw(nn1:end,n2:nn2,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,5,3),imagesc(imrotate(PRINLMfim1(nn1:end,n2:nn2,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,5,4),imagesc(imrotate(PRINLMfim11(nn1:end,n2:nn2,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,5,5),imagesc(imrotate(dpim1(nn1:end,n2:nn2,n),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,5,6),imagesc(imrotate(pdpim1(nn1:end,n2:nn2,n),90));set(gca,'xtick',[],'ytick',[]);



%%
% %IXI 0.5% noise/threshold, where threshold=300 
% psnrnix=zeros(1,5);ssimnix=zeros(1,5);
% psnrpix=zeros(1,5);ssimpix=zeros(1,5);
% psnrbix=zeros(1,5);ssimbix=zeros(1,5);
% psnrPpix=zeros(1,5);ssimPpix=zeros(1,5);
% psnrppix=zeros(1,5);ssimppix=zeros(1,5);
% psnrbwpix=zeros(1,5);ssimbwpix=zeros(1,5);
% psnrPbhix=zeros(1,5);ssimPbhix=zeros(1,5);
% psnrpbhix=zeros(1,5);ssimpbhix=zeros(1,5);
% 
% ma4=400;
% noiselr=ma4*0.01*[1 3 5 7 9];
% Imix=round(imix/ma4*255);
% Imix(Imix>255)=255;
% rt=ma4/255;
% ma=max(imix,[],'all');
% % noiselr=ma*0.01*[1 3 5 7 9];
% % Imix=round(imix/ma*255);
% % rt=ma/255;

% for i=1:5
% 
%         nimix=round(ricernd(imix, noiselr(i)*ones(281,277,10)));
%         %round(ricernd(Imix, noisesm(i)*ones(260,311,10)));
%         %psnrn3(i) = 20*log10(4096/sqrt(mean((imix(index3)-nim3(index3)).^2)));
%         psnrnix(i) = 20*log10(255/sqrt(mean((Imix(indexix)-dis255(nimix(indexix)/rt)).^2)));
%         ssimnix(i)=ssim_index3d( dis255(nimix/rt) , Imix , ker , indexix );
% 
%         dnim=NLPCAnp( (nimix),noiselr(i),2.2);%dnnim=NLPCAnp( (nnimix),noisesm(i),2.2);
%         dnim = RiC(dnim,noiselr(i));%dnnim=RiC(dnnim,noisesm(i));
% 
%         %psnrpix(i) = 20*log10(4096/sqrt(mean((imix(indexix)-round(dnim(indexix))).^2)));
%         psnrpix(i) = 20*log10(255/sqrt(mean((Imix(indexix)-dis255(dnim(indexix)/rt)).^2)));
%         ssimpix(i)=ssim_index3d( dis255(dnim/rt) , Imix , ker , indexix );
% 
%         %[~,bim,~] = bm4d(1, nim3/4096, noiselr(i)/4096, 1,1,1,0,0);bim=bim*4096;
%         [~,bim,~] = bm4d(1, nimix/ma4, noiselr(i)/ma4, 1,1,1,0,0);bim=bim*ma4;
%         %psnrb3(i) = 20*log10(4096/sqrt(mean((imix(index3)-round(bim(index3))).^2)));
%         psnrbix(i) = 20*log10(255/sqrt(mean((Imix(indexix)-dis255(bim(indexix)/rt)).^2)));
%         ssimbix(i)=ssim_index3d( dis255(bim/rt) , Imix , ker , indexix );
% 
%         Ppim = cPRI_NL_PCA( nimix, 9,1,noiselr(i)*ones(281,277,10),dnim,1);
%         nnimix=(nimix/rt);dnnim=dnim/rt;nnimix(nnimix>255)=255;dnnim(dnnim>255)=255;
%         ppim = RINLMmy(dnnim, nnimix, noisesm(i),alpha,8);%ma3
%         %psnrPp3(i) = 20*log10(4096/sqrt(mean((imix(index3)-round(Ppim(index3))).^2)));
%         psnrPpix(i) = 20*log10(255/sqrt(mean((Imix(indexix)-dis255(Ppim(indexix)/rt)).^2)));
%         ssimPpix(i)=ssim_index3d( dis255(Ppim/rt) , Imix , ker , indexix );
%         %psnrppix(i) = 20*log10(4096/sqrt(mean((imix(indexix)-round(ppim(indexix))).^2)));
%         psnrppix(i)= 20*log10(255/sqrt(mean((Imix(indexix)-dis255(ppim(indexix))).^2)));
%         ssimppix(i)=ssim_index3d( dis255(ppim) , Imix , ker , indexix );
% 
%         % bwpim=bm4dw(nimix,dnim,noiselr(i));
%         %bwpim = bwp(1,nimix/4096,noiselr(i)/4096,1,1,1);bwpim=bwpim*4096;
%         %bwpim = bwp(1,nimix/ma4,noiselr(i)/ma4,1,1,1);bwpim=bwpim*ma4;
%         %psnrbwpix(i) = 20*log10(4096/sqrt(mean((imix(indexix)-round(bwpim(indexix))).^2)));
%         %psnrbwpix(i) = 20*log10(255/sqrt(mean((Imix(indexix)-dis255(bwpim(indexix)/rt)).^2)));
%         %ssimbwpix(i)=ssim_index3d( dis255(bwpim/rt) , Imix , ker , indexix );
% 
%         % bhim=bm4dh(nimix,noiselr(i));
%         % bhim=RiC(bhim,noiselr(i));
%         % %[~,bhim,~]=bm4d(1, nimix/4096, noiselr(i)/4096, 1,1,0,0,0);bhim=bhim*4096;
%         % 
%         % % [~,bhim,~]=bm4d(1, nimix/ma3, noiselr(i)/ma3, 1,1,0,0,0);bhim=bhim*ma3;
%         % Pbhim=cPRI_NL_PCA(nimix,9,1,noiselr(i)*ones(260,311,10),bhim,1);
%         % bhhim=bhim/rt;bhhim(bhhim>255)=255;%nnim3=(nim3/rt);
%         % pbhim=RINLMmy(bhhim,nnimix,noisesm(i),alpha,8);
%         % 
%         % %psnrPbhix(i) = 20*log10(4096/sqrt(mean((imix(indexix)-round(Pbixim(indexix))).^2)));
%         % 
%         % psnrPbhix(i) = 20*log10(255/sqrt(mean((Imix(indexix)-dis255(Pbhim(indexix)/rt)).^2)));
%         % ssimPbhix(i)=ssim_index3d( dis255(Pbhim/rt) , Imix , ker , indexix );
%         % 
%         % %psnrpbhix(i) = 20*log10(4096/sqrt(mean((imix(indexix)-round(pbhim(indexix))).^2)));
%         % psnrpbhix(i) = 20*log10(255/sqrt(mean((Imix(indexix)-dis255(pbhim(indexix)/rt)).^2)));
%         % ssimpbhix(i)=ssim_index3d( dis255(pbhim/rt) , Imix , ker , indexix );
% 
% 
% end
% %residualPp=Imix-dis255(Ppim/rt);
% residualpp=Imix-dis255(ppim);
% residualb=Imix-dis255(bim/rt);
% residualbwp=Imix-dis255(bwpim/rt);
% %Ppim=dis255(Ppim/rt);
% ppim=dis255(ppim);
% bim=dis255(bim/rt);
% bwpim=dis255(bwpim/rt);
% colormap(gray);
% 
% subplot(4,10,31),imagesc(imrotate(Imix(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,10,32),imagesc(imrotate(nnimix(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
% %subplot(4,10,33),imagesc(imrotate(Ppim(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,10,34),imagesc(imrotate(ppim(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,10,35),imagesc(imrotate(bim(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,10,36),imagesc(imrotate(bwpim(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
% %subplot(4,10,37),imagesc(imrotate(residualPp(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,10,38),imagesc(imrotate(residualpp(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,10,39),imagesc(imrotate(residualb(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
% subplot(4,10,40),imagesc(imrotate(residualbwp(:,:,3),90));set(gca,'xtick',[],'ytick',[]);

%%
%real clinical data
%HCP 1%
%real clinical data


%  Often there are only a few high-value pixels, so the whole image looks 
% very dark. We can improve this by letting all pixel values above a
% certain level be displayed at maximum brightness(McRobbie et al., 2017). 
% On the other hand it is also helpful to make the background noise as dark 
% as possible, and this is done by setting all pixel values below the noise 
% level to have minimum brightness. 
% However, remember that you are only changing the displayed voxel
% intensities, not the values in the underlying MR images.

        nim3=imh;


        dnim=NLPCAnp( (nim3),noimch,2.2);
        dnim = RiC(dnim,noimch);
        

        
        [~,bim,~] =bm4d(1, nim3/1024, noimch/1024, 1,1,1,0,0);bim=bim*1024;


        Ppim = cPRI_NL_PCA( nim3, 9,1,noimch*ones(260,311,260),dnim,1);
        rt=1024/255;
        nnim3=(nim3/rt);dnnim=dnim/rt;nnim3(nnim3>255)=255;dnnim(dnnim>255)=255;
        ppim = RINLMmy(dnnim, nnim3, noimch/rt,alpha,8);


        bwpim=bm4dw(nim3,dnim,noimch);


        % bhim=bm4dh(nim4,noimch);
        % bhim=RiC(bhim,noimch);
        % Pbhim=cPRI_NL_PCA(nim4,9,1,noimch*ones(256,256,128),bhim,1);
        % pbhim=RINLMmy(bhim,nim4,noimch,alpha,12,4095);



% residualPp=imch-Ppim;
% residualpp=imch-ppim;
% residualb=imch-bim;

imh=dis255(imh/rt);
residualPp=imh-dis255(Ppim/rt);
residualpp=imh-dis255(ppim);
residualb=imh-dis255(bim/rt);
residualbwp=imh-dis255(bwpim/rt);
Ppim=dis255(Ppim/rt);
ppim=dis255(ppim);
bim=dis255(bim/rt);
bwpim=dis255(bwpim/rt);


colormap(gray);

subplot(4,10,11),imagesc(imrotate(imh(:,:,130),90));set(gca,'xtick',[],'ytick',[]);
%subplot(4,10,22),imagesc(imrotate(nim2(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,13),imagesc(imrotate(Ppim(:,:,130),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,14),imagesc(imrotate(ppim(:,:,130),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,15),imagesc(imrotate(bim(:,:,130),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,16),imagesc(imrotate(bwpim(:,:,130),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,17),imagesc(imrotate(residualPp(:,:,130),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,18),imagesc(imrotate(residualpp(:,:,130),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,19),imagesc(imrotate(residualb(:,:,130),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,20),imagesc(imrotate(residualbwp(:,:,130),90));set(gca,'xtick',[],'ytick',[]);



%%
%real clinical data

%OASIS 9.44%
%  Often there are only a few high-value pixels, so the whole image looks 
% very dark. We can improve this by letting all pixel values above a
% certain level be displayed at maximum brightness(McRobbie et al., 2017). 
% On the other hand it is also helpful to make the background noise as dark 
% as possible, and this is done by setting all pixel values below the noise 
% level to have minimum brightness. 
% However, remember that you are only changing the displayed voxel
% intensities, not the values in the underlying MR images.


        nim4=imo;


        dnim=NLPCAnp( (nim4),noimo,2.2);
        dnim = RiC(dnim,noimo);
        

        
        %[~,bim,~] =bm4d(1, nim4/1024, noimo/1024, 1,1,1,0,0);bim=bim*1024;
        [~,bim,~] =bm4d(1, nim4/3072, noimo/3072, 1,1,1,0,0);bim=bim*3072;

        Ppim = cPRI_NL_PCA( nim4, 9,1,noimo*ones(256,256,128),dnim,1);
        rt=3*1024/255;
        nnim4=(nim4/rt);dnnim=dnim/rt;nnim4(nnim4>255)=255;dnnim(dnnim>255)=255;
        ppim = RINLMmy(dnnim, nnim4, noimo/rt,alpha,8);


        %bwpim=bm4dw(nim4,dnim,noimo);


        % bhim=bm4dh(nim4,noimo);
        % bhim=RiC(bhim,noimo);
        % Pbhim=cPRI_NL_PCA(nim4,9,1,noimo*ones(256,256,128),bhim,1);
        % pbhim=RINLMmy(bhim,nim4,noimo,alpha,12,4095);



% residualPp=imo-Ppim;
% residualpp=imo-ppim;
% residualb=imo-bim;

imo=dis255(imo/rt);
residualPp=imo-dis255(Ppim/rt);
residualpp=imo-dis255(ppim);
residualb=imo-dis255(bim/rt);
residualbwp=imo-dis255(bwpim/rt);
Ppim=dis255(Ppim/rt);
ppim=dis255(ppim);
bim=dis255(bim/rt);
bwpim=dis255(bwpim/rt);


colormap(gray);






subplot(4,10,21),imagesc(imrotate(imo(:,:,64),90));set(gca,'xtick',[],'ytick',[]);
%subplot(4,10,22),imagesc(imrotate(nim2(:,:,3),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,23),imagesc(imrotate(Ppim(:,:,64),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,24),imagesc(imrotate(ppim(:,:,64),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,25),imagesc(imrotate(bim(:,:,64),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,26),imagesc(imrotate(bwpim(:,:,64),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,27),imagesc(imrotate(residualPp(:,:,64),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,28),imagesc(imrotate(residualpp(:,:,64),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,29),imagesc(imrotate(residualb(:,:,64),90));set(gca,'xtick',[],'ytick',[]);
subplot(4,10,30),imagesc(imrotate(residualbwp(:,:,64),90));set(gca,'xtick',[],'ytick',[]);

