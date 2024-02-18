%optimization of RINLMmy
noisesm=255*0.01*[1 3 5 7 9];
noiselr=4096*0.01*[1 3 5 7 9];
%x1=[2.12 1.16 2.46];
alpha=0.4;
imsizes=[181,217,10];
imsize=[181,217,181];

ker=[1 1 1];
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

nii=load_nii('T1w_acpc_dc_restore.nii.gz');
imh=nii.img;
nii=load_nii('OAS1_0001_MR1_mpr-1_anon.img');
imo=nii.img;
imo=single(imo);


noimch=RicianSTD(imh);

noimo = RicianSTD(imo);
%%

%
%change β, fix α
%im2
im2=truncateslice(im2,5);index2=find(im2>0);
betaPSNR=zeros(21,5);
betaSSIM=zeros(21,5);
for i=1:5
    level=2*i-1;
    noise=2.55*level;
    nim=ricernd(im2,noise*ones(181,217,10));
    dnim=NLPCAnp(nim,noise,2.2);dnim=RiC(dnim,noise);

    for beta=1:0.1:3
        pim=RINLMmy(dnim,nim,noise,alpha,8,beta);
        %Comment out the code for calculating beta in RINLMmy, and add a
        %independent variable
        bb=round(10*beta-9);
        betaPSNR(bb,i)=20*log10(255/sqrt(mean((im2(index2)-pim(index2)).^2)));
        betaSSIM(bb,i)=ssim_index3d(im2,pim,[1 1 1],index2);

    end
end



%%
%full image time（层厚怎么定）
%10
index1=find(im1>0);

for i=1
    level=2*i-1;
    noise=2.55*level;
    nim=ricernd(im1,noise*ones(181,217,181));
    dnim=NLPCAnp(nim,noise,2.2);dnim=RiC(dnim,noise);
    tic;pim=RINLMmy(dnim,nim,noise,alpha,8);toc;
    ans=20*log10(255/sqrt(mean((im1(index1)-pim(index1)).^2)));
end

%%
% %change α, fix β，alpha = 0.57*level^(-0.31);
% %
% afiPSNR=zeros(1,5) ;afiSSIM=zeros(1,5);
% aflPSNR=zeros(1,5) ;aflSSIM=zeros(1,5);
% for i=1:5
%     noise=2.55*(2*i-1);
%     nim=ricernd(im2,noise*ones(181,217,10));
%     dnim=NLPCAnp(nim,noise,2.2);dnim=RiC(dnim,noise);
%     aflim=RINLMmy(dnim,nim,noise,0,8,3);%
%     afim=RINLMmy(dnim,nim,noise,0.4,8,3);
%     afiPSNR(i)=20*log10(255/sqrt(mean((im2(index2)-afim(index2)).^2)));
%     afiSSIM(i)=ssim_index3d(im2,afim,[1 1 1],index2);
%     aflPSNR(i)=20*log10(255/sqrt(mean((im2(index2)-aflim(index2)).^2)));
%     aflSSIM(i)=ssim_index3d(im2,aflim,[1 1 1],index2);
% 
% end
%%
%formula(2) vs exact noise level input
myPSNR=zeros(1,5) ;mySSIM=zeros(1,5);
majPSNR=zeros(1,5) ;majSSIM=zeros(1,5);
for i=1:5
    noise=2.55*(2*i-1);
    nim=ricernd(im1,noise*ones(181,217,181));
    dnim=NLPCAnp(nim,noise,2.2);dnim=RiC(dnim,noise);
    pim1=RINLMmy(dnim,nim,noise,alpha,8);
    pim2=cPRI_NL_PCA(nim,9,1,noise*ones(181,217,181),dnim,1); 
    myPSNR(i)=20*log10(255/sqrt(mean((im1(index1)-pim1(index1)).^2)));
    mySSIM(i)=ssim_index3d(im1,pim1,[1 1 1],index1);
    majPSNR(i)=20*log10(255/sqrt(mean((im1(index1)-pim2(index1)).^2)));
    majSSIM(i)=ssim_index3d(im1,pim2,[1 1 1],index1);

end



