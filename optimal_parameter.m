%%
% read the file
file ='t1_icbm_normal_1mm_pn0_rf0.rawb';
fid = fopen(file,'r');    
imsize=[181,217,181];
global im

im=zeros(imsize(1:3));
for z=1:imsize(3)    
    im(:,:,z) = fread(fid,imsize(1:2));
end
fclose(fid);
sigma = 0.01*max(im(:));
%%
%load precomputed table , the values of etta and fai shouldn't be changed.
% global etta
% global fai
% load('precomputation.mat')

%%
%we add rician noise but with no median filter
%global Nim
global nnim
%Nim=ricernd(im, sigma*ones(imsize(1:3)));%im+randn(size(im))*sigma;
nnim = ricernd(im, sigma*ones(imsize(1:3)));%Nim;
%nnim =normrnd(im, sigma*ones(imsize(1:3)));%Nim;
%%
%finding the optimal parameter d, M, w, tau, beta for original algorithm

%lb = [2 8 2 0.1 1]; ub = [5 216 3 4 5];
lb = [2 8 2 0.1 1]; ub = [4 216 3 3 3];
%M=[3 27 3 2.4002 2.4576;4,64,3,2, 2];

%options = optimoptions('particleswarm','SwarmSize',50,'MaxIterations',50,...
%     'FunctionTolerance', 1e-3,'Display','iter','InitialSwarmMatrix',M ,...
%     'UseVectorized',true,'PlotFcn','pswplotbestf');
options = optimoptions('particleswarm','SwarmSize',50,'MaxIterations',50,...
     'FunctionTolerance', 1e-3,'Display','iter',...
     'UseVectorized',true,'PlotFcn','pswplotbestf');
[xit ,fval,exitflag,output] = particleswarm(@psnrall, 5 ,  lb,ub ,options );
save('result.mat',xit,fval,exitflag,output);
%%
%finding the optimal parameter d, M, w, tau, beta for original algorithm
nnim=normrnd(im, sigma*ones(imsize(1:3)));
lb=[0.1 1]; ub=[3 3];
%lb = [2 8 2 0.1 1]; ub = [4 216 3 3 3];
%M=[3 27 3 2.4002 2.4576;4,64,3,2, 2];
M=[2.46 2.46; 2.71 2];
options = optimoptions('particleswarm','SwarmSize',20,'MaxIterations',50,...
    'FunctionTolerance', 1e-3,'Display','iter','InitialSwarmMatrix',M ,...
    'UseVectorized',true,'PlotFcn','pswplotbestf');%

[xit ,fval,exitflag,output] = particleswarm(@psnrallM, 2 ,  lb,ub ,options );
save('result2.mat',xit,fval,exitflag,output);


%%
nnim =ricernd(im, sigma*ones(imsize(1:3)));%Nim;
%change tb, but fix T=2.4576
psnrcf=zeros(1,5);
for i=1:5
    tb=2.4-(i-1)*0.1;
    dnim=NLPCApso(single(nnim),3, 27, 3, tb, 2.4576);
    index = find(im>0);
    psnrcf(i) = 20*log10(255/sqrt(mean((im(index)-dnim(index)).^2)));
end

%change tb, but fix T=2.4576, this time plus
psnrcfplus=zeros(1,5);
for i=1:5
    tb=2.4+(i-1)*0.1;
    dnim=NLPCApso(single(nnim),3, 27, 3, tb, 2.4576);
    index = find(im>0);
    psnrcfplus(i) = 20*log10(255/sqrt(mean((im(index)-dnim(index)).^2)));
end

%let tb=T and change them
psnrcc=zeros(1,5);
for i=1:5
    bian=2.4-(i-1)*0.1;
    dnim=NLPCApso(single(nnim),3, 27, 3, bian, bian);
    index = find(im>0);
    psnrcc(i) = 20*log10(255/sqrt(mean((im(index)-dnim(index)).^2)));    
end
psnrccinv=zeros(1,5);
for i=1:5
    bian=2.4+(i)*0.02;
    dnim=NLPCApso(single(nnim),3, 27, 3, bian, bian);
    index = find(im>0);
    psnrccinv(i) = 20*log10(255/sqrt(mean((im(index)-dnim(index)).^2)));    
end

%%
%precomputation
ptheta=linspace(0,50,1000);
effSNR=sqrt((2+ptheta.^2)./kexi(ptheta)-1);
global etta
global fai
load('precomputation.mat')

%%
%Boundary checking
[dnim,~]=NLPCApso(single(nnim),4, 64, 3, 2.46, 2.46);
psnryj = 20*log10(255/sqrt(mean((im(index)-dnim(index)).^2)));  
[dnim,~]=NLPCApso(single(nnim),2, 8, 3, 2.46, 2.46);
psnrzj = 20*log10(255/sqrt(mean((im(index)-dnim(index)).^2))); 






