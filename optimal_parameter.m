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
im = truncateslice(im, 8);
sigma = 0.01*max(im(:));



%we add rician noise but with no median filter
%global Nim
global nnim
%Nim=ricernd(im, sigma*ones(imsize(1:3)));%im+randn(size(im))*sigma;
%nnim = ricernd(im, sigma*ones(imsize(1:3)));%Nim;
nnim =normrnd(im, sigma*ones(181,217,16));%Nim;

%im = truncateslice(im, 3);


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

lb=[0.1 1]; ub=[4 4];

M=[2.46 2.46; 2.83 2];
options = optimoptions('particleswarm','SwarmSize',20,'MaxIterations',50,...
    'FunctionTolerance', 1e-3,'Display','iter','InitialSwarmMatrix',M ,...
    'UseVectorized',true,'PlotFcn','pswplotbestf');%

[xit ,fval,exitflag,output] = particleswarm(@psnrallM, 2 ,  lb,ub ,options );
save('result2.mat',xit,fval,exitflag,output);




%%
%Boundary checking
%dnim=NLPCApso(single(nnim),4, 64, 3, 2.46, 2.46);
%psnryj = 20*log10(255/sqrt(mean((im(index)-dnim(index)).^2)));  
%dnim=NLPCApso(single(nnim),2, 8, 3, 2.46, 2.46);
%psnrzj = 20*log10(255/sqrt(mean((im(index)-dnim(index)).^2))); 

%%
%for T2w and PDw normal, find taubeta in 1~3 iteratively
%preparation
file2 ='t2_icbm_normal_1mm_pn0_rf0.rawb';
fid2 = fopen(file2,'r');    
imsize=[181,217,181];

im2=zeros(imsize(1:3));
for z=1:imsize(3)    
    im2(:,:,z) = fread(fid2,imsize(1:2));
end
fclose(fid2);

file3 ='pd_icbm_normal_1mm_pn0_rf0.rawb';
fid3 = fopen(file3,'r');    
imsize=[181,217,181];

im3=zeros(imsize(1:3));
for z=1:imsize(3)    
    im3(:,:,z) = fread(fid3,imsize(1:2));
end
fclose(fid3);

nnim2 =normrnd(im2, sigma*ones(imsize(1:3)));
nnim3 =normrnd(im3, sigma*ones(imsize(1:3)));
psnr2=zeros(1,21);
psnr3=zeros(1,21);
%%
%for T2w and PDw normal, find taubeta in 1~3 iteratively
%start
index2 = find(im2>0); index3 = find(im3>0);

for i=1:0.1:3
    T = double(i);
    [dnim2,~] = NLPCA(nnim2,1,T, T);
    [dnim3,~] = NLPCA(nnim3,1,T, T);
    psnr2(round(10*i-9)) = 20*log10(255/sqrt(mean((im2(index2)-dnim2(index2)).^2)));
    psnr3(round(10*i-9)) = 20*log10(255/sqrt(mean((im3(index3)-dnim3(index3)).^2)));
end

%%
file4 ='t2_ai_msles2_1mm_pn0_rf0.rawb';
fid4 = fopen(file4,'r');    
imsize=[181,217,181];

im4=zeros(imsize(1:3));
for z=1:imsize(3)    
    im4(:,:,z) = fread(fid4,imsize(1:2));
end
fclose(fid4);



nnim4 =normrnd(im4, sigma*ones(imsize(1:3)));
psnr4=zeros(1,21);

%for T2w msles, find taubeta in 3~1 iteratively
%start
index4 = find(im4>0); 
for i=3:-0.1:1
    T = double(i);
    [dnim4,~] = NLPCA(nnim4,1,T, T);
    psnr4(round(10*i-9)) = 20*log10(255/sqrt(mean((im4(index4)-dnim4(index4)).^2)));
end







