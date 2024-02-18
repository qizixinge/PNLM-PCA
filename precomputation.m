%precomputation
%precomputed table
fai1 = linspace (0,50 , 1000) ; 
etta1=exp(-fai1.^2/4).*sqrt(pi/2).*((fai1.^ 2 / 2+ 1 ).* ...
    besseli(0,fai1.^2/4)+fai1.^2/2.*besseli(1,fai1.^2/4));
fai2 = linspace (50.01,500,13500) ; 
%100对应1%
%etta2=0.9863*fai2.^1.003+0.1241;
etta2=sqrt(fai2.^2+1);
fai = [ fai1 , fai2] ; etta = [ etta1 , etta2 ] ;
save ( ' precomputation.mat'  , 'fai'   ,  'etta' )







