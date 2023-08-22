function [B] = ricernd(v,s)
%there must be size(v)=size(s)
%this function is applicable to both stationary noise and spatially varying
%noise
%Exactly right
n1=normrnd(zeros(size(v)),s);
n2=normrnd(zeros(size(v)),s);
B=sqrt((v+n1).^2+(n2).^2);
%noise modulation function beta is not considered here, see 2010 Manjon in details.
end