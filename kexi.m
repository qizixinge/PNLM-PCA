function [correction] = kexi(SNR_real)
%简写SNR_real为theta

theta=SNR_real;
%matlab有bug,SNR超过37就会变成NAN
%if theta<37
    correction=2+theta.^2-pi/8.*(exp(-theta.^2/2)).*((2+theta.^2).*besseli(0,theta.^2/4)+theta.^2.*besseli(1,theta.^2/4)).^2;
    %correction=2+theta^2-pi/8*(exp(-theta^2/2))*((2+theta^2)*besseli(0,theta^2/4)+theta^2*besseli(1,theta^2/4))^2;
%else
    %correction=1;
%end
id0= isnan(correction);  correction(id0) =  1;
id1=(abs(correction)>1);correction(id1)=1;
%注意是修正贝塞尔函数
end