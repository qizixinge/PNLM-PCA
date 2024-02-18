function [gim] = RiC(dnim,noise)
global etta
global fai
%should be sqrt(fai^2+1)
A=size(dnim,1);
B=size(dnim,2);
C=size(dnim,3);
noise_map=noise*ones(A,B,C);
eetta = dnim ./ noise_map;
            eetta=reshape(eetta,1,A*B*C);
            idxe = (eetta<min(etta)); eetta(idxe) = min(etta);%avoid outliers

            temp = interp1(etta, fai, eetta, 'linear');%arrayfun(@ettainv, eetta);
            temp=reshape(temp,A,B,C);
            gim = noise_map.* temp;
end