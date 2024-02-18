function [PRINLMfim] = RINLMmy( gim, nim, noise, alpha,dtr)

%0~255

%Round only after the pre filter is processed

if dtr==8
    ma=1;%255
    level = noise/(255*ma)*100;
    noise = noise/ma;
% elseif dtr == 12
%     %ma=4096/255;
%     ma=ma0/255;
%     level = noise/ma0*100;
%     noise = 2.55*level;
end
%ma=1;

%level = noise/(max(gim,[],'all'))*100;


%level = noise/(255*ma)*100;

beta = 3.2-(level+1)/2*0.3;
%beta=3;
%gim(gim<0)=0;
gim = gim/ma;




miu = imboxfilt3( gim  , 'padding' , 'symmetric' ) ;
% kernel = ones(3,3,3) / 27;
%     miu = convn(gim, kernel, 'same');

ggim = ceil(gim);
%ggim = round(gim);

A0=size(ggim,1);
B0=size(ggim,2);
C0=size(ggim,3);
ABC=A0*B0*C0;

%sg=1 represents a case of slice
sg=ceil(C0/10);




mmiu = ceil(miu);
%mmiu = round(miu);

ma1 = (max(ggim,[],'all')); mi1 = min(ggim,[],'all');
ma2 = (max(mmiu,[],'all')); mi2 = min(mmiu,[],'all');



% gggim = reshape(ggim,A0*B0*C0,1)+1;
% mmmiu = reshape(mmiu,A0*B0*C0,1)+1;
ggim = reshape(ggim,ABC,1);
mmiu = reshape(mmiu,ABC,1);
gim = reshape(gim,ABC,1);
miu = reshape(miu,ABC,1);


%Save with sparse matrix?
% gm_hist = accumarray([gggim,mmmiu], 1, [ma1+1, ma2+1]); 



w = zeros(ma1-mi1+1,ma2-mi2+1);
uw = zeros(ma1-mi1+1,ma2-mi2+1);

% if alph==0
%     alpha = 0.57*level^(-0.31);
% else
%     alpha=alph;
% end
h=alpha*noise;





%pri = zeros(ma2+1,ma1+1);



n2 = reshape(nim.^2,ABC,1)/ma^2;




%[w,uw] = RINLMmyin_mex(gim,ggim,miu,mmiu,n2,h,level,ma2,ma1);



%Further improvements can be made here to increase speed, such as 
% considering that some μ will not appear, so there is no need to iterate 
% from 0 to the maximum value
%for i=0:ma1
for j=mi2:ma2
    listi = (unique(ggim(mmiu==j)))';
    %idxz1 = find(abs(ggim-i)<3*h);

    idxz1 = find(abs(miu-j)<beta*h);
    idxz1 = randsample(idxz1,round(length(idxz1)/(sg*level)));

    %g2 = (ggim(idxz1) - i).^2;

    %u2 = 3*(miu(idxz1) - j).^2;
    %u2 = 3*(miu - j).^2;

    %g2 = (ggim - i).^2;

    ni2=n2(idxz1);

    for idi=1:length(listi)
    %for j=listj
        i=listi(idi);
        %idxz2 = find((g2<2.5*h)&(abs(mmiu-j)<2.5*h));
        %idxz = find((abs(ggim-i)<2.5*h)&(abs(mmiu-j)<2.5*h));

        ggm = gim(idxz1);
        idxz2 = find((abs(ggm-i)<beta*h));

        %idxz2 = randsample(idxz2,round(length(idxz2)/5));
        %u2 = 3*(mmu(idxz2) - j).^2;
        %u2 = 3*(mmiu - j).^2;

        g2 = (ggm(idxz2) - i).^2;
        %g2= (gim - i).^2;

        %temp1 = exp(-(g2+u2(idxz2))/(4*h^2));
        temp1 = exp(-(g2)/(4*h^2));
        %w(j+1,i+1) = sum(temp1(:));
        w(i+1-mi1,j+1-mi2) = sum(temp1(:));
        temp2 = temp1.*(ni2(idxz2));
        %temp2 = temp1.*(ni1(idxz2));
        %temp2 = temp1.*n2;
        uw(i+1-mi1,j+1-mi2) = sum(temp2(:));
        %PRI(i+1,j+1)=sum(temp2(:))/sum(temp1(:))-;
    end
end



% for i=mi1:ma1
%     listj = (unique(mmiu(ggim==i)))';
%     %idxz1 = find(abs(ggim-i)<3*h);
% 
%     idxz1 = find(abs(gim-i)<3*h);
%     idxz1 = randsample(idxz1,round(length(idxz1)/level));
% 
%     %g2 = (ggim(idxz1) - i).^2;
%     g2 = (gim(idxz1) - i).^2;
% 
% 
%     %u2 = 3*(miu(idxz1) - j).^2;
%     %u2 = 3*(miu - j).^2;
% 
%     %g2 = (ggim - i).^2;
% 
%     ni2=n2(idxz1);
% 
%     %ni1=n1(idxz1);
%     for idj=1:length(listj)
%     %for j=listj
%         j=listj(idj);
%         %idxz2 = find((g2<2.5*h)&(abs(mmiu-j)<2.5*h));
%         %idxz = find((abs(ggim-i)<2.5*h)&(abs(mmiu-j)<2.5*h));
%         mmu = miu(idxz1);
%         idxz2 = find((abs(mmu-j)<3*h));
% 
%         %idxz2 = randsample(idxz2,round(length(idxz2)/5));
%         %u2 = 3*(mmu(idxz2) - j).^2;
%         %u2 = 3*(mmiu - j).^2;
% 
% 
%         %g2= (gim - i).^2;
% 
%         %temp1 = exp(-(g2+u2(idxz2))/(4*h^2));
%         temp1 = exp(-(g2(idxz2))/(4*h^2));
%         %w(j+1,i+1) = sum(temp1(:));
%         w(i+1-mi1,j+1-mi2) = sum(temp1(:));
%         temp2 = temp1.*(ni2(idxz2));
%         %temp2 = temp1.*(ni1(idxz2));
%         %temp2 = temp1.*n2;
%         uw(i+1-mi1,j+1-mi2) = sum(temp2(:));
%         %PRI(i+1,j+1)=sum(temp2(:))/sum(temp1(:))-;
%     end
% end


PRI=sqrt(max(uw./w-2*noise^2,0));

%look up table or interpolate

%ind = sub2ind(size(PRI),mmiu+1-mi2,ggim+1-mi1);
% %ind = sub2ind(size(w),ggggim+1,mmmmiu+1);
%ind = sub2ind(size(PRI),ggim+1,mmiu+1);
%PRINLMfim = PRI(ind);
% ww = w(ind);
% uuww = uw(ind);


[pri1,pri2] = find(PRI~=0);
pri=PRI(PRI~=0);
F = scatteredInterpolant(pri1-1, pri2-1, pri, 'nearest');



% maiv=1/ma;
% [xq,yq] = meshgrid(mi1:maiv:ma1, mi2:maiv:ma2);
% vq = griddata(pri1-1,pri2-1,pri,xq,yq,'cubic');
% vq(isnan(vq))=0;
% PRINLMfim = interp2(xq,yq,vq,gim,miu,'nearest',0);
PRINLMfim = F([gim,miu]);
PRINLMfim=reshape(PRINLMfim,A0,B0,C0)*ma;



end



