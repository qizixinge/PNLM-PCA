function [dnim] = NLPCApso(nnim,d,M,w,tb, T)
xw=2*w+1; W=xw^3;
XW=zeros(1,xw);
for i=1:xw
    XW(i)=mod(i,xw);
end
knim=padarray(nnim , [ XW(xw-mod(size(nnim , 1) ,xw)) ,XW( xw-mod(size(nnim , 2) ,xw)) ,XW(xw-mod(size(nnim , 3),xw))  ] , 'symmetric' , 'post');
%numofwindow=(size(knim, 1)-xw+1)*(size(knim, 2)-xw+1)*(size(knim, 3)-xw+1);
numofwindow=(size(knim, 1)/xw)*(size(knim, 2)/xw)*(size(knim, 3)/xw);
%A=size(knim, 1)-xw+1; B=size(knim, 2)-xw+1; C=size(knim, 3)-xw+1;
A=size(knim, 1)/xw; B=size(knim, 2)/xw; C=size(knim, 3)/xw;
a=size(knim, 1); b=size(knim, 2); c=size(knim, 3);
D=xw-d+1; D3=d^3 ; D2=d^2;
numofpatchw=(floor(D))^3; Ds2=(floor(D))^2; Ds = (floor(D));
%initialiazation
pt = zeros(numofpatchw , d^3  ); dist2 =zeros(1 , numofpatchw ); 
XTO = zeros(M ,D3  ); %order = zeros(1 , D3 );
sm = zeros(a, b, c  );  COUNT = zeros( a, b, c );
for j=1:numofwindow
    jj = j-1;
    for cc=1:d
        for rr=1:d
            for zz=1:d
                for k=1:numofpatchw
                    kk = k-1;
                    %temp0 = [xw*floor(mod(jj,(B*A))/B)+s*floor(mod(kk ,Ds2)/Ds) , xw*mod(jj,B)+s*mod(kk,Ds) ,  xw*floor(jj/(B*A))+s*floor(kk/Ds2)];
                    pt(k, (zz-1)*d^2+(rr-1)*d+cc)=knim( xw*floor(mod(jj,(B*A))/B)+floor(mod(kk ,Ds2)/Ds)+rr , xw*mod(jj,B)+mod(kk,Ds)+cc ,  xw*floor(jj/(B*A))+floor(kk/Ds2)+zz );
                    %pt(k , :) = reshape( knim( temp0(1)+1:temp0(1)+d , temp0(2)+1:temp0(2)+d , temp0(3)+1:temp0(3)+d) , 1 , D3 );
                end
            end
        end
    end
    for k=1:numofpatchw
        for r=1:numofpatchw
            dist2(r) = sum((pt(k , : )-pt(r , : )).^2);
        end
        [~ , distsx] = mink( dist2 ,M);
        for rr=1:M
            XTO(rr , :) = pt(distsx(rr) , : );
        end
        XT = XTO-ones(size(XTO,1),1)*mean(XTO) ;
        XX = XT' *XT/(M-1);
        [Me , Ve] = eig(XX);%
        Ve = diag(Ve);
        [Ve, order]=sort(Ve , 'descend'); 
        Me = Me(: , order); 
        Ve = sqrt(abs(Ve));
        %
        medlmd = median(Ve );%linear median
        temp3 = sum(Ve<=T*medlmd); medlmdt = median(Ve((D3-temp3+1):end ));
        Tau = tb* (medlmdt);
        %Tau = tau*sigma;
        Mes = Me( : , 1 : sum(Ve>=Tau));%linear threshold
        Y=XT*Mes;
        XFT=Y*Mes'+ones(size(XTO,1),1)*mean(XTO);
        for cc=1:d
            for rr=1:d
                for zz=1:d
                    for i=1:M
                        temp=distsx(i)-1;
                        temp2=[ xw*floor(mod(jj,(B*A))/B)+floor(mod(temp ,Ds2)/Ds)+rr , xw*mod(jj,B)+mod(temp,Ds)+cc ,  xw*floor(jj/(B*A))+floor(temp/Ds2)+zz];
                        %temp2=[ xw*floor(mod(jj,(B*A))/B)+s*floor(mod(temp ,Ds2)/Ds) , xw*mod(jj,B)+s*mod(temp,Ds) ,  xw*floor(jj/(B*A))+s*floor(temp/Ds2)];
                        %sm(temp2(1)+1:temp2(1)+d , temp2(2)+1:temp2(2)+d , temp2(3)+1:temp2(3)+d)= sm(temp2(1)+1:temp2(1)+d , temp2(2)+1:temp2(2)+d , temp2(3)+1:temp2(3)+d)+ reshape(XFT(i , :), d, d, d);
                        %smg(temp2(1)+1:temp2(1)+d , temp2(2)+1:temp2(2)+d , temp2(3)+1:temp2(3)+d) =smg(temp2(1)+1:temp2(1)+d , temp2(2)+1:temp2(2)+d , temp2(3)+1:temp2(3)+d)+(sigma)*ones(d,d,d);
                        %COUNT(temp2(1)+1:temp2(1)+d , temp2(2)+1:temp2(2)+d , temp2(3)+1:temp2(3)+d)=COUNT(temp2(1)+1:temp2(1)+d , temp2(2)+1:temp2(2)+d , temp2(3)+1:temp2(3)+d)+ones(d,d,d);
                        sm(temp2(1) , temp2(2) , temp2(3))= sm(temp2(1) , temp2(2) , temp2(3))+ (XFT(i , (zz-1)*D2+(rr-1)*d+cc));
                        %smg(temp2(1) , temp2(2) , temp2(3)) =smg(temp2(1) , temp2(2) , temp2(3))+(sigma);
                        COUNT(temp2(1) , temp2(2) , temp2(3))=COUNT(temp2(1) , temp2(2) , temp2(3))+1;
                    end
                end
            end
        end
    end    
end    
kdnim = sm./COUNT;%biased estimation
%knoise_mapr = smg./COUNT;%linear average
dnim = kdnim( 1:size(nnim,1) , 1:size(nnim,2) , 1:size(nnim,3));
end

