function [PRINLMfim] = RI_NLM( gim, nim, noise, miu,alpha)
A=size(noise,1);
B=size(noise,2);
C=size(noise,3);
kPRINLMfim=zeros(A,B,C);
%record1=zeros(181,217,181);
%record2=zeros(181,217,181);
for j=6:B-5
	for i=6:A-5
		for k=6:C-5
			qiuheweigh = 0;
			qiuhepatchweigh = 0; 
			h = alpha*noise(i,j,k);
            if h==0
                kPRINLMfim(i,j,k)=gim(i,j,k);
            else
                for n=0:10
                    for m=0:10
                        for p=0:10
                            g2 = (gim(i,j,k)-gim(i+m-5,j+n-5,k+p-5))^2; 
                            u2 = 3*(miu(i,j,k)-miu(i+m-5,j+n-5,k+p-5))^2;
                            if(abs(miu(i,j,k)-miu(i+m-5,j+n-5,k+p-5))<h)
                                weigh = exp(-0.5*((g2+u2)/(2*h*h)));	
                            else
                                weigh=0;
                            end
%                             if sqrt(g2 )<0.1
%                                 record1( i-5, j-5, k-5 )=record1(i-5, j-5, k-5 )+ 1;
%                             else
%                                 record2(i-5, j-5, k-5 )=record2(i-5, j-5, k-5 )+ 1;
%                             end
                            qiuheweigh = qiuheweigh+weigh;
                            patch_weigh = weigh*nim(i+m-5,j+n-5,k+p-5)^2;
                            %patch_weigh = weigh*nim(i+m-5,j+n-5,k+p-5);
                            qiuhepatchweigh = qiuhepatchweigh+patch_weigh;
%                             if vim(i+m-5,j+n-5,k+p-5)==vim(i,j,k)
%                                 record1(i-5,j-5,k-5)=record1(i-5,j-5,k-5)+weigh;%sumweight
%                             end
                            
                        end
                    end
                end
%                 record1(i-5,j-5,k-5)=record1(i-5,j-5,k-5)/qiuheweigh;
                temp = qiuhepatchweigh/qiuheweigh-2*noise(i,j,k)^2 ;
                kPRINLMfim(i,j,k) = sqrt( max(temp ,0.0) );
                %kPRINLMfim(i,j,k) = sig_patch_weigh/sig_weigh;
%                 if temp<0
%                     record1( i-5, j-5, k-5 )=1;
%                 end    
            end
        end
    end
end
PRINLMfim=kPRINLMfim(6:(A-5),6:(B-5),6:(C-5));

end

