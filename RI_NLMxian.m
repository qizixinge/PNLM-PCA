function [PRINLMfim] = RI_NLMxian(gim, nim, noise, miu,alpha)
A=size(noise,1);
B=size(noise,2);
C=size(noise,3);
kPRINLMfim=zeros(191,227,191);
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
                            qiuheweigh = qiuheweigh+weigh;
                            %patch_weigh = weigh*nim(i+m-5,j+n-5,k+p-5)^2;
                            patch_weigh = weigh*nim(i+m-5,j+n-5,k+p-5);
                            qiuhepatchweigh = qiuhepatchweigh+patch_weigh;
                        end
                    end
                end
                %temp = qiuhepatchweigh/qiuheweigh-2*noise(i,j,k)^2 ;
                %kPRINLMfim(i,j,k) = sqrt( max(temp ,0.0) );
                kPRINLMfim(i,j,k) = qiuhepatchweigh/qiuheweigh;
            end
        end
    end
end
PRINLMfim=kPRINLMfim(6:186,6:222,6:186);
end

