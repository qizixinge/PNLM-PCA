function [out] = dis255(img)
img=round(img);
img(img>255)=255;
out=img;


end