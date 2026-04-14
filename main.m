clc;close all;clear all

%% 加解密图像展示
im = double(imread('date/color_image_512/avion.ppm'));
cover = double(imread('date/color_image_512/baboon.ppm'));

key='8d5ab8ba5340fce4420829ad5d12a0e45dacb0858544163d04c1d02b73e3697d';
kt='haar';

[m,n,k]=size(im);

[cip,NC,dnkey,max_x1,min_x1,index] = encryption(im,cover,key,kt);
[rim] = dencryption(cip,cover,dnkey,kt,max_x1,min_x1,index,[m,n,k]);


figure(11)
imshow(uint8(im))
figure(12)
imshow(uint8(cover))
figure(133)
imshow(NC/(max(NC(:))))
figure(13)
imshow(uint8(cip))
figure(14)
imshow(uint8(20*abs(cip-cover)))
figure(15)
imshow(uint8(rim))
figure(16)
imshow(uint8(50*abs(im-rim)))


