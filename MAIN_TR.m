clc;close all;clear all

%% 加解密图像展示
% 灰度

im = imread('E:\date\ccia_CVG_image\color_image_512/avion.ppm');
im = double(rgb2gray(im));
cover = imread('E:\date\ccia_CVG_image\color_image_512/sailboat.ppm');
cover = double(rgb2gray(cover));

% im = imread('date/color_image_512/windmill.ppm');
% im = double(rgb2gray(im));
% cover = imread('date/color_image_512/utahmtn.ppm');
% cover = double(rgb2gray(cover));

% im = imread('date/color_image_512/elephant.ppm');
% im = double(rgb2gray(im));
% cover = imread('date/color_image_512/house.ppm');
% cover = double(rgb2gray(cover));
% 
% 彩色
% im = double(imread('date/color_image_512/avion.ppm'));
% cover = double(imread('date/color_image_512/baboon.ppm'));
% 
% im = double(imread('date/color_image_512/frog.ppm'));
% cover = double(imread('date/color_image_512/beeflowr.ppm'));
% 
% 封面彩色，明文灰度
% im = imread('date/color_image_512/peppers.ppm');
% im = double(rgb2gray(im));
% cover = double(imread('date/color_image_512/portofino.ppm'));
% 
% 封面灰度，明文彩色
% im = imread('date/color_image_512/raiz1.ppm');
% im = double(imresize(im,[256,256]));
% cover = imread('date/color_image_512/sailboat.ppm');
% cover = double(rgb2gray(cover));


% im = imread('E:\date\ccia_CVG_image\color_image_512/lena.ppm');
% cover = imread('E:\date\ccia_CVG_image\color_image_512/peppers.ppm');
% im = imread('E:\date\ccia_CVG_image\color_image_512/baboon.ppm');
% im = rgb2gray(im);
% cover = rgb2gray(cover);
% im = imresize(im, 0.5);
% cover = imresize(cover, 0.5);

im = double(im);
cover=double(cover);

key='8d5ab8ba5340fce4420829ad5d12a0e45dacb0858544163d04c1d02b73e3697d';
% kt='spa';   % spa空间域其余为变换域
kt='haar';   % spa空间域其余为变换域

[m,n,k]=size(im);

[cip,NC,dnkey,max_x1,min_x1,index,X1,X2] = encryption(im,cover,key,kt);
[rim,RNC,RX2,RX1] = dencryption(cip,cover,dnkey,kt,max_x1,min_x1,index,[m,n,k]);

psnr_cip = psnr(double(uint8(cip)),double(uint8(cover)),255)
ssim_cip = ssim(double(uint8(cip)),double(uint8(cover)))

psnr_rim = psnr(double(uint8(im)),double(uint8(rim)),255)
ssim_rim = ssim(double(uint8(im)),double(uint8(rim)))

[c_NC1,c_NC2,c_NC3] = image_coor(NC);
en_nc = double_entropy(NC(:), 11*5*5*5);

[c_cip1,c_cip2,c_cip3] = image_coor(cip);
en_cip = double_entropy(cip(:), 11*5*5*5);


figure(11)
imshow(uint8(im))
figure(12)
imshow(uint8(cover))
figure(133)
imshow(NC/(max(NC(:))))
figure(13)
imshow(uint8(cip))
figure(14)
imshow(uint8(50*abs(cip-cover)))
figure(15)
imshow(uint8(rim))
figure(16)
imshow(uint8(30*abs(im-rim)))