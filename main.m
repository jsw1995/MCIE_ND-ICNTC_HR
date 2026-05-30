clc;close all;clear all

%% 加解密图像展示
% 灰度

% % im = imread('E:\date\ccia_CVG_image\color_image_512/avion.ppm');
% % im = double(rgb2gray(im));
% % cover = imread('E:\date\ccia_CVG_image\color_image_512/sailboat.ppm');
% % cover = double(rgb2gray(cover));
% 
% im = imread('date/color_image_512/windmill.ppm');
% im = double(rgb2gray(im));
% cover = imread('date/color_image_512/utahmtn.ppm');
% cover = double(rgb2gray(cover));
% 
% % im = imread('date/color_image_512/elephant.ppm');
% % im = double(rgb2gray(im));
% % cover = imread('date/color_image_512/house.ppm');
% % cover = double(rgb2gray(cover));

% 彩色
% im = double(imread('date/color_image_512/avion.ppm'));
% cover = double(imread('date/color_image_512/baboon.ppm'));

% im = double(imread('date/color_image_512/frog.ppm'));
% cover = double(imread('date/color_image_512/beeflowr.ppm'));
% 
% 封面彩色，明文灰度
% im = imread('date/color_image_512/peppers.ppm');
% im = double(rgb2gray(im));
% cover = double(imread('date/color_image_512/portofino.ppm'));

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
kt='spa';   % spa空间域其余为变换域
% kt='haar';   % spa空间域其余为变换域

[m,n,k]=size(im);
tem = load(strcat('tem/',kt , '_tem.mat'));
tem = getfield(tem, 'tem2');

[cip,NC,dnkey,max_x1,min_x1,index,X1,X2] = encryption(im,cover,key,kt,tem);
[rim,RNC,RX2,RX1] = dencryption(cip,cover,dnkey,kt,max_x1,min_x1,index,tem,[m,n,k]);

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

%% 不同变换核下密文对比
% im = imread('date/color_image_512/avion.ppm');
% im = double(rgb2gray(im));
% cover = imread('date/color_image_512/baboon.ppm');
% cover = double(rgb2gray(cover));
% 
% 
% key='8d5ab8ba5340fce4420829ad5d12a0e45dacb0858544163d04c1d02b73e3697d';
% % kt='spa';
% kt='haar';
% 
% [m,n,k]=size(im);
% 
% [cip,NC,dnkey,max_x1,min_x1,index] = encryption(im,cover,key,kt);
% 
% psnr_cip = psnr(double(uint8(cip)),double(uint8(cover)),255)
% ssim_cip = ssim(double(uint8(cip)),double(uint8(cover)))
% 
% figure(11)
% imshow(uint8(im))
% figure(12)
% imshow(uint8(cover))
% figure(13)
% imshow(uint8(cip))
% figure(14)
% imshow(uint8(30*abs(cip-cover)))

%% 不同类高斯噪声下密文对比
% im = imread('date/color_image_512/avion.ppm');
% im = double(rgb2gray(im));
% cover = imread('date/color_image_512/baboon.ppm');
% cover = double(rgb2gray(cover));
% 
% 
% key='8d5ab8ba5340fce4420829ad5d12a0e45dacb0858544163d04c1d02b73e3697d';
% kt='haar';
% nc_ty='cs';  % 'cs','dc', '4f','co'
% 
% [m,n,k]=size(im);
% 
% [cip,NC,dnkey,max_x1,min_x1,index] = encryption3(im,cover,key,kt,nc_ty);
% 
% psnr_cip = psnr(double(uint8(cip)),double(uint8(cover)),255)
% ssim_cip = ssim(double(uint8(cip)),double(uint8(cover)))
% 
% figure(11)
% imshow(uint8(im))
% figure(12)
% imshow(uint8(cover))
% figure(13)
% imshow(uint8(cip))
% figure(14)
% imshow(uint8(30*abs(cip-cover)))

%% 类高斯噪声对密文影响
% % plain_file = 'E:\date\\ccia_CVG_image\color_image_512/';
% % fileExt = '*.ppm';  %待读取图像的后缀名
% % img_size = 512;
% 
% plain_file = 'E:\date\BSDS100/';
% fileExt = '*.png';  %待读取图像的后缀名
% img_size = 256;
% 
% % plain_file = 'E:\date\DIV2K\DIV2K_valid_HR/';
% % fileExt = '*.png';  %待读取图像的后缀名
% % img_size = 1024;
% 
% plain_files = dir(fullfile(plain_file,fileExt)); 
% plain_len = size(plain_files,1);
% 
% psnr11=ones(plain_len-1,2);
% ssim11=ones(plain_len-1,2);
% 
% key='8d5ab8ba5340fce4420829ad5d12a0e45dacb0858544163d04c1d02b73e3697d';
% kt='db4';   % 依次改变
% nc_ty='dc';  % 'cs','dc', '4f','co'
% img_size_half = int16(img_size/2);
% 
% for i=1:plain_len-1
%     i
%     tim = double(rgb2gray(imread(strcat(plain_file,plain_files(i).name))));
%     [m,n]=size(tim);
%     if m>img_size && n>img_size
%         im = tim(floor(m/2)-img_size_half:floor(m/2)+img_size_half-1,floor(n/2)-img_size_half:floor(n/2)+img_size_half-1,:);
%     elseif m<img_size && n>img_size
%         im = tim(:,floor(n/2)-img_size_half:floor(n/2)+img_size_half-1,:);
%         im = double(imresize(uint8(im),[img_size,img_size]));
%     elseif m>img_size && n<img_size
%         im = tim(floor(m/2)-img_size_half:floor(m/2)+img_size_half-1,:,:);
%         im = double(imresize(uint8(im),[img_size,img_size]));
%     else
%         im = double(imresize(uint8(tim),[img_size,img_size]));
%     end
% 
%     tcover = double(rgb2gray(imread(strcat(plain_file,plain_files(i+1).name))));
%     [m,n]=size(tcover);
%     if m>img_size && n>img_size
%         cover = tcover(floor(m/2)-img_size_half:floor(m/2)+img_size_half-1,floor(n/2)-img_size_half:floor(n/2)+img_size_half-1,:);
%     elseif m<img_size && n>img_size
%         cover = tcover(:,floor(n/2)-img_size_half:floor(n/2)+img_size_half-1,:);
%         cover = double(imresize(uint8(cover),[img_size,img_size]));
%     elseif m>img_size && n<img_size
%         cover = tcover(floor(m/2)-img_size_half:floor(m/2)+img_size_half-1,:,:);
%         cover = double(imresize(uint8(cover),[img_size,img_size]));
%     else
%         cover = double(imresize(uint8(tcover),[img_size,img_size]));
%     end
%     
%     
%     [cip2,NC,dnkey,max_x1,min_x1,index] = encryption3(im,cover,key,kt,nc_ty);
%     
%     
%     psnr1 = psnr(double(uint8(cip2)),cover,255);
%     ssim1 = ssim(uint8(cip2),uint8(cover));
%     psnr11(i,2) = psnr1;
%     ssim11(i,2) = ssim1;
%     
% end
% 
% psnr2 = mean(psnr11(:,2))
% ssim2 = mean(ssim11(:,2))


%% 变换核对密文影响
% plain_file = 'E:\date\\ccia_CVG_image\color_image_512/';
% fileExt = '*.ppm';  %待读取图像的后缀名
% img_size = 512;
% 
% % plain_file = 'E:\date\BSDS100/';
% % fileExt = '*.png';  %待读取图像的后缀名
% % img_size = 256;
% 
% % plain_file = 'E:\date\DIV2K\DIV2K_valid_HR/';
% % fileExt = '*.png';  %待读取图像的后缀名
% % img_size = 1024;
% 
% plain_files = dir(fullfile(plain_file,fileExt)); 
% plain_len = size(plain_files,1);
% 
% psnr11=ones(plain_len-1,2);
% ssim11=ones(plain_len-1,2);
% 
% key='8d5ab8ba5340fce4420829ad5d12a0e45dacb0858544163d04c1d02b73e3697d';
% kt='db4';   % 依次改变
% img_size_half = int16(img_size/2);
% 
% for i=1:plain_len-1
%     i
%     tim = double(rgb2gray(imread(strcat(plain_file,plain_files(i).name))));
%     [m,n]=size(tim);
%     if m>img_size && n>img_size
%         im = tim(floor(m/2)-img_size_half:floor(m/2)+img_size_half-1,floor(n/2)-img_size_half:floor(n/2)+img_size_half-1,:);
%     elseif m<img_size && n>img_size
%         im = tim(:,floor(n/2)-img_size_half:floor(n/2)+img_size_half-1,:);
%         im = double(imresize(uint8(im),[img_size,img_size]));
%     elseif m>img_size && n<img_size
%         im = tim(floor(m/2)-img_size_half:floor(m/2)+img_size_half-1,:,:);
%         im = double(imresize(uint8(im),[img_size,img_size]));
%     else
%         im = double(imresize(uint8(tim),[img_size,img_size]));
%     end
% 
%     tcover = double(rgb2gray(imread(strcat(plain_file,plain_files(i+1).name))));
%     [m,n]=size(tcover);
%     if m>img_size && n>img_size
%         cover = tcover(floor(m/2)-img_size_half:floor(m/2)+img_size_half-1,floor(n/2)-img_size_half:floor(n/2)+img_size_half-1,:);
%     elseif m<img_size && n>img_size
%         cover = tcover(:,floor(n/2)-img_size_half:floor(n/2)+img_size_half-1,:);
%         cover = double(imresize(uint8(cover),[img_size,img_size]));
%     elseif m>img_size && n<img_size
%         cover = tcover(floor(m/2)-img_size_half:floor(m/2)+img_size_half-1,:,:);
%         cover = double(imresize(uint8(cover),[img_size,img_size]));
%     else
%         cover = double(imresize(uint8(tcover),[img_size,img_size]));
%     end
%     
%     
%     [cip2,NC,dnkey,max_x1,min_x1,index] = encryption(im,cover,key,kt);
%     
%     
%     psnr1 = psnr(double(uint8(cip2)),cover,255);
%     ssim1 = ssim(uint8(cip2),uint8(cover));
%     psnr11(i,2) = psnr1;
%     ssim11(i,2) = ssim1;
%     
% end
% 
% psnr2 = mean(psnr11(:,2));
% ssim2 = mean(ssim11(:,2));


%% 算法常规安全性分析
% 统计特性分析 
% im = double(imread('date/color_image_512/avion.ppm'));
% cover = double(imread('date/color_image_512/baboon.ppm'));
% 
% key = '8d5ab8ba5340fce4420829ad5d12a0e45dacb0858544163d04c1d02b73e3697d';
% kt='haar';
% 
% [m,n,k]=size(im);
% 
% [cip,tcip,dnkey,max_x1,min_x1,index] = encryption(im,cover,key,kt);
% [rim] = dencryption(cip,cover,dnkey,kt,max_x1,min_x1,index,[m,n,k]);
% 
% psnr_cip = psnr(double(uint8(cip)),double(uint8(cover)),255)
% ssim_cip = ssim(double(uint8(cip)),double(uint8(cover)))
% 
% figure(1)
% imshow(uint8(im))
% figure(2)
% imshow(uint8(cover))
% figure(3)
% imshow(uint8(tcip))
% figure(4)
% imshow(uint8(cip))
% figure(5)
% imshow(uint8(rim))
% % 直方图
% H_IM = plothist( im,11 );
% H_CO = plothist( cover,12 );
% H_NC = plothist( tcip,13 );
% H_C = plothist( cip,14 );
% % 相邻像素相关性
% [ch_im,cv_im,cd_im]=plotcoor( im,21 );
% [ch_co,cv_co,cd_co]=plotcoor( cover,22 );
% [ch_nc,cv_nc,cd_nc]=plotcoor( tcip,23 );
% [ch_c,cv_c,cd_c]=plotcoor( cip,24 );
% % 三通道相关性
% [crg_im, crb_im, cgb_im]=plotchcoor( im,31 );
% [crg_co, crb_co, cgb_co]=plotchcoor( cover,32 );
% [crg_nc, crb_nc, cgb_nc]=plotchcoor( tcip,33 );
% [crg_c, crb_c, cgb_c]=plotchcoor( cip,34 );



% 密钥敏感性

% im = double(imread('date/color_image_512/avion.ppm'));
% cover = double(imread('date/color_image_512/baboon.ppm'));
% 
% key = '8d5ab8ba5340fce4420829ad5d12a0e45dacb0858544163d04c1d02b73e3697d';
% kt='haar';
% 
% [m,n,k]=size(im);
% 
% [cip,tcip,dnkey,max_x1,min_x1,index] = encryption(im,cover,key,kt);
% dnkey(20)=dnkey(20)+10^-14;  % 1,10,11,20
% [rim] = dencryption(cip,cover,dnkey,kt,max_x1,min_x1,index,[m,n,k]);
% 
% figure(1)
% imshow(uint8(im))
% figure(2)
% imshow(uint8(cover))
% figure(3)
% imshow(uint8(tcip))
% figure(4)
% imshow(uint8(cip))
% figure(5)
% imshow(uint8(rim))

%% 鲁棒性

% im = double(imread('date/color_image_512/avion.ppm'));
% cover = double(imread('date/color_image_512/baboon.ppm'));
% 
% 
% key = '8d5ab8ba5340fce4420829ad5d12a0e45dacb0858544163d04c1d02b73e3697d';
% kt='spa';
% % kt='haar';
% 
% [m,n,k]=size(im);
% 
% [cip,tcip,dnkey,max_x1,min_x1,index] = encryption(im,cover,key,kt);
% [rim] = dencryption(cip,cover,dnkey,kt,max_x1,min_x1,index,[m,n,k]);
% 
% % 噪声
% psnr_all = zeros(1, 11);
% ssim_all = zeros(1, 11);
% for i=0:10
%     noise_density = 0.0005 * 2*i;
%     cip1 = double(imnoise(uint8(cip),'salt & pepper',noise_density));
%     [rim1] = dencryption(cip1,cover,dnkey,kt,max_x1,min_x1,index,[m,n,k]);
%     psnr_all(i+1) = psnr(double(uint8(im)),double(uint8(rim1)),255);
%     ssim_all(i+1) = ssim(double(uint8(im)),double(uint8(rim1)));
%     
%     save_rim1_path = sprintf('lu/sp/rim1_noise_%d.png', i);
%     imwrite(uint8(rim1), save_rim1_path);
%     save_cip1_path = sprintf('lu/sp/cip1_noise_%d.png', i);
%     imwrite(uint8(cip1), save_cip1_path);
% end
% 
% xlable = (0:10)*0.0005*2;
% figure('Position', [100, 100, 700, 400]);
% 
% yyaxis left
% h1 = plot(xlable, psnr_all, '-o', 'LineWidth', 2, 'MarkerSize', 6, ...
%      'MarkerFaceColor', 'b', 'Color', [0, 0.4470, 0.7410]);
% ylabel('PSNR (dB)', 'FontSize', 14, 'FontName', 'Times New Roman');
% xlabel('S&P noise density', 'FontSize', 14, 'FontName', 'Times New Roman');
% ylim([min(psnr_all)-5, max(psnr_all)+5]);
% set(gca, 'YColor', [0, 0.4470, 0.7410]);
% 
% % 右Y轴 - SSIM
% yyaxis right
% h2 = plot(xlable, ssim_all, '-s', 'LineWidth', 2, 'MarkerSize', 6, ...
%      'MarkerFaceColor', 'r', 'Color', [0.8500, 0.3250, 0.0980]);
% ylabel('SSIM', 'FontSize', 14, 'FontName', 'Times New Roman');
% ylim([min(ssim_all)-0.05, max(ssim_all)+0.05]);
% set(gca, 'YColor', [0.8500, 0.3250, 0.0980]);
% 
% % 设置坐标轴字体
% set(gca, 'FontSize', 14, 'FontName', 'Times New Roman');
% 
% % 添加网格和图例
% legend([h1, h2], {'PSNR', 'SSIM'}, 'Location', 'best', 'FontSize', 12, 'FontName', 'Times New Roman');

% % 剪切
% psnr_all = zeros(1, 11);
% ssim_all = zeros(1, 11);
% for i=0:10
%     cut_size = uint16(8 * i);
%     cip1 = cip;
%     if cut_size~=0
%         cip1(1:cut_size,1:cut_size,:)=0;
%     end
%     [rim1] = dencryption(cip1,cover,dnkey,kt,max_x1,min_x1,index,[m,n,k]);
%     psnr_all(i+1) = psnr(double(uint8(im)),double(uint8(rim1)),255);
%     ssim_all(i+1) = ssim(double(uint8(im)),double(uint8(rim1)));
%     
%     save_cip1_path = sprintf('lu/cut/cip1_cut_%d.png', i);
%     imwrite(uint8(cip1), save_cip1_path);
%     save_rim_path = sprintf('lu/cut/rim1_cut_%d.png', i);
%     imwrite(uint8(rim1), save_rim_path);
% end
% 
% xlable = (0:10)*8;
% figure('Position', [100, 100, 700, 400]);
% 
% yyaxis left
% h1 = plot(xlable, psnr_all, '-o', 'LineWidth', 2, 'MarkerSize', 6, ...
%      'MarkerFaceColor', 'b', 'Color', [0, 0.4470, 0.7410]);
% ylabel('PSNR (dB)', 'FontSize', 14, 'FontName', 'Times New Roman');
% xlabel('Cut size', 'FontSize', 14, 'FontName', 'Times New Roman');
% ylim([min(psnr_all)-5, max(psnr_all)+5]);
% set(gca, 'YColor', [0, 0.4470, 0.7410]);
% 
% % 右Y轴 - SSIM
% yyaxis right
% h2 = plot(xlable, ssim_all, '-s', 'LineWidth', 2, 'MarkerSize', 6, ...
%      'MarkerFaceColor', 'r', 'Color', [0.8500, 0.3250, 0.0980]);
% ylabel('SSIM', 'FontSize', 14, 'FontName', 'Times New Roman');
% ylim([min(ssim_all)-0.05, max(ssim_all)+0.05]);
% set(gca, 'YColor', [0.8500, 0.3250, 0.0980]);
% 
% % 设置坐标轴字体
% set(gca, 'FontSize', 14, 'FontName', 'Times New Roman');
% 
% % 添加网格和图例
% legend([h1, h2], {'PSNR', 'SSIM'}, 'Location', 'best', 'FontSize', 12, 'FontName', 'Times New Roman');


% cip1 = cip;cip1(256:272,256:272,:)=0;
% cip2 = cip;cip2(256:288,256:288,:)=0;
% dcc1 = 30*abs(cip1-cip);
% dcc2 = 30*abs(cip2-cip);

% cip1 = double(imnoise(uint8(cip),'salt & pepper',0.001));
% cip2 = double(imnoise(uint8(cip),'salt & pepper',0.005));
% dcc1 = abs(cip1-cip);
% for i=2:512
%     for j=2:512
%         if sum(dcc1(i,j,:))~=0
%             dcc1(i-1:i,j-1:j,:)=255;
%         end
%     end
% end
% dcc2 = abs(cip2-cip);
% for i=2:512
%     for j=2:512
%         if sum(dcc2(i,j,:))~=0
%             dcc2(i-1:i,j-1:j,:)=255;
%         end
%     end
% end
% 
% [rim1] = dencryption(cip1,cover,dnkey,kt,max_x1,min_x1,index,[m,n,k]);
% [rim2] = dencryption(cip2,cover,dnkey,kt,max_x1,min_x1,index,[m,n,k]);
% 
% [psnr_rim1, ~] = psnr(im,rim1,255)
% [psnr_rim2, ~] = psnr(im,rim2,255)

% figure(1)
% imshow(uint8(im))
% figure(2)
% imshow(uint8(cover))
% figure(3)
% imshow(uint8(tcip))
% figure(4)
% imshow(uint8(cip))
% figure(41)
% imshow(uint8(30*abs(cip-cip)))
% figure(5)
% imshow(uint8(rim))
% figure(6)
% imshow(uint8(cip1))
% figure(61)
% imshow(uint8(dcc1))
% figure(7)
% imshow(uint8(rim1))
% figure(8)
% imshow(uint8(cip2))
% figure(81)
% imshow(uint8(dcc2))
% figure(9)
% imshow(uint8(rim2))


%% 差分攻击和选择明密文攻击分析
% im = double(imread('date/color_image_512/avion.ppm'));
% cover = double(imread('date/color_image_512/baboon.ppm'));

% im = double(imread('E:\date\ccia_CVG_image\color_image_512/toucan.ppm'));
% cover = double(imread('E:\date\ccia_CVG_image\color_image_512/portofino.ppm'));
% 
% im1 = im;im1(1,1,:)=mod(im1(1,1,:)+1,256);
% im2 = im;im2(128,512,:)=mod(im2(128,512,:)-1,256);
% im3 = im;im3(1,1,:)=im(128,512,:);im3(128,512,:)=im(1,1,:);
% 
% 
% key = '8d5ab8ba5340fce4420829ad5d12a0e45dacb0858544163d04c1d02b73e3697d';
% kt='haar';
% 
% [m,n,k]=size(im);
% 
% 
% [cip,tcip,dnkey,max_x1,min_x1,index] = encryption(im,cover,key,kt);
% [cip1,tcip1,dnkey1,max_x11,min_x11,index1] = encryption(im1,cover,key,kt);
% [cip2,tcip2,dnkey2,max_x12,min_x12,index2] = encryption(im2,cover,key,kt);
% [cip3,tcip3,dnkey3,max_x13,min_x13,index3] = encryption(im3,cover,key,kt);
% 
% [rim] = dencryption(cip,cover,dnkey,kt,max_x1,min_x1,index,[m,n,k]);
% 
% 
% [mssim_tcip1] = ssim(tcip1,tcip);
% [mssim_tcip2] = ssim(tcip2,tcip);
% [mssim_tcip3] = ssim(tcip3,tcip);
% 
% [ NP_tcip1 ] = NPCR( tcip1,tcip );
% [ NP_tcip2 ] = NPCR( tcip2,tcip );
% [ NP_tcip3 ] = NPCR( tcip3,tcip );
% 
% [ UA_tcip1 ] = UACI( tcip1,tcip,674 );
% [ UA_tcip2 ] = UACI( tcip2,tcip,674 );
% [ UA_tcip3 ] = UACI( tcip3,tcip,674 );
% 
% [mssim_cip1] = ssim(cip1,cip);
% [mssim_cip2] = ssim(cip2,cip);
% [mssim_cip3] = ssim(cip3,cip);
% 
% [ NP_cip1 ] = NPCR( cip1,cip );
% [ NP_cip2 ] = NPCR( cip2,cip );
% [ NP_cip3 ] = NPCR( cip3,cip );
% 
% [ UA_cip1 ] = UACI( cip1,cip,255 );
% [ UA_cip2 ] = UACI( cip2,cip,255 );
% [ UA_cip3 ] = UACI( cip3,cip,255 );
% 
% [rim1] = dencryption(cip,cover,dnkey1,kt,max_x1,min_x1,index,[m,n,k]);
% 
% figure(1)
% imshow(uint8(im))
% figure(11)
% imshow(uint8(im1))
% figure(2)
% imshow(uint8(cover))
% figure(3)
% imshow(uint8(tcip))
% figure(4)
% imshow(uint8(cip))
% figure(5)
% imshow(uint8(rim))
% figure(6)
% imshow(uint8(rim1))




