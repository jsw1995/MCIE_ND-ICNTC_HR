function [rim] = dencryption(cip,cover,dnkey,kt,max_x1,min_x1,index,im_size)
%DENCRYPTION 此处显示有关此函数的摘要
%   此处显示详细说明

m=im_size(1);n=im_size(2);k=im_size(3);

% 根据kt提取出分解基与index—kt
tem = load(strcat('tem/',kt , '_tem.mat'));
tem = getfield(tem, 'tem2');
b1=2*max(abs(tem(:,1)))+1;
b2=2*max(abs(tem(:,2)))+1;
b3=2*max(abs(tem(:,3)))+1;
b4=2*max(abs(tem(:,4)))+1;
deba = [b1,b2,b3,b4];  % 分解基

tem(tem(:,1)<0,1) = tem(tem(:,1)<0,1) + b1;
tem(tem(:,2)<0,2) = tem(tem(:,2)<0,2) + b2;
tem(tem(:,3)<0,3) = tem(tem(:,3)<0,3) + b3;
tem(tem(:,4)<0,4) = tem(tem(:,4)<0,4) + b4;
index_kt = tem(:,1)*(b2*b3*b4) + tem(:,2)*(b3*b4) + tem(:,3)*b4 + tem(:,4);

%  密钥
a=dnkey(1:10);init=dnkey(11:20);
b=dnkey(21:29);c=dnkey(30:38);
d=dnkey(39:47);e=dnkey(48:56);

%  混沌序列生成
AA = CCP( a,b,c,d,e );
[ Y ] = CM_10D( init,AA,m*n*k );

tic
% 提取
RNC = extract( cip,cover,Y(7:10,:),[m/2,n/2,k],deba );
toc
tic
%反向bit置乱
RX2 = bit_dscram( RNC,Y(3:6,:),1,deba );

% 直方图反向重组
RX1 = HRRE( RX2,index,index_kt );

% 重建
r1 = Y(1,1:0.5*m*n*k);
r2 = Y(2,1:0.5*m*n*k);
rim=refact( RX1,r1,r2,0.5,max_x1, min_x1 );
toc
end


function [ A ] = CCP( a,b,c,d,e )
%   混沌控制参数
%   a为对角线元素N
%  b其他元素N-1
%  c控制b位置参数N-1个（1~I-1）
%  d控制A1放到哪里（0，1）
%  e控制b放到B中还是C中（0，1）

N=length(a);
A(1,1)=a(1);
for i=2:N
    if d(i-1)<0.5
        A1=A; A2=a(i);
        B=zeros(i-1,1); C=zeros(1,i-1);
    else
        A1=a(i); A2=A;
        B=zeros(1,i-1); C=zeros(i-1,1);
    end
    if e(i-1)<0.5
        C(c(i-1))=b(i-1);
    else
        B(c(i-1))=b(i-1);
    end
    A = [A1,B;
         C,A2];
end

end
function [ y ] = CM_10D( init,A,L )
%   10维混沌映射
%   init初值，A参数矩阵，L迭代次数

y = [];
y(:,1) = init;

for i=1:L+1000
    y(:,i+1)=mod(A * y(:,i) + pi,1);
end

y = y(:,1001:L+1000);

end


function [ rim ] = refact( image,x1,y1,cr,max_x3, min_x3 )
%   REFACT 此处显示有关此函数的摘要
%   此处显示详细说明

x = 1-2*mod(x1 * 10^2 + y1 * 10^2, 1);
y = 1-2*mod(x1 * pi - y1 * exp(1), 1);

X3 = image / 255 * (max_x3 - min_x3) + min_x3;

[m, n, hight] = size(image);
rows = m/cr;
columns = n/cr;
if hight==1
   rim = refact_gry( X3,x(1:int32(rows*cr*rows)), y(1:int32(columns*cr*columns)),cr );
else
   rimr = refact_gry( X3(:,:,1), x(1:int32(rows*cr*rows)), y(1:int32(columns*cr*columns)),cr );
   rimg = refact_gry( X3(:,:,2), x(int32(rows*cr*rows)+1:int32(2*rows*cr*rows)), y(int32(columns*cr*columns)+1:int32(2*columns*cr*columns)),cr );
   rimb = refact_gry( X3(:,:,3), x(int32(2*rows*cr*rows)+1:int32(3*rows*cr*rows)), y(int32(2*columns*cr*columns)+1:int32(3*columns*cr*columns)),cr );
   rim = cat(3,rimr,rimg,rimb);
end

end
function [ rim ] = refact_gry( image,x,y,cr )
%   REFACT 此处显示有关此函数的摘要
%   此处显示详细说明

X3 = image;

[m,n] = size(image);
rows = m/cr;
columns = n/cr;

R1 = reshape(x, [int16(rows * cr), rows]);
R2 = reshape(y, [int16(columns * cr), columns]);

% 优化测量矩阵
R1(:, 1:int16(cr * rows)) = R1(:, 1:int16(cr * rows)) * 500;
R1 = orth(R1')';
R2(:, 1:int16(cr * columns)) = R2(:, 1:int16(cr * columns)) * 500;
R2 = orth(R2');

% 重构
[ rec ] = nsl0_2d(X3, R1, R2);

rim = idct2(rec);
rim = double(uint8(rim));

end
function [ rim ] = nsl0_2d( y,A,B )
%   NSL0_2D 此处显示有关此函数的摘要
%   此处显示详细说明

sigma_min = 0.0001;
sigma_decrease_factor = 0.1;  
ksai = 0.01;

A_pinv = pinv(A);
B_pinv = pinv(B);

s = A_pinv * y * B_pinv;

sigma = 4 * max(max(abs(s)));
r = 0;
r0 = y - A * s * B;

while (sigma>sigma_min)

    if sum(sum((r-r0).^2)) < ksai  
%         sum(sum((r-r0).^2))
        d = -(sigma^2 * s) ./ (s.*s + sigma^2);
        s = s + d;
        s = s - A_pinv * (A * s * B - y) * B_pinv;
        r0 = y - A * s * B;
        
    end

    sigma = sigma * sigma_decrease_factor;
 
end

rim = s;

end


function [ rim ] = HRRE( X,index,index_kt )
%   反向直方图重组

rim = zeros(size(X));
for i=1:256
    rim(X==index_kt(i))=index(i);
end

end


function [ cip ] = scram( im,r,type )
%   SCRAM 乱序循环移位
%   此处显示详细说明

[m1,n1,k1]=size(im);
m=m1; n=k1*n1;
im=reshape(im,[m,n]);
cip = im;
if type==1
    [~,r1] = sort(r(1:m+n));
    [~,r2] = sort(r(m+n+1:2*(m+n)));
    v1 = floor(mod(r(1:2*m)*10^12,n));
    v2 = floor(mod(r(2*m+1:2*(m+n))*10^12,m));
    ii=0;jj=0;

    for i=1:m+n
        if r1(i)<=m
           ii=ii+1;
           cip(r1(i),:)=circshift(cip(r1(i),:),[0,v1(ii)]); 
        else
            jj=jj+1;
            cip(:,r1(i)-m)=circshift(cip(:,r1(i)-m),[v2(jj),0]);
        end
    end

    for i=1:m+n
        if r2(i)<=m
           ii=ii+1;
           cip(r2(i),:)=circshift(cip(r2(i),:),[0,v1(ii)]); 
        else
            jj=jj+1;
            cip(:,r2(i)-m)=circshift(cip(:,r2(i)-m),[v2(jj),0]);
        end
    end
else
    [~,r1] = sort(r(1:m+n));
    v1 = floor(mod(r(1:m)*10^12,n));
    v2 = floor(mod(r(m+1:(m+n))*10^12,m));
    ii=0;jj=0;

    for i=1:m+n
        if r1(i)<=m
           ii=ii+1;
           cip(r1(i),:)=circshift(cip(r1(i),:),[0,v1(ii)]); 
        else
            jj=jj+1;
            cip(:,r1(i)-m)=circshift(cip(:,r1(i)-m),[v2(jj),0]);
        end
    end
end

cip = reshape(cip,[m1,n1,k1]);

end
function [ rim ] = dscram( cip,r,type )
%   DSCRAM 解密
%   此处显示详细说明

[m1,n1,k1]=size(cip);
m=m1; n=k1*n1;
rim = cip;
if type==1
    [~,r1] = sort(r(1:m+n));
    [~,r2] = sort(r(m+n+1:2*(m+n)));
    v1 = -floor(mod(r(1:2*m)*10^12,n));
    v2 = -floor(mod(r(2*m+1:2*(m+n))*10^12,m));
    ii=2*m+1; jj=2*n+1;

    for i=(m+n):-1:1
        if r2(i)<=m
            ii=ii-1;
            rim(r2(i),:)=circshift(rim(r2(i),:),[0,v1(ii)]); 
        else
            jj=jj-1;
            rim(:,r2(i)-m)=circshift(rim(:,r2(i)-m),[v2(jj),0]);
        end
    end

    for i=(m+n):-1:1
        if r1(i)<=m
            ii=ii-1;
            rim(r1(i),:)=circshift(rim(r1(i),:),[0,v1(ii)]); 
        else
            jj=jj-1;
            rim(:,r1(i)-m)=circshift(rim(:,r1(i)-m),[v2(jj),0]);
        end
    end
else
    [~,r1] = sort(r(1:m+n));
    v1 = -floor(mod(r(1:m)*10^12,n));
    v2 = -floor(mod(r(m+1:(m+n))*10^12,m));
    ii=m+1; jj=n+1;

    for i=(m+n):-1:1
        if r1(i)<=m
            ii=ii-1;
            rim(r1(i),:)=circshift(rim(r1(i),:),[0,v1(ii)]); 
        else
            jj=jj-1;
            rim(:,r1(i)-m)=circshift(rim(:,r1(i)-m),[v2(jj),0]);
        end
    end
end
rim = reshape(rim,[m1,n1,k1]);
end
function [ rim ] = bit_dscram( NC,r,type,deba )
%   SCRAM 乱序循环移位
%   此处显示详细说明
a=deba(1);b=deba(2);c=deba(3);d=deba(4); 

CA11 = floor(NC/(b*c*d));
CH11 = floor(mod(NC,(b*c*d))/(c*d));
CV11 = floor(mod(NC,(c*d))/d);
CD11 = mod(NC,d);

r1 = r(1,:);
r2 = r(2,:);
r3 = r(3,:);
r4 = r(4,:);

CA11 = dscram( CA11,r1,type );
CH11 = dscram( CH11,r2,type );
CV11 = dscram( CV11,r3,type );
CD11 = dscram( CD11,r4,type );

rim = CA11*(b*c*d) + CH11*(c*d) + CV11*d + CD11;

end


function [ rim ] = extract( cip,cover,r,im_shape,deba )
%   EXTRACT 自己的嵌入操作(小波变换域)
%   此处显示详细说明

m1=im_shape(1); n1=im_shape(2);k1=im_shape(3);
[m,n,k] = size(cip);

a=deba(1);b=deba(2);c=deba(3);d=deba(4); 
x = r(1,:);  y = r(2,:); z= r(3,:); w= r(4,:);

% 修改封面 (整数小波才需要修改)
max1=252;min1=3;
cover(cover<min1)=min1;
cover(cover>max1)=max1;

% 提升的小波变换
LS=liftwave('haar','Int2Int');
cip=double(cip);
cover=double(cover);
if k==1
    [CA3,CH3,CV3,CD3]=lwt2(cip,LS);
    [CA31,CH31,CV31,CD31]=lwt2(cover,LS);
else
    [CA3r,CH3r,CV3r,CD3r]=lwt2(cip(:,:,1),LS);
    [CA3g,CH3g,CV3g,CD3g]=lwt2(cip(:,:,2),LS);
    [CA3b,CH3b,CV3b,CD3b]=lwt2(cip(:,:,3),LS);
    CA3 = cat(3,CA3r,CA3g,CA3b);
    CH3 = cat(3,CH3r,CH3g,CH3b);
    CV3 = cat(3,CV3r,CV3g,CV3b);
    CD3 = cat(3,CD3r,CD3g,CD3b);
    
    [CA31r,CH31r,CV31r,CD31r]=lwt2(cover(:,:,1),LS);
    [CA31g,CH31g,CV31g,CD31g]=lwt2(cover(:,:,2),LS);
    [CA31b,CH31b,CV31b,CD31b]=lwt2(cover(:,:,3),LS);
    CA31 = cat(3,CA31r,CA31g,CA31b);
    CH31 = cat(3,CH31r,CH31g,CH31b);
    CV31 = cat(3,CV31r,CV31g,CV31b);
    CD31 = cat(3,CD31r,CD31g,CD31b);
end
% 置乱
CA3 = scram( CA3,x,2 );
CH3 = scram( CH3,y,2 );
CV3 = scram( CV3,z,2 );
CD3 = scram( CD3,w,2 );

CA31 = scram( CA31,x,2 );
CH31 = scram( CH31,y,2 );
CV31 = scram( CV31,z,2 );
CD31 = scram( CD31,w,2 );

CA4 = mod(mod(CA3(1:m1*n1*k1),a) - mod(CA31(1:m1*n1*k1),a), a);
CH4 = mod(mod(CH3(1:m1*n1*k1),b) - mod(CH31(1:m1*n1*k1),b), b);
CV4 = mod(mod(CV3(1:m1*n1*k1),c) - mod(CV31(1:m1*n1*k1),c), c);
CD4 = mod(mod(CD3(1:m1*n1*k1),d) - mod(CD31(1:m1*n1*k1),d), d);


rim = CA4*(b*c*d) + CH4*(c*d) + CV4*d + CD4;
rim = reshape(rim,[m1,n1,k1]);

end


