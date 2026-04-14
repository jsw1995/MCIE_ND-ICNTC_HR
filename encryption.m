function [cip,NC,dnkey,max_x1,min_x1,index] = encryption(im,cover,key,kt)
%ENCRYPTION2 此处显示有关此函数的摘要
%   此处显示详细说明

[m,n,k]=size(im);

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

%   动态密钥生成
dnkey = dkey( im,key );
a=dnkey(1:10);init=dnkey(11:20);
b=dnkey(21:29);c=dnkey(30:38);
d=dnkey(39:47);e=dnkey(48:56);

%  混沌序列生成
AA = CCP( a,b,c,d,e );
[ Y ] = CM_10D( init,AA,m*n*k );

% tic
%   压缩
r1 = Y(1,1:0.5*m*n*k);
r2 = Y(2,1:0.5*m*n*k);
[ X1,max_x1,min_x1 ] = compre( im,r1,r2,0.5 );

% figure(1)
% h=histogram(X1,256);
% h.FaceColor = [0 0 0];
% h.EdgeColor = 'g';
% % title('重组前')

% 直方图重组
[ X2,index ] = HRE( X1,index_kt );

% figure(2)
% h=histogram(X2,675);
% h.FaceColor = [0 0 0];
% h.EdgeColor = 'g';

% bit置乱
[ NC ] = bit_scram( X2,Y(3:6,:), 1, deba );

% figure(3)
% h=histogram(NC,675);
% h.FaceColor = [0 0 0];
% h.EdgeColor = 'g';
% toc
% tic
% 嵌入
[ cip ] = embed( NC,cover,Y(7:10,:), deba );
% toc
end


function [ dnkey ] = dkey( p,key )
%   DKEY 动态密钥生成
%   p明文，key动态密钥 

x = ones(1,256);
sha_sum = 0;
time = clock;
sha_time = SHA(time,'SHA-256');
% sha_time='825ab8b35340fce4440825ad5d62a0e75dacb0858944163d04c1d02b73e36971';
sha_p = SHA(p,'SHA-256');
for i=1:32  % hex2dec只能到2^52，所以运用循环每8位来一次，也可以其他位数
    tem = ones(1,8);
    tem2 = ones(1,8);
    sn = dec2bin(hex2dec(sha_p((i-1)*2+1:(i-1)*2+2)),8);
    tn = dec2bin(hex2dec(sha_time((i-1)*2+1:(i-1)*2+2)),8);
    kn = dec2bin(hex2dec(key((i-1)*2+1:(i-1)*2+2)),8);
    tem(sn==kn) = 0;
    tem2(char(tem+'0')==tn) = 0;
    x((i-1)*8+1:(i-1)*8+8) = tem2;
    sha_sum = sha_sum + bin2dec(num2str(tem2));
end

a=[];b=[];c=[];init=[];
d=x(1:9);e=x(248:256);
for i = 1:10
    a(i)=exp( (sha_sum/(256*32)+ bin2dec(num2str(x((i-1)*8+1:i*8)))/256) * 5 );
    init(i)= (sha_sum/(256*32)+ bin2dec(num2str(x((i-1)*8+1:i*8)))/256) * sqrt(2);
    if i<=9
        b(i)= (sha_sum/(256*32)+ bin2dec(num2str(x((i-1)*8+1:i*8)))/256) ;
        c(i)=int8( mod( bin2dec(num2str(x((i+22)*8+1:(i+23)*8))), i)+1 );
    end
end

dnkey = [a,init,b,c,d,e];

end

function h = SHA(inp,meth)
% HASH - Convert an input variable into a message digest using any of
%        several common hash algorithms
%
% USAGE: h = hash(inp,'meth')
%
% inp  = input variable, of any of the following classes:
%        char, uint8, logical, double, single, int8, uint8,
%        int16, uint16, int32, uint32, int64, uint64
% h    = hash digest output, in hexadecimal notation
% meth = hash algorithm, which is one of the following:
%        MD2, MD5, SHA-1, SHA-256, SHA-384, or SHA-512
%
% NOTES: (1) If the input is a string or uint8 variable, it is hashed
%            as usual for a byte stream. Other classes are converted into
%            their byte-stream values. In other words, the hash of the
%            following will be identical:
%                     'abc'
%                     uint8('abc')
%                     char([97 98 99])
%            The hash of the follwing will be different from the above,
%            because class "double" uses eight byte elements:
%                     double('abc')
%                     [97 98 99]
%            You can avoid this issue by making sure that your inputs
%            are strings or uint8 arrays.
%        (2) The name of the hash algorithm may be specified in lowercase
%            and/or without the hyphen, if desired. For example,
%            h=hash('my text to hash','sha256');
%        (3) Carefully tested, but no warranty. Use at your own risk.
%        (4) Michael Kleder, Nov 2005
%
% EXAMPLE:
%
% algs={'MD2','MD5','SHA-1','SHA-256','SHA-384','SHA-512'};
% for n=1:6
%     h=hash('my sample text',algs{n});
%     disp([algs{n} ' (' num2str(length(h)*4) ' bits):'])
%     disp(h)
% end

inp=inp(:);
% convert strings and logicals into uint8 format
if ischar(inp) || islogical(inp)
    inp=uint8(inp);
else % convert everything else into uint8 format without loss of data
    inp=typecast(inp,'uint8');
end

% verify hash method, with some syntactical forgiveness:
meth=upper(meth);
switch meth
    case 'SHA1'
        meth='SHA-1';
    case 'SHA256'
        meth='SHA-256';
    case 'SHA384'
        meth='SHA-384';
    case 'SHA512'
        meth='SHA-512';
    otherwise
end
algs={'MD2','MD5','SHA-1','SHA-256','SHA-384','SHA-512'};
if isempty(strcmp(algs,meth))
    error(['Hash algorithm must be ' ...
        'MD2, MD5, SHA-1, SHA-256, SHA-384, or SHA-512']);
end

% create hash
x=java.security.MessageDigest.getInstance(meth);
x.update(inp);
h=typecast(x.digest,'uint8');
h=dec2hex(h)';
if(size(h,1))==1 % remote possibility: all hash bytes  128, so pad:
    h=[repmat('0',[1 size(h,2)]);h];
end
h=lower(h(:)');

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


function [ X3,max_x3,min_x3 ] = compre( image,x1,y1,cr )
%   COMPRE 此处显示有关此函数的摘要
%   此处显示详细说明

x = 1-2*mod(x1 * 10^2 + y1 * 10^2, 1);
y = 1-2*mod(x1 * pi - y1 * exp(1), 1);

[rows, columns, hight] = size(image);
if hight == 1
    X3=compre_gry( image,x(1:int32(rows * cr * rows)),y(1:int32(columns * cr *columns)),cr );
else
    X3R=compre_gry( image(:,:,1), x(1:int32(rows*cr*rows)), y(1:int32(columns*cr*columns)), cr );
    X3G=compre_gry( image(:,:,2), x(int32(rows*cr*rows)+1:int32(2*rows*cr*rows)), y(int32(columns*cr*columns)+1:int32(2*columns*cr*columns)),cr );
    X3B=compre_gry( image(:,:,3), x(int32(2*rows*cr*rows)+1:int32(3*rows*cr*rows)),y(int32(2*columns*cr*columns)+1:int32(3*columns*cr*columns)),cr );
    X3=cat(3,X3R,X3G,X3B);
end

max_x3 = max(X3(:));
min_x3 = min(X3(:));

X3 = round((X3 - min_x3) / (max_x3 - min_x3) * 255);

end
function [ X3 ] = compre_gry( image,x,y,cr )
%   COMPRE 此处显示有关此函数的摘要
%   此处显示详细说明

[rows, columns] = size(image);
X1 = dct2(image);

x = x(1:int32(rows * cr * rows));
y = y(1:int32(columns * cr *columns));
R1 = reshape(x, [int16(rows * cr), rows]);
R2 = reshape(y, [int16(columns * cr), columns]);

% 优化测量矩阵
R1(:, 1:uint16(cr * rows)) = R1(:, 1:uint16(cr * rows)) * 500;
R1 = orth(R1')';
R2(:, 1:uint16(cr * columns)) = R2(:, 1:uint16(cr * columns)) * 500;
R2 = orth(R2')';

X3 = R1 * X1 * R2';

end


function [ X,index ] = HRE( im,index_kt )
%   直方图重组

[hist,~] = imhist(uint8(im));

[~,index]= sort(hist,'descend');
X = im;
for i=1:256
    X(im==index(i))=index_kt(i);
end
% figure(2)
% histogram(X,525)
% title('重组后')

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
function [ NC ] = bit_scram( im,r,type,deba )
%   SCRAM 乱序循环移位
%   此处显示详细说明
a=deba(1);b=deba(2);c=deba(3);d=deba(4); 
CA11 = floor(im/(b*c*d));
CH11 = floor(mod(im,(b*c*d))/(c*d));
CV11 = floor(mod(im,(c*d))/d);
CD11 = mod(im,d);

r1 = r(1,:);
r2 = r(2,:);
r3 = r(3,:);
r4 = r(4,:);

CA11 = scram( CA11,r1,type );
CH11 = scram( CH11,r2,type );
CV11 = scram( CV11,r3,type );
CD11 = scram( CD11,r4,type );

NC = CA11*(b*c*d) + CH11*(c*d) + CV11*d + CD11;

end


function [ cip ] = embed( im,cover,r,deba )
%   embed 自己的嵌入操作(小波变换域)
%   此处显示详细说明

[m, n, k] = size(cover);
[m1,n1,k1] = size(im);
im = reshape(im,[1,m1*n1*k1]);

a=deba(1);b=deba(2);c=deba(3);d=deba(4);
x = r(1,:);  y = r(2,:); z= r(3,:); w= r(4,:);

% 修改封面 (整数小波才需要修改)
max1=252;min1=3;
cover(cover<min1)=min1;
cover(cover>max1)=max1;

% 提升的小波变换
cover=double(cover);
im=double(im);
LS=liftwave('haar','Int2Int');
if k==1
    [CA,CH,CV,CD]=lwt2(cover,LS);
else
    [CAr,CHr,CVr,CDr]=lwt2(cover(:,:,1),LS);
    [CAg,CHg,CVg,CDg]=lwt2(cover(:,:,2),LS);
    [CAb,CHb,CVb,CDb]=lwt2(cover(:,:,3),LS);
    CA = [CAr,CAg,CAb];
    CH = [CHr,CHg,CHb];
    CV = [CVr,CVg,CVb];
    CD = [CDr,CDg,CDb];
end

% 置乱
CA = scram( CA,x,2 );
CH = scram( CH,y,2 );
CV = scram( CV,z,2 );
CD = scram( CD,w,2 );


% 嵌入
CA1 = reshape(CA,[1,int32(m/2*n/2*k)]);
CH1 = reshape(CH,[1,int32(m/2*n/2*k)]);
CV1 = reshape(CV,[1,int32(m/2*n/2*k)]);
CD1 = reshape(CD,[1,int32(m/2*n/2*k)]);


CA11 = floor(im/(b*c*d));
CH11 = floor(mod(im,(b*c*d))/(c*d));
CV11 = floor(mod(im,(c*d))/d);
CD11 = mod(im,d);


CA1(1:m1*n1*k1) = mod(CA11 + mod(CA1(1:m1*n1*k1),a), a) + floor(CA1(1:m1*n1*k1)/a)*a;
CH1(1:m1*n1*k1) = mod(CH11 + mod(CH1(1:m1*n1*k1),b), b) + floor(CH1(1:m1*n1*k1)/b)*b;
CV1(1:m1*n1*k1) = mod(CV11 + mod(CV1(1:m1*n1*k1),c), c) + floor(CV1(1:m1*n1*k1)/c)*c;
CD1(1:m1*n1*k1) = mod(CD11 + mod(CD1(1:m1*n1*k1),d), d) + floor(CD1(1:m1*n1*k1)/d)*d;


CA1 = reshape(CA1,[m/2,k*n/2]);
CH1 = reshape(CH1,[m/2,k*n/2]);
CV1 = reshape(CV1,[m/2,k*n/2]);
CD1 = reshape(CD1,[m/2,k*n/2]);

% 修改
CA2 = correction(CA,CA1,a);
CH2 = correction(CH,CH1,b);
CV2 = correction(CV,CV1,c);
CD2 = correction(CD,CD1,d);

% 反向置乱
CA2 = dscram( CA2,x,2 );
CH2 = dscram( CH2,y,2 );
CV2 = dscram( CV2,z,2 );
CD2 = dscram( CD2,w,2 );

% 合并
if k==1
    cip = ilwt2(CA2,CH2,CV2,CD2,LS);
else
    cipr = ilwt2(CA2(:,1:n/2),CH2(:,1:n/2),CV2(:,1:n/2),CD2(:,1:n/2),LS);
    cipg = ilwt2(CA2(:,n/2+1:2*n/2),CH2(:,n/2+1:2*n/2),CV2(:,n/2+1:2*n/2),CD2(:,n/2+1:2*n/2),LS);
    cipb = ilwt2(CA2(:,2*n/2+1:3*n/2),CH2(:,2*n/2+1:3*n/2),CV2(:,2*n/2+1:3*n/2),CD2(:,2*n/2+1:3*n/2),LS);
    cip = cat(3,cipr,cipg,cipb);
end

end
function [ cip2 ] = correction( cover,cip,e )
%   CORRECTION 修正
%   此处显示详细说明

tem = cip-cover;
cip2 = cip;
cip2(tem < -e/2) = cip(tem < -e/2) + e;
cip2(tem > e/2) = cip(tem > e/2) - e;
end




