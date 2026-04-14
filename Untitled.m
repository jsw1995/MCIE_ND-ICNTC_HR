%% 进行最优直方图移位参数选择
clc;clear all;close all;clear;

% 变换域
dwt_typ='haar';
LS=liftwave(dwt_typ,'Int2Int');
ca = round(255*rand(1024,1024)); ch=round(40*rand(1024,1024))-20; cv=round(40*rand(1024,1024))-20; cd=round(40*rand(1024,1024))-20;
[ im ] = ilwt2( ca,ch,cv,cd,LS );

num=1;
sum_tem=[];
tem = [];
for i1=-2:+2
    for i2=+2:-1:-2
        for i3=-2:+2
            for i4=+3:-1:-3
                tem(num,:)=[i1,i2,i3,i4];
                [ rim ] = ilwt2( ca+i1,ch+i2,cv+i3,cd+i4,LS );
                sum_tem(num) = psnr(double(uint8(im)),double(uint8(rim)),255);
                num=num+1;
            end  
        end
    end
end
[~,index]= sort(sum_tem,'descend');
tem = tem(index(1:256),:);
tem2 = tem;

b1=max(tem(:,1))-min(tem(:,1))+1;
b2=max(tem(:,2))-min(tem(:,2))+1;
b3=max(tem(:,3))-min(tem(:,3))+1;
b4=max(tem(:,4))-min(tem(:,4))+1;
deba = [b1,b2,b3,b4];  % 分解基

tem(tem(:,1)<0,1) = tem(tem(:,1)<0,1) + b1;
tem(tem(:,2)<0,2) = tem(tem(:,2)<0,2) + b2;
tem(tem(:,3)<0,3) = tem(tem(:,3)<0,3) + b3;
tem(tem(:,4)<0,4) = tem(tem(:,4)<0,4) + b4;
index_kt = tem(:,1)*(b2*b3*b4) + tem(:,2)*(b3*b4) + tem(:,3)*b4 + tem(:,4);





% 空域
% dwt_typ='haar';
% LS=liftwave(dwt_typ,'Int2Int');
% ca = round(255*rand(1024,1024)); ch=round(255*rand(1024,1024)); cv=round(255*rand(1024,1024)); cd=round(255*rand(1024,1024));
% im = [ca,ch;cv,cd]; 
% 
% num=1;
% sum_tem=[];
% tem = [];
% for i1=-1:+2
%     for i2=-1:+2
%         for i3=-1:+2
%             for i4=-1:+2
%                 tem(num,:)=[i1,i2,i3,i4];
%                 rim = [ca+i1,ch+i2;cv+i3,cd+i4];
%                 sum_tem(num) = psnr(double(uint8(im)),double(uint8(rim)),255);
%                 num=num+1;
%             end  
%         end
%     end
% end
% [~,index]= sort(sum_tem,'descend');
% tem = tem(index(1:256),:);
% tem2 = tem;
% 
% b1=max(tem(:,1))-min(tem(:,1))+1;
% b2=max(tem(:,2))-min(tem(:,2))+1;
% b3=max(tem(:,3))-min(tem(:,3))+1;
% b4=max(tem(:,4))-min(tem(:,4))+1;
% deba = [b1,b2,b3,b4];  % 分解基
% 
% tem(tem(:,1)<0,1) = tem(tem(:,1)<0,1) + b1;
% tem(tem(:,2)<0,2) = tem(tem(:,2)<0,2) + b2;
% tem(tem(:,3)<0,3) = tem(tem(:,3)<0,3) + b3;
% tem(tem(:,4)<0,4) = tem(tem(:,4)<0,4) + b4;
% index_kt = tem(:,1)*(b2*b3*b4) + tem(:,2)*(b3*b4) + tem(:,3)*b4 + tem(:,4);


