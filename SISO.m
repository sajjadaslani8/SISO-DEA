clc;
clear;
close all;


% x = [ 20  19  25  27  22  55  33  31  30  50  53  38];
% y = [100 150 160 180 94 230 220 152 190 250 260 250];

% x = [ 20  19  25  27  22  55  33  31  30  50  35  38];
% y = [100 150 160 180 94 230 220 152 190 250 260 250];


% x = [ 20  19  25  27  22  55  33  31  30  50  53  83];
% y = [100 150 220 180 94 230 240 152 190 250 260 250];


% x = [20  190  250  28  24  55  330  31  30  50  53  38];
% y = [184 150 160 230 210 230 220 152 190 250 260 250];


x = [ 50	60	20	22	23	30	25	29	24	33	47	53	50	60	80];
y = [130   140	150	100	90	150	170	190	220	240	260	260	250	240	250];

% P0 = [35 250];

% P0 = [29 255];
P0 = [24 218];

figure(1)
data = [x;y];
subplot(1,2,1)
plot(x,y,'b*');
hold on

Nprim = 6;

xmin = min(x);
ymax = max(y);
alpha(1) = 0;
beta(1) = 1/ymax;

for i=2:Nprim
    alpha(i) = 1/xmin;
    beta(i) = 0;
end

i = 3;
Mprim = [x;y];

[ymax, index] = max(Mprim(2,:));
index = find(Mprim(2,:)==ymax);
xmintemp = min(Mprim(1,index));
idx = find(Mprim(1,:)==xmintemp);
Pmax = [Mprim(1,idx);Mprim(2,idx)];

DE(1,1) = Pmax(1);
DE(2,1) = Pmax(2);

while(i<(Nprim+1))&(size(Mprim,2) > 1)
    
    minslope = 1000;
    for k=1:size(Mprim,2)
        slope = (Pmax(2)-Mprim(2,k))/(Pmax(1)-Mprim(1,k));
        if slope > 0
            if slope < minslope
                minslope = slope;
                temp = inv([Mprim(1,k) Mprim(2,k);Pmax(1) Pmax(2)])*[1;1];
                alpha(i) = temp(1);
                beta(i) = temp(2);
                DE(1,i+1) = Mprim(1,k);
                DE(2,i+1) = Mprim(2,k);
                idx2 = k;
            end
        end
    end
    pstar = Mprim(:,idx2);
    
    tempp = find(Mprim(1,:) > pstar(1));
    Mprim(:,tempp)=[];
    
    Pmax = pstar;
    i = i+1;
end

DE(:,find(DE(1,:) == 0))=[];
alphabeta = [alpha;beta];

x = [0:0.1:100];
for i=1:Nprim
    if alphabeta(2,i) == 0
        y = [0:0.1:450];
        plot((1/alphabeta(1,i))*ones(1,size(y,2)),y);
        hold on
    else
        y = (1/alphabeta(2,i))*(1-alphabeta(1,i)*x);
        plot(x,y)
        hold on
    end
end
xlim([0 100])
ylim([0 400])
alphabetapre = alphabeta;


plot(P0(1),P0(2),'r*')
hold on

x0 = P0(1);
y0 = P0(2);
Tp = [];

for i=1:size(DE,2)-1
    
    x1 = DE(1,i);
    y1 = DE(2,i);
    
    x2 = DE(1,i+1);
    y2 = DE(2,i+1);
    
    if y0 > (((y2-y1)/(x2-x1))*(x0-x1) + y1)
        Tp = [Tp DE(:,i) DE(:,i+1)];
    end
end

if x0 < min(DE(1,:))
    idx = find(DE(1,:)==min(DE(1,:)));
    Tp = [Tp DE(:,idx)];
end
if y0 > max(DE(2,:))
    idx = find(DE(2,:)==max(DE(2,:)));
    Tp = [Tp DE(:,idx)];
end

if size(Tp,1) ~= 0    
    [UniXY,Index]=unique(Tp','rows');
    DupIndex=setdiff(1:size(Tp',1),Index);
    Tp(:,DupIndex)=[]
end
Mprim = [Tp P0'];
Mprimpre = Mprim;

if min(Mprim(1,:)) < xmin
    xmin = min(Mprim(1,:));
    alphabeta(1,find(alphabeta(2,:)==0))=1/xmin;
end
alphabeta
NumTalphabeta = 0;
for k = 1 : size(Tp,2)-1
    temp = inv([Tp(1,k) Tp(2,k);Tp(1,k+1) Tp(2,k+1)])*[1;1];
    indxx = find(alphabeta(1,:) == temp(1));
    indxy = find(alphabeta(2,:) == temp(2));
    if indxx == indxy
        alphabeta(1,indxx) = 1/xmin;
        alphabeta(2,indxx) = 0;
        NumTalphabeta = NumTalphabeta + 1;
    end
end
alphabeta

alpha = alphabeta(1,:);
beta = alphabeta(2,:);
sortindex = (-alpha./beta)';
sortindex(find(sortindex==-Inf))=Inf;
temp = [alpha' beta' sortindex];
temp = sortrows(temp,3);
alpha = temp(:,1)';
beta = temp(:,2)';
alphabeta(1,:) = alpha;
alphabeta(2,:) = beta;



[ymax, index] = max(Mprim(2,:));
index = find(Mprim(2,:)==ymax);
xmintemp = min(Mprim(1,index));
idx = find(Mprim(1,:)==xmintemp);
Pmax = [Mprim(1,idx);Mprim(2,idx)];

% n = size(DE,2) + 1;
% i = n - NumTalphabeta + 1;

i = min(find(alphabeta(2,:) == 0));
newalphabetaindex = i;

while(i<(Nprim+1))&(size(Mprim,2) > 1)
    pstar = Mprim(:,find(Mprim(2,:)==max(Mprim(2,:))));
    minslope = 1000;
    idx2 = 0;
    for k=1:size(Mprim,2)
        slope = (Pmax(2)-Mprim(2,k))/(Pmax(1)-Mprim(1,k));
        if slope > 0
            if slope < minslope
                minslope = slope;
                temp = inv([Mprim(1,k) Mprim(2,k);Pmax(1) Pmax(2)])*[1;1];
                alphabeta(1,i) = temp(1);
                alphabeta(2,i) = temp(2);
                %                 DE(1,i+1) = Mprim(1,k);
                %                 DE(2,i+1) = Mprim(2,k);
                idx2 = k;
            end
        end
    end

    if (idx2 ~= 0)

        pstar = Mprim(:,idx2);
        tempp = find(Mprim(1,:) > pstar(1));
        Mprim(:,tempp)=[];
        Pmax = pstar;
    end

    i = i+1;
end
alphabeta

if P0(2)> max(data(2,:))
    alphabeta(1,1) = 0;
    alphabeta(2,1) = 1/P0(2);
end


% end

subplot(1,2,2)
x = [0:0.1:100];
for i=1:Nprim
    if alphabeta(2,i) == 0
        y = [0:0.1:450];
        plot((1/alphabeta(1,i))*ones(1,size(y,2)),y);
        hold on
    else
        y = (1/alphabeta(2,i))*(1-alphabeta(1,i)*x);
        plot(x,y)
        hold on
    end
    
end
% data = [x;y];

plot(data(1,:),data(2,:),'b*');
hold on
plot(P0(1),P0(2),'r*');

% alphabeta
xlim([0 100])
ylim([0 400])


figure(2)
plot(data(1,:),data(2,:),'b*');
hold on

x = [0:0.1:100];
for i=1:Nprim
    if alphabetapre(2,i) == 0
        y = [0:0.1:450];
        plot((1/alphabetapre(1,i))*ones(1,size(y,2)),y);
        hold on
    else
        y = (1/alphabetapre(2,i))*(1-alphabetapre(1,i)*x);
        plot(x,y)
        hold on
    end
end
xlim([0 100])
ylim([0 400])

ii = 1;
title('Algorithm #2')
print('test2','-r600','-dpng');
I2 = imread('test2.png');
frame = im2frame(I2);
F(ii) = frame;
ii = ii+1;

pause(2)

plot(P0(1),P0(2),'r*');
hold on

title('Algorithm #2')
print('test2','-r600','-dpng');
I2 = imread('test2.png');
frame = im2frame(I2);
F(ii) = frame;
ii = ii+1;

pause(2)
% alphabetapretemp = alphabetapre;


for i=1:Nprim
    if (P0(2) > ((1/alphabetapre(2,i))*(1-alphabetapre(1,i)*P0(1))))&(alphabetapre(2,i)~=0)
        indexalbe(i) = 1;
    else
        indexalbe(i) = 0;
    end
       
end

alphabetapre(:,find(indexalbe==1))=[];



figure(3)
plot(data(1,:),data(2,:),'b*');
hold on

xlim([0 100])
ylim([0 400])

x = [0:0.1:100];
for i=1:size(alphabetapre,2)
    if alphabetapre(2,i) == 0
        y = [0:0.1:450];
        plot((1/alphabetapre(1,i))*ones(1,size(y,2)),y);  
        hold on
    else
        y = (1/alphabetapre(2,i))*(1-alphabetapre(1,i)*x);
        plot(x,y)     
        hold on
    end
end

plot(P0(1),P0(2),'r*');
hold on
title('Algorithm #2')
print('test2','-r600','-dpng');
I2 = imread('test2.png');
frame = im2frame(I2);
F(ii) = frame;
ii = ii+1;

pause(2)
x = [0:0.1:100];
for i=newalphabetaindex:newalphabetaindex +(size(Mprimpre,2)-1)-1
    y = (1/alphabeta(2,i))*(1-alphabeta(1,i)*x);
    plot(x,y,'k--')
    title('Algorithm #2')
    print('test2','-r600','-dpng');
    I2 = imread('test2.png');
    frame = im2frame(I2);
    F(ii) = frame;
    ii = ii+1;
    hold on
end

writerObj = VideoWriter('Algorithm2.mp4');
writerObj.FrameRate = 1;
open(writerObj)
for i= 1:ii-1

    writeVideo(writerObj,F(i))

end

close(writerObj)


    
    

