%% Abel Inverse filtered flame images
% Abel inverse for 550 nm filtered image
clear UI5 UT5 UI7 UT7 m5 m7 p5 p7 M5 M7
figure
imagesc(F5);
title '550 nm Filtered Image';
axis(gca, 'image');
hold on
% select the center line of the flame
disp( 'Pick a point along the central line. ;Press "Enter" to confirm.');
[CenF5x,CenF5y]=getpts;
CenF5x=round(CenF5x);
CenF5y=round(CenF5y);
% trimming or padding image
[f5height,f5width]=size(F5);
dif=f5width-2*CenF5x;
if (dif > 0)
    F5abel=zeros(f5height,f5width-dif); % if the left side is shorter, trim the right side
    [f5abelheight,f5abelwidth]=size(F5abel);
    for row = 1:f5abelheight
        for col = 1:f5abelwidth
            F5abel(row,col)=F5(row,col);         
        end
    end
elseif (dif < 0)    % if the right side is shorter, trim the left side
    F5abel=zeros(f5height,f5width+dif);
    [f5abelheight,f5abelwidth]=size(F5abel);
    for row = 1:f5abelheight
        for col = 1:f5abelwidth
            F5abel(row,col)=F5(row,col-dif);            
        end
    end
else
    F5abel=F5;
    [f5abelheight,f5abelwidth]=size(F5abel);
end
close
% show the trimmed image with the selected center line
figure
imagesc(F5abel);
title 'Trimmed 550 nm Filtered Image';
axis(gca, 'image');
% define the new center line coordinates
if (dif >= 0)
    CenFabel5x=CenF5x;
else
    CenFabel5x=CenF5x-abs(dif);
end
line([CenFabel5x CenFabel5x],[1 f5height],'Color','w','LineStyle','-.','LineWidth',Lw);
% define parameters
promptAbel = {'\fontsize{12} \Delta r:'};
dlgtitleAbel = 'Abel Inversion';
defaultInputAbel = {'0.15'};
dimsAbel = [1 70];
optsAbel.Interpreter = 'tex';
AnsAbel = inputdlg(promptAbel,dlgtitleAbel,dimsAbel,defaultInputAbel,optsAbel);
deltar= str2double(AnsAbel{1}); % Set deltar to be the spacing between lines-of-sight.
Nproj5=CenFabel5x;     % the x component of the image/flame center
D550= tabel(Nproj5);    % Create the D matrix based on the number of projections
F550= zeros(Nproj5,1);  % F is a vector with the discrete values of the function being invserted.

v5=double(F5abel(:,1:f5abelwidth));
va5 = v5(:,((end/2)+1):end);
vb5 = v5(:,((end/2)):-1:1);

vc5 = va5';

vd5 = vb5';

for k=1:f5abelheight
for i=1:Nproj5
for j=1:Nproj5
a5=D550(i,j);

z5 = vc5(j,k); 
w5 = vd5(j,k); 

% b is a column vector with your line-of-sight projections.
m5(j,k)=a5*double(z5);
p5(j,k)=a5*double(w5);
% Take the contraction of D with P.
end
F550(i,k)=(1/deltar)*sum(m5(:,k));
M5(i,k)=(1/deltar)*sum(p5(:,k));

end
end
UI5(CenFabel5x:-1:1,:) = M5(1:CenFabel5x,:);
UI5(CenFabel5x+1:f5abelwidth,:) = F550(1:CenFabel5x,:);

figure
subplot(1,2,1)
imagesc(v5);
colormap(jet);
colorbar;
daspect([1 1 1]);   %sets the aspect ratio of the figure as 1:1 rather than automatic
title('550 nm Filtered Image');
subplot(1,2,2)
UT5=UI5';
imagesc(UT5); 
colormap(jet);
colorbar;
daspect([1 1 1]);
title('Abel Inverted Image');
set(gcf, 'Position',  [100, 100, 1000, 390])


% Abel inverse for 750 nm filtered image
figure
imagesc(F7);
title '750 nm Filtered Image';
axis(gca, 'image');
hold on
% select the center line of the flame
%disp( 'Pick a point along the central line. ;Press "Enter" to confirm.');
%[CenF7x,CenF7y]=getpts;
CenF7x=CenF5x;
CenF7y=CenF5y;
% trimming or padding image
[f7height,f7width]=size(F7);
dif=f7width-2*CenF7x;
if (dif > 0)
    F7abel=zeros(f7height,f7width-dif); % if the left side is shorter, trim the right side
    [f7abelheight,f7abelwidth]=size(F7abel);
    for row = 1:f7abelheight
        for col = 1:f7abelwidth
            F7abel(row,col)=F7(row,col);         
        end
    end
elseif (dif < 0)    % if the right side is shorter, trim the left side
    F7abel=zeros(f7height,f7width+dif);
    [f7abelheight,f7abelwidth]=size(F7abel);
    for row = 1:f7abelheight
        for col = 1:f7abelwidth
            F7abel(row,col)=F7(row,col-dif);            
        end
    end
else
    F7abel=F7;
    [f7abelheight,f7abelwidth]=size(F7abel);
end
close
% show the trimmed image with the selected center line
figure
imagesc(F7abel);
title 'Trimmed 750 nm Filtered Image';
axis(gca, 'image');
% define the new center line coordinates
if (dif >= 0)
    CenFabel7x=CenF7x;
else
    CenFabel7x=CenF7x-abs(dif);
end
line([CenFabel7x CenFabel7x],[1 f7height],'Color','w','LineStyle','-.','LineWidth',Lw);
% define parameters
Nproj7=CenFabel7x;     % the x component of the image/flame center
D750= tabel(Nproj7);    % Create the D matrix based on the number of projections
F750= zeros(Nproj7,1);  % F is a vector with the discrete values of the function being invserted.

v7=double(F7abel(:,1:f7abelwidth));
va7 = v7(:,((end/2)+1):end);
vb7 = v7(:,((end/2)):-1:1);

vc7 = va7';

vd7 = vb7';

for k=1:f7abelheight
for i=1:Nproj7
for j=1:Nproj7
a7=D750(i,j);

z7 = vc7(j,k); 
w7 = vd7(j,k); 

% b is a column vector with your line-of-sight projections.
m7(j,k)=a7*double(z7);
p7(j,k)=a7*double(w7);
% Take the contraction of D with P.
end
F750(i,k)=(1/deltar)*sum(m7(:,k));
M7(i,k)=(1/deltar)*sum(p7(:,k));

end
end
UI7(CenFabel7x:-1:1,:) = M7(1:CenFabel5x,:);
UI7(CenFabel7x+1:f7abelwidth,:) = F750(1:CenFabel7x,:);

figure
subplot(1,2,1)
imagesc(v7);
colormap(jet);
colorbar;
daspect([1 1 1]);   %sets the aspect ratio of the figure as 1:1 rather than automatic
title('750 nm Filtered Image');
subplot(1,2,2)
UT7=UI7';
imagesc(UT7); 
colormap(jet);
colorbar;
daspect([1 1 1]);
title('Abel Inverted Image');
set(gcf, 'Position',  [100, 100, 1000, 390])

% calculate appearing temperature T_a pixel by pixel
T_a550=zeros(size(f5abelheight,f5abelwidth));
T_a750=zeros(size(f5abelheight,f5abelwidth));
for row = 1:f5abelheight
    for col = 1:f5abelwidth
        T_a550(row,col)=B_550/(real(log(UT5(row,col))-A_550));
        T_a750(row,col)=B_750/(real(log(UT7(row,col))-A_750));
    end
end
figure
imagesc(T_a550,[300 1600]);
colormap('jet')
colorbar
title 'Abel Inverted Apparent Temperature of 550 nm (K)';
axis(gca,'image');
figure
imagesc(T_a750,[300 1600]);
colormap('jet')
colorbar
title 'Abel Inverted Apparent Temperature of 750 nm (K)';
axis(gca,'image');

%% Flame temperature calculation: step 2. Newton's method
% calcualte temperature pixel by pixel by Newton's method
alph_absor = str2double(AnsStp23{8});
T2D=zeros(size(f5abelheight,f5abelwidth),'double');
for row = 1:f5abelheight
    for col = 1:f5abelwidth
        T_a1=T_a550(row,col);
        T_a2=T_a750(row,col);
        TOL=10; % define tolerence
        N_MaxIt=50;  % define maximum mumber of iterations
        T_ap0=(T_a1+T_a2)/2; % first guess: appealling to intermediate value theorem
        T2D(row,col)=NewtonMethod(Lamda_eqv1,Lamda_eqv2,T_ap0,T_a1,T_a2,TOL,N_MaxIt,alph_absor);
    end
end
figAbelTemp =figure('PaperType','<custom>','PaperSize',[13 11]);  %set paper size
imagesc(real(T2D),[800 1200])
colormap(jet)
title 'Abel Inverted Temperature (K)'
colorbar;
axis(gca, 'image');



