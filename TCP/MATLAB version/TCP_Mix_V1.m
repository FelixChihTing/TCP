%% Two-color pyrometer (TCP): Method selectable version (14 bit camera)
%=========================================================================
% step 1: Method option
% step 2: Input information
% step 3.1: Load calibration images and subtract the dark current
% step 3.2: Extract data from calibration images
% step 3.3: ln(I) chart creation
% step 3.4: ln(R) chart creation (skip for KL method)
% step 4.1: Split image method selection
% step 4.2: Load flame image
% step 4.3: Flame analysis region selection
% step 4.4: Flame image crop
% step 5: Temperature calculation
% ========================================================================
clear
%% STEP 1: PYROMETRY METHOD OPTIONS
MethodTCP = questdlg('Select one pyrometric method.',...
    'Pyrometry Method','Intensity Ratio','KL method','Cancel');
switch MethodTCP
    case 'Intensity Ratio'
        TCPm = 1;
        promptStp1 ={'\fontsize{10} \epsilon_1 of the target:','\fontsize{10} \epsilon_2 of the target:'};
        dlgtitleStp1 ='Emissivity';
        defaultInput1 = {'0.9','0.9'};
        dims = [1 50; 1 50];
        optStp1.Interpreter = 'tex';
        AnsStp1 = inputdlg(promptStp1,dlgtitleStp1,dims,defaultInput1,optStp1);
    case 'KL method'
        TCPm = 2;
        promptStp1 ={'\fontsize{10} \alpha of the target:'};
        dlgtitleStp1 ='Emissivity factor';
        defaultInput1 = {'1.39'};
        dims = [1 50];
        optStp1.Interpreter = 'tex';
        AnsStp1 = inputdlg(promptStp1,dlgtitleStp1,dims,defaultInput1,optStp1);
end
%% STEP 2: INPUT INFORMATION
promptStp23 = {'\fontsize{12} File name of the dark current image:',...
    '\fontsize{12} Maximum temperature of calibration:',...
    '\fontsize{12} Minimum temperature of calibration',...
    '\fontsize{12} Temperature interval:','\fontsize{12} Image format:',...
    '\fontsize{12} Bit depth','\fontsize{12} File name of the flame image:'};
dlgtitleStp23 = 'Calibration Data';
defaultInput = {'darkCurrent','1000','840','10','.png','16','f2'};
dims = [1 70; 1 70; 1 70; 1 70; 1 70; 1 70; 1 70];
optsStp23.Interpreter = 'tex';
AnsStp23 = inputdlg(promptStp23,dlgtitleStp23,dims,defaultInput,optsStp23);
darCFile=strcat(AnsStp23{1},AnsStp23{5});
% define the conversion factor
if str2double(AnsStp23{6}) == 8
    convFt = 255;
elseif str2double(AnsStp23{6}) == 16
    convFt = 16383; % 14bit camera
end
BGNoise=double(imread(darCFile))/convFt; % convert the 16-/8-bit image to a double image
%% STEP 3.1: LOAD CALIBRATION IMAGES 
ImgNum = 1+((str2double(AnsStp23{2})-str2double(AnsStp23{3}))/str2double(AnsStp23{4})); % calculate the total number of images to be loaded
IT=cell([1 ImgNum]);
for n = 1 : ImgNum
    % specify the file to be loaded
    CalifileName = strcat(num2str(str2double(AnsStp23{3})+(n-1)*str2double(AnsStp23{4})),AnsStp23{5});    
    IT(1,n) = {double(imread(CalifileName))/convFt-BGNoise};    % convert to double image and remove dark current
end
%%  STEP 3.2: EXTRACT DATA FROM CALIBRATION IMAGES
% REGION OF INTEREST SELECTION
figure
imagesc(IT{1,ImgNum})   % the last image is supposed to have the brightest illumination
axis(gca, 'image');
colormap('parula')
zoom on
ZoCfdlg=warndlg('Close the dialog to continue.','Zoom In');
waitfor(ZoCfdlg);
hold on
[qheight,qwidth]=size(IT{1,17});
Lw=1;
% Define point of reference for each blackbody target
disp( 'Pick a reference point for 550 nm (left part) filtered image. ;Press "Enter" to confirm.');
[BBr5x,BBr5y]=getpts;
BBr5x=round(BBr5x);
BBr5y=round(BBr5y);
disp( 'Pick a reference point for 750 nm (right part) filtered image. ;Press "Enter" to confirm.');
[BBr7x,BBr7y]=getpts;
BBr7x=round(BBr7x);
BBr7y=round(BBr7y);
% calculate the distance between the two reference points
ReBBL=abs(BBr7x-BBr5x);
%selecting three points on the 550 nm filtered image (left part) to define the rectangular region to be analyzed for
disp( 'Pick a point as the upper left bound of the rectangular analyzing window. ;Press "Enter" to confirm.');
[L5bx,L5by]=getpts;
L5bx=round(L5bx);
L5by=round(L5by);
disp( 'Pick a point as the right bound of the rectangular analyzing window. ;Press "Enter" to confirm.');
[R5bx,R5by]=getpts;
R5bx=round(R5bx);
R5by=round(R5by);
disp( 'Pick a point as the lower bound of the rectangular analyzing window. ;Press "Enter" to confirm.');
[B5bx,B5by]=getpts;
B5bx=round(B5bx);
B5by=round(B5by);
line([L5bx R5bx],[L5by L5by],'Color','w','LineStyle','--','LineWidth',Lw);
line([R5bx R5bx],[L5by B5by],'Color','w','LineStyle','--','LineWidth',Lw);
line([R5bx L5bx],[B5by B5by],'Color','w','LineStyle','--','LineWidth',Lw);
line([L5bx L5bx],[B5by L5by],'Color','w','LineStyle','--','LineWidth',Lw);
% create a rectangular region of interest in the 750nm filtered image (right part)
L7bx=L5bx+ReBBL;
L7by=L5by;
R7bx=R5bx+ReBBL;
B7by=B5by;
line([L7bx R7bx],[L7by L7by],'Color','w','LineStyle','--','LineWidth',Lw);
line([R7bx R7bx],[L7by B7by],'Color','w','LineStyle','--','LineWidth',Lw);
line([R7bx L7bx],[B7by B7by],'Color','w','LineStyle','--','LineWidth',Lw);
line([L7bx L7bx],[B7by L7by],'Color','w','LineStyle','--','LineWidth',Lw);
% Calculate the total number of pixels to be added up in the rectangular region
CaliPix=(abs(R5bx-L5bx)+1)*(abs(B5by-L5by)+1);
%--------------------------------------------------------------------------
% Averaging the indices within the region
IndexT = cell([2,ImgNum]);  % column 1 for 550nm and column 2 for 750nm
IndexT(:,:) = {0};  % making the entries numerical
for filterNum = 1 : 2
    for n = 1: ImgNum        
        if filterNum == 1
            % averaging the index for the 550nm filtered region
            for col = L5bx:R5bx
                for row = L5by:B5by
                    IndexT{filterNum,n} = IndexT{filterNum,n}+ IT{1,n}(row,col)/CaliPix;
                end
            end
        else
            % averaging the index for the 750nm filtered region
            for col = L7bx:R7bx
                for row = L7by:B7by
                    IndexT{filterNum,n}=IndexT{filterNum,n}+IT{1,n}(row,col)/CaliPix;
                end
            end       
        end
    end    
end

%% STEP 3.3: lnI CHARTS CREATION 
% create the temperature array T
T=zeros([1 ImgNum],'double'); 
for n = 1:ImgNum
    T(1,n) = 273+str2double(AnsStp23{2})-(n-1)*str2double(AnsStp23{4}); % in the order of [highest temperature lowest temperature]
end
% create the lnI matrix for 550nm and 750nm
lnI=zeros([2 ImgNum],'double'); % column 1 for 550nm and column 2 for 750nm
for n = 1:ImgNum
    lnI(1,n) = log(IndexT{1,ImgNum-n+1});   % 550 nm, Index for higest temperature is placed at the left
    lnI(2,n) = log(IndexT{2,ImgNum-n+1});   % 750 nm
end
% create lnI-1/T graph
figCalib = figure('PaperType','<custom>','PaperSize',[13 11]);  %set paper size
plot(1./T,lnI(1,:),'o',1./T,lnI(2,:),'o')
hold on
% curve fit for each wavelength
LinearFit_550=fit(1./T.',lnI(1,:).','poly1');  % column vector is assumed for the fit fuction
LinearFit_750=fit(1./T.',lnI(2,:).','poly1');
plot(LinearFit_550,'g-');
plot(LinearFit_750,'r-');
P_550=coeffvalues(LinearFit_550);
P_750=coeffvalues(LinearFit_750);
A_550=P_550(2);
B_550=P_550(1);
A_750=P_750(2);
B_750=P_750(1);
txt_550 = ['ln(I_6_9_0)=',num2str(round(B_550)),'*(1/T) +(',num2str(round(A_550)),')\rightarrow'];
txt_750 = ['ln(I_7_5_0)=',num2str(round(B_750)),'*(1/T) +(',num2str(round(A_750)),')\downarrow'];
text(7.85*10^-4,-4,txt_550,'fontSize',14)
text(8.0*10^-4,-2,txt_750,'fontSize',14)
legend('690 nm','750 nm')
xlabel('1/T (1/K)');
ylabel('ln(I)');
% calculate for effective wavelengths
c2_planck=(1.98644*10^-25)/(1.3806*10^-23); % c2=hc/k
Lamda_eqv1=(-c2_planck/B_550)*10^9;    % in nm
Lamda_eqv2=(-c2_planck/B_750)*10^9;    % in nm
txt_550_wavelength = ['\lambda_e_q_v_1 = ',num2str(round(Lamda_eqv1)),' nm'];
txt_750_wavelength = ['\lambda_e_q_v_2 = ',num2str(round(Lamda_eqv2)),' nm'];
text(7.85*10^-4,-3.7,txt_550_wavelength,'fontSize',14)
text(8.0*10^-4,-1.7,txt_750_wavelength,'fontSize',14)
%%  STEP 3.4: lnR CHART CREATION
if TCPm ==1
    % create the lnR matrix
    lnR=zeros([1 ImgNum],'double'); % column 1 for 550nm and column 2 for 750nm
    for n = 1:ImgNum
        lnR(1,n) = log(IndexT{1,ImgNum-n+1}/IndexT{2,ImgNum-n+1});   % Index for higest temperature is placed at the left
    end
    % create lnR graph
    figCalibR = figure('PaperType','<custom>','PaperSize',[13 11]);  %set paper size
    plot(1./T,lnR(1,:),'o')
    hold on
    % curve fit for each wavelength
    LinearFit_R=fit(1./T.',lnR(1,:).','poly1');  % column vector is assumed for the fit fuction
    plot(LinearFit_R,'g-');
    P_R=coeffvalues(LinearFit_R);
    A_R=P_R(2);
    B_R=P_R(1);
    txt_R = ['ln(R)=',num2str(round(B_R)),'*(1/T) +(',num2str(round(A_R)),')\rightarrow'];
    text(7.85*10^-4,-0.68,txt_R,'fontSize',14)
    legend('Intensity Ratio')
    xlabel('1/T (1/K)');
    ylabel('ln(R)');
end
%% STEP 4.1: SPLIT IMAGE METHOD OPTIONS
MethodSPLT = questdlg('Select one image split method.',...
    'Image Split Method','Central line','Manual selection','Cancel');
switch MethodSPLT
    case 'Central line'
        SPTm = 1;
    case 'Manual selection'
        SPTm = 2;
end
%% STEP 4.2: LOAD FLAME IMAGE
%load the image with flame
    flamFile=strcat(AnsStp23{7},AnsStp23{5});
    IMflame=double(imread(flamFile))/convFt-BGNoise;
    figure
    imagesc(IMflame);
    axis(gca, 'image');
    colormap('parula')
    zoom on
    ZoCfdlg=warndlg('Close the dialog to continue.','Zoom In');
    waitfor(ZoCfdlg);
    hold on
%% STEP 4.3: FLAME ANALYSIS REGION SELECTION
    %selecting three points to define the rectangular region to be cropped for
    %550 nm filtered image (left part)
    disp( 'Pick a point as the upper left bound of the rectangular analyzing window. ;Press "Enter" to confirm.');
    [Lf5bx,Lf5by]=getpts;
    Lf5bx=round(Lf5bx);
    Lf5by=round(Lf5by);
    disp( 'Pick a point as the right bound of the rectangular analyzing window. ;Press "Enter" to confirm.');
    [Rf5bx,Rf5by]=getpts;
    Rf5bx=round(Rf5bx);
    Rf5by=round(Rf5by);
    disp( 'Pick a point as the lower bound of the rectangular analyzing window. ;Press "Enter" to confirm.');
    [Bf5bx,Bf5by]=getpts;
    Bf5bx=round(Bf5bx);
    Bf5by=round(Bf5by);
    line([Lf5bx Rf5bx],[Lf5by Lf5by],'Color','w','LineStyle','--','LineWidth',Lw);
    line([Rf5bx Rf5bx],[Lf5by Bf5by],'Color','w','LineStyle','--','LineWidth',Lw);
    line([Rf5bx Lf5bx],[Bf5by Bf5by],'Color','w','LineStyle','--','LineWidth',Lw);
    line([Lf5bx Lf5bx],[Bf5by Lf5by],'Color','w','LineStyle','--','LineWidth',Lw);

%% STEP 4.4: FLMAE IMAGE CROP
if SPTm == 1
    % Define point of reference for the left field by the bottom left corner
    Fr5x=1;
    Fr5y=qheight;
    % Define point of reference for the right field by the bottom left corner
    Fr7x=round(qwidth/2)-2;
    Fr7y=qheight;
    % calculate the distance between the two reference points
    RefL=abs(Fr7x-Fr5x);
    RefW=Fr7y-Fr5y;
    % create a rectangular region of interest in the 750nm filtered image
    % (right part)
    Lf7bx=Lf5bx+RefL;
    Lf7by=Lf5by+RefW;
    Rf7bx=Rf5bx+RefL;
    Bf7by=Bf5by+RefW;
    line([Lf7bx Rf7bx],[Lf7by Lf7by],'Color','w','LineStyle','--','LineWidth',Lw);
    line([Rf7bx Rf7bx],[Lf7by Bf7by],'Color','w','LineStyle','--','LineWidth',Lw);
    line([Rf7bx Lf7bx],[Bf7by Bf7by],'Color','w','LineStyle','--','LineWidth',Lw);
    line([Lf7bx Lf7bx],[Bf7by Lf7by],'Color','w','LineStyle','--','LineWidth',Lw);
    % CROP IMAGE ----------------------------------------------------------
    pheight=abs(Lf5by-Bf5by)+1;
    pwidth=abs(Rf5bx-Lf5bx)+1;

    F5=zeros(size(pheight,pwidth));
    F7=zeros(size(pheight,pwidth));
    for col = 1:pwidth
        for row = 1:pheight
            F5(row,col)=IMflame(Lf5by-1+row,Lf5bx-1+col);
            F7(row,col)=IMflame(Lf7by-1+row,Lf7bx-1+col);
        end
    end
    figure
    imagesc(F5);
    title '690 nm Filtered Image';
    axis(gca,'image');
    figure
    imagesc(F7);
    title '750 nm Filtered Image';
    axis(gca,'image');
elseif SPTm ==2
    % Define point of reference for each flame, e.g., the position of flame tip
    disp( 'Pick a reference point for 550 nm filtered image (left part). ;Press "Enter" to confirm.');
    [Fr5x,Fr5y]=getpts;
    Fr5x=round(Fr5x);
    Fr5y=round(Fr5y);
    disp( 'Pick a reference point for 750 nm filtered image (right part). ;Press "Enter" to confirm.');
    [Fr7x,Fr7y]=getpts;
    Fr7x=round(Fr7x);
    Fr7y=round(Fr7y);
    % calculate the distance between the two reference points
    RefL=abs(Fr7x-Fr5x);
    RefW=Fr7y-Fr5y;
    % create a rectangular region of interest in the 750nm filtered image (right part)
    Lf7bx=Lf5bx+RefL;
    Lf7by=Lf5by+RefW;
    Rf7bx=Rf5bx+RefL;
    Bf7by=Bf5by+RefW;
    line([Lf7bx Rf7bx],[Lf7by Lf7by],'Color','w','LineStyle','--','LineWidth',Lw);
    line([Rf7bx Rf7bx],[Lf7by Bf7by],'Color','w','LineStyle','--','LineWidth',Lw);
    line([Rf7bx Lf7bx],[Bf7by Bf7by],'Color','w','LineStyle','--','LineWidth',Lw);
    line([Lf7bx Lf7bx],[Bf7by Lf7by],'Color','w','LineStyle','--','LineWidth',Lw);
    % CROP IMAGE-----------------------------------------------------------
    pheight=abs(Lf5by-Bf5by)+1;
    pwidth=abs(Rf5bx-Lf5bx)+1;

    F5=zeros(size(pheight,pwidth));
    F7=zeros(size(pheight,pwidth));
    for col = 1:pwidth
        for row = 1:pheight
            F5(row,col)=IMflame(Lf5by-1+row,Lf5bx-1+col);
            F7(row,col)=IMflame(Lf7by-1+row,Lf7bx-1+col);
        end
    end
    figure
    imagesc(F5);
    title '550 nm Filtered Image';
    axis(gca,'image');
    figure
    imagesc(F7);
    title '750 nm Filtered Image';
    axis(gca,'image');
end
%% STEP 5: TEMPERATURE CALCULATION
if TCPm == 1
    % Flame temperature calculation: step 2. modify the A_R constant by considering emissivity of the target
    A_Rmod =log(exp(A_R)*(str2double(AnsStp1{1})/str2double(AnsStp1{2})));
    % Flame temperature calculation: step 3. T = B_R/(ln(F5/F7)-A_Rmod)
    T2D = zeros(size(pwidth,pheight),'double');
    for row = 1:pheight
        for col = 1:pwidth
            T2D(row,col) = B_R/(log(F5(row,col)/F7(row,col))-A_Rmod);
        end
    end
    figTemp =figure('PaperType','<custom>','PaperSize',[13 11]);  %set paper size
    imagesc(real(T2D),[1000 2800])
    colormap(jet)
    title 'Temperature (K)'
    colorbar;
    axis(gca,'image');
elseif TCPm == 2
    % calculate appearing temperature T_a pixel by pixel
    T_a550=zeros(size(pwidth,pheight));
    T_a750=zeros(size(pwidth,pheight));
    for row = 1:pheight
        for col = 1:pwidth
            T_a550(row,col)=B_550/(real(log(F5(row,col))-A_550));
            T_a750(row,col)=B_750/(real(log(F7(row,col))-A_750));
        end
    end
    figure
    imagesc(T_a550);
    colormap('jet')
    colorbar
    title 'Apparent Temperature of 690 nm (K)';
    axis(gca,'image');
    figure
    imagesc(T_a750);
    colormap('jet')
    colorbar
    title 'Apparent Temperature of 750 nm (K)';
    axis(gca,'image');
    % Flame temperature calculation: step 2. Newton's method
    % define alpha_absor value
    alph_absor = str2double(AnsStp1{1});
    % setting for Newton's method
    promptNewton = {'\fontsize{12} Tolerence:',...
    '\fontsize{12} Maximum number of iteration:'};
    dlgtitle_Newton = 'Calibration Data';
    defaultInput_Newton = {'10','50'};
    dims_Newton = [1 50; 1 50];
    opts_Newton.Interpreter = 'tex';
    Ans_Newton = inputdlg(promptNewton,dlgtitle_Newton,dims_Newton,defaultInput_Newton,opts_Newton);
    % calcualte temperature pixel by pixel by Newton's method
    T2D=zeros(size(pwidth,pheight),'double');
    for row = 1:pheight
        for col = 1:pwidth
            T_a1=T_a550(row,col);
            T_a2=T_a750(row,col);
            TOL=str2double(Ans_Newton{1}); % define tolerence
            N_MaxIt=str2double(Ans_Newton{2});  % define maximum mumber of iterations
            T_ap0=(T_a1+T_a2)/2; % first guess: appealling to intermediate value theorem
            T2D(row,col)=NewtonMethod(Lamda_eqv1,Lamda_eqv2,T_ap0,T_a1,T_a2,TOL,N_MaxIt,alph_absor);
        end
    end
    figTemp =figure('PaperType','<custom>','PaperSize',[13 11]);  %set paper size
    imagesc(real(T2D),[1000 2800])
    colormap(jet)
    title 'Temperature (K)'
    colorbar;
    axis(gca,'image');
end
