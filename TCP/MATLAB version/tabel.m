function [D]=tabel(N);
% Tabel returns D for the Abel inversion given N, the # of projections.
% The final version, the one that agrees with Varghese's code
% Matlab M-file to calculate Abel inversion
% i,j=0..M, M=N-1, where N is number of projections
M=N-1;
% Defining I0 Matrix:
% I0=zeros(N)+10 %debugging statement; remove later
for i=0:M
    for j=0:M
        if j==i & j==0 %first condition
            I0(i+1,j+1)=0;
        end
        if j<i | j==M %first condition also
            I0(i+1,j+1)=0;
        end
        if j==i & j<M & j>0 %third condition
            if j==M-1
            I0(i+1,j+1)=double((1/(2*pi))*log((sqrt(((2*j+2)*(2*j+2)-4*i*i)) + 2*j +...
                2)/(2*j)));
            else
            I0(i+1,j+1)=double((1/(2*pi))*log((sqrt(((2*j+1)*(2*j+1)-4*i*i)) + 2*j +...
                1)/(2*j)));
            end
        end
        if j>i & j<M %fourth condition
            if j==M-1
            I0(i+1,j+1)=double((1/(2*pi))*log((sqrt(((2*j+2)*(2*j+2)-4*i*i)) + 2*j +...
                2)/(sqrt(((2*j-1)*(2*j-1)-4*i*i)) + 2*j - 1)));
            else
            I0(i+1,j+1)=double((1/(2*pi))*log((sqrt(((2*j+1)*(2*j+1)-4*i*i)) + 2*j +...
                1)/(sqrt(((2*j-1)*(2*j-1)-4*i*i)) + 2*j - 1)));
            end
        end
    end
end
showmeI0=I0;

%Defining I1 Matrix:
% I1=zeros(N)+10 %debugging statement; remove later
for i=0:M
    for j=0:M
        if j<i | j==M %first condition also
            I1(i+1,j+1)=0;
        end
        if j==i & j~=M %second condition
            if j==M-1
                I1(i+1,j+1)=double((1/(2*pi))*(sqrt((2*j+2)*(2*j+2)-4*i*i)) - 2*j*...
                    I0(i+1,j+1));
            else
                I1(i+1,j+1)=double((1/(2*pi))*(sqrt((2*j+1)*(2*j+1)-4*i*i)) - 2*j*...
                    I0(i+1,j+1));
            end
        end
        if j>i & j~=M %fourth condition
            if j==M-1
            I1(i+1,j+1)=double((1/(2*pi))*(sqrt((2*j+2)*(2*j+2)-4*i*i) - sqrt((2*j-1)*(2*j-1) -...
                4*i*i)) - 2*j*I0(i+1,j+1));
            else
            I1(i+1,j+1)=double((1/(2*pi))*(sqrt((2*j+1)*(2*j+1)-4*i*i) - sqrt((2*j-1)*(2*j-1) -...
                4*i*i)) - 2*j*I0(i+1,j+1));
            end
        end
    end
end
showmeI1=I1;


%Defining D Matrix:
D=zeros(N)+25;
for i=0:M
    for j=0:M
        if j<i-1 %first condition
            D(i+1,j+1)=0;
        end
        if j==i-1 %second condition
            D(i+1,j+1)=double(I0(i+1,j+2)-I1(i+1,j+2));
        end
        if j==i
            if j<M %third condition
                D(i+1,j+1)=double(I0(i+1,j+2)-I1(i+1,j+2)+2*I1(i+1,j+1));
            else %j=M
                D(i+1,j+1)=double(2*I1(i+1,j+1));
            end
        end
        if j>i
            if j==M %fourth condition
                D(i+1,j+1)=double(2*I1(i+1,j+1)-I0(i+1,j)-I1(i+1,j));
            else
                D(i+1,j+1)=double(I0(i+1,j+2)-I1(i+1,j+2)+2*I1(i+1,j+1)-I0(i+1,j)-...
                    I1(i+1,j));
            end
        end
    end
end
D(1,2)=double(I0(1,3)-I1(1,3)+2*I1(1,2)-2*I1(1,1)); %redefine D(1,2)
%k=5/0

%fid=fopen('8_29co2.dat')
%DATA=fscanf(fid,'%g %g %g %g',[4 inf])
%A = fscanf(fid,format) reads all the data from the file specified by
%fid, converts it according to the specified format string, and returns
%it in matrix A. Argument fid is an integer file identifier obtained
%from fopen. format is a string specifying the format of the data to be
%read.

%fid = fopen('exp.txt');
%a = fscanf(fid,'%g %g',[2 inf]) % It has two rows now.
%a = a';
%fclose(fid)
%DATA = dlmread('8_29co2.dat',',',2,1)
%M = DLMREAD(FILENAME,DLM,R,C) reads data from the DLM-delimited
%file format FILENAME. R and C specify the row R and column C
%where the upper-left corner of the data lies in the file. R and C
%are zero-based so that R=0 and C=0 specifies the first value in the
%file.