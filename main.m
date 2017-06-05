clear;close all;
%% Truncated SVD Image Compression
global A;

img = imread('images/1.jpg');
if size(size(img),2) == 3
   %if the image is RGB, convert to the gray-scale 
   img = rgb2gray(img); 
end


%Convert image into double precision
A = im2double(img);

%Image size
[m,n]=size(A);

%Display  the input image
figure();imshow(A,[]);

%Reshape the matrix (image)
s=1;
A = reshape(A,[m/s,n*s]);

%Define the number of dominant terms
r = 100;

%Truncated SVD for small-sized images
repeat = 100;
t1=zeros(repeat,1);
t2=zeros(repeat,1);
for i=1:repeat
    tic,
    [U,S,V] = trunc_svd1a(A,r);
    t1(i)=toc;
    tic,
    [U1,S1,V1] = trunc_svd1b(A,r);
    t2(i)=toc;
end

%Mean CPU times
t1 = mean(t1) %CPU time for trunc_svd1a
t2 = mean(t2) %CPU time for trunc_svd1b
 

%Truncated SVD for large-sized images
repeat = 5;
t1=zeros(repeat,1);
t2=zeros(repeat,1);
for i=1:repeat
    tic,
    %A^TA (or AA^T) is not used
    [U2,S2,V2] = trunc_svd2(A,r);
    t1(i)=toc;
    
    tic,
    %A^TA (or AA^T) is used
    [U3,S3,V3] = trunc_svd2b(A,r);
    t2(i)=toc;
end

%Mean CPU times
t1 = mean(t1)
t2 = mean(t2)

%Reconstructed image
A_reconst = U2*S2*V2';

%Calculate MSE (Mean square error)
MSE = sum(sum((A_reconst-A)).^2)/(m*n)

%Reshape reconstructed image
A_reconst = reshape(A_reconst,[m,n]);

%Display reconstructed image
figure();imshow(A_reconst,[]);

%Calculate memory requirement (MR)
MR = (m/s+n*s+1)*r/(m*n)
