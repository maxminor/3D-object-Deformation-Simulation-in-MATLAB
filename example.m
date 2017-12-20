% example code to demonstrate strel3d

%% create a 3D dataset
sw=7;
[y,x]=meshgrid(-sw:sw,-sw:sw);
m=sqrt(x.^2+y.^2);
m=imcomplement(m);
m=m-min(m(:));
ses2=ceil(sw/2);            
m(m<m(ses2,ses2))=0;
m(m>=m(ses2,ses2))=1;

tmp=zeros(15,15,15);
tmp(:,:,6)=m;
tmp(:,:,8)=m;

%% strel vs strel3d comparison for dilation
img1=tmp;
img1=imdilate(img1,strel('disk',5)); % dilation with "strel"

img2=tmp;
img2=imdilate(img2,strel3d(5)); % image dilation with "strel3d"

figure(1); hold on
subplot(3,1,1)
imagesc(squeeze(tmp(7,:,:)))
xlabel('original image','fontsize',14)
axis image
subplot(3,1,2)
imagesc(squeeze(img1(7,:,:)))
xlabel('strel (disk) does NOT connect objects in z-dimension','fontsize',14)
axis image
subplot(3,1,3)
imagesc(squeeze(img2(7,:,:)))
xlabel('strel3d (sphere) does connect objects in z-dimension','fontsize',14)
axis image
colormap(gray)
set(gcf,'color','white')

