clear all; close all; clc;
load Testdata
L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x; 
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);

[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

Ut_av = zeros(n,n,n);
for j=1:20
  Un(:,:,:)=reshape(Undata(j,:),n,n,n);
  Ut_av = Ut_av + fftn(Un);
end

Ut_av = Ut_av/20;
Ut_av = fftshift(Ut_av);
Ut_abs = abs(Ut_av);

[M,I] = max(Ut_abs(:));
%normalize
Ut_abs = Ut_abs/M;

%Plot frequency domain
close all, isosurface(Kx,Ky,Kz,abs(Ut_abs),0.6) 
axis([-20 20 -20 20 -20 20]), grid on, drawnow 
pause(1)


%Get coordinate of peak signal
[a,b,c]=ind2sub(size(Ut_abs),I);
kp = [ks(a),ks(b),ks(c)];


filter = exp(-2*((Ky-kp(1)).^2+(Kx-kp(2)).^2+(Kz-kp(3)).^2));

movement = zeros(20,3);
for j=1:20
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    Ut = fftn(Un);
    Ut = fftshift(Ut);
    unft = Ut.*filter;
    
    unft = ifftshift(unft);
    Unf = ifftn(unft);
    %Unf = ifftshift(Unf);
    [M,I] = max(Unf(:));
    [a,b,c]=ind2sub(size(Unf),I);
    location = [x(a),x(b),x(c)];
    movement(j,:) = location;
    isosurface(X,Y,Z,abs(Unf),0.1), 
    axis([-20 20 -20 20 -20 20]), grid on, drawnow, hold on
end
