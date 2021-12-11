%% Program Rekonstruksi Circle Phantom 128 Pixel%%
%%Fadilla Sofa Amatullah%%
%%10217012%%

clc  
clear all

%% Inisiasi phantom uji

%Menentukan parameter citra : Jumlah pixel, ukuran pixel, posisi titik 
%pusat, radius dan koefisien atenuasi
nx = 128;                       
ny = 128;                       
ps = 2;                         
cx = [0; 0; 25; -25];           
cy = [0; 0; 0; 0];
radius = [100; 95; 40; 25];     
at = [2; -1.5; -0.2; 0.3];            

x = zeros(1,nx);
y = zeros(ny,1);
xx = zeros(nx,ny);
yy = zeros(nx,ny);
phantom = zeros(nx,ny);

x(1,1) = -nx+(ps/2);
y(1,1) = ny-(ps/2);

for i = 2 : nx
    x(1,i)= x(1,1)+(i-1)*ps;
    y(i,1)= y(1,1)-(i-1)*ps;
end

for i = 1 : nx
    for j = 1 : ny
        xx(i,j) = x(1,j);
        yy(i,j) = y(i,1);
    end
end

% Menentukan nilai gray level/ nilai koefisien atenuasi phantom
for i = 1:size(cx,1)
    sx = cx(i,1);
    sy = cy(i,1);
    rad = radius(i,1);
    a = at(i,1);
    for j=1:nx
        for k=1:ny
            if ((((xx(j,k)-sx)/rad).^2 + ((yy(j,k)-sy)/rad).^2) <= 1);
                if phantom(j,k)==0
                    phantom(j,k) = a;
                else
                    phantom(j,k)=a+phantom(j,k);
                end
            end
        end
    end
end

%Menampilkan gambar phantom
figure(1)
   imagesc(x, y, phantom)               
    colormap('gray')
    title('Circle Phantom')
    xlabel('Position')
    ylabel('Position')
    colorbar

%% Akuisisi data proyeksi (sinogram)
% Penentuan jumlah berkas sinar dan jumlah sudut proyeksi

disp(sprintf('Akuisisi data proyeksi...'));

nr = 128;
dr = 2;
na = nr*2;
r = zeros(1,nr);
rr = zeros(na,nr);
angle = zeros(na,1);
sg = zeros(na,nr);
r(1,1) = -nr + (dr/2);

for i = 2:nr
    r(1,i) = r(1,1)+(i-1)*dr;
end 

for i = 1:na
    angle(i,1) = (i-1)/na*pi;
end

for i = 1:na
    for j = 1:nr
        rr(i,j) = r(1,j);
    end
end

% Penentuan data nilai tiap pixel pada sinogram
for i = 1:size(cx,1)
    sx = cx(i,1);
    sy = cy(i,1);
    rad = radius(i,1);
    a = at(i,1);
    tau = zeros(nr,na);
    for j = 1:na
        for k = 1:nr
            tau(j,k) = sx*cos(angle(j,1)) + sy*sin(angle(j,1));
            if((rr(j,k)-tau(j,k))^2 <= rad^2)
                sg(j,k)=sg(j,k)+a*2*sqrt(rad^2-(rr(j,k)-tau(j,k))^2);
            end
        end
    end
end

% Menampilkan gambar hasil sinogram
figure(2)
   imagesc(angle/pi*180,r, sg')   
    colormap('gray')                 
    title('Sinogram')
    xlabel('Angle')
    ylabel('Ray Position')
    colorbar

%% Rekonstruksi Citra dengan Proyeksi Balik Tanpa Filter
disp(sprintf('\nProyeksi Balik Citra...'));
% Proses proyeksi balik dari data sinogram
sinogram = sg'; 
lamin = zeros(nx,ny);  
  for ia = 1:na
    projection_ia=sinogram(:,ia); 
    projection_smear=repmat (projection_ia,1,128); 
    %fungsi repmat = menyalin baris matrix yang sama untuk semua kolom            
    rot= imrotate(projection_smear', ia*180/256, 'bicubic','crop');
    %fungsi imrotate = memproyeksikan balik data sesuai sudut proyeksi
    lamin=lamin+rot';   
  end

%Kalibrasi nilai gray level
m = max(max(lamin))+abs(min(min(lamin)));
lamin = ((lamin + abs(min(min(lamin))))/m)*2;
  
% Menampilkan hasil citra dengan proyeksi balik
figure(3)
imagesc(x, y, lamin'); colormap('gray'); axis('image')
title('Simple Backprojection Image')
xlabel('mm')  
ylabel('mm')
colorbar

%Hitung Nilai PSNR
input = phantom;
output = lamin';
M = 128;
N = 128;
peakval = 2;
Fungsi(input,output,M,N,peakval);
%% Rekonstruksi Citra dengan Penerapan Filter Ram-Lak
disp(sprintf('\nFilter Ram-Lak...'));
% Perkalian sinogram dalam domain frekuensi dengan nilai frekuensi pada
% filter Ram-Lak

sinogramfiltered=fftshift(fft(sinogram));     
a = length(sinogram);
freq=linspace(-1, 1, a/2).';
Filter = abs(freq);
Filter = repmat(Filter,1,256);
sinogramfilt0 = (real(sinogramfiltered.*Filter));
sinogramfilt=abs(ifft(ifftshift(sinogramfiltered.*Filter))); 

figure(4)
plot(r, real(sinogramfiltered(:,45)));
title('Respon Frekuensi Data Proyeksi Sebelum Difilter')
xlabel('Ray Position')
ylabel('Freq (Hz)')

figure(5)
    plot(r, sinogramfilt0(:,45));
    title('Respon Frekuensi Data Proyeksi Setelah Difilter Ram-Lak')
    xlabel('Ray Position')
    ylabel('Freq (Hz)')

figure(6)
   imagesc(angle/pi*180,r, sinogramfilt)   
    colormap('gray')                 
    title('Sinogram Filtered Ram-Lak')
    xlabel('Angle')
    ylabel('Ray Position')
    colorbar

% Proyeksi balik citra
  bpf_recon = zeros(nx,ny);
  for ia = 1:na
    bpf_ia=sinogramfilt(:,ia);
    bpf_smear=repmat(bpf_ia,1,128);
    rot1= imrotate(bpf_smear', ia*180/256, 'bicubic','crop');   
    bpf_recon=bpf_recon+rot1';
  end
  
  m = max(max(bpf_recon))+abs(min(min(bpf_recon)));
  bpf_recon = ((bpf_recon + abs(min(min(bpf_recon))))/m)*2;
  
% Plot respon frekuensi
figure(7)
plot(freq, Filter);  
title('Respon Filter Ram-Lak')
xlabel('Freq(w)')
ylabel('H(w)')
  
% Menampilkan rekonstruksi citra dengan penerapan Filter Ram-Lak
  figure(8)
  imagesc(x, y, bpf_recon'); colormap('gray'); 
  axis('image')  
  title('Citra dengan Filter Ram-Lak')
  xlabel('Position')
  ylabel('Position')
  colorbar

%Hitung Nilai PSNR
input = phantom;
output = bpf_recon';
M = 128;
N = 128;
peakval = 2;
Fungsi(input,output,M,N,peakval);
%% Rekonstruksi Citra dengan Penerapan Filter Hamming
disp(sprintf('\nFilter Hamming...'));
% Perkalian sinogram dalam domain frekuensi dengan nilai frekuensi pada
% filter Hamming
sinogramfiltered2=fftshift(fft(sinogram));

b = length(sinogram);
hamm = zeros(b/2,1);
freq=linspace(-1, 1, b/2).';
absfreq = abs(freq);
for i=1:b/2
    hamm(i,1)=0.54+0.46*cos(pi*freq(i,1));
end
hamm = hamm.*absfreq;
hamm2 = repmat(hamm,1,256);
sinogramfilt0 = (real(sinogramfiltered2.*hamm2));
sinogramfilt2=abs(ifft(ifftshift(sinogramfiltered2.*hamm2)));  

figure(9)
plot(r, sinogramfilt0(:,45));
title('Respon Frekuensi Data Proyeksi Setelah Difilter Hamming')
xlabel('Ray Position')
ylabel('Frekuensi(Hz)')

% Proyeksi balik citra
bpf_recon2 = zeros(nx,ny);
  for ia = 1:na
    bpf_ia2=sinogramfilt2(:,ia);
    bpf_smear2=repmat(bpf_ia2,1,128);
    rot2= imrotate(bpf_smear2', ia*180/256, 'bicubic','crop');   
    bpf_recon2=bpf_recon2+rot2';
  end
  
% Kalibrasi nilai gray level   
m = max(max(bpf_recon2))+abs(min(min(bpf_recon2)));
bpf_recon2 = ((bpf_recon2 + abs(min(min(bpf_recon2))))/m)*2;
  
% Plot respon frekuensi
figure(10)
plot(freq, hamm);  
title('Respon Filter Hamming')
xlabel('Freq(w)')
ylabel('H(w)')
  
% Menampilkan rekonstruksi citra dengan penerapan Filter Hamming
figure(11)
imagesc(x, y, bpf_recon2'); colormap('gray'); 
axis('image')  
title('Citra dengan Filter Hamming')
xlabel('Position')
ylabel('Position')
colorbar

%Hitung Nilai PSNR
input = phantom;
output = bpf_recon2';
M = 128;
N = 128;
peakval = 2;
Fungsi(input,output,M,N,peakval);
%% Rekonstruksi Citra dengan Penerapan Filter Hanning
disp(sprintf('\nFilter Hanning...'));
% Perkalian sinogram dalam domain frekuensi dengan nilai frekuensi pada
% filter Hanning
sinogramfiltered3=fftshift(fft(sinogram));
c = length(sinogram);
hann = zeros(c/2,1);
freq=linspace(-1, 1, c/2).';
absfreq = abs(freq);
for i=1:c/2
    hann(i,1)=0.5*(1+cos(pi*freq(i,1)));
end
hann = hann.*absfreq;
hann2 = repmat(hann,1,256);
sinogramfilt0 = (real(sinogramfiltered3.*hann2));
sinogramfilt3=abs(ifft(ifftshift(sinogramfiltered3.*hann2)));    

figure(12)
plot(r, sinogramfilt0(:,45));
title('Respon Frekuensi Data Proyeksi Setelah Difilter Hanning')
xlabel('Ray Position')
ylabel('Frekuensi(Hz)')

%Proyeksi balik data sinogram
bpf_recon3 = zeros(nx,ny);
  for ia = 1:na
    bpf_ia3=sinogramfilt3(:,ia);
    bpf_smear3=repmat(bpf_ia3,1,128);
    rot3= imrotate(bpf_smear3', ia*180/256, 'bicubic','crop');   
    bpf_recon3=bpf_recon3+rot3';
  end

% Kalibrasi nilai gray level
m = max(max(bpf_recon3))+abs(min(min(bpf_recon3)));
bpf_recon3 = ((bpf_recon3 + abs(min(min(bpf_recon3))))/m)*2;
  
% Plot respon frekuensi Hanning
figure(13)
plot(freq, hann);  
title('Respon Filter Hanning')
xlabel('Freq(w)')
ylabel('H(w)')

% Menampilkan rekonstruksi citra dengan penerapan Filter Hanning
figure(14)
imagesc(x, y, bpf_recon3'); colormap('gray'); 
axis('image')  
title('Citra dengan Filter Hanning')
xlabel('Position')
ylabel('Position')
colorbar

%Hitung Nilai PSNR
input = phantom;
output = bpf_recon3';
M = 128;
N = 128;
peakval = 2;
Fungsi(input,output,M,N,peakval);
 
function [MSE,PSNR] = Fungsi(input,output,M,N,peakval)
    MSE = sum(sum((input-output).^2))/(M*N);
    PSNR = 10*log10((peakval^2)/MSE);
    fprintf('MSE : %7.2f ', MSE);
    fprintf('\nPSNR : %9.7f dB \n', PSNR);
end