clear all
close all
clc

%parameteres defined 
M = 256;
K = log2(M);
nFFT = 512;
cyc_pre_len = 128;
snr =25;
data_length = 384;
limit_temp = sqrt(M);


%read the data- image
img = imread('ictc.jpg');
reshape_img = imresize(img, [512 512]);
img_bin_temp = dec2bin(reshape_img(:))';
img_bin = img_bin_temp(:);

%convert to symbols 
reshape_img_bin = reshape(img_bin,K,length(img_bin)/K)';
symbols = bin2dec(reshape_img_bin);

%modulation qam
dataMod = qammod(symbols,M,'bin');

%zero padding
reshape_dataMod_blocks = reshape(dataMod,data_length,length(dataMod)/data_length);

null_3 = zeros(20,length(reshape_dataMod_blocks));
null_11 = zeros(88,length(reshape_dataMod_blocks));
demux_first = reshape_dataMod_blocks(1:data_length/2,:);
demux_second = reshape_dataMod_blocks(data_length/2+1:data_length,:);

reshape_blocks_zero = [null_3;demux_first;null_11;demux_second;null_3 ];

%ifft
ifft_out = ifft(reshape_blocks_zero);

%cyclic prefix
cyc_temp = ifft_out(end-cyc_pre_len+1:end,:);
cyc_pre = [cyc_temp;ifft_out];
cyc_serial = cyc_pre(:);

%awgn and rayleigh fading

% rayChan = comm.RayleighChannel( 'SampleRate',1e4,'MaximumDopplerShift',100,'PathDelays',[0 2e-4],'AveragePathGains',[0 -9]);
% fadedsig = rayChan(cyc_serial);

awgn_sig = awgn(cyc_serial,snr);
g = exp(-(0:3));
g = g/norm(g);
fad_sig = conv(awgn_sig,g,'same');


%serial to parallel
rx_sig_par = reshape(fad_sig,nFFT+cyc_pre_len,length(fad_sig)/(nFFT+cyc_pre_len));

%remove cyclic prefix
rx_sig = rx_sig_par(cyc_pre_len+1:end,:);

%fft convert
fft_out = fft(rx_sig);

%channel estimation 
G = fft_out(:,1)./reshape_blocks_zero(:,1);
fft_out = fft_out./repmat(G,1,size(fft_out,2));


%zero padding removed
zero_padding_remove = [fft_out(21:212,:);fft_out(301:492,:)];


%demodulation symbols
x_out = zero_padding_remove(:);   

dataDemod = qamdemod(x_out,M,'bin');

desymbols  = dec2bin(dataDemod);   
de_bin = reshape(desymbols',numel(desymbols),1);
gmi_bin = reshape(de_bin,8,numel(de_bin)/8);
gmi = uint8(bin2dec(gmi_bin'));
gmi_im = reshape(gmi,size(reshape_img)); 



%figparts
%input and output images

f1 = figure;
figure(f1);
subplot(1,2,1)
imshow(reshape_img);
title('Input Image')

subplot(1,2,2)
imshow(gmi_im);
title('Output Image')

set(gcf,'position',[10,10,550,440])

%input and output image part ends
%constellation diagram part


constDiagram_modulation = comm.ConstellationDiagram('ShowReferenceConstellation',false,'XLimits',[-limit_temp limit_temp],'YLimits',[-limit_temp limit_temp]);
% constDiagram_modulation(dataMod)
scatterplot(dataMod)

constDiagram_demodulation = comm.ConstellationDiagram('ReferenceConstellation',dataMod,'ShowReferenceConstellation',true);
% constDiagram_demodulation(x_out)    %x_out(1:1280,:)
scatterplot(x_out)

%constellation diagram end
%ber vs snr part 
f2 = figure;
figure(f2);

a = 0;
bit_ratio = ones(1,41);
for snr = 0:40
    a = a+1;
    awgn_sig = awgn(cyc_serial,snr);
    %serial to parallel
    rx_sig_par = reshape(awgn_sig,nFFT+cyc_pre_len,length(fad_sig)/(nFFT+cyc_pre_len));

    %remove cyclic prefix
    rx_sig = rx_sig_par(cyc_pre_len+1:end,:);

    %fft convert
    fft_out = fft(rx_sig);
    
    %zero padding removed
    zero_padding_remove = [fft_out(21:212,:);fft_out(301:492,:)];

    %demodulation symbols
    x_out = zero_padding_remove(:);   
    dataDemod = qamdemod(x_out,M,'bin');
    
    [bit_err,bit_ratio(a)]=biterr(symbols,dataDemod);
end
snr=0:40;
semilogy(snr,bit_ratio,'-ok');
grid;
title('OFDM Bit Error Rate VS Signal To Noise Ratio');
ylabel('BER');
xlabel('SNR [dB]');



