% Use FFT on image and its strain field
close all;
% make fftim
fftim = double(dydxRotated(20:120, 50:300));
smallest = min(min(fftim));
largest = max(max(fftim));
fftim = (fftim+abs(smallest))./abs(largest-smallest); % fftim is of range 0-1.


% subtract out mean to avoid peak in one corner
fftim = fftim - mean(mean(fftim));

range = [min(fftim(:)) max(fftim(:))];
f1 = figure; imshow(fftim, range);

% Perform FFT calculation
ft = fftshift(fft2(fftim));

figure; imshow(abs(ft), []); colormap('gray'); %title('FFT');
figure; imshow(log(abs(ft)), []); colormap('gray'); %title('FFT');

% Perform 1D FFT on a sample row across image

% row = round(size(fftim,1)/2);
% % Show the selected row
% %figure(f1); hold on; plot([1 size(fftim,2)], [row row], 'r-'); hold off
% % Extract the row
% x = fftim((row-20):1:(row+20), :);
% x = mean(x);
% % Display the grey-level profile
% figure; plot(x); title('Horizontal strain variance');
% 
% X = fft(x);
% figure; plot((1:length(X)), abs(X)); %title('Amplitudes as a function of frequency');
% axis([0 40 0 15]);
% xlabel('frequency in image');

% Perform the same test in the y direction
% column = round(size(fftim,2)/2);
% Show the selected column
%figure(f1); hold on; plot([column column], [1, size(fftim,1)], 'r-'); hold off
% Extract the row
% y = fftim(:,(column-50):2:(column+50))';
% y = mean(y);
% Display the grey-level profile
% figure; plot(y); title('Grey-level profile');
% 
% Y = fft(y);
% figure; plot((1:length(Y)), abs(Y)); title('Amplitude as a function of frequency');
% axis([0 40 0 15]);
% xlabel('frequency across image');
% ylabel('Amplitude');