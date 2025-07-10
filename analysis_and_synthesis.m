clearvars
close all


% Setup the constants------------------------------------------------------


global f samples frame_size frame_loc frame_count fft_size freq_bin_size;
f = 16000;
samples = 32000;
frame_size = 640;
fft_size = 2^nextpow2(frame_size);
freq_bin_size = f/fft_size;
frame_loc = (frame_size/2 : frame_size/2 : samples - frame_size/2)';
frame_count = length(frame_loc);


% Load and prune data------------------------------------------------------


[x_1, f_rec_1] = audioread('Audiofiles/haiti_jamaica_belize.mp3');
% sound(x_1, f_rec_1);
x_1 = resample(x_1, f, f_rec_1);
audio_1 = x_1(1:samples);

[x_2, f_rec_2] = audioread('Audiofiles/hello_everyone_im.mp3');
% sound(x_2, f_rec_2);
x_2 = resample(x_2, f, f_rec_2);
audio_2 = x_2(1:samples);

audio_mix = audio_1/(max(abs(audio_1))) + 0.5*audio_2/max(abs(audio_2));


% Analysis setup-----------------------------------------------------------


global filter_count f_p;

filter_count = 512;
filter_num = (filter_count:-1:1)';
f_p = (f/2)*(10.^(-0.004.*filter_num));
q_p = (5:5/(filter_count - 1):10)';
bw_p = f_p./q_p;

theta_p = (2*pi/f).*f_p;
r_p = exp(-(pi/f).*bw_p);

b1_p = -2.*r_p.*cos(theta_p);
b2_p = r_p.^2;

num = [1, 0, -1];
den = [ones(filter_count, 1), b1_p, b2_p];


% Filter coefficient approximations----------------------------------------


imp_len = frame_size/2;
impulse = [1; zeros(imp_len-1, 1)];

h_a = zeros(filter_count, imp_len);
h_hat_a = zeros(filter_count, imp_len);
for j = 1:filter_count
    h_a(j, :) = filter(num, den(j, :), impulse);
    h_hat_a(j, :) = fft(h_a(j, :), imp_len);
    h_a(j, :) = h_a(j, :)/max(abs(h_hat_a(j, :)));
end
% h_a = h_a ./ sqrt(sum(h_a.^2, 2));
% h_a = h_a ./ max(abs(sum(h_a, 1)));

h_s = flip(h_a, 2);
% h_s = flip(h_a, 2) ./ sqrt(sum(h_a.^2, 2));

h_as = zeros(filter_count, 2*imp_len-1);
% h_hat_as = zeros(filter_count, 2*imp_len-1);
for j = 1:filter_count
    h_as(j, :) = conv(h_a(j, :), h_s(j, :));
    % h_hat_as(j, :) = fft(h_as(j, :), 2*imp_len-1);
    % h_as(j, :) = h_as(j, :)/max(abs(h_hat_as(j, :)));
end
h_as = h_as ./ max(abs(sum(h_as, 1)));


% Constructing the filter--------------------------------------------------


function y = transfer_func(num, den, freqs)
    n_exp = (0:-1:-numel(num)+1);
    d_exp = (0:-1:-numel(den)+1);
    y = sum(num.*exp(log(freqs(:))*n_exp), 2)./sum(den.*exp(log(freqs(:))*d_exp), 2);
end


% Evaluating impulse responses---------------------------------------------


h_plot = [5, 100, 250, 400, 512];

figure
for j = 1:numel(h_plot)
    h = h_as(h_plot(j), :);
    subplot(numel(h_plot), 1, j);
    plot(h);
    title(['$f_p = ', num2str(f_p(h_plot(j))), '$'], Interpreter="latex");
end

figure
reconstructed_impulse = sum(h_as, 1);
plot([zeros(imp_len, 1); impulse]);
hold on
plot(reconstructed_impulse);
hold off
legend("Original Impulse", "Reconstructed Impulse");


% Magnitude and phase response plots---------------------------------------


filter_plot_res = fft_size;
frequencies = (-filter_plot_res/2:filter_plot_res/2-1)*(f/filter_plot_res);
arg = exp(1i*(2/f)*pi*frequencies);

h = transfer_func(h_a(372, :), 1, arg);
% h = fftshift(fft(h_a(372, :), filter_plot_res));
mag_arr_a = abs(h);
pha_arr_a = angle(h);

h = transfer_func(h_s(372, :), 1, arg);
% h = fftshift(fft(h_s(372, :), filter_plot_res));
mag_arr_s = abs(h);
pha_arr_s = angle(h);

h = transfer_func(h_as(372, :), 1, arg);
% h = fftshift(fft(h_as(372, :), filter_plot_res));
mag_arr_as = abs(h);
pha_arr_as = angle(h);

mag_arr_a = mag2db(mag_arr_a);
mag_arr_s = mag2db(mag_arr_s);
mag_arr_as = mag2db(mag_arr_as);
pha_arr_s = unwrap(pha_arr_s);
pha_arr_as = unwrap(pha_arr_as);

figure
subplot(2, 3, 1)
plot(frequencies(end/2+1:end), mag_arr_a(end/2+1:end));
subplot(2, 3, 4)
plot(frequencies(end/2+1:end), pha_arr_a(end/2+1:end)); 

subplot(2, 3, 2)
plot(frequencies(end/2+1:end), mag_arr_s(end/2+1:end));
subplot(2, 3, 5)
plot(frequencies(end/2+1:end), pha_arr_s(end/2+1:end)); 

subplot(2, 3, 3)
plot(frequencies(end/2+1:end), mag_arr_as(end/2+1:end));
subplot(2, 3, 6)
plot(frequencies(end/2+1:end), pha_arr_as(end/2+1:end)); 


% Testing reconstruction of filter bank------------------------------------


analysis = zeros(length(audio_mix), filter_count);

for j = 1:filter_count
    analysis(:, j) = filter(h_a(j, :), 1, audio_mix);
end

reconstructed_analysis = sum(analysis, 2);
% playblocking(audioplayer(reconstructed_analysis, f));
% audiowrite('Audiofiles/analysis_output.mp3', reconstructed_analysis, f);

s_after_a = zeros(length(audio_mix), filter_count);

for j = 1:filter_count
    s_after_a(:, j) = filter(h_s(j, :), 1, analysis(:, j));
end

reconstructed_s_after_a = sum(s_after_a, 2);
% playblocking(audioplayer(reconstructed_s_after_a, f));


a_and_s = zeros(length(audio_mix), filter_count);

for j = 1:filter_count
    % h_as_shift = circshift(h_as(j, :), -floor(size(h_as, 2)/2));
    h_as_shift = h_as(j, :);
    a_and_s(:, j) = filter(h_as_shift, 1, audio_mix);
end

reconstructed_a_and_s = sum(a_and_s, 2);
% playblocking(audioplayer(audio_mix, f));
% playblocking(audioplayer(reconstructed_a_and_s, f));
% audiowrite('Audiofiles/analysis_synthesis_output.mp3', reconstructed_a_and_s, f);


% Reomving frequencies through filter bank---------------------------------


function removed = filter_bank_removal(filter_bank, wav, freqs, wav_2)
    global filter_count f_p

    bank_select = true(1, filter_count);
    filter_idx = 2;
    for freq = freqs
        while (f_p(filter_idx) < freq) || (f_p(filter_idx-1) > freq)
            filter_idx = filter_idx + 1;
            if filter_idx > 512
                break
            end
        end
        bank_select(max(filter_idx-3, 1)) = false;
        bank_select(max(filter_idx-2, 1)) = false;
        bank_select(max(filter_idx-1, 1)) = false;
        bank_select(min(filter_idx, filter_count)) = false;
        bank_select(min(filter_idx+1, filter_count)) = false;
        bank_select(min(filter_idx+2, filter_count)) = false;
        filter_idx = filter_idx + 1;
        if filter_idx > 512
            break
        end
    end

    filters = 1:filter_count;
    filter_iter = filters(bank_select);

    removed = zeros(length(wav), 1);
    for j = filters
        if ismember(j, filter_iter)
            removed = removed + filter(filter_bank(j, :), 1, wav);
        else
            removed = removed + 2 .* filter(filter_bank(j, :), 1, wav_2);
        end
    end
end


% Frame-wise pitch detection-----------------------------------------------


audio_mix_1 = zeros(samples, 1);

for i = 1:frame_count
    frame_mix = audio_mix(frame_loc(i)-frame_size/2+1 : frame_loc(i)+frame_size/2);
    frame_1 = audio_1(frame_loc(i)-frame_size/2+1 : frame_loc(i)+frame_size/2);
    frame_2 = audio_2(frame_loc(i)-frame_size/2+1 : frame_loc(i)+frame_size/2);
    % freq = amdf(frame_1, 50, 500);
    freq = pitch(frame_2, f, WindowLength=frame_size, OverlapLength=0, Range=[51, 500]);
    % freq = acr(frame_1, 500);

    if freq == 0
        frame_new = frame_mix;
    else
        harmonics = (1:floor(f/(2*freq)))*freq;
        frame_new = filter_bank_removal(h_as, frame_mix, harmonics, frame_1);
    end

    audio_mix_1(frame_loc(i)-(frame_size/2)+1 : frame_loc(i)+frame_size/2)...
            = audio_mix_1(frame_loc(i)-(frame_size/2)+1 : frame_loc(i)+frame_size/2)...
            + frame_new(1:frame_size);
end
% audio_mix_1 = audio_mix_1 ./ max(abs(audio_mix_1));

% playblocking(audioplayer(audio_mix, f));
playblocking(audioplayer(audio_mix_1, f));
% audiowrite('Audiofiles/filter_bank_female_removed.mp3', audio_mix_1, f);

figure
subplot(2, 1, 1)
plot(audio_mix);
subplot(2, 1, 2)
plot(audio_mix_1);