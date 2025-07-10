clearvars
close all


% Setup the constants------------------------------------------------------


global f samples frame_size frame_loc frame_count fft_size freq_bin_size;
f = 16000;
samples = 32000;
frame_size = 320;
% fft_size = 2^nextpow2(frame_size);
fft_size = 512;
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

audio_mix = audio_1/(max(abs(audio_1))) + 0.3*audio_2/max(abs(audio_2));


% Pitch estimation---------------------------------------------------------


function pitch = amdf(x, f_min, f_max)
    global f;
    m_min = round(f/f_max);
    m_max = round(f/f_min);

    delta_m = zeros(1, m_max-m_min+1);
    for m = m_min:m_max
        delta_m(m-m_min+1) = mean(abs(x(1:end-m) - x(m+1:end)));
    end
    minimas = find(delta_m < 0.6*(max(delta_m) + min(delta_m)));
    if (isempty(minimas))
        pitch = f_min;
        return;
    end

    fundamental = zeros(1, length(minimas));
    j = 2;
    while (j < length(minimas)) && (minimas(j) - minimas(j-1) == 1)
        fundamental(j-1) = minimas(j-1);
        j = j+1;
    end
    fundamental(j-1) = minimas(j-1);
    fundamental = fundamental(fundamental ~= 0);

    [~, index] = min(delta_m(fundamental));
    index = index + fundamental(1) - 1;
    pitch = f/(index+m_min-1);
end


function pitch = acr(x, f_max)
    global f;
    m_min = round(f/f_max);
    [acf, lags] = autocorr(x, NumLags= length(x) - 1);
    [m, max_idx] = max(acf(m_min:end));
    pitch = f/lags(max_idx + m_min - 1);
end


% Spectrum modifications---------------------------------------------------


function frame_new = remove_freq(frame_1, frame_2, frame_mix, freq)
    global f freq_bin_size fft_size;
    harmonics = (1:floor(f/(2*freq)))*freq;
    remove_idx = [floor(harmonics/freq_bin_size) - 1, ...
                    floor(harmonics/freq_bin_size), ...
                    ceil(harmonics/freq_bin_size), ...
                    ceil(harmonics/freq_bin_size) + 1];

    remove_idx = remove_idx((remove_idx >= 1) & (remove_idx <= fft_size/2));

    spec_1 = fft(frame_1, fft_size);
    spec_2 = fft(frame_2, fft_size);
    spec_mix = fft(frame_mix, fft_size);

    for i = remove_idx
        spec_mix(i) = spec_mix(i) * 0;
        % spec_mix(i) = max((abs(spec_mix(i)) - 1.2*abs(spec_1(i))), 0).*exp(1i*angle(spec_mix(i)));
        % spec_mix(i) = abs(spec_1(i)).*exp(1i*angle(spec_mix(i)));

        spec_mix(fft_size+2-i) = conj(spec_mix(i));
    end
    frame_new = ifft(spec_mix);
end


audio_mix_1 = zeros(samples, 1);

for i = 1:frame_count
    frame_mix = audio_mix(frame_loc(i)-frame_size/2+1 : frame_loc(i)+frame_size/2);
    frame_1 = audio_1(frame_loc(i)-frame_size/2+1 : frame_loc(i)+frame_size/2);
    frame_2 = audio_2(frame_loc(i)-frame_size/2+1 : frame_loc(i)+frame_size/2);
    % freq = amdf(frame_1, 50, 500);
    freq = pitch(frame_1, f, WindowLength=frame_size, OverlapLength=0, Range=[51, 500]);
    % freq = acr(frame_1, 500);

    if freq == 0
        frame_new = frame_mix;
    else
        frame_new = remove_freq(frame_1, frame_2, frame_mix, freq);
    end

    audio_mix_1(frame_loc(i)-(frame_size/2)+1 : frame_loc(i)+frame_size/2)...
            = audio_mix_1(frame_loc(i)-(frame_size/2)+1 : frame_loc(i)+frame_size/2)...
            + frame_new(1:frame_size);
end

audio_mix_1 = real(audio_mix_1) / max(abs(audio_mix_1));


% Data presentation--------------------------------------------------------


% playblocking(audioplayer(audio_mix, f));
playblocking(audioplayer(audio_mix_1, f));

% figure
% subplot(4, 1, 1);
% plot((1:samples), audio_1, Color='#00c0c0');
% grid on
% subplot(4, 1, 2);
% plot((1:samples), audio_2, Color='#c000c0');
% grid on
% subplot(4, 1, 3);
% plot((1:samples), audio_mix, Color='#8080ff');
% grid on
% subplot(4, 1, 4);
% plot((1:samples), audio_mix_1, Color='#ffa000');
% grid off

% audiowrite('Audiofiles/mix.mp3', audio_mix, f);
% audiowrite('Audiofiles/female_removed_new.mp3', audio_mix_1, f);