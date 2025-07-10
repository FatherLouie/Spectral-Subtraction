clear all
close all

% Constants----------------------------------------------------------------

global f samples frame_size frame_loc frame_count;
f = 16000;
samples = 32000;
frame_size = 320;
frame_loc = [frame_size/2 : frame_size/2 : samples - frame_size/2];
frame_count = length(frame_loc);

[wave, f_rec] = audioread('Audiofiles/ecuador.mp3');
% wave = resample(wave, f, f_rec);
wave = wave(1:samples);
% wave = wave(601:920);


% AMDF---------------------------------------------------------------------


function pitch = amdf(x, f_min, f_max)
    global f;
    m_min = round(f/f_max);
    m_max = round(f/f_min);

    delta_m = zeros(1, m_max-m_min+1);
    for m = [m_min:m_max]
        delta_m(m-m_min+1) = mean(abs(x(1:end-m) - x(m+1:end)));
    end
   
    minimas = find(delta_m < 0.22*(max(delta_m) + min(delta_m)));
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


function pitch = amdf_new(x, f_min, f_max)
    global f;
    m_min = round(f/f_max);
    m_max = round(f/f_min);
    num_samples = numel(x);

    delta_m = zeros(1, m_max-m_min+1);
    for m = [m_min:m_max]
        % delta_m(m-m_min+1) = mean(abs(x(1:end-m) - x(m+1:end)));
        x_shift = [x(m:num_samples); x(1:m-1)];
        delta_m(m-m_min+1) = mean(abs(x_shift - x));
    end
   
    minimas = [];
    for m = 2:(numel(delta_m)-1)
        if (delta_m(m) <= delta_m(m-1)) && (delta_m(m) <= delta_m(m+1))
            minimas(end+1) = m;
        end
    end

    index = 0;
    for j = 1:numel(minimas)
        if delta_m(minimas(j)) < min(delta_m) + 0.75*(max(delta_m) - min(delta_m))
            index = minimas(j);
            break
        end
    end
    % disp(index)
    % disp(delta_m(minimas))
    % [~, indices] = sort(delta_m(minimas));
    % index = min(indices(1:1));
    pitch = f/(m_min+index-1);
end


% Autocorrelation----------------------------------------------------------


function pitch = acr(x, f_max)
    global f;
    m_min = round(f/f_max);
    [acf, lags] = autocorr(x, NumLags= length(x) - 1);
    [m, max_idx] = max(acf(m_min:end));
    pitch = f/lags(max_idx + m_min - 1);
end


% Cepstrum-----------------------------------------------------------------


function pitch = cepstrum(x, f_max)
    global f;
    m_min = round(f/f_max);
    ceps = ifft(log(abs(fft(x))));
    [~, max_idx] = max(abs(ceps(m_min:end - 10)));
    pitch = f/(max_idx + m_min - 1);
end

% All methods--------------------------------------------------------------


f_amdf = zeros(1, frame_count);
f_ncf = zeros(1, frame_count);
f_acr = zeros(1, frame_count);
f_pef = zeros(1, frame_count);
f_cep = zeros(1, frame_count);
f_lhs = zeros(1, frame_count);
f_srh = zeros(1, frame_count);

% freq = acr(wave(frame_loc(100)-frame_size/2+1 : frame_loc(100)+frame_size/2))

for i = 1:frame_count
    frame = wave(frame_loc(i)-frame_size/2+1 : frame_loc(i)+frame_size/2);
    f_amdf(i) = amdf_new(frame, 50, 500);
    f_ncf(i) = pitch(frame, f, WindowLength=frame_size, OverlapLength=0, Method='NCF', Range=[51, 500]);
    f_pef(i) = pitch(frame, f, WindowLength=frame_size, OverlapLength=0, Method='PEF', Range=[51, 500]);
    f_cep(i) = pitch(frame, f, WindowLength=frame_size, OverlapLength=0, Method='CEP', Range=[51, 500]);
    f_lhs(i) = pitch(frame, f, WindowLength=frame_size, OverlapLength=0, Method='LHS', Range=[51, 500]);
    f_srh(i) = pitch(frame, f, WindowLength=frame_size, OverlapLength=0, Method='SRH', Range=[51, 500]);
    f_acr(i) = acr(frame, 500);
end

figure
subplot(5, 1, 1)
plot(frame_loc, f_ncf, Color='#c00000');

subplot(5, 1, 2)
plot(frame_loc, f_pef, Color='#00c000');

subplot(5, 1, 3)
plot(frame_loc, f_cep, Color='#0000c0');

subplot(5, 1, 4)
plot(frame_loc, f_lhs, Color='#ffa000');

subplot(5, 1, 5)
plot(frame_loc, f_srh, Color='#a000ff');

figure
subplot(3, 1, 1)
plot(wave, Color='#787878');
title('Arbitrary speech signal', FontName='Palatino Linotype');
subplot(3, 1, 2)
plot(frame_loc, f_amdf, Color='#a07000');
title('Average magintude difference function', FontName='Palatino Linotype');
subplot(3, 1, 3)
plot(frame_loc, f_acr, Color='#00a088');
title('Autocorrelation function', FontName='Palatino Linotype');