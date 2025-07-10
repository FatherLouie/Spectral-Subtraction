clearvars
close all

f = 16000;
imp_len = 1024;

function y = transfer_func(num, den, z)
    n_exp = (0:-1:-numel(num)+1);
    d_exp = (0:-1:-numel(den)+1);
    y = sum(num.*(z.^n_exp))/sum(den.*(z.^d_exp));
end

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
num_rev = [zeros(1, imp_len), flip(num, 2)];
den_rev = flip(den, 2);
num_as = conv(num, num_rev);
% num_as = num_as./num_as(1);
% num_as = [zeros(1, imp_len), num_as];
den_as = zeros(filter_count, 5);
for j = 1:filter_count
    den_as(j, :) = conv(den(j, :), den_rev(j, :));
    % den_as(j, :) = den_as(j, :)*(den_as(j, end)/den_as(j, 1));
    den_as(j, :) = den_as(j, :)*(num_as(1)/den_as(j, 1));
end

impulse = [1; zeros(imp_len-1, 1)];

h_fir = filter(num, den(367, :), impulse);
h_fir_rev = h_fir(numel(h_fir):-1:1);
h_fir_as = conv(h_fir, h_fir_rev);

figure
subplot(3, 1, 1)
plot(h_fir);
subplot(3, 1, 2)
plot(h_fir_rev);
subplot(3, 1, 3)
plot(h_fir_as);

filter_plot_res = 5000;
frequencies = (1:filter_plot_res/2)*(f/filter_plot_res);

mag_arr_a = zeros(filter_plot_res/2, 1);
pha_arr_a = zeros(filter_plot_res/2, 1);
mag_arr_s = zeros(filter_plot_res/2, 1);
pha_arr_s = zeros(filter_plot_res/2, 1);
mag_arr_as = zeros(filter_plot_res/2, 1);
pha_arr_as = zeros(filter_plot_res/2, 1);

for j = 1:filter_plot_res/2
    arg = exp(1i*(2/f)*pi*frequencies(j));

    h_z_a = transfer_func(h_fir', 1, arg);
    mag_arr_a(j) = abs(h_z_a);
    pha_arr_a(j) = angle(h_z_a);

    h_z_s = transfer_func(h_fir_rev', 1, arg);
    mag_arr_s(j) = abs(h_z_s);
    pha_arr_s(j) = angle(h_z_s);

    h_z_as = transfer_func(h_fir_as', 1, arg);
    mag_arr_as(j) = abs(h_z_as);
    pha_arr_as(j) = angle(h_z_as);
end

mag_arr_a = mag2db(mag_arr_a);
mag_arr_s = mag2db(mag_arr_s);
mag_arr_as = mag2db(mag_arr_as);
pha_arr_s = unwrap(pha_arr_s);
pha_arr_as = unwrap(pha_arr_as);

figure
subplot(2, 1, 1)
plot(frequencies, mag_arr_a);
subplot(2, 1, 2)
plot(frequencies, pha_arr_a); 

figure
subplot(2, 1, 1)
plot(frequencies, mag_arr_s);
subplot(2, 1, 2)
plot(frequencies, pha_arr_s); 

figure
subplot(2, 1, 1)
plot(frequencies, mag_arr_as);
subplot(2, 1, 2)
plot(frequencies, pha_arr_as); 