function [rad_vox]= vesicle_rad_function(distribution_type, ...
    min_vesicle_rad_vox, max_vesicle_rad_vox, mat_length)
% Christopher Zahasky
% 12/20/2017
% This function provides a random vesicle radius based on a given type of
% distribution, max_radius, min_radius, and mean_radius

% close all,
% clear all
% set(0,'DefaultAxesFontSize',14)
% set(0,'defaultlinelinewidth',2)
%
% % input
% min_vesicle_rad_vox =  2;
% max_vesicle_rad_vox = 9;

% optional input
histogram_bin_width = 0.2;
% mat_length = 5000;
x = linspace(min_vesicle_rad_vox, max_vesicle_rad_vox, mat_length);
rand_mat = rand(mat_length, 1);
mean_rad = mean([max_vesicle_rad_vox,  min_vesicle_rad_vox]);

% option 1, random values across vesicle size range
if distribution_type == 1
    rad_vox = (rand_mat.*(max_vesicle_rad_vox - min_vesicle_rad_vox))+ min_vesicle_rad_vox;
    
%     figure(1)
%     subplot(1,4,1)
%     hist(rad_vox, [min_vesicle_rad_vox: histogram_bin_width: max_vesicle_rad_vox])
%     title('Random radius')
    
% option 2, normal distribution
elseif distribution_type == 2
    sig = std([mean_rad, max_vesicle_rad_vox])/2;
    f =  1/sqrt(2*pi*sig^2).*exp(-(x-mean_rad).^2./(2*sig^2));
    tot = trapz(x,f);
    cdf_f = zeros(length(f)-1,1);
    for i=2:length(f)
        cdf_f(i) = trapz(x(1:i), f(1:i))/tot;
    end
    
    rad_vox = zeros(size(rand_mat));
    
    for i=1:length(cdf_f)-1
        cind = find(rand_mat >= cdf_f(i) & rand_mat < cdf_f(i+1));
        rad_vox(cind) = x(i);
    end
    
    % % optional plot of normal distribution
    % figure
    % subplot(1,2,1)
    % plot(x,f)
    % title('PDF')
    % axis([min_vesicle_rad_vox, max_vesicle_rad_vox, 0, 1])
    % subplot(1,2,2)
    % plot(x,cdf_f, '-ob')
    % title('CDF')
    % axis([min_vesicle_rad_vox, max_vesicle_rad_vox, 0, 1])
    
%     figure(1)
%     subplot(1,4,2)
%     hist(rad_vox, [min_vesicle_rad_vox: histogram_bin_width: max_vesicle_rad_vox])
%     title('Normal distribution radius')
    
% option 3, logrithmic distribution
elseif distribution_type == 3
    sig = 1;
    mu = 0;
    x = linspace(mean_rad-exp(sig/2), max_vesicle_rad_vox, mat_length);
    
    % mu = min_vesicle_rad_vox;
    % sig = std([mean_rad, max_vesicle_rad_vox])/2;
    f =  1./(x-x(1)+eps) .*(1/(sig*sqrt(2*pi))).*exp(-((log(x-x(1)+eps)-mu).^2)./(2.*sig^2));
    tot = trapz(x,f);
    cdf_f = zeros(length(f)-1,1);
    for i=2:length(f)
        cdf_f(i) = trapz(x(1:i), f(1:i))/tot;
    end
    
    rad_vox = zeros(size(rand_mat));
    
    for i=1:length(cdf_f)-1
        cind = find(rand_mat >= cdf_f(i) & rand_mat < cdf_f(i+1));
        rad_vox(cind) = x(i);
    end
    
    % figure
    % subplot(1,2,1)
    % plot(x,f)
    % title('PDF')
    % axis([min_vesicle_rad_vox, max_vesicle_rad_vox, 0, 1])
    % subplot(1,2,2)
    % plot(x,cdf_f, '-ob')
    % title('CDF')
    % axis([min_vesicle_rad_vox, max_vesicle_rad_vox, 0, 1])
%     log_mean = mean(rad_vox)
    
%     figure(1)
%     subplot(1,4,3)
%     hist(rad_vox, [min_vesicle_rad_vox: histogram_bin_width: max_vesicle_rad_vox])
%     title('Lognormal distribution radius')
    
% option 4, bimodal distribution
elseif distribution_type == 4
    x = linspace(min_vesicle_rad_vox, max_vesicle_rad_vox, mat_length);
    mean1 = (max_vesicle_rad_vox-min_vesicle_rad_vox)/4 + min_vesicle_rad_vox;
    mean2 = 3*(max_vesicle_rad_vox-min_vesicle_rad_vox)/4 + min_vesicle_rad_vox;
    sig = std([mean_rad, max_vesicle_rad_vox])/4;
    
    f1 =  1/sqrt(2*pi*sig^2).*exp(-(x-mean1).^2./(2*sig^2));
    f2 = 1/sqrt(2*pi*sig^2).*exp(-(x-mean2).^2./(2*sig^2));
    f = f1 + f2;
    tot = trapz(x,f);
    cdf_f = zeros(length(f)-1,1);
    for i=2:length(f)
        cdf_f(i) = trapz(x(1:i), f(1:i))/tot;
    end
    
    rad_vox = zeros(size(rand_mat));
    
    for i=1:length(cdf_f)-1
        cind = find(rand_mat >= cdf_f(i) & rand_mat < cdf_f(i+1));
        rad_vox(cind) = x(i);
    end
    
%     figure(1)
%     subplot(1,4,4)
%     hist(rad_vox, [min_vesicle_rad_vox: histogram_bin_width: max_vesicle_rad_vox])
%     title('Bimodal distribution radius')
else
    disp(['First input into "vesicle_rad_function" must be ',...
        'the distribution type (1 = random, 2 = normal ', ...
        '3 = lognormal, 4 = bimodal). Distribution input out of range, ',...
        'returned distribution is random.'])
    rad_vox = (rand_mat.*(max_vesicle_rad_vox - min_vesicle_rad_vox))+ min_vesicle_rad_vox;
end