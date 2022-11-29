%% Experiment setup
clc; clear; close all;
addpath(genpath(cd));

%% Load origin figure
ori_img = imread('../assets/mountain.jpeg');
ori_img_d = double(ori_img);
sigma_length = min(size(ori_img_d));
color_list = {'red','green','blue'};

%% calculate and plot figure
figure(5);
subplot(2,2,1)
imshow(ori_img,[]);
TitleH = title('(a) an image example');
set(gcf,'color','none'); % background color -> none
set(gca,'color','none'); % axis color -> none
% set(TitleH, 'Position', [0.5, 0.02], ...
%   'VerticalAlignment', 'bottom', ...
%   'HorizontalAlignment', 'center')

for idx = 1:3
    subplot(2,2,idx+1)
    svd_sigma = svd(ori_img_d(:,:,idx));
    plot(1:100,svd_sigma(1:100));
    TitleH = title(['(',char(96+idx),...
            ') ',color_list{idx},' channel']);
    xlabel('the i-th singluar value')
    ylabel('Magnitude')
    set(gcf,'color','none'); % background color -> none
    set(gca,'color','none','xtick',0:10:100); % axis color -> none
    % set(TitleH, 'Position', [0.5, 0.02], ...
    %     'VerticalAlignment', 'bottom', ...
    %     'HorizontalAlignment', 'center')
end

%% saving the results
save_dir = '../output/fig1/';
% setting save dir
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

print(gcf, '-dpdf',[save_dir,'/fig1.pdf'])
print(gcf,[save_dir,'/fig1.jpeg'] ,'-djpeg','-r300')
close all;