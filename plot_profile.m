function [outputArg1,outputArg2] = plot_profile(model,pat_profile)
% This function will plot patient profiles and the TLE fingerprint on a
% radar plot as shown in Morgan et al. Brain Communications 2022 Figure 5
% subplots
%
% inputs - 
% model = matrix of the TLE fingerprint in 40 x 6 x 2
%   40 are 40 permutations of seizure free patients
%   6 are lobes = prefrontal, parietal, occipital, temporal,
%   motor/somatosensory and subcortical
%   2 are FC and SC
% pat_profile = matrix of the invidual patient profile output from
%   TLE_fingerprint_distance.m matrix in 6 x 2
%   6 are lobes = prefrontal, parietal, occipital, temporal,
%   motor/somatosensory and subcortical
%   2 types of connectivity are FC and SC

% created by VL Morgan April 2022

% set some colors
% col_epb=[0 118 192]./256; % blue
col_epb=[0 138 255]./256; % blue
% col_epr=[163 1 52]./256;  % red
col_epr=[161 0 64]./256;  % red
% col_epg=[103 119 24]./256; % green
col_epg=[84 115 6]./256; % green
% col_epo=[227 124 29]./256; %orange
col_epo=[235 99 0]./256; %orange
col_epy=[255 227 89]./256; % new yellow
col_epgb=[0 89 64]./256; % new green_blue
col_w=[1 1 1];
col_b=[0 0 0];

% create plot model data and close the gap
clear mean_model_permFCSC_lobe  model_perm_FCSC_lobe std_model_permFCSC_lobe
plot_model(:,1:6)=model(:,1:6,1);
plot_model(:,7:12)=model(:,1:6,2);;
plot_model(:,13)=model(:,1,1);

mean_plot_model=mean(plot_model,1)';

% create plot pat using pat_profile from TLE_fingerprint_distance and close
% the gap
plot_pat(1:6)=pat_profile(1:6,1);
plot_pat(7:12)=pat_profile(1:6,2);
plot_pat(13)=pat_profile(1,1);


figure; 
% radar plot
rad_max=40;
rad_min=-40;

clear tht
tht=[0:2*pi/12:2*pi-(2*pi/12) 0]+deg2rad(15);%tht_rot=tht;

% format image
clear tht2
tht2=[0:2*pi/360:2*pi-(2*pi/360) 0]+deg2rad(15);%tht_rot=tht;
polarplot([tht2],(ones(length(tht2),1).*rad_max),'linew',3,'Color',col_w); % edge line
rlim([rad_min rad_max]); 
rticks([rad_min:20:rad_max])
hold
box off
set(gcf,'color','w');
set(findobj(gca,'type','line'),'linew',2)

ax=gca;
ax.ThetaZeroLocation='top';
deg=tht./pi*180;
ax.ThetaTickLabel={'FC pref', 'FC par','FC occ','FC temp','FC mot/som','FC sub','SC sub','SC mot/som',...
    'SC temp','SC occ','SC par','SC pref'};
ax.ThetaTick=[deg(1:12)];
ax.RAxisLocation=270; %
ax.FontName='Arial';
ax.FontSize=14;
ax.FontWeight='bold';
ax.GridColor='k';
ax.LineWidth=2;

% plot zero line
polarplot([tht2],(zeros(length(tht2),1)),'linew',2,'Color',col_w);
polarplot([tht2],(zeros(length(tht2),1)),'--','linew',1,'Color',col_b);

% plot model
polarplot(tht,mean_plot_model,'linew',4,'Color',col_b); % 

% plot pat
polarplot(tht,plot_pat,'linew',4,'Color',col_epg); % 


end

