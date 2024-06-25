%% Supplementary Analysis for manuscript "Shifting Phenology as a Key Driver of Shelf Zooplankton Population Variability"
% MATLAB code for the Lotka-Volterra model (Discussion and Figure S10)
% Adapted from model provided by Greg L. Britten
% Isabel Honda
% May 2024

%% Lotka-Volterra Model (Appendix S1, Figure S10)

clc; clear; close all;
dt = 0.1;
bloom_time1 = 90/dt; % bloom timing differs by 20 days
bloom_time2 = 70/dt;
diap_init = 250;

% prescribe phytoplankton as instant bloom and exponential decay
p1 = [zeros(1, bloom_time1-1), 1*exp(-0.002*(0:(365/dt)-bloom_time1))+0.2];
p2 = [zeros(1, bloom_time2-1), 1*exp(-0.002*(0:(365/dt)-bloom_time2))+0.2];

% initial conditions
z1 = zeros(1, 365/dt);
c1 = zeros(1, 365/dt);
z2 = zeros(1, 365/dt);
c2 = zeros(1, 365/dt);
z1(1) = 0.1;
c1(1) = 0.1;
z2(1) = 0.1;
c2(1) = 0.1;

z_g = 0.03; % zooplankton grazing rate on phytoplankton
c_g = 0.03; % consumer grazing rate on zooplankton
c_m = 0.001; % consumer mortality
c_emerge = 120/dt; % predators appear at day 120

% biomass evolution for type 1 functional response
% assumes same grazing rate for calanus and top predator
for i = 2:length(p1)
    if i > c_emerge
        z1(i) = z1(i-1) + z_g*z1(i-1)*p1(i-1)*dt - c_g*c1(i-1)*z1(i-1)*dt;
        z2(i) = z2(i-1) + z_g*z2(i-1)*p2(i-1)*dt - c_g*c2(i-1)*z2(i-1)*dt;
        c1(i) = c1(i-1) + c_g*c1(i-1)*z1(i-1)*dt - c_m*c1(i-1)*dt;
        c2(i) = c2(i-1) + c_g*c2(i-1)*z2(i-1)*dt - c_m*c2(i-1)*dt;
    else
        z1(i) = z1(i-1) + z_g*z1(i-1)*p1(i-1)*dt;
        z2(i) = z2(i-1) + z_g*z2(i-1)*p2(i-1)*dt;
        c1(i) = c1(i-1);
        c2(i) = c2(i-1);
    end
end

figure(1)
clf
set(gcf,'color','w');

time = (1:(365/dt))*dt;
subplot(3,1,1)
plot(time, p1,'g','LineWidth',2)
hold on
plot(time, p2,'g--','LineWidth',2)
ylabel('P')
legend('Reference', 'Earlier Bloom', 'Location', 'northeast')
xlim([1 365])

subplot(3,1,2)
plot(time, z1, 'r','LineWidth',2)
hold on
plot(time, z2, 'r--','LineWidth',2)
ylabel('Z')
xline(diap_init, '--');
y_limits = ylim(gca);
max_y = y_limits(2);
text(diap_init + 2, max_y - max_y/10, 'Diapause Date')
xlim([1 365])

subplot(3,1,3)
plot(time, c1, 'k','LineWidth',2)
hold on
plot(time, c2, 'k--','LineWidth',2)
ylabel('C')
xlabel('Day of Year')
xlim([1 365])
