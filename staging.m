% Find best dv distribution for different criteria
clear, clc, close all
%% Case 1: 3 stages
M_pl = 25000; % kg
dV_total = 9300; % m/s
% translunar insertion burn for third stage
dv_TLI = 3150; % m/s

% % PROPELLANT % %
ue_H2 = 4273;
% ue_RP1 = 3136;

ue_1 = ue_H2; % m/s % LOX/LH2 exit velocity L2 p16
ue_2 = ue_H2;
ue_TLI = ue_H2;

% increased by 10% for reusibility
inert_fraction_H2 = 0.075;
% inert_fraction_RP1 = 0.063;
delta_1 = inert_fraction_H2 * 1.1;  % From L2 p 16 for LH2
delta_2 = inert_fraction_H2;
delta_3 = inert_fraction_H2;

% calculate m0 for each value in the ratio
% m_pl, ue_1, ue_2, delta1, and delta2 stay the same.
dV_ratio = linspace(0,1);
% get all the possible dV1 and dV2
dV1s = dV_total .* dV_ratio;
dV2s = dV_total .* (1 - dV_ratio);

% mass ratio of final stage
r_3s = exp(-dv_TLI./ue_TLI);
% mass ratio of second stage
r_2s = exp(-dV2s./ue_2);
% mass ratio of first stage
r_1s = exp(-dV1s./ue_1);

% Pl fraction of final stage
l_TLIs = r_3s - delta_3;
% pl fraction of second stage
l_2s = r_2s - delta_2;
% pl fraction of first stage
l_1s = r_1s - delta_1;

% initial mass of final stage
m0_TLIs = M_pl ./ l_TLIs;
% initial mass of second stage
m0_2s = m0_TLIs ./ l_2s;
% initial mass
m0s = m0_2s ./ l_1s;

% find the minimum value
[m0_min, index] = min(m0s);
fprintf("Minimum initial mass = %4.3f kg.\n", m0_min)

% dV at this ratio
dV1 = dV1s(index);
dV2 = dV2s(index);
fprintf("Optimal delta dV for second stage = %4.3f and first stage = %4.3f.\n", dV2, dV1)

%% plotting
plot( dV_ratio, m0s)
title("Total initial mass vs ratio of delta V distribution")
xlabel("Ratio (1 == all dV1)")
hold on, grid on
% plot ideal points
plot(dV_ratio(index), m0_min, '*')
% find point where dV1 = 3000
tol = 100;
[ii,jj] = find(abs(dV1s-3000)<tol);
q = jj(1);
r = dV_ratio(q);
M0_3000 = m0s(q);
plot(r, M0_3000, 'o')
fprintf("When ratio = %4.3f, dV1 = %4.3f m/s", r, dV1s(q))
fprintf("Total mass = %4.3f kg.\n", M0_3000)

% % find the other values
% m_in_1 = delta1*m0
% m_pr_1 = m0 - m_in_1 - m_pl_1s(index)
% m0_2 = m_pl_1s(index);
% m_in_2 = delta2*m0_2
% m_pr_2 = m0_2 - m_in_2 - m_pl
