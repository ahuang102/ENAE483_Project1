% inert mass heuristics
% Lecture 6
clear, clc, close all

% gravity on earth
g = 9.81; %m/s^2

%% Project Requirements
% required payload
m_pl = 25000; % kg
% total delta V required to reach 300 km altitude
dV_total = 9300; % m/s
% translunar insertion burn for third stage
dv_TLI = 3150; % m/s
% radius of the final stage that does the burn
radius_TLI_stage = 2.5; % m % Upper stage diameter at payload interface (NOT LESS THAN 5 meters)

%% Selected Variables
% TODO: Determine number of stages. Assuming 2
% % STAGE SIZING % %
% First stage to get to 3000 m/s
r_stage1 = 3.5; % m 
% Second stage for remaining 6300 m/s
r_stage2 = 2.5; % m

% % PROPELLANT % %
ue_H2 = 4273;
ue_1 = ue_H2; % m/s % LOX/LH2 exit velocity L2 p16
% ue_2 = 3136; % m/s % LOX/RP-1 Exit velocity 
ue_2 = ue_H2;
ue_TLI = ue_H2;

% initial inert mass fraction estimation determined by used propellant
% increased by 10% for reusibility for stage 1
inert_fraction_H2 = 0.075;
% inert_fraction_RP1 = 0.063;
delta_1 = inert_fraction_H2 * 1.1;  % From L2 p 16 for LH2
delta_2 = inert_fraction_H2;
delta_3 = inert_fraction_H2;

% dV_1 optimal near 4500 but caps at 3000/3500 depending on method of reuse
dV_1 = 3000;
dV_2 = dV_total - dV_1;
% dV for first stage must be less than or equal to 3000 m/s for reusibility
% as per requirements if returning to landing site

% need 15% of initial proplelant load must remain -> 15% of propellant lost
% from increasing 15%
% so need 15% more initial propellant 


%% Calculations
% mass ratio between final and initial of first stage
r_1 = exp(-dV_1/ue_1)
r_2 = exp(-dV_2/ue_2)
r_TLI = exp(-dv_TLI/ue_TLI)

% payload fraction // lambda
y_1 = r_1 - delta_1;
y_2 = r_2 - delta_2;
y_3 = r_TLI - delta_3;

% total mass at stages 
M0_TLI = m_pl / y_3;

M0_2 = M0_TLI / y_2;
M0_1 = M0_2 / y_1; % inaccurate since Mp_1 was increased

% inert mass of stages
% M_inert_tli = M0_TLI * delta_3;
M_inert_2 = M0_2 * delta_2;
M_inert_1 = M0_1 * delta_1;


%% Inert Mass Heursitics
% propellant mass estimates
% stage 1
m_pr_1 = M0_1 - M_inert_1 - M0_2;
% stage 2
m_pr_2 = M0_2 - M_inert_2 - M0_TLI;
% stage 3
% m_pr_TLI = M0_TLI - M_inert_tli - M_pl;

% the equation on page 19 for stage mass
M_stage1 = M_inert_1 + m_pr_1;
M_stage2 = M_inert_2 + m_pr_2;

% M_stageTLI = M_inert_tli + m_pr_TLI;

% stage 1 - LOX/LH2 
e_1 = 0.987 * M_stage1^(-0.183);

% stage 2 - LOX/LH2
e_2 = 0.987 * M_stage2^(-0.183);
% stage 3
% e_tli = 0.987 * M_stageTLI^(-0.183);

%% With these values, iterate through to find the new masses necessary
% to reach delta V of 9300

disp("Heuristics")
% Keep M_pl constant but vary masses of stage


% new values of inert and propellant mass
% m_in_3 = e_tli * M_stageTLI;
m_in_2 = e_2 * M_stage2
m_in_1 = e_1 * M_stage1

% new values of propellant
% m_pr_3 = M_stageTLI - m_in_3
m_pr_2 = M_stage2 - m_in_2
m_pr_1 = M_stage1 - m_in_1

% new initial masses
% M0_3_2 = M_stageTLI + M_pl;
M0_2_2 = M_stage2 + M0_TLI;
M0_1_2 = M_stage1 + M0_2_2;
% New values for ratio
% r_3 = (m_in_3 + M_pl)/(M0_TLI)
r_2 = (m_in_2 + M0_TLI)/(M0_2_2)
r_1 = (m_in_1 + M0_2_2)/(M0_1_2)

% New delta V
dV1 = -ue_1 * log(r_1)
dV2 = -ue_2 * log(r_2)
% dV3 = -ue_TLI * log(r_3)
fprintf("New dV total = %4.3f vs %4.3f\n", dV1+dV2, dV_total)

% calculate new delta values
d1 = m_in_1 / M0_1_2;
d2 = m_in_2 / M0_2_2;
fprintf("New inert mass fractions are stage 1 = %4.3f, stage 2 = %4.3f.\n", d1, d2)


% 
% 
% % Check from lecture
% m_st1 = 483700;
% m_st2 = 115000;
% 
% m_in2 = 13460;
% m_in1 = 21240;
% 
% m_pr2 = 101500;
% m_pr1 = 462500;
% 
% 
% m_pr2 + m_in2
% m_pr1 + m_in1
% 
% m_pl = 25000;
% M0_2 = m_st2 + m_pl;
% M0_1 = m_st1 + M0_2;
% 
% 
% r2 = (m_in2 + m_pl)/ ( m_st2 + m_pl)
% r1 = (m_in1 + M0_2)/ (m_st1 + M0_2)
% 
% d1 = m_in1 / M0_1
% d2 = m_in2 / M0_2
% 
% 
% 
% return

% % m_in = E * m_st
% 
% % calculate m0 for each value in the ratio
% % m_pl, ue_1, ue_2, delta1, and delta2 stay the same.
% dV_ratio = linspace(0,1);
% % get all the possible dV1 and dV2
% dV1s = dV_total .* dV_ratio;
% dV2s = dV_total .* (1 - dV_ratio);
% 
% % mass ratio of final stage (a number)
% r_3 = exp(-dv_TLI./ue_TLI);
% % mass ratio of second stage
% r_2s = exp(-dV2s./ue_2);
% % mass ratio of first stage
% r_1s = exp(-dV1s./ue_1);
% 
% % calculate total mass at each stage
% M_total = @(r, e) ((1-e)./(r-e)) .* m_pl;
% M0_1 = M_total(r_1s, e_1);
% M0_2 = M_total(r_2s, e_2);
% M0_3 = M_total(r_3, e_tli);
% M_total = M0_1 + M0_2 + M0_3 + m_pl;
% [M_total_min, bi] = min(M_total);
% fprintf("New total mass: M0 = %4.3f, M0_2 = %4.3f, M0_3 = %4.3f kg.\n", M0_1(bi), M0_2(bi), M0_3)
% fprintf("Sum : %4.3f.\n", M_total_min)
% 
% % find inert mass fraction
% 
% l1 = m_pl / M0_1(bi)
% l2 = m_pl / M0_2(bi)
% l3 = m_pl / M0_3

% % find the minimum value
% [m0_min, index] = min(m0s);
% fprintf("Minimum initial mass = %4.3f kg.\n", m0_min)

