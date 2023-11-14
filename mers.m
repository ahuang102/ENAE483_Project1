% Mass estimation relations
% Started by Luke Brauch
clear, clc, close all

% gravity on earth
g = 9.81; %m/s^2

%% Project Requirements
% required payload
M_pl = 25000; % kg
% total delta V required to reach 300 km altitude
dV_total = 9300; % m/s
% translunar insertion burn for third stage
dv_TLI = 3150; % m/s
% radius of the upper stage that does the burn
radius_TLI_stage = 2.5; % m % Upper stage diameter at payload interface (NOT LESS THAN 5 meters)

%% Used values
% TODO: Determine number of stages. Assuming 2
r_stage1 = 3.5; % m 
r_stage2 = 3.5; % m

% specific impulse (TODO: based on chosen rockets)
Isp_1 = 380; % s
Isp_2 = 450; % s
Isp_TLI = 450; % s

% oxygen/fuel ratio
OF_ratio_1 = 3.2;
OF_ratio_2 = 6.1;
OF_ratio_TLI = 6; % random value

% initial inert mass fraction estimation
delta_1 = 0.07; 
delta_2 = 0.07;
delta_3 = 0.07;

% assume 50/50 split for the dV between stages sTODO
r = 0.5; % r = 1 means all the dV comes from the first stage

%% Constants
% densities of fuels
rho_LOX = 1140;
rho_LCH4 = 423;
rho_LH2 = 71;

%% Calculations
% exit velocity
ue_1 = Isp_1 * g; % m/s
ue_2 = Isp_2 * g; 
ue_TLI = Isp_TLI * g;
% delta V for each stage
dV_1 = dV_total * r;
dV_2 = dV_total * (1-r);

% mass ratio between final and initial of first stage
r_1 = exp(-dV_1/ue_1);
r_2 = exp(-dV_2/ue_2);
r_TLI = exp(-dv_TLI/ue_TLI);

% payload fraction
y_1 = r_1 - delta_1;
y_2 = r_2 - delta_2;
y_3 = r_TLI - delta_3;

% Initial mass of stages
M0_TLI = M_pl / y_3;
M0_2 = M0_TLI / y_2;
M0_1 = M0_2 / y_1;

% inert mass of stages
M_inert_tli = delta_3 * M0_TLI; 
M_inert_2 = delta_2*M0_2;
M_inert_1 = delta_1*M0_1;

% total inert mass estimate
Min_total = M_inert_1 + M_inert_2 + M_inert_tli;
M_total_pool = 0;
%% Propellant masses
Mp_1 = M0_1 * (1-r_1);
Mp_2 = M0_2 * (1-r_2);
Mp_3 = M0_TLI * (1-r_TLI);

%% Mass of propellant components
% mass of LOX component
M_LOX_1 = OF_ratio_1/(OF_ratio_1 + 1) * Mp_1;
% mass of LCH4 Component
M_LCH4_1 = 1/(OF_ratio_1 + 1) * Mp_1;
% Stage 2
M_LOX_2 = OF_ratio_2/(OF_ratio_2 + 1) * Mp_2;
M_LH2_2 = 1/(OF_ratio_2 + 1) * Mp_2;
% Stage 3 (for the TLI)
M_LOX_3 = OF_ratio_TLI/(OF_ratio_TLI + 1) * Mp_3;
M_LH2_3 = 1/(OF_ratio_TLI + 1) * Mp_3;


%% Mass of tanks
% mass of LOX tank
M_LOX_tank_1 = 0.0107 * M_LOX_1;
% mass of LCH4 tank
M_LCH4_tank_1 = 0.0287 * M_LCH4_1;
% (stage 2)
M_LOX_tank_2 = 0.0107 * M_LOX_2;
% mass of LH2 tank for second stage
M_LH2_tank_2 = 0.128 * M_LH2_2;
% stage 3
M_LOX_tank_TLI = 0.0107 * M_LOX_3;
M_LH2_tank_TLI = 0.128 * M_LH2_3;

% sum of tank masses
M_tank_total = M_LOX_tank_1 +M_LCH4_tank_1 + M_LOX_tank_2 + M_LH2_tank_2 + M_LOX_tank_TLI + M_LH2_tank_TLI;

%% Stage 1 Insulation Mass
% Volume of LOX tank

% total volume of tank needed
V_LOX_tank_1 = M_LOX_1 / rho_LOX;

% returns total height (including spherical caps) and surface area
[A_LOX_tank_1, h_LOX_total, r_LOX] = GetTankDimensions(r_stage1, V_LOX_tank_1);

% mass of insulation for tank
M_LOX_insul_1 = 1.123 * A_LOX_tank_1;

% Volume of LCH4 tank
V_LCH4_tank = M_LCH4_1 / rho_LCH4;
[A_LCH4_tank, h_LCH4_total, r_LCH4] = GetTankDimensions(r_stage1, V_LCH4_tank);

% mass of insulation for LCH4 tank
M_LCH4_insul_1 = 1.123 * A_LCH4_tank;

% height of stage 1
h_stage1 = h_LOX_total + h_LOX_total;
 
M_ins_1 = M_LOX_insul_1 + M_LCH4_insul_1;

%% Stage 2 Insulation Mass
% Volume of LOX tank

V_LOX_tank_2 = M_LOX_2 / rho_LOX;
[A_LOX_tank_2, h_LOX_tank_2, r_LOX_2] = GetTankDimensions(r_stage2, V_LOX_tank_2);

M_LOX_insul_2 = 1.123 * A_LOX_tank_2

% Volume of LH2 tank 
V_LH2_tank = M_LH2_2/ rho_LH2;

[A_LH2_tank, h_LH2_tank, r_LH2] = GetTankDimensions(r_stage2, V_LH2_tank);

% mass of insulation for LCH4 tank
M_LH2_insul_2 = 2.88 * A_LH2_tank;

h_stage2 = h_LH2_tank + h_LOX_tank_2;

M_ins_2 = M_LOX_insul_2 + M_LH2_insul_2;

%% Stage 3 Insulation Mass
% LOX tank volume/mass
V_LOX_tank_3 = M_LOX_3 / rho_LOX;
[A_LOX_tank_3, h_LOX_tank_3, r_LOX_3] = GetTankDimensions(radius_TLI_stage, V_LOX_tank_3);

M_LOX_insul_3 = 1.123 * A_LOX_tank_3;

% LH2 tank volume/mass
V_LH2_tank_3 = M_LH2_3 / rho_LH2;
[A_LH2_tank_3, h_LH2_tank_3, r_LH2_3] = GetTankDimensions(radius_TLI_stage, V_LH2_tank_3);

% mass of insulation for LCH4 tank
M_LH2_insul_3 = 2.88 * A_LH2_tank_3;

h_stage3 = h_LH2_tank_3 + h_LOX_tank_3;

M_ins_3 = M_LOX_insul_3 + M_LH2_insul_3;

fprintf("Total height of rocket = %4.3f m\n", h_stage1 + h_stage2 + h_stage3)

% total mass of tanks
M_ins_total = M_ins_3 + M_ins_2 + M_ins_1;

%% Fairings
% aft_faring of stage 1
% add 3 meters for engine and stuff
A_cylinder_1 = 2*pi*(h_stage1 + 3) * r_stage1;
M_aft_fairing_1 = 4.95*(A_cylinder_1)^(1.15);

% Fairing for stage 1 to 2
A_frustrum = pi*(r_stage1 + r_stage2)*sqrt((r_stage1-r_stage2)^2+(r_stage1+r_stage2)^2);
M_intertank = 4.95*(A_frustrum)^(1.15)

% Stage 2
% aft_faring of stage 2
% add 3 meters for engine and stuff
A_cylinder_2 = 2*pi*(h_stage2 + 3) * r_stage2;
M_aft_fairing_2 = 4.95*(A_cylinder_2)^(1.15)

% Fairing for stage 2 - 3
A_frustrum_2 = pi*(radius_TLI_stage + r_stage2)*sqrt((radius_TLI_stage-r_stage2)^2+(radius_TLI_stage+r_stage2)^2);
M_intertank_2 = 4.95*(A_frustrum_2)^(1.15)

% Stage 3
% aft_faring of stage 3
A_cylinder_3 = 2*pi*(h_stage3 + 3) * radius_TLI_stage;
M_aft_fairing_3 = 4.95*(A_cylinder_3)^(1.15)

% playload fairing
A_pl_fairing = pi*radius_TLI_stage*(sqrt(radius_TLI_stage^2 + 7^2))
M_pl_fairing = 4.95*(A_pl_fairing)^(1.15)

% total fairing mass
M_aft_fairing = M_aft_fairing_1 + M_aft_fairing_2 + M_aft_fairing_3;
M_intertank_Total = M_intertank_2 + M_intertank;

M_total_fairing = M_aft_fairing + M_pl_fairing + M_intertank_Total; 

%% Mass of avionics
M_avionics_1 = 10*M0_1^(0.361);
M_avionics_2 = 10*M0_2^(0.361);
M_avionics_3 = 10*M0_TLI^(0.361);

M_avionics = M_avionics_3 + M_avionics_2 + M_avionics_1;

%% mass of wiring
M_wiring_1 = 1.058*sqrt(M0_1)*h_stage1^(.25) ;
M_wiring_2 = 1.058*sqrt(M0_2)*h_stage2^(.25);
M_wiring_3 = 1.058*sqrt(M0_TLI)*h_stage3^(.25);

M_wiring = M_wiring_3 + M_wiring_2 + M_wiring_1;

%% Engine mass
% number of engines for each stage
n_engines = 2;

% stage 1
TW_ratio = 1.3;
% chamber pressure
P_1 = 13.4e6; % Pa
Ae_At = 14.5; % exit area over throat area
 
% engine thrust (each)
T = M0_1*g*TW_ratio/(n_engines); % N
% engine mass (each)
M_engine_1 = 7.81e-4*T + 3.37e-5*T*sqrt(Ae_At) + 59;
% mass of thrust structure (each)
M_thrust_struc_1 = 2.55e-4 * T;

M_engines_1 = (M_engine_1 + M_thrust_struc_1) * n_engines;

% Stage 2
TW_ratio_2 = 0.76;
P_2 = 7.66e6; % Pa
Ae_At = 80;
T_2 = M0_2*g*TW_ratio_2/(n_engines); % N
M_engine_2 = 7.81e-4*T_2 + 3.37e-5*T_2*sqrt(Ae_At) + 59;
M_thrust_struc_2 = 2.55e-4 * T_2;

M_engines_2 = (M_engine_2 + M_thrust_struc_2) * n_engines;

% Stage 3
TW_ratio_3 = 0.76; % TODO: Find values
P_3 = 7.66e6; % Pa
Ae_At = 80; 
T_3 = M0_2*g*TW_ratio_3/(n_engines); % N
M_engine_3 = 7.81e-4*T_3 + 3.37e-5*T_3*sqrt(Ae_At) + 59;
M_thrust_struc_3 = 2.55e-4 * T_3;

M_engines_3 = (M_engine_3 + M_thrust_struc_3) * n_engines;

% Total engine mass
M_engines = (M_engine_3 + M_engine_2 + M_engine_1) * n_engines;
M_thrust_structures = (M_thrust_struc_3 + M_thrust_struc_2 + M_thrust_struc_1) * n_engines;

%% gimbal mass
M_gimbals_1 = 237.8 * (T/P_1)^(0.9375);
% stage 2
M_gimbals_2 = 237.8 * (T_3/P_2)^(0.9375);
% stage 3
M_gimbals_3 = 237.8 * (T_3/P_3)^(0.9375);

M_gimbals = M_gimbals_1 + M_gimbals_2 + M_gimbals_3;

%% Total Inert Mass

M_stage1 = M_LOX_tank_1 + M_LCH4_tank_1 + M_ins_1 + M_aft_fairing_1 + M_intertank + M_avionics_1 + M_wiring_1 + M_engines_1 + M_gimbals_1 ;
Margin_1 = (M_inert_1 - M_stage1)/M_stage1 * 100;
fprintf("Inert mass of stage 1 = %4.3f kg. This has a margin of %4.2f %%.\n", M_stage1, Margin_1)

M_stage2 = M_LOX_tank_2 + M_LH2_tank_2 + M_ins_2 + M_aft_fairing_2 + M_intertank_2 + M_avionics_2 + M_wiring_2 + M_engines_2 + M_gimbals_2;
Margin_2 = (M_inert_2 - M_stage2)/M_stage2 * 100;
fprintf("Inert mass of stage 2 = %4.3f kg. This has a margin of %4.2f %%.\n", M_stage2, Margin_2)

M_stage3 = M_LOX_tank_TLI + M_LH2_tank_TLI + M_ins_3 + M_aft_fairing_3 + M_pl_fairing + M_avionics_3 + M_wiring_3 + M_engines_3 + M_gimbals_3;
Margin_3 = (M_inert_tli - M_stage3)/M_stage3 * 100;
fprintf("Inert mass of stage 3 = %4.3f kg. This has a margin of %4.2f %%.\n", M_stage3, Margin_3)

M_total = M_stage1 + M_stage2 + M_stage3;
fprintf("Total inert mass = %4.3f kg\n", M_total)

Margin = (Min_total - M_total)/M_total * 100;
fprintf("Margin = %4.2f%% \n", Margin)

% Compare alternative way of summing to check answers
M_total2 = M_gimbals + M_thrust_structures + M_engines + M_wiring + M_avionics + M_total_fairing + M_ins_total + M_tank_total;

function [A, h, r] = GetTankDimensions(stage_radius, volume)
    % volume of the semispherical heads to the tank
    V_ends = 4/3*pi*stage_radius^3;
    % volume remaining
    V_cyl = volume - V_ends;
    % use smaller sphere if this happens
    if V_cyl < 0 
        % just estimate a smaller sphere
        r = (volume/(4*pi/3))^(1/3);
        % surface area of sphere
        A = 4*pi*r^2;
        fprintf("radius = %4.3f m sphere.\n", r)
        h = r*2;
    else
        % height of cylindrical portion of tank
        h_cyl = V_cyl / (pi*stage_radius^2);
        % total height of this tank
        h = h_cyl + stage_radius*2;
        fprintf("height = %4.3f m with spherical caps\n", h)
        % surface area of tank (area of sphere ends + area of cylinder)
        A = 4*pi*stage_radius^2 + 2*pi*stage_radius*h_cyl;
        r = stage_radius;
    end
end
