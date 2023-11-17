% Mass estimation relations
% Started by Luke Brauch
% Lecture 6
% clear, clc, close all

% gravity on earth
g = 9.81; %m/s^2

%% Project Requirements
% required payload
M_pl = 25000; % kg
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
ue_RP1 = 3136;

ue_1 = ue_H2; % m/s % LOX/LH2 exit velocity L2 p16
% ue_2 = 3136; % m/s % LOX/RP-1 Exit velocity 
ue_2 = ue_H2;
ue_TLI = ue_H2;
% ue_TLI = 3136; % m/s % needs to be liquid / controllable since it is a maneuver

% oxygen/fuel ratio of H2
OF_ratio_LH2 = 3.2;
% https://kinetics.nist.gov/RealFuels/macccr/macccr2008/Bruno2.pdf
% typical  OF ratio for RP1
% OF_ratio_RP1 = 2.56;

% Uses H2
OF_ratio_TLI = OF_ratio_LH2; 

% initial inert mass fraction estimation
% determined by used propellant
% increased by 10% for reusibility
inert_fraction_H2 = 0.075;
% inert_fraction_RP1 = 0.063;
% delta_1 = inert_fraction_H2 * 1.1;  % From L2 p 16 for LH2
% delta_2 = inert_fraction_H2;
delta_3 = inert_fraction_H2;

% data from inert mass heuristic calculations
delta_1 = 0.0512;
delta_2 = 0.0812; 

% % ENGINES % % 
n_engines_1 = 1;
n_engines_2 = 1;
n_engines_3 = 1;

TW_ratio_1 = 1.3;
TW_ratio_2 = 0.76;
TW_ratio_3 = 0.76;

% Exit area to throat area ratio
area_ratio_1 = 80;
area_ratio_2 = 80;
area_ratio_3 = 80;

% chamber pressure
P_1 = 7.66e6; % Pa
P_2 = 7.66e6;
P_3 = 7.66e6;

% % STAGING % %
% Assume 
dV_1 = 3000;
dV_2 = dV_total - dV_1;
% dV for first stage must be less than or equal to 3000 m/s for reusibility
% as per requirements if returning to landing site

%% Constants
% densities of fuels (kg/m^3)
rho_LOX = 1140;
rho_LCH4 = 423;
rho_RP1 = 820;
rho_LH2 = 71;

% establish functions for convenience for tank calculations
m_LH2_tank = @(M) 0.128 * M;
m_LOX_tank = @(M) 0.0107 * M;
m_RP1_tank = @(M) 0.0148 * M;

%% Calculations
M0_2=     4.194658337752676e+05

%% STAGE 1 %%
% mass ratio between final and initial of first stage
r_1 = exp(-dV_1/ue_1);

% payload fraction // lambda
y_1 = r_1 - delta_1;

% Initial mass of stages
M0_1 = M0_2 / y_1; % inaccurate since Mp_1 was increased

% inert mass of stages
M_inert_1 = M0_1 * delta_1;

% Propellant Mass

Mp_1 = M0_1 * (1-r_1);

% need 15% of initial proplelant load must remain -> 15% of propellant lost
% from increasing 15%
% so need 15% more initial propellant 
% stage mass
E_LH2 = 0.987 * (M_inert_1 + Mp_1)^(-.275);
% Find new stage 1 mass accounting for 15% more unused propellant
Mp_1 = Mp_1 * 1.15; % increase propellant mass
% find new inert mass

Mass1.Inert = Mp_1 * E_LH2 / (1-E_LH2);

% Mass of propellant
M_LOX_1 = OF_ratio_LH2/(OF_ratio_LH2 + 1) * Mp_1;
M_LH2_1 = 1/(OF_ratio_LH2 + 1) * Mp_1;

% Mass of Tanks

% mass of LOX tank
Mass1.LOX_tank = m_LOX_tank(M_LOX_1);
% mass of H2 tank 
Mass1.LH2_tank = m_LH2_tank(M_LH2_1);

% Insulation Mass

% Volume of LOX tank
V_LOX_tank_1 = M_LOX_1 / rho_LOX;

% returns total height (including spherical caps) and surface area
[A_LOX_tank_1, h_LOX_total, r_LOX] = GetTankDimensions(r_stage1, V_LOX_tank_1);

% mass of insulation for tank
Mass1.LOX_insulation = 1.123 * A_LOX_tank_1;

% Volume of LH2 tank
V_LH2_tank_1 = M_LH2_1 / rho_LH2;
[A_LH2_tank_1, h_LH2_tank_1, r_LH2_1] = GetTankDimensions(r_stage1, V_LH2_tank_1);

% mass of insulation for LH2 tank
Mass1.LH2_insulation = 2.88 * A_LH2_tank_1;

% height of stage 1
h_stage1 = h_LOX_total + h_LH2_tank_1;

% Fairings

% aft_faring of stage 1
% add 3 meters for engine and misc
A_cylinder_1 = 2*pi*(h_stage1 + 3) * r_stage1;
Mass1.aft_fairing = 4.95*(A_cylinder_1)^(1.15);

% Fairing for stage 1 to 2 (or payload of stage 1)
A_frustrum = pi*(r_stage1 + r_stage2)*sqrt((r_stage1-r_stage2)^2+(r_stage1+r_stage2)^2);
Mass1.intertank_fairing = 4.95*(A_frustrum)^(1.15);

% Mass of avionics
Mass1.avionics = 10*M0_1^(0.361);

% Mass of Wiring
Mass1.wiring = 1.058*sqrt(M0_1)*h_stage1^(.25);

% Engine mass

% engine thrust 
T = M0_1*g*TW_ratio_1/(n_engines_1); % N
% engine mass 
Mass1.engines = (7.81e-4*T + 3.37e-5*T*sqrt(area_ratio_1) + 59) * n_engines_1;
% mass of thrust structure 
Mass1.thrust_structure = 2.55e-4 * T * n_engines_1;

% gimbal mass
Mass1.gimbals = 237.8 * (T/P_1)^(0.9375);


%% Total Inert Mass

Masses.Stage1 = Mass1;
Table = struct2table(Mass1);
q = sum(Table(1,2:end), 'all').Variables;
Mass1.Total_inert_mass = q;

Margin_1 = (M_inert_1 - q)./q .* 100;
fprintf("Inert mass of stage 1 = %4.3f kg. This has a margin of %4.2f %%.\n", q, Margin_1)

return

%% STAGE 2 %%
% exit velocity from Isp
% ue_1 = Isp_1 * g; % m/s
% ue_2 = Isp_2 * g; 
% ue_TLI = Isp_TLI * g;

% delta V for each stage using ratio
% r = 0.5; % r = 1 means all the dV comes from the first stage
% dV_1 = dV_total * r;
% dV_2 = dV_total * (1-r);

% mass ratio between final and initial of first stage
r_1 = exp(-dV_1/ue_1);
r_2 = exp(-dV_2/ue_2);
r_TLI = exp(-dv_TLI/ue_TLI);

% payload fraction // lambda
y_1 = r_1 - delta_1;
y_2 = r_2 - delta_2;
y_3 = r_TLI - delta_3;

% Initial mass of stages
M0_TLI = M_pl / y_3;
M0_2 = M0_TLI / y_2;
M0_1 = M0_2 / y_1; % inaccurate since Mp_1 was increased

% inert mass of stages
M_inert_tli = M0_TLI * delta_3; 
M_inert_2 = M0_2 * delta_2;
M_inert_1 = M0_1 * delta_1;

M_total_pool = 0;
%% Propellant masses
Mp_1 = M0_1 * (1-r_1);
Mp_2 = M0_2 * (1-r_2);
Mp_3 = M0_TLI * (1-r_TLI);

% Stage inert mass for stage 1 (
% E_RP1 = 1.6062 * (M_inert_1 + Mp_1)^(-.275);
E_LH2 = 0.987 * (M_inert_1 + Mp_1)^(-.275);
% Find new stage 1 mass accounting for 15% more unused propellant
Mp_1 = Mp_1 * 1.15; % increase propellant mass
% find new inert mass
% M_inert_1 = Mp_1 * E_RP1 / (1-E_RP1)
M_inert_1 = Mp_1 * E_LH2 / (1-E_LH2);

% total inert mass estimate
Min_total = M_inert_1 + M_inert_2 + M_inert_tli;
Mass1.inert = Min_total;
fprintf("Total estimated inert mass: %4.3f kg.\n", Min_total)

%% Mass of propellant components
% mass of LOX component
M_LOX_1 = OF_ratio_LH2/(OF_ratio_LH2 + 1) * Mp_1;
M_LH2_1 = 1/(OF_ratio_LH2 + 1) * Mp_1;
% Stage 2
M_LOX_2 = OF_ratio_LH2/(OF_ratio_LH2 + 1) * Mp_2;
M_LH2_2 = 1/(OF_ratio_LH2 + 1) * Mp_2;

% Stage 3 (for the TLI)
M_LOX_3 = OF_ratio_TLI/(OF_ratio_TLI + 1) * Mp_3;
M_LH2_3 = 1/(OF_ratio_TLI + 1) * Mp_3;


%% Mass of tanks
% establish functions for convenience
m_LH2_tank = @(M) 0.128 * M;
m_LOX_tank = @(M) 0.0107 * M;
m_RP1_tank = @(M) 0.0148 * M;

% mass of LOX tank
M_LOX_tank_1 = m_LOX_tank(M_LOX_1);
% mass of H2 tank for s1
M_LH2_tank_1 = m_LH2_tank(M_LH2_1);

% (stage 2)
M_LOX_tank_2 = m_LOX_tank(M_LOX_2);
% mass of LH2 tank for second stage
M_LH2_tank_2 = m_LH2_tank(M_LH2_2);

% stage 3
M_LOX_tank_TLI = m_LOX_tank(M_LOX_3);
M_LH2_tank_TLI = m_LH2_tank(M_LH2_3);

% sum of tank masses
M_tank_total = M_LOX_tank_1 +M_LH2_tank_1 + M_LOX_tank_2 + M_LH2_tank_2 + M_LOX_tank_TLI + M_LH2_tank_TLI;
% Mass.tanks = M_tank_total;
Mass1.LOX_tank = M_LOX_tank_TLI + M_LOX_tank_2 + M_LOX_tank_1;
Mass1.LH2_tank = M_LH2_tank_TLI + M_LH2_tank_2 + M_LH2_tank_1;

%% Stage 1 Insulation Mass
% Volume of LOX tank

% total volume of tank needed
V_LOX_tank_1 = M_LOX_1 / rho_LOX;

% returns total height (including spherical caps) and surface area
[A_LOX_tank_1, h_LOX_total, r_LOX] = GetTankDimensions(r_stage1, V_LOX_tank_1);

% mass of insulation for tank
M_LOX_insul_1 = 1.123 * A_LOX_tank_1;

% Volume of LCH4 tank
V_LH2_tank_1 = M_LH2_1 / rho_LH2;
[A_LH2_tank_1, h_LH2_tank_1, r_LH2_1] = GetTankDimensions(r_stage1, V_LH2_tank_1);

% mass of insulation for LH2 tank
M_LH2_insul_1 = 2.88 * A_LH2_tank_1;

% height of stage 1
h_stage1 = h_LOX_total + h_LH2_tank_1;

% mass of stage 1 insulation
M_ins_1 = M_LOX_insul_1 + M_LH2_insul_1;

%% Stage 2 Insulation Mass
% Volume of LOX tank

V_LOX_tank_2 = M_LOX_2 / rho_LOX;
[A_LOX_tank_2, h_LOX_tank_2, r_LOX_2] = GetTankDimensions(r_stage2, V_LOX_tank_2);
% mass of insulation
M_LOX_insul_2 = 1.123 * A_LOX_tank_2;

% Volume of LH2 tank 
V_LH2_tank = M_LH2_2/ rho_LH2;

% dimensions of tank
[A_LH2_tank_2, h_LH2_tank_2, r_LH2_2] = GetTankDimensions(r_stage2, V_LH2_tank);

% mass of insulation for LCH4 tank
M_LH2_insul_2 = 2.88 * A_LH2_tank_2;

% height of stage
h_stage2 = h_LH2_tank_2 + h_LOX_tank_2;
% mass of stage 2 insulation
M_ins_2 = M_LOX_insul_2 + M_LH2_insul_2;

%% Stage 3 Insulation Mass
% LOX tank volume/mass
V_LOX_tank_3 = M_LOX_3 / rho_LOX;
[A_LOX_tank_3, h_LOX_tank_3, r_LOX_3] = GetTankDimensions(radius_TLI_stage, V_LOX_tank_3);

M_LOX_insul_3 = 1.123 * A_LOX_tank_3;

% LH2 tank volume/mass
V_LH2_tank_3 = M_LH2_3 / rho_LH2;
[A_LH2_tank_3, h_LH2_tank_3, r_LH2_3] = GetTankDimensions(radius_TLI_stage, V_LH2_tank_3);

M_LH2_insul_3 = 2.88 * A_LH2_tank_3;

h_stage3 = h_LH2_tank_3 + h_LOX_tank_3;

M_ins_3 = M_LOX_insul_3 + M_LH2_insul_3;

fprintf("Total height of rocket = %4.3f m\n", h_stage1 + h_stage2 + h_stage3)

% total mass of tanks
M_ins_total = M_ins_3 + M_ins_2 + M_ins_1;
% Mass.insulation = M_ins_total;
Mass1.LOX_insulation = M_LOX_insul_3 + M_LOX_insul_2 + M_LOX_insul_1;
Mass1.LH2_insulation = M_LH2_insul_3 + M_LH2_insul_2 + M_LH2_insul_1;

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

Mass1.payload_fairing = M_pl_fairing;
Mass1.intertank_fairing = M_intertank_Total;
Mass1.aft_fairing = M_aft_fairing;

%% Mass of avionics
M_avionics_1 = 10*M0_1^(0.361);
M_avionics_2 = 10*M0_2^(0.361);
M_avionics_3 = 10*M0_TLI^(0.361);

M_avionics = M_avionics_3 + M_avionics_2 + M_avionics_1;
Mass1.avionics = M_avionics;

%% mass of wiring
M_wiring_1 = 1.058*sqrt(M0_1)*h_stage1^(.25) ;
M_wiring_2 = 1.058*sqrt(M0_2)*h_stage2^(.25);
M_wiring_3 = 1.058*sqrt(M0_TLI)*h_stage3^(.25);

M_wiring = M_wiring_3 + M_wiring_2 + M_wiring_1;
Mass1.wiring = M_wiring;

%% Engine mass
% number of engines for each stage
n_engines_1 = 1;

% stage 1
TW_ratio_1 = 1.3;
% chamber pressure
P_1 = 13.4e6; % Pa
area_ratio_1 = 14.5; % exit area over throat area
 
% engine thrust (each)
Table = M0_1*g*TW_ratio_1/(n_engines_1); % N
% engine mass (each)
M_engine_1 = 7.81e-4*Table + 3.37e-5*Table*sqrt(area_ratio_1) + 59;
% mass of thrust structure (each)
M_thrust_struc_1 = 2.55e-4 * Table;

M_engines_1 = (M_engine_1 + M_thrust_struc_1) * n_engines_1;

% Stage 2
TW_ratio_2 = 0.76;
P_2 = 7.66e6; % Pa
area_ratio_1 = 80;
T_2 = M0_2*g*TW_ratio_2/(n_engines_1); % N
M_engine_2 = 7.81e-4*T_2 + 3.37e-5*T_2*sqrt(area_ratio_1) + 59;
M_thrust_struc_2 = 2.55e-4 * T_2;

M_engines_2 = (M_engine_2 + M_thrust_struc_2) * n_engines_1;

% Stage 3
TW_ratio_3 = 0.76; % TODO: Find values
P_3 = 7.66e6; % Pa
area_ratio_1 = 80; 
T_3 = M0_2*g*TW_ratio_3/(n_engines_1); % N
M_engine_3 = 7.81e-4*T_3 + 3.37e-5*T_3*sqrt(area_ratio_1) + 59;
M_thrust_struc_3 = 2.55e-4 * T_3;

M_engines_3 = (M_engine_3 + M_thrust_struc_3) * n_engines_1;

% Total engine mass
M_engines = (M_engine_3 + M_engine_2 + M_engine_1) * n_engines_1;
M_thrust_structures = (M_thrust_struc_3 + M_thrust_struc_2 + M_thrust_struc_1) * n_engines_1;

Mass1.engines = M_engines;
Mass1.thrust_structure = M_thrust_structures;

%% gimbal mass
M_gimbals_1 = 237.8 * (Table/P_1)^(0.9375);
% stage 2
M_gimbals_2 = 237.8 * (T_3/P_2)^(0.9375);
% stage 3
M_gimbals_3 = 237.8 * (T_3/P_3)^(0.9375);

M_gimbals = M_gimbals_1 + M_gimbals_2 + M_gimbals_3;
Mass1.gimbals = M_gimbals;

Table = struct2table(Mass1)
%% Total Inert Mass

M_stage1 = M_LOX_tank_1 + M_LH2_tank_1 + M_ins_1 + M_aft_fairing_1 + M_intertank + M_avionics_1 + M_wiring_1 + M_engines_1 + M_gimbals_1 ;
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

%% FINAL STAGE %%






% Returns: A = surface area of tank, h = total height of tank, r = radius
% of tank (if a sphere is smaller than the stage radius)
function [A, h, r] = GetTankDimensions(stage_radius, volume)
    % volume of the 2 semispherical heads to the tank
    V_ends = 4/3*pi*stage_radius^3;
    % volume remaining for the inner cylinder
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
