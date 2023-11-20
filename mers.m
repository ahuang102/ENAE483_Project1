% Mass estimation relations
clear, clc, close all

% gravity on earth
g = 9.81; %m/s^2

%% Project Requirements
% required payload
M_pl = 25000; % kg
% total delta V required to reach 300 km altitude
dV_total = 9300; % m/s
% translunar insertion burn for third stage
dV_TLI = 3150; % m/s
% radius at payload interface greater than or equal to this
r_PL = 2.5; % Upper stage diameter at payload interface (NOT LESS THAN 5 meters)
% radius of the final stage that does the burn

%% Selected Variables
radius_TLI_stage = 1.6; % m % adjust to fit rest of ship 

% % STAGE SIZING % %
% First stage to get to 3000 m/s
r_stage1 = 4.4; % m 
% Second stage for remaining 6300 m/s
r_stage2 = 4; % m

% % PROPELLANT % %
ue_H2 = 4273;
ue_RP1 = 3136;

ue_1 = ue_H2; % m/s % LOX/LH2 exit velocity L2 p16
% ue_2 = 3136; % m/s % LOX/RP-1 Exit velocity 
ue_2 = ue_H2;
ue_3 = ue_RP1;
% ue_TLI = 3136; % m/s % needs to be liquid / controllable since it is a maneuver

% oxygen/fuel ratio of H2
OF_ratio_LH2 = 6;
% https://kinetics.nist.gov/RealFuels/macccr/macccr2008/Bruno2.pdf
% typical  OF ratio for RP1
OF_ratio_RP1 = 2.56;

% first two stages use H2
OF_ratio_1 = OF_ratio_LH2;
OF_ratio_2 = OF_ratio_LH2;
% Uses RP1
OF_ratio_3 = OF_ratio_RP1;

% initial inert mass fraction estimation
% determined by used propellant
delta_3 = 0.06;
delta_2 = 0.0726;
delta_1 = .0526;

% % ENGINES % % 
n_engines_1 = 6;
n_engines_2 = 6;
n_engines_3 = 2; % backup engine for burn

TW_ratio_1 = 1.3;
TW_ratio_2 = 1.2;
TW_ratio_3 = 0.76;

% Exit area to throat area ratio
area_ratio_1 = 30;
area_ratio_2 = 30;
area_ratio_3 = 80;

% chamber pressure
P_1 = 7.66e6; % Pa
P_2 = 7.66e6;
P_3 = 7.66e6;

% % STAGING % %
% Assume 
dV_1 = 3500;
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

D_TLI = 2*radius_TLI_stage;
D_2 = 2*r_stage2;
D_1 = 2*r_stage1;


%% STAGE 3 %% 
% mass ratio
r_3 = exp(-dV_TLI/ue_3);
% payload fraction // lambda
y_3 = r_3 - delta_3;

% Initial mass of stages
M0_3 = M_pl / y_3;

% inert mass of stages
Mass3.Inert = M0_3 * delta_3;

% Propellant Mass

Mp_3 = M0_3 * (1-r_3);

% Mass of individual propellant
M_LOX_3 = OF_ratio_3/(OF_ratio_3 + 1) * Mp_3;
M_RP1_3 = 1/(OF_ratio_3 + 1) * Mp_3;

% Mass of Tanks

% mass of LOX tank
Mass3.LOX_tank = m_LOX_tank(M_LOX_3);
% mass of H2 tank 
Mass3.Fuel_tnak = m_RP1_tank(M_RP1_3);

% Insulation Mass

% Volume of LOX tank
V_LOX_tank_3 = M_LOX_3 / rho_LOX ;

% returns total height (including spherical caps) and surface area 
[A_LOX_tank_3, h_LOX_total_3, r_LOX_3] = GetTankDimensions(radius_TLI_stage, V_LOX_tank_3);

% mass of insulation for tank 
Mass3.LOX_insulation = 1.123 * A_LOX_tank_3;

% Volume of LH2 tank (Larger tank)
V_fuel_tank_3 = M_RP1_3 / rho_RP1; 
[A_fuel_tank_3, h_fuel_tank_3, r_fuel_3] = GetTankDimensions(radius_TLI_stage, V_fuel_tank_3);

% mass of insulation for LH2 tank % No insulation for RP1
% Mass3.LH2_insulation = 2.88 * A_fuel_tank_3;
Mass3.LH2_insulation = 0;

% height
h_stage3 = h_LOX_total_3 + h_fuel_tank_3;

%% Fairings

% playload fairing
h_pl = 6;
A_pl_fairing = pi*r_PL*(sqrt(r_PL^2 + h_pl^2));
Mass3.playload_fairing = 4.95*(A_pl_fairing)^(1.15);

% Intertank fairing
% height between center of tanks
h_intertank = r_fuel_3 + r_LOX_3;
A_frustrum_3 = pi*(r_fuel_3 + r_LOX_3)*sqrt((r_fuel_3 - r_LOX_3)^2+(h_intertank)^2);
mi1 = 4.95*(A_frustrum_3)^(1.15);
% save total intertank
Mass3.intertank_fairing  = mi1;

% aft_faring for between stages
% add a few meters for thrusters
h_aft = r_LOX_3 + 2;
A_cylinder_3 = 2*pi * r_LOX_3 * h_aft;
m_aft = 4.95*(A_cylinder_3)^(1.15);
Mass3.aft_fairing = m_aft;

% Mass of avionics
Mass3.avionics = 10*M0_3^(0.361);

% Mass of Wiring
Mass3.wiring = 1.058*sqrt(M0_3)*h_stage3^(.25);

% Engine mass

% engine thrust 
T = M0_3*g*TW_ratio_3/(n_engines_3); % N
% engine mass 
Mass3.engines = (7.81e-4*T + 3.37e-5*T*sqrt(area_ratio_3) + 59) * n_engines_3;
% mass of thrust structure 
Mass3.thrust_structure = 2.55e-4 * T * n_engines_3;

% gimbal mass
Mass3.gimbals = 237.8 * (T/P_3)^(0.9375) * n_engines_3;

%% Total Inert Mass
Masses.Stage3 = Mass3;
Table = struct2table(Mass3);
q = sum(Table(1,2:end), 'all').Variables;
% Mass3.Total_inert_mass = q;

Margin_3 = (Mass3.Inert - q)./q .* 100;
fprintf("Inert mass of stage 3 = %4.3f kg. This has a margin of %4.2f %%.\n", q, Margin_3)
% Mass3.Margin = Margin_3;

% Length over diameter not to exceed 10
LDr_3 = h_stage3/D_TLI;

fprintf("Length/radius ratio: %4.3f.\n", LDr_3)

Mass3 = struct2table(Mass3)

%% STAGE 2 %% 
% mass ratio between final and initial of first stage
r_2 = exp(-dV_2/ue_2);

% payload fraction // lambda
y_2 = r_2 - delta_2;

% Initial mass of stages
M0_2 = M0_3 / y_2;

% inert mass of stages
Mass2.Inert = M0_2 * delta_2;

% Propellant Mass

Mp_2 = M0_2 * (1-r_2);
% Mass of individual propellant
M_LOX_2 = OF_ratio_2/(OF_ratio_2 + 1) * Mp_2;
M_LH2_2 = 1/(OF_ratio_2 + 1) * Mp_2;

% Mass of Tanks

% mass of LOX tank
Mass2.LOX_tank = m_LOX_tank(M_LOX_2);
% mass of H2 tank 
Mass2.LH2_tank = m_LH2_tank(M_LH2_2);

% Insulation Mass

% Volume of LOX tank
V_LOX_tank_2 = M_LOX_2 / rho_LOX;

% returns total height (including spherical caps) and surface area
[A_LOX_tank_2, h_LOX_total_2, r_LOX_2] = GetTankDimensions(r_stage2, V_LOX_tank_2);

% mass of insulation for tank
Mass2.LOX_insulation = 1.123 * A_LOX_tank_2;

% Volume of LH2 tank
V_LH2_tank_2 = M_LH2_2 / rho_LH2;
[A_LH2_tank_2, h_LH2_tank_2, r_LH2_2] = GetTankDimensions(r_stage2, V_LH2_tank_2);

% mass of insulation for LH2 tank
Mass2.LH2_insulation = 2.88 * A_LH2_tank_2;

% height
h_stage2 = h_LOX_total_2 + h_LH2_tank_2;

%% Fairings

% Fairing for stage 2-3
A_frustrum_2 = pi*(radius_TLI_stage + r_stage2)*sqrt((radius_TLI_stage-r_stage2)^2+(radius_TLI_stage+r_stage2)^2);
mi1 = 4.95*(A_frustrum_2)^(1.15);

% intertank fairing between tanks
A_frustrum_2 = pi*(r_LH2_2 + r_LOX_2)*sqrt((r_LH2_2 - r_LOX_2)^2+(r_LH2_2 + r_LOX_2)^2);
mi2 = 4.95*(A_frustrum_2)^(1.15);
Mass2.intertank_fairing = mi1 + mi2;

% aft fairing for bottom tank and thrusters
% add 2 meters for engine and misc
A_cylinder_2 = 2*pi*(r_LH2_2 + 2) * r_stage2;
Mass2.aft_fairing = 4.95*(A_cylinder_2)^(1.15);

% Mass of avionics
Mass2.avionics = 10*M0_2^(0.361);

% Mass of Wiring
Mass2.wiring = 1.058*sqrt(M0_2)*h_stage2^(.25);

% Engine mass

% engine thrust 
T = M0_2*g*TW_ratio_2/(n_engines_2); % N
% engine mass 
Mass2.engines = (7.81e-4*T + 3.37e-5*T*sqrt(area_ratio_2) + 59) * n_engines_2;
% mass of thrust structure 
Mass2.thrust_structure = 2.55e-4 * T * n_engines_2;

% gimbal mass
Mass2.gimbals = 237.8 * (T/P_2)^(0.9375) * n_engines_2;

%% Total Inert Mass
Masses.Stage2 = Mass2;
Table = struct2table(Mass2);
q = sum(Table(1,2:end), 'all').Variables;
% Mass2.Total_inert_mass = q;

Margin_2 = (Mass2.Inert - q)./q .* 100;
fprintf("Inert mass of stage 2 = %4.3f kg. This has a margin of %4.2f %%.\n", q, Margin_2)
% Mass2.Margin = Margin_2;

LDr_2 = h_stage2/D_2;
fprintf("Length/radius ratio: %4.3f.\n", LDr_2)
Mass2 = struct2table(Mass2)

%% STAGE 1 %%
% mass ratio between final and initial of first stage
r_1 = exp(-dV_1/ue_1);

% payload fraction // lambda
y_1 = r_1 - delta_1;

% Initial mass of stages
M0_1 = M0_2 / y_1;

% inert mass of stages
M_inert_1 = M0_1 * delta_1;

% Propellant Mass

Mp_1 = M0_1 * (1-r_1);

% need 15% of initial proplelant load must remain -> 15% of propellant lost
M_inert_1 = M_inert_1 * 1.15;
Mp_1 = Mp_1 * 1.15;
% find new inert mass
Mass1.Inert = M_inert_1;

% New initial mass
M0_1 = Mass1.Inert / delta_1;

% Mass of propellant
M_LOX_1 = OF_ratio_1/(OF_ratio_1 + 1) * Mp_1;
M_LH2_1 = 1/(OF_ratio_1 + 1) * Mp_1;

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

%% Fairings

% Fairing for stage 1 to 2 (or payload of stage 1)
A_frustrum = pi*(r_stage1 + r_stage2)*sqrt((r_stage1-r_stage2)^2+(r_stage1+r_stage2)^2);
Mass1.intertank_fairing = 4.95*(A_frustrum)^(1.15);

% aft_faring of stage 1
% add 2 meters for engine and misc
A_cylinder_1 = 2*pi*(r_LH2_1 + 2) * r_stage1;
Mass1.aft_fairing = 4.95*(A_cylinder_1)^(1.15);

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
Mass1.gimbals = 237.8 * (T/P_1)^(0.9375) * n_engines_1;

%% Total Inert Mass

Masses.Stage1 = Mass1;
Table = struct2table(Mass1);
q = sum(Table(1,2:end), 'all').Variables;
% Mass1.Total_inert_mass = q;

Margin_1 = (Mass1.Inert - q)./q .* 100;
fprintf("Inert mass of stage 1 = %4.3f kg. This has a margin of %4.2f %%.\n", q, Margin_1)
% Mass1.Margin = Margin_1;
LDr_1 = h_stage1/D_1;
fprintf("Length/radius ratio: %4.3f.\n", LDr_1)

Mass1 = struct2table(Mass1)
%% Totals

fprintf("Total L/D ratio: %4.3f.\n", LDr_1 + LDr_2 + LDr_3)
% find using weighted average
r_total = r_stage1 + r_stage2 + radius_TLI_stage;
h_total = h_stage1 + h_stage2 + h_stage3;
r_ave = (r_stage1*h_stage1 + r_stage2*h_stage2 + radius_TLI_stage*h_stage3) / (h_total);
fprintf("Total L/D ratio (weighted by height): %4.3f\n", h_total/(r_ave*2))

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
        fprintf("height of cylinder:  %4.3f \n", h_cyl)
        % surface area of tank (area of sphere ends + area of cylinder)
        A = 4*pi*stage_radius^2 + 2*pi*stage_radius*h_cyl;
        r = stage_radius;
    end
end
