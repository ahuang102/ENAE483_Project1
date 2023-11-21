% Gavin Bramble, Term 1 Project
%% 1. Staging Trade Study - 2 Stage
clear
clc

mpl = 25000; % [kg]
ve = 4273; % [m/s]
imf = .075;
dv1 = 9300; % [m/s]
dv2 = 3150; % [m/s]
cost_LOXLH2 = 1.102857143; % [$/kg]

syms mpr1 min1 mpr2 min2

eq1 = imf == min1/(min1+mpr1+min2+mpr2+mpl);
eq2 = imf == min2/(min2+mpr2+mpl);
eq3 = dv1 == -ve*log((min1+min2+mpr2+mpl)/(min1+mpr1+min2+mpr2+mpl));
eq4 = dv2 == -ve*log((min2+mpl)/(min2+mpr2+mpl));

arr = [eq1 eq2 eq3 eq4];
unknown = [mpr1 min1 mpr2 min2];
ans = vpasolve(arr,unknown);
mo = mpl + ans.mpr1 + ans.min1 + ans.mpr2 + ans.min2
mpr_tot = ans.mpr1 + ans.mpr2
min_tot = ans.min1 + ans.min2
cost_pr = mpr_tot*cost_LOXLH2
cost_in_nonrec = (12.73*min_tot^.55)*(10^6);
cost_in_unit1 = (.3024*min_tot^.662)*(10^6);
cost_in_total = cost_in_nonrec + cost_in_unit1
cost_total = cost_pr + cost_in_total
cost_per_kgmpl = cost_total/mpl
min1 = ans.min1
min2 = ans.min2
mpr1 = ans.mpr1
mpr2 = ans.mpr2


%% 2. Staging Trade Study - 3 Stage Non-Reusable (Split Stage1&2 dv in half)
clear
clc

mpl = 25000; % [kg]
ve = 4273; % [m/s]
imf = .075;
dv1 = 9300/2; % [m/s]
dv2 = 9300/2; % [m/s]
dv3 = 3150; % [m/s]
cost_LOXLH2 = 1.102857143; % [$/kg]

syms mpr1 min1 mpr2 min2 mpr3 min3

eq1 = imf == min1/(min1+mpr1+min2+mpr2+min3+mpr3+mpl);
eq2 = imf == min2/(min2+mpr2+min3+mpr3+mpl);
eq3 = imf == min3/(min3+mpr3+mpl);
eq4 = dv1 == -ve*log((mpl+min1+mpr2+min2+mpr3+min3)/(mpl+mpr1+min1+mpr2+min2+mpr3+min3));
eq5 = dv2 == -ve*log((mpl+min2+mpr3+min3)/(mpl+mpr2+min2+mpr3+min3));
eq6 = dv3 == -ve*log((mpl+min3)/(mpl+mpr3+min3));

arr = [eq1 eq2 eq3 eq4 eq5 eq6];
unknown = [mpr1 min1 mpr2 min2 mpr3 min3];
ans = vpasolve(arr,unknown);
mo = mpl + ans.mpr1 + ans.min1 + ans.mpr2 + ans.min2 + ans.mpr3 + ans.min3
mpr_tot = ans.mpr1 + ans.mpr2 + ans.mpr3
min_tot = ans.min1 + ans.min2 + ans.min3
cost_pr = mpr_tot*cost_LOXLH2
cost_in_nonrec = (12.73*min_tot^.55)*(10^6);
cost_in_unit1 = (.3024*min_tot^.662)*(10^6);
cost_in_total = cost_in_nonrec + cost_in_unit1
cost_total = cost_pr + cost_in_total
cost_per_kgmpl = cost_total/mpl
min1 = ans.min1
min2 = ans.min2
min3 = ans.min3
mpr1 = ans.mpr1
mpr2 = ans.mpr2
mpr3 = ans.mpr3



%% 3. Staging Trade Study - 4 Stage Non-Reusable (Stage1,2,3 dv equal split)
clear
clc

mpl = 25000; % [kg]
ve = 4273; % [m/s]
imf = .075;
dv1 = 9300/3; % [m/s]
dv2 = 9300/3; % [m/s]
dv3 = 9300/3; % [m/s]
dv4 = 3150; % [m/s]
cost_LOXLH2 = 1.102857143; % [$/kg]

syms mpr1 min1 mpr2 min2 mpr3 min3 mpr4 min4

eq1 = imf == min1/(min1+mpr1+min2+mpr2+min3+mpr3+min4+mpr4+mpl);
eq2 = imf == min2/(min2+mpr2+min3+mpr3+min4+mpr4+mpl);
eq3 = imf == min3/(min3+mpr3+min4+mpr4+mpl);
eq4 = imf == min4/(min4+mpr4+mpl);
eq5 = dv1 == -ve*log((mpl+min1+mpr2+min2+mpr3+min3+mpr4+min4)/(mpl+mpr1+min1+mpr2+min2+mpr3+min3+mpr4+min4));
eq6 = dv2 == -ve*log((mpl+min2+mpr3+min3+mpr4+min4)/(mpl+mpr2+min2+mpr3+min3+mpr4+min4));
eq7 = dv3 == -ve*log((mpl+min3+mpr4+min4)/(mpl+mpr3+min3+mpr4+min4));
eq8 = dv4 == -ve*log((mpl+min4)/(mpl+mpr4+min4));

arr = [eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8];
unknown = [mpr1 min1 mpr2 min2 mpr3 min3 mpr4 min4];
ans = vpasolve(arr,unknown);
mo = mpl + ans.mpr1 + ans.min1 + ans.mpr2 + ans.min2 + ans.mpr3 + ans.min3 + ans.mpr4 + ans.min4
mpr_tot = ans.mpr1 + ans.mpr2 + ans.mpr3 + ans.mpr4
min_tot = ans.min1 + ans.min2 + ans.min3 + ans.min4
cost_pr = mpr_tot*cost_LOXLH2
cost_in_nonrec = (12.73*min_tot^.55)*(10^6);
cost_in_unit1 = (.3024*min_tot^.662)*(10^6);
cost_in_total = cost_in_nonrec + cost_in_unit1
cost_total = cost_pr + cost_in_total
cost_per_kgmpl = cost_total/mpl
min1 = ans.min1
min2 = ans.min2
min3 = ans.min3
min4 = ans.min4
mpr1 = ans.mpr1
mpr2 = ans.mpr2
mpr3 = ans.mpr3
mpr4 = ans.mpr4


%% 4. Staging Trade Study - 3 Stage Reusable RTL (No dv optimization split)
clear
clc

mpl = 25000; % [kg]
ve = 4273; % [m/s]
imf = .075;
dv1 = 3000; % [m/s]
dv2 = 6300; % [m/s]
dv3 = 3150; % [m/s]
cost_LOXLH2 = 1.102857143; % [$/kg]

syms mpr1 min1 mpr2 min2 mpr3 min3

eq1 = imf == 1.1*min1/(1.1*min1+1.15*mpr1+min2+mpr2+min3+mpr3+mpl);
eq2 = imf == min2/(min2+mpr2+min3+mpr3+mpl);
eq3 = imf == min3/(min3+mpr3+mpl);
eq4 = dv1 == -ve*log((mpl+1.1*min1+mpr2+min2+mpr3+min3)/(mpl+1.15*mpr1+1.1*min1+mpr2+min2+mpr3+min3));
eq5 = dv2 == -ve*log((mpl+min2+mpr3+min3)/(mpl+mpr2+min2+mpr3+min3));
eq6 = dv3 == -ve*log((mpl+min3)/(mpl+mpr3+min3));

arr = [eq1 eq2 eq3 eq4 eq5 eq6];
unknown = [mpr1 min1 mpr2 min2 mpr3 min3];
ans = vpasolve(arr,unknown);
mo = mpl + ans.mpr1 + ans.min1 + ans.mpr2 + ans.min2 + ans.mpr3 + ans.min3
mpr_tot = ans.mpr1 + ans.mpr2 + ans.mpr3
min_tot = ans.min1 + ans.min2 + ans.min3
cost_pr = mpr_tot*cost_LOXLH2
cost_in_nonrec = (12.73*min_tot^.55)*(10^6);
cost_in_nonrec = (12.73*min_tot^.55)*(10^6);
cost_in_unit1 = (.3024*min_tot^.662)*(10^6);
cost_in_total = cost_in_nonrec + cost_in_unit1
cost_total = cost_pr + cost_in_total
cost_per_kgmpl = cost_total/mpl
min1 = ans.min1
min2 = ans.min2
min3 = ans.min3
mpr1 = ans.mpr1
mpr2 = ans.mpr2
mpr3 = ans.mpr3


%% 5. Staging Trade Study - 4 Stage Reusable RTL (Half dv split)
clear
clc

mpl = 25000; % [kg]
ve = 4273; % [m/s]
imf = .075;
dv1 = 3000; % [m/s]
dv2 = (9300-3000)/2; % [m/s]
dv3 = (9300-3000)/2; % [m/s]
dv4 = 3150; % [m/s]
cost_LOXLH2 = 1.102857143; % [$/kg]

syms mpr1 min1 mpr2 min2 mpr3 min3 mpr4 min4

eq1 = imf == 1.1*min1/(1.1*min1+1.15*mpr1+min2+mpr2+min3+mpr3+min4+mpr4+mpl);
eq2 = imf == min2/(min2+mpr2+min3+mpr3+min4+mpr4+mpl);
eq3 = imf == min3/(min3+mpr3+min4+mpr4+mpl);
eq4 = imf == min4/(min4+mpr4+mpl);
eq5 = dv1 == -ve*log((mpl+1.1*min1+mpr2+min2+mpr3+min3+mpr4+min4)/(mpl+1.15*mpr1+1.1*min1+mpr2+min2+mpr3+min3+mpr4+min4));
eq6 = dv2 == -ve*log((mpl+min2+mpr3+min3+mpr4+min4)/(mpl+mpr2+min2+mpr3+min3+mpr4+min4));
eq7 = dv3 == -ve*log((mpl+min3+mpr4+min4)/(mpl+mpr3+min3+mpr4+min4));
eq8 = dv4 == -ve*log((mpl+min4)/(mpl+mpr4+min4));

arr = [eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8];
unknown = [mpr1 min1 mpr2 min2 mpr3 min3 mpr4 min4];
ans = vpasolve(arr,unknown);
mo = mpl + ans.mpr1 + ans.min1 + ans.mpr2 + ans.min2 + ans.mpr3 + ans.min3 + ans.mpr4 + ans.min4
mpr_tot = ans.mpr1 + ans.mpr2 + ans.mpr3 + ans.mpr4
min_tot = ans.min1 + ans.min2 + ans.min3 + ans.min4
cost_pr = mpr_tot*cost_LOXLH2
cost_in_nonrec = (12.73*min_tot^.55)*(10^6);
cost_in_unit1 = (.3024*min_tot^.662)*(10^6);
cost_in_total = cost_in_nonrec + cost_in_unit1
cost_total = cost_pr + cost_in_total
cost_per_kgmpl = cost_total/mpl
min1 = ans.min1
min2 = ans.min2
min3 = ans.min3
min4 = ans.min4
mpr1 = ans.mpr1
mpr2 = ans.mpr2
mpr3 = ans.mpr3
mpr4 = ans.mpr4


%% 6. Staging Trade Study - 3 Stage Non-Reusable (Optimized dv split)
clear
clc

mpl = 25000; % [kg]
ve = 4273; % [m/s]
imf = .075;
dv = 9300; % [m/s]
dv3 = 3150; % [m/s]
cost_LOXLH2 = 1.102857143; % [$/kg]

syms mpr1 min1 mpr2 min2 mpr3 min3

min_cost_per_kgmpl = 5e6;
min_dv1 = 1;

for i = 1:dv
    dv1 = i;
    dv2 = dv - dv1;
    eq1 = imf == min1/(min1+mpr1+min2+mpr2+min3+mpr3+mpl);
    eq2 = imf == min2/(min2+mpr2+min3+mpr3+mpl);
    eq3 = imf == min3/(min3+mpr3+mpl);
    eq4 = dv1 == -ve*log((mpl+min1+mpr2+min2+mpr3+min3)/(mpl+mpr1+min1+mpr2+min2+mpr3+min3));
    eq5 = dv2 == -ve*log((mpl+min2+mpr3+min3)/(mpl+mpr2+min2+mpr3+min3));
    eq6 = dv3 == -ve*log((mpl+min3)/(mpl+mpr3+min3));
    arr = [eq1 eq2 eq3 eq4 eq5 eq6];
    unknown = [mpr1 min1 mpr2 min2 mpr3 min3];
    ans = vpasolve(arr,unknown);
    mo = mpl + ans.mpr1 + ans.min1 + ans.mpr2 + ans.min2 + ans.mpr3 + ans.min3;
    mpr_tot = ans.mpr1 + ans.mpr2 + ans.mpr3;
    min_tot = ans.min1 + ans.min2 + ans.min3;
    cost_pr = mpr_tot*cost_LOXLH2;
    cost_in_nonrec = (12.73*min_tot^.55)*(10^6);
    cost_in_unit1 = (.3024*min_tot^.662)*(10^6);
    cost_in_total = cost_in_nonrec + cost_in_unit1;
    cost_total = cost_pr + cost_in_total;
    cost_per_kgmpl = cost_total/mpl;
    if cost_per_kgmpl < min_cost_per_kgmpl
        min_mo = mo;
        min_mpr_tot = mpr_tot;
        min_min_tot = min_tot;
        min_cost_pr = cost_pr;
        min_cost_in_total = cost_in_total;
        min_cost_total = cost_total;
        min_cost_per_kgmpl = cost_per_kgmpl;
        min_dv1 = i;
        min_dv2 = dv2;
        min_min1 = ans.min1;
        min_min2 = ans.min2;
        min_min3 = ans.min3;
        min_mpr1 = ans.mpr1;
        min_mpr2 = ans.mpr2;
        min_mpr3 = ans.mpr3;
    end
    i
end

min_mo
min_mpr_tot
min_min_tot
min_cost_pr
min_cost_in_total
min_cost_total
min_cost_per_kgmpl
min_dv1
min_dv2
min_min1
min_min2
min_min3
min_mpr1
min_mpr2
min_mpr3


%% 7. Staging Trade Study - 3 Stage Reusable RTL (dv split optimized)
clear
clc

mpl = 25000; % [kg]
ve = 4273; % [m/s]
imf = .075;
dv = 3000; % [m/s]
dv_ = 9300; % [m/s]
dv3 = 3150; % [m/s]
cost_LOXLH2 = 1.102857143; % [$/kg]

syms mpr1 min1 mpr2 min2 mpr3 min3

min_cost_per_kgmpl = 5e6;
min_dv1 = 1;

for i = 1:dv
    dv1 = i;
    dv2 = dv_ - dv1;
    eq1 = imf == 1.1*min1/(1.1*min1+1.15*mpr1+min2+mpr2+min3+mpr3+mpl);
    eq2 = imf == min2/(min2+mpr2+min3+mpr3+mpl);
    eq3 = imf == min3/(min3+mpr3+mpl);
    eq4 = dv1 == -ve*log((mpl+1.1*min1+mpr2+min2+mpr3+min3)/(mpl+1.15*mpr1+1.1*min1+mpr2+min2+mpr3+min3));
    eq5 = dv2 == -ve*log((mpl+min2+mpr3+min3)/(mpl+mpr2+min2+mpr3+min3));
    eq6 = dv3 == -ve*log((mpl+min3)/(mpl+mpr3+min3));
    arr = [eq1 eq2 eq3 eq4 eq5 eq6];
    unknown = [mpr1 min1 mpr2 min2 mpr3 min3];
    ans = vpasolve(arr,unknown);
    mo = mpl + ans.mpr1 + ans.min1 + ans.mpr2 + ans.min2 + ans.mpr3 + ans.min3;
    mpr_tot = ans.mpr1 + ans.mpr2 + ans.mpr3;
    min_tot = ans.min1 + ans.min2 + ans.min3;
    cost_pr = mpr_tot*cost_LOXLH2;
    cost_in_nonrec = (12.73*min_tot^.55)*(10^6);
    cost_in_nonrec = (12.73*min_tot^.55)*(10^6);
    cost_in_unit1 = (.3024*min_tot^.662)*(10^6);
    cost_in_total = cost_in_nonrec + cost_in_unit1;
    cost_total = cost_pr + cost_in_total;
    cost_per_kgmpl = cost_total/mpl;
    if cost_per_kgmpl < min_cost_per_kgmpl
        min_mo = mo;
        min_mpr_tot = mpr_tot;
        min_min_tot = min_tot;
        min_cost_pr = cost_pr;
        min_cost_in_total = cost_in_total;
        min_cost_total = cost_total;
        min_cost_per_kgmpl = cost_per_kgmpl;
        min_dv1 = i;
        min_dv2 = dv2;
        min_min1 = ans.min1;
        min_min2 = ans.min2;
        min_min3 = ans.min3;
        min_mpr1 = ans.mpr1;
        min_mpr2 = ans.mpr2;
        min_mpr3 = ans.mpr3;
    end
    i
end

min_mo
min_mpr_tot
min_min_tot
min_cost_pr
min_cost_in_total
min_cost_total
min_cost_per_kgmpl
min_dv1
min_dv2
min_min1
min_min2
min_min3
min_mpr1
min_mpr2
min_mpr3


%% 8. Staging Trade Study - 4 Stage Reusable RTL (dv split optimized)
clear
clc

mpl = 25000; % [kg]
ve = 4273; % [m/s]
imf = .075;
dv1 = 3000; % [m/s]
dv = 9300 - dv1; % [m/s]
dv4 = 3150; % [m/s]
cost_LOXLH2 = 1.102857143; % [$/kg]

syms mpr1 min1 mpr2 min2 mpr3 min3 mpr4 min4

min_cost_per_kgmpl = 5e6;
min_dv1 = 1;

for i = 1:dv
    dv2 = i;
    dv3 = dv - dv2;
    eq1 = imf == 1.1*min1/(1.1*min1+1.15*mpr1+min2+mpr2+min3+mpr3+min4+mpr4+mpl);
    eq2 = imf == min2/(min2+mpr2+min3+mpr3+min4+mpr4+mpl);
    eq3 = imf == min3/(min3+mpr3+min4+mpr4+mpl);
    eq4 = imf == min4/(min4+mpr4+mpl);
    eq5 = dv1 == -ve*log((mpl+1.1*min1+mpr2+min2+mpr3+min3+mpr4+min4)/(mpl+1.15*mpr1+1.1*min1+mpr2+min2+mpr3+min3+mpr4+min4));
    eq6 = dv2 == -ve*log((mpl+min2+mpr3+min3+mpr4+min4)/(mpl+mpr2+min2+mpr3+min3+mpr4+min4));
    eq7 = dv3 == -ve*log((mpl+min3+mpr4+min4)/(mpl+mpr3+min3+mpr4+min4));
    eq8 = dv4 == -ve*log((mpl+min4)/(mpl+mpr4+min4));
    arr = [eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8];
    unknown = [mpr1 min1 mpr2 min2 mpr3 min3 mpr4 min4];
    ans = vpasolve(arr,unknown);
    mo = mpl + ans.mpr1 + ans.min1 + ans.mpr2 + ans.min2 + ans.mpr3 + ans.min3 + ans.mpr4 + ans.min4;
    mpr_tot = ans.mpr1 + ans.mpr2 + ans.mpr3 + ans.mpr4;
    min_tot = ans.min1 + ans.min2 + ans.min3 + ans.min4;
    cost_pr = mpr_tot*cost_LOXLH2;
    cost_in_nonrec = (12.73*min_tot^.55)*(10^6);
    cost_in_unit1 = (.3024*min_tot^.662)*(10^6);
    cost_in_total = cost_in_nonrec + cost_in_unit1;
    cost_total = cost_pr + cost_in_total;
    cost_per_kgmpl = cost_total/mpl;
    if cost_per_kgmpl < min_cost_per_kgmpl
        min_mo = mo;
        min_mpr_tot = mpr_tot;
        min_min_tot = min_tot;
        min_cost_pr = cost_pr;
        min_cost_in_total = cost_in_total;
        min_cost_total = cost_total;
        min_cost_per_kgmpl = cost_per_kgmpl;
        min_dv2 = i;
        min_dv3 = dv3;
        min_min1 = ans.min1;
        min_min2 = ans.min2;
        min_min3 = ans.min3;
        min_min4 = ans.min4;
        min_mpr1 = ans.mpr1;
        min_mpr2 = ans.mpr2;
        min_mpr3 = ans.mpr3;
        min_mpr4 = ans.mpr4;
    end
    i
end

min_mo
min_mpr_tot
min_min_tot
min_cost_pr
min_cost_in_total
min_cost_total
min_cost_per_kgmpl
min_dv2
min_dv3
min_min1
min_min2
min_min3
min_min4
min_mpr1
min_mpr2
min_mpr3
min_mpr4


%% 9. Staging Trade Study - 3 Stage Reusable Drone Ship (No dv optimization split)
clear
clc

mpl = 25000; % [kg]
ve = 4273; % [m/s]
imf = .075;
dv1 = 3500; % [m/s]
dv2 = 5800; % [m/s]
dv3 = 3150; % [m/s]
cost_LOXLH2 = 1.102857143; % [$/kg]

syms mpr1 min1 mpr2 min2 mpr3 min3

eq1 = imf == 1.1*min1/(1.1*min1+1.05*mpr1+min2+mpr2+min3+mpr3+mpl);
eq2 = imf == min2/(min2+mpr2+min3+mpr3+mpl);
eq3 = imf == min3/(min3+mpr3+mpl);
eq4 = dv1 == -ve*log((mpl+1.1*min1+mpr2+min2+mpr3+min3)/(mpl+1.05*mpr1+1.1*min1+mpr2+min2+mpr3+min3));
eq5 = dv2 == -ve*log((mpl+min2+mpr3+min3)/(mpl+mpr2+min2+mpr3+min3));
eq6 = dv3 == -ve*log((mpl+min3)/(mpl+mpr3+min3));

arr = [eq1 eq2 eq3 eq4 eq5 eq6];
unknown = [mpr1 min1 mpr2 min2 mpr3 min3];
ans = vpasolve(arr,unknown);
mo = mpl + ans.mpr1 + ans.min1 + ans.mpr2 + ans.min2 + ans.mpr3 + ans.min3
mpr_tot = ans.mpr1 + ans.mpr2 + ans.mpr3
min_tot = ans.min1 + ans.min2 + ans.min3
cost_pr = mpr_tot*cost_LOXLH2
cost_in_nonrec = (12.73*min_tot^.55)*(10^6);
cost_in_nonrec = (12.73*min_tot^.55)*(10^6);
cost_in_unit1 = (.3024*min_tot^.662)*(10^6);
cost_in_total = cost_in_nonrec + cost_in_unit1
cost_total = cost_pr + cost_in_total
cost_per_kgmpl = cost_total/mpl
min1 = ans.min1
min2 = ans.min2
min3 = ans.min3
mpr1 = ans.mpr1
mpr2 = ans.mpr2
mpr3 = ans.mpr3


%% 10. Staging Trade Study - 4 Stage Reusable Drone Ship (Half dv split)
clear
clc

mpl = 25000; % [kg]
ve = 4273; % [m/s]
imf = .075;
dv1 = 3500; % [m/s]
dv2 = (9300-3500)/2; % [m/s]
dv3 = (9300-3500)/2; % [m/s]
dv4 = 3150; % [m/s]
cost_LOXLH2 = 1.102857143; % [$/kg]

syms mpr1 min1 mpr2 min2 mpr3 min3 mpr4 min4

eq1 = imf == 1.1*min1/(1.1*min1+1.05*mpr1+min2+mpr2+min3+mpr3+min4+mpr4+mpl);
eq2 = imf == min2/(min2+mpr2+min3+mpr3+min4+mpr4+mpl);
eq3 = imf == min3/(min3+mpr3+min4+mpr4+mpl);
eq4 = imf == min4/(min4+mpr4+mpl);
eq5 = dv1 == -ve*log((mpl+1.1*min1+mpr2+min2+mpr3+min3+mpr4+min4)/(mpl+1.05*mpr1+1.1*min1+mpr2+min2+mpr3+min3+mpr4+min4));
eq6 = dv2 == -ve*log((mpl+min2+mpr3+min3+mpr4+min4)/(mpl+mpr2+min2+mpr3+min3+mpr4+min4));
eq7 = dv3 == -ve*log((mpl+min3+mpr4+min4)/(mpl+mpr3+min3+mpr4+min4));
eq8 = dv4 == -ve*log((mpl+min4)/(mpl+mpr4+min4));

arr = [eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8];
unknown = [mpr1 min1 mpr2 min2 mpr3 min3 mpr4 min4];
ans = vpasolve(arr,unknown);
mo = mpl + ans.mpr1 + ans.min1 + ans.mpr2 + ans.min2 + ans.mpr3 + ans.min3 + ans.mpr4 + ans.min4
mpr_tot = ans.mpr1 + ans.mpr2 + ans.mpr3 + ans.mpr4
min_tot = ans.min1 + ans.min2 + ans.min3 + ans.min4
cost_pr = mpr_tot*cost_LOXLH2
cost_in_nonrec = (12.73*min_tot^.55)*(10^6);
cost_in_unit1 = (.3024*min_tot^.662)*(10^6);
cost_in_total = cost_in_nonrec + cost_in_unit1
cost_total = cost_pr + cost_in_total
cost_per_kgmpl = cost_total/mpl
min1 = ans.min1
min2 = ans.min2
min3 = ans.min3
min4 = ans.min4
mpr1 = ans.mpr1
mpr2 = ans.mpr2
mpr3 = ans.mpr3
mpr4 = ans.mpr4


%% 11. Staging Trade Study - 3 Stage Reusable Ship (dv split optimized)
clear
clc

mpl = 25000; % [kg]
ve = 4273; % [m/s]
imf = .075;
dv = 3500; % [m/s]
dv_ = 9300; % [m/s]
dv3 = 3150; % [m/s]
cost_LOXLH2 = 1.102857143; % [$/kg]

syms mpr1 min1 mpr2 min2 mpr3 min3

min_cost_per_kgmpl = 5e6;
min_dv1 = 1;

for i = 1:dv
    dv1 = i;
    dv2 = dv_ - dv1;
    eq1 = imf == 1.1*min1/(1.1*min1+1.05*mpr1+min2+mpr2+min3+mpr3+mpl);
    eq2 = imf == min2/(min2+mpr2+min3+mpr3+mpl);
    eq3 = imf == min3/(min3+mpr3+mpl);
    eq4 = dv1 == -ve*log((mpl+1.1*min1+mpr2+min2+mpr3+min3)/(mpl+1.05*mpr1+1.1*min1+mpr2+min2+mpr3+min3));
    eq5 = dv2 == -ve*log((mpl+min2+mpr3+min3)/(mpl+mpr2+min2+mpr3+min3));
    eq6 = dv3 == -ve*log((mpl+min3)/(mpl+mpr3+min3));
    arr = [eq1 eq2 eq3 eq4 eq5 eq6];
    unknown = [mpr1 min1 mpr2 min2 mpr3 min3];
    ans = vpasolve(arr,unknown);
    mo = mpl + ans.mpr1 + ans.min1 + ans.mpr2 + ans.min2 + ans.mpr3 + ans.min3;
    mpr_tot = ans.mpr1 + ans.mpr2 + ans.mpr3;
    min_tot = ans.min1 + ans.min2 + ans.min3;
    cost_pr = mpr_tot*cost_LOXLH2;
    cost_in_nonrec = (12.73*min_tot^.55)*(10^6);
    cost_in_nonrec = (12.73*min_tot^.55)*(10^6);
    cost_in_unit1 = (.3024*min_tot^.662)*(10^6);
    cost_in_total = cost_in_nonrec + cost_in_unit1;
    cost_total = cost_pr + cost_in_total;
    cost_per_kgmpl = cost_total/mpl;
    if cost_per_kgmpl < min_cost_per_kgmpl
        min_mo = mo;
        min_mpr_tot = mpr_tot;
        min_min_tot = min_tot;
        min_cost_pr = cost_pr;
        min_cost_in_total = cost_in_total;
        min_cost_total = cost_total;
        min_cost_per_kgmpl = cost_per_kgmpl;
        min_dv1 = i;
        min_dv2 = dv2;
        min_min1 = ans.min1;
        min_min2 = ans.min2;
        min_min3 = ans.min3;
        min_mpr1 = ans.mpr1;
        min_mpr2 = ans.mpr2;
        min_mpr3 = ans.mpr3;
    end
    i
end

min_mo
min_mpr_tot
min_min_tot
min_cost_pr
min_cost_in_total
min_cost_total
min_cost_per_kgmpl
min_dv1
min_dv2
min_min1
min_min2
min_min3
min_mpr1
min_mpr2
min_mpr3


%% 12. Staging Trade Study - 4 Stage Reusable Drone Ship (dv split optimized)
clear
clc

mpl = 25000; % [kg]
ve = 4273; % [m/s]
imf = .075;
dv1 = 3500; % [m/s]
dv = 9300 - dv1; % [m/s]
dv4 = 3150; % [m/s]
cost_LOXLH2 = 1.102857143; % [$/kg]

syms mpr1 min1 mpr2 min2 mpr3 min3 mpr4 min4

min_cost_per_kgmpl = 5e6;
min_dv1 = 1;

for i = 1:dv
    dv2 = i;
    dv3 = dv - dv2;
    eq1 = imf == 1.1*min1/(1.1*min1+1.05*mpr1+min2+mpr2+min3+mpr3+min4+mpr4+mpl);
    eq2 = imf == min2/(min2+mpr2+min3+mpr3+min4+mpr4+mpl);
    eq3 = imf == min3/(min3+mpr3+min4+mpr4+mpl);
    eq4 = imf == min4/(min4+mpr4+mpl);
    eq5 = dv1 == -ve*log((mpl+1.1*min1+mpr2+min2+mpr3+min3+mpr4+min4)/(mpl+1.05*mpr1+1.1*min1+mpr2+min2+mpr3+min3+mpr4+min4));
    eq6 = dv2 == -ve*log((mpl+min2+mpr3+min3+mpr4+min4)/(mpl+mpr2+min2+mpr3+min3+mpr4+min4));
    eq7 = dv3 == -ve*log((mpl+min3+mpr4+min4)/(mpl+mpr3+min3+mpr4+min4));
    eq8 = dv4 == -ve*log((mpl+min4)/(mpl+mpr4+min4));
    arr = [eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8];
    unknown = [mpr1 min1 mpr2 min2 mpr3 min3 mpr4 min4];
    ans = vpasolve(arr,unknown);
    mo = mpl + ans.mpr1 + ans.min1 + ans.mpr2 + ans.min2 + ans.mpr3 + ans.min3 + ans.mpr4 + ans.min4;
    mpr_tot = ans.mpr1 + ans.mpr2 + ans.mpr3 + ans.mpr4;
    min_tot = ans.min1 + ans.min2 + ans.min3 + ans.min4;
    cost_pr = mpr_tot*cost_LOXLH2;
    cost_in_nonrec = (12.73*min_tot^.55)*(10^6);
    cost_in_unit1 = (.3024*min_tot^.662)*(10^6);
    cost_in_total = cost_in_nonrec + cost_in_unit1;
    cost_total = cost_pr + cost_in_total;
    cost_per_kgmpl = cost_total/mpl;
    if cost_per_kgmpl < min_cost_per_kgmpl
        min_mo = mo;
        min_mpr_tot = mpr_tot;
        min_min_tot = min_tot;
        min_cost_pr = cost_pr;
        min_cost_in_total = cost_in_total;
        min_cost_total = cost_total;
        min_cost_per_kgmpl = cost_per_kgmpl;
        min_dv2 = i;
        min_dv3 = dv3;
        min_min1 = ans.min1;
        min_min2 = ans.min2;
        min_min3 = ans.min3;
        min_min4 = ans.min4;
        min_mpr1 = ans.mpr1;
        min_mpr2 = ans.mpr2;
        min_mpr3 = ans.mpr3;
        min_mpr4 = ans.mpr4;
    end
    i
end

min_mo
min_mpr_tot
min_min_tot
min_cost_pr
min_cost_in_total
min_cost_total
min_cost_per_kgmpl
min_dv2
min_dv3
min_min1
min_min2
min_min3
min_min4
min_mpr1
min_mpr2
min_mpr3
min_mpr4


%% Resiliency
clear
clc

r = 6;          % [missions/year]
d1 = 6/12;      % [years]
d2 = 9/12;      % [years]
k = 0.8;        % [unitless]
surge = 9;      % [missions/year]
S = surge/r;    % [unitless]

m1 = (S*r*k*d1)/(S-1) % [missions]
m2 = (S*r*k*d2)/(S-1) % [missions]