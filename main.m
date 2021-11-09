clear; clc; close all;

index  = [0 0 0 0 0; 0 1 1 1 0; 0 1 2 1 0; 0 1 1 1 0; 0 0 0 0 0];
index2 = [0 0 0 0 0; 0 1 1 1 0; 0 1 1 1 0; 0 1 1 1 0; 0 0 0 0 0];
order  = [0 0 0 0 0; 0 1 2 3 0; 0 4 5 6 0; 0 7 8 9 0; 0 0 0 0 0];
num_act = max(max(order));                % total number of active blocks
[u, v] = size(order);

dx = 500 * index2; dy = 500 * index2; h = 50 * index2; Vb = dx.*dy.*h;
k_x = 100 * index2; k_y = 100 * index2; rw = 0.25; beta_c = 1.127e-3;
phi = 0.2*index2; skin = 0; q_sp = -20; dt = 5;

Po_0 = 4800*index2;
sat_w_0 = 0.3 * index2; 
sat_g_0 = 0.1 * index2; 
Bo_0 = zeros(u,v); Bw_0 = zeros(u,v); Bg_0 = zeros(u,v); Rso_0 = zeros(u,v);
Pw_0 = zeros(u,v); Pg_0 = zeros(u,v); Pcow_0 = zeros(u,v); Pcgo_0 = zeros(u,v);

for i = 1 : u
    for j = 1 :v
        if index2(i,j) ~= 0
            [~, Bo_0(i,j), ~, Rso_0(i,j)] = Oil_table(Po_0(i,j));
            [~, ~, Pcow_0(i,j)] = Sw_table(sat_w_0(i,j));
            [~, ~, Pcgo_0(i,j)] = Sg_table(sat_g_0(i,j));
            Pw_0(i,j) = Po_0(i,j) - Pcow_0(i,j);
            Pg_0(i,j) = Po_0(i,j) + Pcgo_0(i,j);
            [~, Bw_0(i,j), ~, ~] = Water_table(Pw_0(i,j));
            [~, Bg_0(i,j), ~] = Gas_table(Pg_0(i,j));
        end
    end
end

Po_old = Po_0;   sat_w_old = sat_w_0;   sat_g_old = sat_g_0;

Bo_old = zeros(u,v); Bw_old = zeros(u,v); Bg_old = zeros(u,v); Rso_old = zeros(u,v);
Pw_old = zeros(u,v); Pg_old = zeros(u,v); Pcow_old = zeros(u,v); Pcgo_old = zeros(u,v);

%% Begin iterations:

qo_sc_record = [];
qw_sc_record = [];
qg_sc_record = [];

for timestep = 1 : 1           % 1 time step = 5 days
    for i = 1 : u
        for j = 1 : v
            if index2(i,j) ~= 0
                [~, Bo_old(i,j), ~, Rso_old(i,j)] = Oil_table(Po_old(i,j));
                [~, ~, Pcow_old(i,j)] = Sw_table(sat_w_old(i,j));
                [~, ~, Pcgo_old(i,j)] = Sg_table(sat_g_old(i,j));
                Pw_old(i,j) = Po_old(i,j) - Pcow_old(i,j);
                Pg_old(i,j) = Po_old(i,j) + Pcgo_old(i,j);
                [~, Bw_old(i,j), ~, ~] = Water_table(Pw_old(i,j));
                [~, Bg_old(i,j), ~] = Gas_table(Pg_old(i,j));
            end
        end
    end

    %% Initialize value of Po, sat_w, sat_g:
    Po = Po_old;
    sat_w = sat_w_old;
    sat_g = sat_g_old;

    %% Iterations:

    Jacob_ = zeros(3*num_act, 3*num_act);

    for n_iter = 1 : 3
        Jacob_1st_columns = diff_p(Po, sat_w, sat_g, Po_old, sat_w_old, sat_g_old);
        Jacob_2nd_columns = diff_sat_w(Po, sat_w, sat_g, Po_old, sat_w_old, sat_g_old);
        Jacob_3rd_columns = diff_sat_g(Po, sat_w, sat_g, Po_old, sat_w_old, sat_g_old);

        for i = 1 : num_act
            Jacob_(:, 3*i-2) = Jacob_1st_columns(:, i);
            Jacob_(:, 3*i-1) = Jacob_2nd_columns(:, i);
            Jacob_(:, 3*i-0) = Jacob_3rd_columns(:, i);
        end

        RHS = Residual(Po, sat_w, sat_g, Po_old, sat_w_old, sat_g_old);

        delta_update = Jacob_ \ RHS;
%         diff = max(abs(delta_update));

        for i = 1 : u
            for j = 1 : v
                if order(i,j) ~= 0
                    Po(i,j) = Po(i,j) - delta_update(3*order(i,j)-2);
                    sat_w(i,j) = sat_w(i,j) - delta_update(3*order(i,j)-1);
                    sat_g(i,j) = sat_g(i,j) - delta_update(3*order(i,j));
                end
            end
        end

    end

    Bo = zeros(u,v); Pcow = zeros(u,v); Pcgo = zeros(u,v);
    Pw = zeros(u,v); Pg = zeros(u,v); Bg = zeros(u,v); Bw = zeros(u,v);
    mu_o = zeros(u,v); mu_w = zeros(u,v); mu_g = zeros(u,v); Rso = zeros(u,v);
    krw = zeros(u,v); krow = zeros(u,v); krg = zeros(u,v); krog = zeros(u,v); kro = zeros(u,v);
    for i = 1 : u
        for j = 1 : v
            if index2(i,j) ~= 0
                [~, Bo(i,j), mu_o(i,j), Rso(i,j)] = Oil_table(Po(i,j));
                [krw(i,j), krow(i,j), Pcow(i,j)] = Sw_table(sat_w(i,j));
                [krg(i,j), krog(i,j), Pcgo(i,j)] = Sg_table(sat_g(i,j));
                kro(i,j) = (krw(i,j) + krow(i,j))*(krg(i,j) + krog(i,j)) - (krw(i,j) + krg(i,j));
                Pw(i,j) = Po(i,j) - Pcow(i,j);
                Pg(i,j) = Po(i,j) + Pcgo(i,j);
                [~, Bw(i,j), mu_w(i,j), ~] = Water_table(Pw(i,j));
                [~, Bg(i,j), mu_g(i,j)] = Gas_table(Pg(i,j));
            end
        end
    end
    
    Po_old = Po; 
    sat_w_old = sat_w;
    sat_g_old = sat_g;          % timestep = timestep + 1
    
    %% Calculate qw_sc and qg_sc for Material Balance Check:
    for i = 1 : u
        for j = 1 : v
            if index(i,j) == 2
                re = 0.198*dx(i,j);
                Gw = 2 * pi * beta_c * sqrt(k_x(i,j)*k_y(i,j))/(log(re/rw)+skin);
                Psf = Po(i,j) + q_sp/(Gw*kro(i,j)/mu_o(i,j)/Bo(i,j));

                qo_sc = q_sp;
                qw_sc = -Gw*krw(i,j)/mu_w(i,j)/Bw(i,j)*(Po(i,j) - Psf);
                qg_sc = -Gw*(krg(i,j)/mu_g(i,j)/Bg(i,j)+ Rso(i,j)*kro(i,j)/mu_o(i,j)/Bo(i,j))*(Po(i,j) - Psf);
            end
        end
    end
    qo_sc_record = [qo_sc_record; qo_sc];
    qw_sc_record = [qw_sc_record; qw_sc];
    qg_sc_record = [qg_sc_record; qg_sc];
end


%% Cumulative Material Balance Check:
CMBC_oil = 0;
CMBC_water = 0;
CMBC_gas = 0;
for i = 1 : u
    for j = 1 : v
        if order(i,j) > 0
            CMBC_oil = CMBC_oil + Vb(i,j)/5.615 * (phi(i,j)*(1-sat_w(i,j)-sat_g(i,j))/Bo(i,j) ...
                - phi(i,j)*(1-sat_w_0(i,j)-sat_g_0(i,j))/Bo_0(i,j))/ sum(qo_sc_record) / dt ;
            CMBC_water = CMBC_water + Vb(i,j)/5.615 * (phi(i,j)*sat_w(i,j)/Bw(i,j) ...
                - phi(i,j)*sat_w_0(i,j)/Bw_0(i,j))/ sum(qw_sc_record) /dt ;
            CMBC_gas = CMBC_gas + Vb(i,j)/5.615 * phi(i,j)*(sat_g(i,j)/Bg(i,j) + Rso(i,j)*(1-sat_w(i,j)-sat_g(i,j))/Bo(i,j))/sum(qg_sc_record)/dt ...
                - Vb(i,j)/5.615 * phi(i,j)*(sat_g_0(i,j)/Bg_0(i,j) + Rso_0(i,j)*(1-sat_w_0(i,j)-sat_g_0(i,j))/Bo_0(i,j))/ sum(qg_sc_record) /dt ;
        end
    end
end