function [Jacob_1st_col] = diff_p(Po, sat_w, sat_g, Po_old, sat_w_old, sat_g_old)
% Calculate Residual of every reservoir blocks
index  = [0 0 0 0 0; 0 1 1 1 0; 0 1 2 1 0; 0 1 1 1 0; 0 0 0 0 0];
index2 = [0 0 0 0 0; 0 1 1 1 0; 0 1 1 1 0; 0 1 1 1 0; 0 0 0 0 0];
order  = [0 0 0 0 0; 0 1 2 3 0; 0 4 5 6 0; 0 7 8 9 0; 0 0 0 0 0];
num_act = max(max(order));      %total number of active blocks
Jacob_1st_col = zeros(3*num_act, num_act);

dx = 500 * index2; dy = 500 * index2; h = 50 * index2; G = 4000 * index2;
Ax = h .* dy; Ay = h .* dx; Vb = h .* dx .* dy;
k_x = 100 * index2; k_y = 100 * index2; [u, v] = size(index);
phi = 0.2 * index2; rw = 0.25; skin = 0; beta_c = 1.127e-3;
q_sp = -20; dt = 5; eps_p = 1e-3;
format long
    
%% Initial reservoir pressure and saturation:
rho_o_old = zeros(u,v); Bo_old = zeros(u,v); mu_o_old = zeros(u,v); Rso_old = zeros(u,v);
for i = 1 : u
    for j = 1 : v
        if index(i,j) ~= 0
            [rho_o_old(i,j), Bo_old(i,j), mu_o_old(i,j), Rso_old(i,j)] = Oil_table(Po_old(i,j));
        end
    end
end

krw_old = zeros(u,v); krow_old = zeros(u,v); Pcow_old = zeros(u,v);
krg_old = zeros(u,v); krog_old = zeros(u,v); Pcgo_old = zeros(u,v);
for i = 1 : u
    for j = 1 : v
        if index(i,j) ~= 0
            [krw_old(i,j), krow_old(i,j), Pcow_old(i,j)] = Sw_table(sat_w_old(i,j));
            [krg_old(i,j), krog_old(i,j), Pcgo_old(i,j)] = Sg_table(sat_g_old(i,j));
        end
    end
end

Pw_old = Po_old - Pcow_old;   Pg_old = Po_old + Pcgo_old;
rho_w_old = zeros(u,v); Bw_old = zeros(u,v); mu_w_old = zeros(u,v); 
rho_g_old = zeros(u,v); Bg_old = zeros(u,v); mu_g_old = zeros(u,v); 
for i = 1 : u 
    for j = 1 : v
        if index(i,j) ~= 0
            [rho_w_old(i,j), Bw_old(i,j), mu_w_old(i,j), ~] = Water_table(Pw_old(i,j));
            [rho_g_old(i,j), Bg_old(i,j), mu_g_old(i,j)] = Gas_table(Pg_old(i,j));
        end
    end
end
        
%% Values of Po, sat_w, sat_g:
rho_o = zeros(u,v); Bo = zeros(u,v); mu_o = zeros(u,v); Rso = zeros(u,v);
krw = zeros(u,v); krow = zeros(u,v); Pcow = zeros(u,v); kro = zeros(u,v);
krg = zeros(u,v); krog = zeros(u,v); Pcgo = zeros(u,v);
for i = 1 : u
    for j = 1 : v
        if index(i,j) ~= 0
            [rho_o(i,j), Bo(i,j), mu_o(i,j), Rso(i,j)] = Oil_table(Po(i,j));
            [krw(i,j), krow(i,j), Pcow(i,j)] = Sw_table(sat_w(i,j));
            [krg(i,j), krog(i,j), Pcgo(i,j)] = Sg_table(sat_g(i,j));
            kro(i,j) = (krow(i,j) + krw(i,j))*(krog(i,j)+krg(i,j)) - (krw(i,j) + krg(i,j));
        end
    end
end

Pw = Po - Pcow;   Pg = Po + Pcgo;
rho_w = zeros(u,v); Bw = zeros(u,v); mu_w = zeros(u,v); 
rho_g = zeros(u,v); Bg = zeros(u,v); mu_g = zeros(u,v); 
for i = 1 : u 
    for j = 1 : v
        if index(i,j) ~= 0
            [rho_w(i,j), Bw(i,j), mu_w(i,j), ~] = Water_table(Pw(i,j));
            [rho_g(i,j), Bg(i,j), mu_g(i,j)] = Gas_table(Pg(i,j));
        end
    end
end

%% Calculate Residual:
Ro = zeros(num_act, 1);
Rw = zeros(num_act, 1);
Rg = zeros(num_act, 1);

for i = 1 : u
    for j = 1 : v
        if index(i,j) ~= 0
            if index(i,j+1) ~= 0
                Po_avg = (Po(i,j) + Po(i,j+1))/2;
                Pw_avg = (Pw(i,j) + Pw(i,j+1))/2;
                Pg_avg = (Pg(i,j) + Pg(i,j+1))/2;
                [rho_o_avg,Bo_avg,mu_o_avg,Rso_avg] = Oil_table(Po_avg);
                [rho_w_avg,Bw_avg,mu_w_avg,~] = Water_table(Pw_avg);
                [rho_g_avg, Bg_avg, mu_g_avg] = Gas_table(Pg_avg);
                
                if Po(i,j)-rho_o(i,j)/144*G(i,j) >= Po(i,j+1)-rho_o(i,j+1)/144*G(i,j+1)
                    kro_up = kro(i,j);
                else; kro_up = kro(i,j+1);
                end
                
                if Pw(i,j)-rho_w(i,j)/144*G(i,j) >= Pw(i,j+1)-rho_w(i,j+1)/144*G(i,j+1)
                    krw_up = krw(i,j);
                else; krw_up = krw(i,j+1);
                end
                
                if Pg(i,j)-rho_g(i,j)/144*G(i,j) >= Pg(i,j+1)-rho_g(i,j+1)/144*G(i,j+1)
                    krg_up = krg(i,j);
                else; krg_up = krg(i,j+1);
                end
                
                Eo = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j) * kro_up/mu_o_avg/Bo_avg;
                Eo_prime = Eo /144 * rho_o_avg;
                Ew = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j) * krw_up/mu_w_avg/Bw_avg;
                Ew_prime = Ew / 144 * rho_w_avg;
                Eg = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j) * krg_up/mu_g_avg/Bg_avg;
                Eg_prime = Eg / 144 * rho_g_avg;
                Ego = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j)* kro_up*Rso_avg/mu_o_avg/Bo_avg;
                Ego_prime = Ego /144 * rho_o_avg;
            else
                Eo = 0; Ew = 0; Eg = 0; Ego = 0; Eo_prime = 0; Ew_prime = 0; Eg_prime = 0; Ego_prime = 0;
            end
            
            if index(i,j-1) ~= 0
                Po_avg = (Po(i,j) + Po(i,j-1))/2;
                Pw_avg = (Pw(i,j) + Pw(i,j-1))/2;
                Pg_avg = (Pg(i,j) + Pg(i,j-1))/2;
                [rho_o_avg,Bo_avg,mu_o_avg,Rso_avg] = Oil_table(Po_avg);
                [rho_w_avg,Bw_avg,mu_w_avg,~] = Water_table(Pw_avg);
                [rho_g_avg, Bg_avg, mu_g_avg] = Gas_table(Pg_avg);
                
                if Po(i,j)-rho_o(i,j)/144*G(i,j) >= Po(i,j-1)-rho_o(i,j-1)/144*G(i,j-1)
                    kro_up = kro(i,j);
                else; kro_up = kro(i,j-1);
                end
                
                if Pw(i,j)-rho_w(i,j)/144*G(i,j) >= Pw(i,j-1)-rho_w(i,j-1)/144*G(i,j-1)
                    krw_up = krw(i,j);
                else; krw_up = krw(i,j-1);
                end
                
                if Pg(i,j)-rho_g(i,j)/144*G(i,j) >= Pg(i,j-1)-rho_g(i,j-1)/144*G(i,j-1)
                    krg_up = krg(i,j);
                else; krg_up = krg(i,j-1);
                end
                
                Wo = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j) * kro_up/mu_o_avg/Bo_avg;
                Wo_prime = Wo /144 * rho_o_avg;
                Ww = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j) * krw_up/mu_w_avg/Bw_avg;
                Ww_prime = Ww / 144 * rho_w_avg;
                Wg = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j) * krg_up/mu_g_avg/Bg_avg;
                Wg_prime = Wg / 144 * rho_g_avg;
                Wgo = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j)* kro_up*Rso_avg/mu_o_avg/Bo_avg;
                Wgo_prime = Wgo /144 * rho_o_avg;
            else
                Wo = 0; Ww = 0; Wg = 0; Wgo = 0; Wo_prime = 0; Ww_prime = 0; Wg_prime = 0; Wgo_prime = 0;
            end
            
            if index(i+1,j) ~= 0
                Po_avg = (Po(i,j) + Po(i+1,j))/2;
                Pw_avg = (Pw(i,j) + Pw(i+1,j))/2;
                Pg_avg = (Pg(i,j) + Pg(i+1,j))/2;
                [rho_o_avg,Bo_avg,mu_o_avg,Rso_avg] = Oil_table(Po_avg);
                [rho_w_avg,Bw_avg,mu_w_avg,~] = Water_table(Pw_avg);
                [rho_g_avg, Bg_avg, mu_g_avg] = Gas_table(Pg_avg);
                
                if Po(i,j)-rho_o(i,j)/144*G(i,j) >= Po(i+1,j)-rho_o(i+1,j)/144*G(i+1,j)
                    kro_up = kro(i,j);
                else; kro_up = kro(i+1,j);
                end
                
                if Pw(i,j)-rho_w(i,j)/144*G(i,j) >= Pw(i+1,j)-rho_w(i+1,j)/144*G(i+1,j)
                    krw_up = krw(i,j);
                else; krw_up = krw(i+1,j);
                end
                
                if Pg(i,j)-rho_g(i,j)/144*G(i,j) >= Pg(i+1,j)-rho_g(i+1,j)/144*G(i+1,j)
                    krg_up = krg(i,j);
                else; krg_up = krg(i+1,j);
                end
                
                So = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j) * kro_up/mu_o_avg/Bo_avg;
                So_prime = So /144 * rho_o_avg;
                Sw = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j) * krw_up/mu_w_avg/Bw_avg;
                Sw_prime = Sw / 144 * rho_w_avg;
                Sg = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j) * krg_up/mu_g_avg/Bg_avg;
                Sg_prime = Sg / 144 * rho_g_avg;
                Sgo = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j)* kro_up*Rso_avg/mu_o_avg/Bo_avg;
                Sgo_prime = Sgo /144 * rho_o_avg;
            else
                So = 0; Sw = 0; Sg = 0; Sgo = 0; So_prime = 0; Sw_prime = 0; Sg_prime = 0; Sgo_prime = 0;
            end
            
            if index(i-1,j) ~= 0
                Po_avg = (Po(i,j) + Po(i-1,j))/2;
                Pw_avg = (Pw(i,j) + Pw(i-1,j))/2;
                Pg_avg = (Pg(i,j) + Pg(i-1,j))/2;
                [rho_o_avg,Bo_avg,mu_o_avg,Rso_avg] = Oil_table(Po_avg);
                [rho_w_avg,Bw_avg,mu_w_avg,~] = Water_table(Pw_avg);
                [rho_g_avg, Bg_avg, mu_g_avg] = Gas_table(Pg_avg);
                
                if Po(i,j)-rho_o(i,j)/144*G(i,j) >= Po(i-1,j)-rho_o(i-1,j)/144*G(i-1,j)
                    kro_up = kro(i,j);
                else; kro_up = kro(i-1,j);
                end
                
                if Pw(i,j)-rho_w(i,j)/144*G(i,j) >= Pw(i-1,j)-rho_w(i-1,j)/144*G(i-1,j)
                    krw_up = krw(i,j);
                else; krw_up = krw(i-1,j);
                end
                
                if Pg(i,j)-rho_g(i,j)/144*G(i,j) >= Pg(i-1,j)-rho_g(i-1,j)/144*G(i-1,j)
                    krg_up = krg(i,j);
                else; krg_up = krg(i-1,j);
                end
                
                No = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j) * kro_up/mu_o_avg/Bo_avg;
                No_prime = No /144 * rho_o_avg;
                Nw = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j) * krw_up/mu_w_avg/Bw_avg;
                Nw_prime = Nw / 144 * rho_w_avg;
                Ng = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j) * krg_up/mu_g_avg/Bg_avg;
                Ng_prime = Ng / 144 * rho_g_avg;
                Ngo = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j)* kro_up*Rso_avg/mu_o_avg/Bo_avg;
                Ngo_prime = Ngo /144 * rho_o_avg;
            else
                No = 0; Nw = 0; Ng = 0; Ngo = 0; No_prime = 0; Nw_prime = 0; Ng_prime = 0; Ngo_prime = 0;
            end
            
            Ro(order(i,j)) = Eo*(Po(i,j+1)-Po(i,j)) + Wo*(Po(i,j-1)-Po(i,j)) + So*(Po(i+1,j)-Po(i,j)) + No*(Po(i-1,j)-Po(i,j)) ...
                -Eo_prime*(G(i,j+1)-G(i,j))-Wo_prime*(G(i,j-1)-G(i,j))-So_prime*(G(i+1,j)-G(i,j))-No_prime*(G(i-1,j)-G(i,j)) ...
                - Vb(i,j)/5.615/dt*phi(i,j)*(1-sat_w(i,j)-sat_g(i,j))/Bo(i,j) ...
                + Vb(i,j)/5.615/dt*phi(i,j)*(1-sat_w_old(i,j)-sat_g_old(i,j))/Bo_old(i,j);
            Rw(order(i,j)) = Ew*(Pw(i,j+1)-Pw(i,j)) + Ww*(Pw(i,j-1)-Pw(i,j)) + Sw*(Pw(i+1,j)-Pw(i,j)) + Nw*(Pw(i-1,j)-Pw(i,j)) ...
                -Ew_prime*(G(i,j+1)-G(i,j))-Ww_prime*(G(i,j-1)-G(i,j))-Sw_prime*(G(i+1,j)-G(i,j))-Nw_prime*(G(i-1,j)-G(i,j)) ...
                - Vb(i,j)/5.615/dt*phi(i,j)*sat_w(i,j)/Bw(i,j) + Vb(i,j)/5.615/dt*phi(i,j)*sat_w_old(i,j)/Bw_old(i,j);
            Rg(order(i,j)) = Eg*(Pg(i,j+1)-Pg(i,j)) + Wg*(Pg(i,j-1)-Pg(i,j)) + Sg*(Pg(i+1,j)-Pg(i,j)) + Ng*(Pg(i-1,j)-Pg(i,j)) ...
                + Ego*(Po(i,j+1)-Po(i,j)) + Wgo*(Po(i,j-1)-Po(i,j)) + Sgo*(Po(i+1,j)-Po(i,j)) + Ngo*(Po(i-1,j)-Po(i,j)) ...
                -Eg_prime*(G(i,j+1)-G(i,j))-Wg_prime*(G(i,j-1)-G(i,j))-Sg_prime*(G(i+1,j)-G(i,j))-Ng_prime*(G(i-1,j)-G(i,j)) ...
                -Ego_prime*(G(i,j+1)-G(i,j))-Wgo_prime*(G(i,j-1)-G(i,j))-Sgo_prime*(G(i+1,j)-G(i,j))-Ngo_prime*(G(i-1,j)-G(i,j)) ...
                - Vb(i,j)/5.615/dt*phi(i,j)*(sat_g(i,j)/Bg(i,j) + Rso(i,j)*(1-sat_w(i,j)-sat_g(i,j))/Bo(i,j)) ...
                + Vb(i,j)/5.615/dt*phi(i,j)*(sat_g_old(i,j)/Bg_old(i,j) + Rso_old(i,j)*(1-sat_w_old(i,j)-sat_g_old(i,j))/Bo_old(i,j));
            
            % Add well specifications:
            if index(i,j) == 2
                re = 0.198*dx(i,j);
                Gw = 2*pi*beta_c*sqrt(k_x(i,j)*k_y(i,j))/(log(re/rw)+skin);
                Psf = Po(i,j) + q_sp/(Gw*kro(i,j)/mu_o(i,j)/Bo(i,j));
                
                qo_sc = q_sp;
                qw_sc = -Gw*krw(i,j)/mu_w(i,j)/Bw(i,j)*(Po(i,j) - Psf);
                qg_sc = -Gw*(krg(i,j)/mu_g(i,j)/Bg(i,j)+ Rso(i,j)*kro(i,j)/mu_o(i,j)/Bo(i,j))*(Po(i,j) - Psf);
            else
                qo_sc = 0;    qw_sc = 0;    qg_sc = 0;
            end
            
            Ro(order(i,j)) = Ro(order(i,j)) + qo_sc;
            Rw(order(i,j)) = Rw(order(i,j)) + qw_sc;
            Rg(order(i,j)) = Rg(order(i,j)) + qg_sc;
            
            %% Derivative respect to Po - E,W,S,N
            if index(i,j+1) ~= 0         % East direction
                Po_avg = (Po(i,j) + Po(i,j+1) + eps_p)/2;
                Pw_avg = (Pw(i,j) + Pw(i,j+1) + eps_p)/2;
                Pg_avg = (Pg(i,j) + Pg(i,j+1) + eps_p)/2;
                
                [rho_o_avg,Bo_avg,mu_o_avg,Rso_avg] = Oil_table(Po_avg);
                [rho_w_avg,Bw_avg,mu_w_avg,~] = Water_table(Pw_avg);
                [rho_g_avg, Bg_avg, mu_g_avg] = Gas_table(Pg_avg);
                
                [rho_o_E, ~, ~, ~] = Oil_table(Po(i,j+1) + eps_p);
                [rho_w_E, ~, ~, ~] = Water_table(Pw(i,j+1) + eps_p);
                [rho_g_E, ~, ~] = Gas_table(Pg(i,j+1) + eps_p);
                
                if Po(i,j)-rho_o(i,j)/144*G(i,j) >= Po(i,j+1)+eps_p-rho_o_E/144*G(i,j+1)
                    kro_up = kro(i,j);
                else; kro_up = kro(i,j+1);
                end
                    
                if Pw(i,j)-rho_w(i,j)/144*G(i,j) >= Pw(i,j+1)+eps_p-rho_w_E/144*G(i,j+1)
                    krw_up = krw(i,j);
                else; krw_up = krw(i,j+1);
                end
                    
                if Pg(i,j)-rho_g(i,j)/144*G(i,j) >= Pg(i,j+1)+eps_p-rho_g_E/144*G(i,j+1)
                    krg_up = krg(i,j);
                else; krg_up = krg(i,j+1);
                end
                    
                Eo_eps = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j) * kro_up/mu_o_avg/Bo_avg;
                Eo_prime_eps = Eo_eps / 144 * rho_o_avg;
                Ew_eps = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j) * krw_up/mu_w_avg/Bw_avg;
                Ew_prime_eps = Ew_eps / 144 * rho_w_avg;
                Eg_eps = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j) * krg_up/mu_g_avg/Bg_avg;
                Eg_prime_eps = Eg_eps / 144 * rho_g_avg;
                Ego_eps = Eo_eps * Rso_avg;
                Ego_prime_eps = Ego_eps / 144 * rho_o_avg;
                
                Ro_E = Ro(order(i,j))-Eo*(Po(i,j+1)-Po(i,j)) + Eo_eps*(Po(i,j+1)+eps_p-Po(i,j)) ...
                    + (Eo_prime - Eo_prime_eps) * (G(i,j+1)-G(i,j));
                Rw_E = Rw(order(i,j))-Ew*(Pw(i,j+1)-Pw(i,j)) + Ew_eps*(Pw(i,j+1)+eps_p-Pw(i,j)) ...
                    + (Ew_prime - Ew_prime_eps) * (G(i,j+1)-G(i,j));
                Rg_E = Rg(order(i,j))-Eg*(Pg(i,j+1)-Pg(i,j)) + Eg_eps*(Pg(i,j+1)+eps_p-Pg(i,j)) ...
                    - Ego*(Po(i,j+1)-Po(i,j)) + Ego_eps*(Po(i,j+1)+eps_p-Po(i,j)) ...
                    + (Eg_prime - Eg_prime_eps + Ego_prime - Ego_prime_eps) * (G(i,j+1)-G(i,j));
                dRo_dp_E = (Ro_E - Ro(order(i,j))) / eps_p;
                dRw_dp_E = (Rw_E - Rw(order(i,j))) / eps_p;
                dRg_dp_E = (Rg_E - Rg(order(i,j))) / eps_p;
                
                Jacob_1st_col(3*order(i,j)-2, order(i,j+1)) = dRo_dp_E;
                Jacob_1st_col(3*order(i,j)-1, order(i,j+1)) = dRw_dp_E;
                Jacob_1st_col(3*order(i,j)-0, order(i,j+1)) = dRg_dp_E;
            end
            
            if index(i,j-1) ~= 0         % West direction
                Po_avg = (Po(i,j) + Po(i,j-1) + eps_p)/2;
                Pw_avg = (Pw(i,j) + Pw(i,j-1) + eps_p)/2;
                Pg_avg = (Pg(i,j) + Pg(i,j-1) + eps_p)/2;
                
                [rho_o_avg,Bo_avg,mu_o_avg,Rso_avg] = Oil_table(Po_avg);
                [rho_w_avg,Bw_avg,mu_w_avg,~] = Water_table(Pw_avg);
                [rho_g_avg, Bg_avg, mu_g_avg] = Gas_table(Pg_avg);
                
                [rho_o_W, ~, ~, ~] = Oil_table(Po(i,j-1) + eps_p);
                [rho_w_W, ~, ~, ~] = Water_table(Pw(i,j-1) + eps_p);
                [rho_g_W, ~, ~] = Gas_table(Pg(i,j-1) + eps_p);
                
                if Po(i,j)-rho_o(i,j)/144*G(i,j) >= Po(i,j-1)+eps_p-rho_o_W/144*G(i,j-1)
                    kro_up = kro(i,j);
                else; kro_up = kro(i,j-1);
                end
                    
                if Pw(i,j)-rho_w(i,j)/144*G(i,j) >= Pw(i,j-1)+eps_p-rho_w_W/144*G(i,j-1)
                    krw_up = krw(i,j);
                else; krw_up = krw(i,j-1);
                end
                    
                if Pg(i,j)-rho_g(i,j)/144*G(i,j) >= Pg(i,j-1)+eps_p-rho_g_W/144*G(i,j-1)
                    krg_up = krg(i,j);
                else; krg_up = krg(i,j-1);
                end
                    
                Wo_eps = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j) * kro_up/mu_o_avg/Bo_avg;
                Wo_prime_eps = Wo_eps / 144 * rho_o_avg;
                Ww_eps = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j) * krw_up/mu_w_avg/Bw_avg;
                Ww_prime_eps = Ww_eps / 144 * rho_w_avg;
                Wg_eps = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j) * krg_up/mu_g_avg/Bg_avg;
                Wg_prime_eps = Wg_eps / 144 * rho_g_avg;
                Wgo_eps = Wo_eps * Rso_avg;
                Wgo_prime_eps = Wgo_eps / 144 * rho_o_avg;
                
                Ro_W = Ro(order(i,j))-Wo*(Po(i,j-1)-Po(i,j)) + Wo_eps*(Po(i,j-1)+eps_p-Po(i,j)) ...
                    + (Wo_prime - Wo_prime_eps) * (G(i,j-1)-G(i,j));
                Rw_W = Rw(order(i,j))-Ww*(Pw(i,j-1)-Pw(i,j)) + Ww_eps*(Pw(i,j-1)+eps_p-Pw(i,j)) ...
                    + (Ww_prime - Ww_prime_eps) * (G(i,j-1)-G(i,j));
                Rg_W = Rg(order(i,j))-Wg*(Pg(i,j-1)-Pg(i,j)) + Wg_eps*(Pg(i,j-1)+eps_p-Pg(i,j)) ...
                    - Wgo*(Po(i,j-1)-Po(i,j)) + Wgo_eps*(Po(i,j-1)+eps_p-Po(i,j)) ...
                    + (Wg_prime - Wg_prime_eps + Wgo_prime - Wgo_prime_eps) * (G(i,j-1)-G(i,j));
                dRo_dp_W = (Ro_W - Ro(order(i,j))) / eps_p;
                dRw_dp_W = (Rw_W - Rw(order(i,j))) / eps_p;
                dRg_dp_W = (Rg_W - Rg(order(i,j))) / eps_p;
                
                Jacob_1st_col(3*order(i,j)-2, order(i,j-1)) = dRo_dp_W;
                Jacob_1st_col(3*order(i,j)-1, order(i,j-1)) = dRw_dp_W;
                Jacob_1st_col(3*order(i,j)-0, order(i,j-1)) = dRg_dp_W;
            end
            
            if index(i+1,j) ~= 0         % South direction
                Po_avg = (Po(i,j) + Po(i+1,j) + eps_p)/2;
                Pw_avg = (Pw(i,j) + Pw(i+1,j) + eps_p)/2;
                Pg_avg = (Pg(i,j) + Pg(i+1,j) + eps_p)/2;
                
                [rho_o_avg,Bo_avg,mu_o_avg,Rso_avg] = Oil_table(Po_avg);
                [rho_w_avg,Bw_avg,mu_w_avg,~] = Water_table(Pw_avg);
                [rho_g_avg, Bg_avg, mu_g_avg] = Gas_table(Pg_avg);
                
                [rho_o_S, ~, ~, ~] = Oil_table(Po(i+1,j) + eps_p);
                [rho_w_S, ~, ~, ~] = Water_table(Pw(i+1,j) + eps_p);
                [rho_g_S, ~, ~] = Gas_table(Pg(i+1,j) + eps_p);
                
                if Po(i,j)-rho_o(i,j)/144*G(i,j) >= Po(i+1,j)+eps_p-rho_o_S/144*G(i+1,j)
                    kro_up = kro(i,j);
                else; kro_up = kro(i+1,j);
                end
                    
                if Pw(i,j)-rho_w(i,j)/144*G(i,j) >= Pw(i+1,j)+eps_p-rho_w_S/144*G(i+1,j)
                    krw_up = krw(i,j);
                else; krw_up = krw(i+1,j);
                end
                    
                if Pg(i,j)-rho_g(i,j)/144*G(i,j) >= Pg(i+1,j)+eps_p-rho_g_S/144*G(i+1,j)
                    krg_up = krg(i,j);
                else; krg_up = krg(i+1,j);
                end
                    
                So_eps = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j) * kro_up/mu_o_avg/Bo_avg;
                So_prime_eps = So_eps / 144 * rho_o_avg;
                Sw_eps = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j) * krw_up/mu_w_avg/Bw_avg;
                Sw_prime_eps = Sw_eps / 144 * rho_w_avg;
                Sg_eps = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j) * krg_up/mu_g_avg/Bg_avg;
                Sg_prime_eps = Sg_eps / 144 * rho_g_avg;
                Sgo_eps = So_eps * Rso_avg;
                Sgo_prime_eps = Sgo_eps / 144 * rho_o_avg;
                
                Ro_S = Ro(order(i,j))-So*(Po(i+1,j)-Po(i,j)) + So_eps*(Po(i+1,j)+eps_p-Po(i,j)) ...
                    + (So_prime - So_prime_eps) * (G(i+1,j)-G(i,j));
                Rw_S = Rw(order(i,j))-Sw*(Pw(i+1,j)-Pw(i,j)) + Sw_eps*(Pw(i+1,j)+eps_p-Pw(i,j)) ...
                    + (Sw_prime - Sw_prime_eps) * (G(i+1,j)-G(i,j));
                Rg_S = Rg(order(i,j))-Sg*(Pg(i+1,j)-Pg(i,j)) + Sg_eps*(Pg(i+1,j)+eps_p-Pg(i,j)) ...
                    - Sgo*(Po(i+1,j)-Po(i,j)) + Sgo_eps*(Po(i+1,j)+eps_p-Po(i,j)) ...
                    + (Sg_prime - Sg_prime_eps + Sgo_prime - Sgo_prime_eps) * (G(i+1,j)-G(i,j));
                dRo_dp_S = (Ro_S - Ro(order(i,j))) / eps_p;
                dRw_dp_S = (Rw_S - Rw(order(i,j))) / eps_p;
                dRg_dp_S = (Rg_S - Rg(order(i,j))) / eps_p;
                
                Jacob_1st_col(3*order(i,j)-2, order(i+1,j)) = dRo_dp_S;
                Jacob_1st_col(3*order(i,j)-1, order(i+1,j)) = dRw_dp_S;
                Jacob_1st_col(3*order(i,j)-0, order(i+1,j)) = dRg_dp_S;
            end
            
            if index(i-1,j) ~= 0         % North direction
                Po_avg = (Po(i,j) + Po(i-1,j) + eps_p)/2;
                Pw_avg = (Pw(i,j) + Pw(i-1,j) + eps_p)/2;
                Pg_avg = (Pg(i,j) + Pg(i-1,j) + eps_p)/2;
                
                [rho_o_avg,Bo_avg,mu_o_avg,Rso_avg] = Oil_table(Po_avg);
                [rho_w_avg,Bw_avg,mu_w_avg,~] = Water_table(Pw_avg);
                [rho_g_avg, Bg_avg, mu_g_avg] = Gas_table(Pg_avg);
                
                [rho_o_N, ~, ~, ~] = Oil_table(Po(i-1,j) + eps_p);
                [rho_w_N, ~, ~, ~] = Water_table(Pw(i-1,j) + eps_p);
                [rho_g_N, ~, ~] = Gas_table(Pg(i-1,j) + eps_p);
                
                if Po(i,j)-rho_o(i,j)/144*G(i,j) >= Po(i-1,j)+eps_p-rho_o_N/144*G(i-1,j)
                    kro_up = kro(i,j);
                else; kro_up = kro(i-1,j);
                end
                    
                if Pw(i,j)-rho_w(i,j)/144*G(i,j) >= Pw(i-1,j)+eps_p-rho_w_N/144*G(i-1,j)
                    krw_up = krw(i,j);
                else; krw_up = krw(i-1,j);
                end
                    
                if Pg(i,j)-rho_g(i,j)/144*G(i,j) >= Pg(i-1,j)+eps_p-rho_g_N/144*G(i-1,j)
                    krg_up = krg(i,j);
                else; krg_up = krg(i-1,j);
                end
                    
                No_eps = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j) * kro_up/mu_o_avg/Bo_avg;
                No_prime_eps = No_eps / 144 * rho_o_avg;
                Nw_eps = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j) * krw_up/mu_w_avg/Bw_avg;
                Nw_prime_eps = Nw_eps / 144 * rho_w_avg;
                Ng_eps = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j) * krg_up/mu_g_avg/Bg_avg;
                Ng_prime_eps = Ng_eps / 144 * rho_g_avg;
                Ngo_eps = No_eps * Rso_avg;
                Ngo_prime_eps = Ngo_eps / 144 * rho_o_avg;
                
                Ro_N = Ro(order(i,j))-No*(Po(i-1,j)-Po(i,j)) + No_eps*(Po(i-1,j)+eps_p-Po(i,j)) ...
                    + (No_prime - No_prime_eps) * (G(i-1,j)-G(i,j));
                Rw_N = Rw(order(i,j))-Nw*(Pw(i-1,j)-Pw(i,j)) + Nw_eps*(Pw(i-1,j)+eps_p-Pw(i,j)) ...
                    + (Nw_prime - Nw_prime_eps) * (G(i-1,j)-G(i,j));
                Rg_N = Rg(order(i,j))-Ng*(Pg(i-1,j)-Pg(i,j)) + Ng_eps*(Pg(i-1,j)+eps_p-Pg(i,j)) ...
                    - Ngo*(Po(i-1,j)-Po(i,j)) + Ngo_eps*(Po(i-1,j)+eps_p-Po(i,j)) ...
                    + (Ng_prime - Ng_prime_eps + Ngo_prime - Ngo_prime_eps) * (G(i-1,j)-G(i,j));
                dRo_dp_N = (Ro_N - Ro(order(i,j))) / eps_p;
                dRw_dp_N = (Rw_N - Rw(order(i,j))) / eps_p;
                dRg_dp_N = (Rg_N - Rg(order(i,j))) / eps_p;
                
                Jacob_1st_col(3*order(i,j)-2, order(i-1,j)) = dRo_dp_N;
                Jacob_1st_col(3*order(i,j)-1, order(i-1,j)) = dRw_dp_N;
                Jacob_1st_col(3*order(i,j)-0, order(i-1,j)) = dRg_dp_N;
            end
            
            %% Central derivative - Po:
            [rho_o_C, Bo_C, mu_o_C, Rso_C] = Oil_table(Po(i,j) + eps_p);
            [rho_w_C, Bw_C, mu_w_C, ~] = Water_table(Pw(i,j) + eps_p);
            [rho_g_C, Bg_C, mu_g_C] = Gas_table(Pg(i,j) + eps_p);
            
            if index(i,j+1) ~= 0
                Po_avg = (Po(i,j) + eps_p + Po(i,j+1))/2;
                Pw_avg = (Pw(i,j) + eps_p + Pw(i,j+1))/2;
                Pg_avg = (Pg(i,j) + eps_p + Pg(i,j+1))/2;
                
                [rho_o_avg,Bo_avg,mu_o_avg,Rso_avg] = Oil_table(Po_avg);
                [rho_w_avg,Bw_avg,mu_w_avg,~] = Water_table(Pw_avg);
                [rho_g_avg, Bg_avg, mu_g_avg] = Gas_table(Pg_avg);
                
                if Po(i,j)+eps_p-rho_o_C/144*G(i,j) >= Po(i,j+1)-rho_o(i,j+1)/144*G(i,j+1)
                    kro_up = kro(i,j);
                else; kro_up = kro(i,j+1); 
                end
                
                if Pw(i,j)+eps_p-rho_w_C/144*G(i,j) >= Pw(i,j+1)-rho_w(i,j+1)/144*G(i,j+1)
                    krw_up = krw(i,j);
                else; krw_up = krw(i,j+1);
                end
                
                if Pg(i,j)+eps_p-rho_g_C/144*G(i,j) >= Pg(i,j+1)-rho_g(i,j+1)/144*G(i,j+1)
                    krg_up = krg(i,j);
                else; krg_up = krg(i,j+1);
                end
                
                Eo_eps = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j) * kro_up/mu_o_avg/Bo_avg;
                Eo_prime_eps = Eo_eps / 144 * rho_o_avg;
                Ew_eps = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j) * krw_up/mu_w_avg/Bw_avg;
                Ew_prime_eps = Ew_eps / 144 * rho_w_avg;
                Eg_eps = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j) * krg_up/mu_g_avg/Bg_avg;
                Eg_prime_eps = Eg_eps / 144 * rho_g_avg;
                Ego_eps = Eo_eps * Rso_avg;
                Ego_prime_eps = Ego_eps / 144 * rho_o_avg;
            else
                Eo_eps = 0; Eo_prime_eps = 0; Ew_eps = 0; Ew_prime_eps = 0;
                Eg_eps = 0; Eg_prime_eps = 0; Ego_eps = 0; Ego_prime_eps = 0;
            end
            
            if index(i,j-1) ~= 0
                Po_avg = (Po(i,j) + eps_p + Po(i,j-1))/2;
                Pw_avg = (Pw(i,j) + eps_p + Pw(i,j-1))/2;
                Pg_avg = (Pg(i,j) + eps_p + Pg(i,j-1))/2;
                
                [rho_o_avg,Bo_avg,mu_o_avg,Rso_avg] = Oil_table(Po_avg);
                [rho_w_avg,Bw_avg,mu_w_avg,~] = Water_table(Pw_avg);
                [rho_g_avg, Bg_avg, mu_g_avg] = Gas_table(Pg_avg);
                
                if Po(i,j)+eps_p-rho_o_C/144*G(i,j) >= Po(i,j-1)-rho_o(i,j-1)/144*G(i,j-1)
                    kro_up = kro(i,j);
                else; kro_up = kro(i,j-1); 
                end
                
                if Pw(i,j)+eps_p-rho_w_C/144*G(i,j) >= Pw(i,j-1)-rho_w(i,j-1)/144*G(i,j-1)
                    krw_up = krw(i,j);
                else; krw_up = krw(i,j-1); 
                end
                
                if Pg(i,j)+eps_p-rho_g_C/144*G(i,j) >= Pg(i,j-1)-rho_g(i,j-1)/144*G(i,j-1)
                    krg_up = krg(i,j);
                else; krg_up = krg(i,j-1); 
                end
                
                Wo_eps = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j) * kro_up/mu_o_avg/Bo_avg;
                Wo_prime_eps = Wo_eps / 144 * rho_o_avg;
                Ww_eps = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j) * krw_up/mu_w_avg/Bw_avg;
                Ww_prime_eps = Ww_eps / 144 * rho_w_avg;
                Wg_eps = beta_c*Ax(i,j)*k_x(i,j)/dx(i,j) * krg_up/mu_g_avg/Bg_avg;
                Wg_prime_eps = Wg_eps / 144 * rho_g_avg;
                Wgo_eps = Wo_eps * Rso_avg;
                Wgo_prime_eps = Wgo_eps / 144 * rho_o_avg;
            else
                Wo_eps = 0; Wo_prime_eps = 0; Ww_eps = 0; Ww_prime_eps = 0;
                Wg_eps = 0; Wg_prime_eps = 0; Wgo_eps = 0; Wgo_prime_eps = 0;
            end
            
            if index(i+1,j) ~= 0
                Po_avg = (Po(i,j) + eps_p + Po(i+1,j))/2;
                Pw_avg = (Pw(i,j) + eps_p + Pw(i+1,j))/2;
                Pg_avg = (Pg(i,j) + eps_p + Pg(i+1,j))/2;
                
                [rho_o_avg,Bo_avg,mu_o_avg,Rso_avg] = Oil_table(Po_avg);
                [rho_w_avg,Bw_avg,mu_w_avg,~] = Water_table(Pw_avg);
                [rho_g_avg, Bg_avg, mu_g_avg] = Gas_table(Pg_avg);
                
                if Po(i,j)+eps_p-rho_o_C/144*G(i,j) >= Po(i+1,j)-rho_o(i+1,j)/144*G(i+1,j)
                    kro_up = kro(i,j);
                else; kro_up = kro(i+1,j); 
                end
                
                if Pw(i,j)+eps_p-rho_w_C/144*G(i,j) >= Pw(i+1,j)-rho_w(i+1,j)/144*G(i+1,j)
                    krw_up = krw(i,j);
                else; krw_up = krw(i+1,j); 
                end
                
                if Pg(i,j)+eps_p-rho_g_C/144*G(i,j) >= Pg(i+1,j)-rho_g(i+1,j)/144*G(i+1,j)
                    krg_up = krg(i,j);
                else; krg_up = krg(i+1,j); 
                end
                
                So_eps = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j) * kro_up/mu_o_avg/Bo_avg;
                So_prime_eps = So_eps / 144 * rho_o_avg;
                Sw_eps = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j) * krw_up/mu_w_avg/Bw_avg;
                Sw_prime_eps = Sw_eps / 144 * rho_w_avg;
                Sg_eps = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j) * krg_up/mu_g_avg/Bg_avg;
                Sg_prime_eps = Sg_eps / 144 * rho_g_avg;
                Sgo_eps = So_eps * Rso_avg;
                Sgo_prime_eps = Sgo_eps / 144 * rho_o_avg;
            else
                So_eps = 0; So_prime_eps = 0; Sw_eps = 0; Sw_prime_eps = 0;
                Sg_eps = 0; Sg_prime_eps = 0; Sgo_eps = 0; Sgo_prime_eps = 0;
            end
            
            if index(i-1,j) ~= 0
                Po_avg = (Po(i,j) + eps_p + Po(i-1,j))/2;
                Pw_avg = (Pw(i,j) + eps_p + Pw(i-1,j))/2;
                Pg_avg = (Pg(i,j) + eps_p + Pg(i-1,j))/2;
                
                [rho_o_avg,Bo_avg,mu_o_avg,Rso_avg] = Oil_table(Po_avg);
                [rho_w_avg,Bw_avg,mu_w_avg,~] = Water_table(Pw_avg);
                [rho_g_avg, Bg_avg, mu_g_avg] = Gas_table(Pg_avg);
                
                if Po(i,j)+eps_p-rho_o_C/144*G(i,j) >= Po(i-1,j)-rho_o(i-1,j)/144*G(i-1,j)
                    kro_up = kro(i,j);
                else; kro_up = kro(i-1,j);
                end
                
                if Pw(i,j)+eps_p-rho_w_C/144*G(i,j) >= Pw(i-1,j)-rho_w(i-1,j)/144*G(i-1,j)
                    krw_up = krw(i,j);
                else; krw_up = krw(i-1,j); 
                end
                
                if Pg(i,j)+eps_p-rho_g_C/144*G(i,j) >= Pg(i-1,j)-rho_g(i-1,j)/144*G(i-1,j)
                    krg_up = krg(i,j);
                else; krg_up = krg(i-1,j); 
                end
                
                No_eps = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j) * kro_up/mu_o_avg/Bo_avg;
                No_prime_eps = No_eps / 144 * rho_o_avg;
                Nw_eps = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j) * krw_up/mu_w_avg/Bw_avg;
                Nw_prime_eps = Nw_eps / 144 * rho_w_avg;
                Ng_eps = beta_c*Ay(i,j)*k_y(i,j)/dy(i,j) * krg_up/mu_g_avg/Bg_avg;
                Ng_prime_eps = Ng_eps / 144 * rho_g_avg;
                Ngo_eps = No_eps * Rso_avg;
                Ngo_prime_eps = Ngo_eps / 144 * rho_o_avg;
            else
                No_eps = 0; No_prime_eps = 0; Nw_eps = 0; Nw_prime_eps = 0;
                Ng_eps = 0; Ng_prime_eps = 0; Ngo_eps = 0; Ngo_prime_eps = 0;
            end
            
            Ro_C = Eo_eps*(Po(i,j+1)-Po(i,j)-eps_p) + Wo_eps*(Po(i,j-1)-Po(i,j)-eps_p) ...
                + So_eps*(Po(i+1,j)-Po(i,j)-eps_p) + No_eps*(Po(i-1,j)-Po(i,j)-eps_p) ...
                - Eo_prime_eps*(G(i,j+1)-G(i,j)) - Wo_prime_eps*(G(i,j-1)-G(i,j)) ...
                - So_prime_eps*(G(i+1,j)-G(i,j)) - No_prime_eps*(G(i-1,j)-G(i,j)) ...
                - Vb(i,j)/5.615/dt*phi(i,j)*(1-sat_w(i,j)-sat_g(i,j))/Bo_C ...
                + Vb(i,j)/5.615/dt*phi(i,j)*(1-sat_w_old(i,j)-sat_g_old(i,j))/Bo_old(i,j);
            
            Rw_C = Ew_eps*(Pw(i,j+1)-Pw(i,j)-eps_p) + Ww_eps*(Pw(i,j-1)-Pw(i,j)-eps_p) ...
                + Sw_eps*(Pw(i+1,j)-Pw(i,j)-eps_p) + Nw_eps*(Pw(i-1,j)-Pw(i,j)-eps_p) ...
                - Ew_prime_eps*(G(i,j+1)-G(i,j)) - Ww_prime_eps*(G(i,j-1)-G(i,j)) ...
                - Sw_prime_eps*(G(i+1,j)-G(i,j)) - Nw_prime_eps*(G(i-1,j)-G(i,j)) ...
                - Vb(i,j)/5.615/dt*phi(i,j)*sat_w(i,j)/Bw_C + Vb(i,j)/5.615/dt*phi(i,j)*sat_w_old(i,j)/Bw_old(i,j);
            Rg_C = Eg_eps*(Pg(i,j+1)-Pg(i,j)-eps_p) + Wg_eps*(Pg(i,j-1)-Pg(i,j)-eps_p) ...
                + Sg_eps*(Pg(i+1,j)-Pg(i,j)-eps_p) + Ng_eps*(Pg(i-1,j)-Pg(i,j)-eps_p) ...
                + Ego_eps*(Po(i,j+1)-Po(i,j)-eps_p) + Wgo_eps*(Po(i,j-1)-Po(i,j)-eps_p) ...
                + Sgo_eps*(Po(i+1,j)-Po(i,j)-eps_p) + Ngo_eps*(Po(i-1,j)-Po(i,j)-eps_p) ...
                - Eg_prime_eps*(G(i,j+1)-G(i,j)) - Wg_prime_eps*(G(i,j-1)-G(i,j)) ...
                - Sg_prime_eps*(G(i+1,j)-G(i,j)) - Ng_prime_eps*(G(i-1,j)-G(i,j)) ...
                - Ego_prime_eps*(G(i,j+1)-G(i,j)) - Wgo_prime_eps*(G(i,j-1)-G(i,j)) ...
                - Sgo_prime_eps*(G(i+1,j)-G(i,j)) - Ngo_prime_eps*(G(i-1,j)-G(i,j)) ...
                - Vb(i,j)/5.615/dt*phi(i,j)*(sat_g(i,j)/Bg_C + Rso_C*(1-sat_w(i,j)-sat_g(i,j))/Bo_C) ...
                + Vb(i,j)/5.615/dt*phi(i,j)*(sat_g_old(i,j)/Bg_old(i,j) + Rso_old(i,j)*(1-sat_w_old(i,j)-sat_g_old(i,j))/Bo_old(i,j));
            
            % Add well specifications:
            if index(i,j) == 2
                re = 0.198*dx(i,j);
                Gw = 2*pi*beta_c*sqrt(k_x(i,j)*k_y(i,j))/(log(re/rw)+skin);
                Psf_eps = Po(i,j)+eps_p + q_sp/(Gw*kro(i,j)/mu_o_C/Bo_C);
                
                qo_sc_eps = q_sp;
                qw_sc_eps = -Gw*krw(i,j)/mu_w_C/Bw_C*(Po(i,j) + eps_p - Psf_eps);
                qg_sc_eps = -Gw*(krg(i,j)/mu_g_C/Bg_C + Rso_C*kro(i,j)/mu_o_C/Bo_C)*(Po(i,j) + eps_p - Psf_eps);
            else
                qo_sc_eps = 0;    qw_sc_eps = 0;    qg_sc_eps = 0;
            end
            
            Ro_C = Ro_C + qo_sc_eps;
            Rw_C = Rw_C + qw_sc_eps;
            Rg_C = Rg_C + qg_sc_eps;
            
            dRo_dp_C = (Ro_C - Ro(order(i,j))) / eps_p;
            dRw_dp_C = (Rw_C - Rw(order(i,j))) / eps_p;
            dRg_dp_C = (Rg_C - Rg(order(i,j))) / eps_p;
            
            Jacob_1st_col(3*order(i,j)-2, order(i,j)) = dRo_dp_C;
            Jacob_1st_col(3*order(i,j)-1, order(i,j)) = dRw_dp_C;
            Jacob_1st_col(3*order(i,j)-0, order(i,j)) = dRg_dp_C;
        end
    end
end
end