function [RHS] = Residual(Po, sat_w, sat_g, Po_old, sat_w_old, sat_g_old)
% Calculate Residual of every reservoir blocks
index  = [0 0 0 0 0; 0 1 1 1 0; 0 1 2 1 0; 0 1 1 1 0; 0 0 0 0 0];
index2 = [0 0 0 0 0; 0 1 1 1 0; 0 1 1 1 0; 0 1 1 1 0; 0 0 0 0 0];
order  = [0 0 0 0 0; 0 1 2 3 0; 0 4 5 6 0; 0 7 8 9 0; 0 0 0 0 0];
num_act = max(max(order));      %total number of active blocks
RHS = zeros(3*num_act, 1);

dx = 500 * index2; dy = 500 * index2; h = 50 * index2; G = 4000 * index2;
Ax = h .* dy; Ay = h .* dx; Vb = h .* dx .* dy;
k_x = 100 * index2; k_y = 100 * index2; [u, v] = size(index);
phi = 0.2 * index2; rw = 0.25; skin = 0; beta_c = 1.127e-3;
q_sp = -20; dt = 5;
format long
    
%% Initial reservoir pressure and saturation:
% Po_old = 4800 * index2;       % oil pressure at t = 0
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
                - Vb(i,j)/5.615/dt*phi(i,j)*sat_w(i,j)/Bw(i,j) ...
                + Vb(i,j)/5.615/dt*phi(i,j)*sat_w_old(i,j)/Bw_old(i,j);
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
            
            RHS(3*order(i,j)-2) = Ro(order(i,j));
            RHS(3*order(i,j)-1) = Rw(order(i,j));
            RHS(3*order(i,j)-0) = Rg(order(i,j));
        end
    end
end
end