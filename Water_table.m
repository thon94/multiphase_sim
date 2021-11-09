function [rho_w, Bw, mu_w, Rsw] = Water_table(P)
% Look up and extrapolate the Water Properties Table
    format long;
    
    Table = [1500	62.228	1.0253	0.52	0
             2000	62.413	1.0222	0.52	0
             2500	62.597	1.0192	0.52	0
             3000	62.782	1.0162	0.52	0
             3500	62.968	1.0132	0.52	0
             4000	63.153	1.0102	0.52	0
             4500	63.337	1.0073	0.52	0
             5000	63.523	1.0051	0.52	0
             5500	63.708	1.0017	0.52	0
             6000	63.893	0.9986	0.52	0];
    for i = 1 : numel(Table(:,1)) - 1
        while P >= Table(i, 1) && P < Table(i+1, 1)
            ratio = (P - Table(i,1))/(Table(i+1,1) - Table(i,1));
            rho_w = ratio*(Table(i+1,2)-Table(i,2)) + Table(i,2);
            Bw = ratio*(Table(i+1,3)-Table(i,3)) + Table(i,3);
            mu_w = ratio*(Table(i+1,4)-Table(i,4)) + Table(i,4);
            Rsw = ratio*(Table(i+1,5)-Table(i,5)) + Table(i,5);
            break
        end
    end
%     output = [rho_w, Bw, mu_w, Rsw];
end