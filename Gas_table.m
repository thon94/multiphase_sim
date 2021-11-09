function [rho_g, Bg, mu_g] = Gas_table(P)
% Look up and extrapolate the Gas Properties Table
    format long;
    
    Table = [1500	5.8267	0.0018	0.015
             2000	8.0573	0.00133	0.0167
             2500	10.228	0.00105	0.0185
             3000	12.208	0.00088	0.0204
             3500	13.942	0.00077	0.0222
             4000	15.431	0.00069	0.0241
             4500	16.705	0.00064	0.026
             5000	17.799	0.0006	0.0278
             5500	18.748	0.00057	0.0296
             6000	19.577	0.00055	0.0313];
%     Bg = 0.005328*exp(-0.001081*P) + 0.000831*exp(-7.244e-5*P);
    for i = 1 : numel(Table(:,1)) - 1
        while P >= Table(i, 1) && P < Table(i+1, 1)
            ratio = (P - Table(i,1))/(Table(i+1,1) - Table(i,1));
            rho_g = ratio*(Table(i+1,2)-Table(i,2)) + Table(i,2);
            Bg = ratio*(Table(i+1,3)-Table(i,3)) + Table(i,3);
            mu_g = ratio*(Table(i+1,4)-Table(i,4)) + Table(i,4);
            break
        end
    end
%     output = [rho_g, Bg, mu_g];
end