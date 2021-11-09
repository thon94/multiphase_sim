function [rho_o, Bo, mu_o, Rso] = Oil_table(P)
% Look up and extrapolate the Oil Properties Table
    format long;
    
    Table = [1500	49.0113	 1.20413	 1.7356	  292.75
             2000	48.5879	 1.2321	     1.5562	  368
             2500	48.1774	 1.26054	 1.4015	  443.75
             3000	47.6939	 1.29208	 1.2516	  522.71
             3500	47.1788	 1.32933	 1.1024	  619
             4000	46.5899	 1.37193	 0.9647	  724.92
             4500	45.5756	 1.42596	 0.918	  818.6
             5000	45.1925	 1.46387	 0.92	  923.12
             5500	45.4413	 1.44983	 0.9243	  965.28
             6000	45.7426	 1.43831	 0.9372	  966.32];

    for i = 1 : numel(Table(:,1)) - 1
        while P >= Table(i, 1) && P < Table(i+1, 1)
            ratio = (P - Table(i,1))/(Table(i+1,1) - Table(i,1));
            rho_o = ratio*(Table(i+1,2)-Table(i,2)) + Table(i,2);
            Bo = ratio*(Table(i+1,3)-Table(i,3)) + Table(i,3);
            mu_o = ratio*(Table(i+1,4)-Table(i,4)) + Table(i,4);
            Rso = ratio*(Table(i+1,5)-Table(i,5)) + Table(i,5);
            break
        end
    end
%     output = [rho_o, Bo, mu_o, Rso];
end

