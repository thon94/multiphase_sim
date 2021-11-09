function [krw, krow, Pcow] = Sw_table(Sw)
% Look up and extrapolate the Water Saturation Table
    format long;
    
    Table = [0.18	0	        1	       9
             0.21	0	        0.92692	   7.26
             0.24	2.00E-05	0.85441	   5.04
             0.27	0.00014	    0.79288	   3.78
             0.3	0.00045	    0.71312	   3
             0.33	0.00111	    0.64526	   2.634
             0.36	0.00232	    0.5798	   2.268
             0.39	0.0043	    0.51709	   1.902
             0.42	0.00733	    0.45744	   1.666
             0.45	0.01175	    0.4011	   1.495
             0.48	0.01791	    0.34831	   1.324
             0.51	0.02623	    0.29924	   1.168
             0.54	0.03714	    0.25403    1.042
             0.57	0.05116	    0.21278	   0.916
             0.6	0.06882	    0.17552	   0.79
             0.63	0.09069	    0.14228	   0.682
             0.66	0.11741	    0.11301	   0.574
             0.69	0.14963	    0.08763	   0.466
             0.72	0.18807	    0.06603	   0.364
             0.75	0.23347	    0.04803	   0.265
             0.78	0.28664	    0.03344	   0.166
             0.81	0.34842	    0.02199	   0.09
             0.84	0.41968	    0.0134	   0.06
             0.87	0.50135	    0.00733	   0.03
             0.9	0.59439	    0.0034	   0.00];
%      krw = 0.004858*exp(10.6*Sw) + 0.00459*exp(10.66*Sw);
%      krow = 614.1*exp(-1.612*Sw) - 612.3*exp(-1.609*Sw);
%      Pcow = 10.05*exp(-4.453*Sw) + 103.5*exp(-17.19*Sw);
    for i = 1 : numel(Table(:,1)) - 1
        while Sw >= Table(i, 1) && Sw < Table(i+1, 1)
            ratio = (Sw - Table(i,1))/(Table(i+1,1) - Table(i,1));
            krw = ratio*(Table(i+1,2)-Table(i,2)) + Table(i,2);
            krow = ratio*(Table(i+1,3)-Table(i,3)) + Table(i,3);
            Pcow = ratio*(Table(i+1,4)-Table(i,4)) + Table(i,4);
            break
        end
    end
%     output = [krw, krow, Pcow];
end