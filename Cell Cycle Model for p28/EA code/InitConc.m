function x0_p =InitConc(x)

    k_cop1_1 = 0.1; % Moscetti et al. 
    k_cop1_2 = 1.1e-3; % Moscetti et al. 
    source_p53 = x(1);
    d_p53 = x(2);
    f1 = x(3); 
    d_Mdm2 = x(4);
    f2 = x(5);
    d_cop1 = x(6);
    k_deg_cop1 = x(7);
    deg_basal = x(9);

    x0_p = zeros(4, 1);
    x0_p(1) = (d_cop1 * d_Mdm2 * (k_cop1_2 + k_deg_cop1) * source_p53)/...
    (d_Mdm2 * f2 * k_cop1_1 * k_deg_cop1 + d_cop1 * (d_Mdm2 * d_p53 + deg_basal * f1) *  (k_cop1_2 + k_deg_cop1)); % p53
    x0_p(2) = f1 / d_Mdm2; % Mdm2
    x0_p(3) = f2 / d_cop1; % Cop1
    x0_p(4) = (d_Mdm2 * f2 * k_cop1_1 * source_p53)/...
    (d_Mdm2 * f2 * k_cop1_1 * k_deg_cop1 + d_cop1 *(d_Mdm2 * d_p53 + deg_basal * f1) * (k_cop1_2 + k_deg_cop1)); % p53/Cop1

end