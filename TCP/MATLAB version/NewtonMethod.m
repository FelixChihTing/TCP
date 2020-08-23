function T_true = NewtonMethod(Lamda_eqv1,Lamda_eqv2,T_ap0,T_a1,T_a2,TOL,N_MaxIt,alph_absor)
% Newton's method: solving the true temperature T given two apparent temperatures Ta1 and Ta2
% define variables
c2_planck=(1.98644*10^-25)/(1.3806*10^-23);
Lamda_eqv1=Lamda_eqv1*10^-9;   % Convert unit from nm to m
Lamda_eqv2=Lamda_eqv2*10^-9;   % Convert unit from nm to m
%T_a1=800;  % Measured parameters (Kelvin)
%T_a2=2500;  % Measured parameters (Kelvin)
%alph_absor=1.39;    % Emperical value

% Input
%T_ap0=1000; % initial approximation
%TOL=0.001;      % tolerence= T_apN-T_apN-1
%N_MaxIt=10;  % maximum mumber of iterations
format compact
% define the problem
n=1;
while n <= N_MaxIt
    fTn_1=real(((1-exp((c2_planck/Lamda_eqv1)*((1/T_ap0)-(1/T_a1))))^(Lamda_eqv1^alph_absor))-((1-exp((c2_planck/Lamda_eqv2)*((1/T_ap0)-(1/T_a2))))^(Lamda_eqv2^alph_absor)));
    derivfTn_1=real((Lamda_eqv1^alph_absor)*((1-exp((c2_planck/Lamda_eqv1)*((1/T_ap0)-(1/T_a1))))^((Lamda_eqv1^alph_absor)-1))*(-exp((c2_planck/Lamda_eqv1)*((1/T_ap0)-(1/T_a1))))*(-c2_planck/(Lamda_eqv1*(T_ap0^2)))-(Lamda_eqv2^alph_absor)*((1-exp((c2_planck/Lamda_eqv2)*((1/T_ap0)-(1/T_a2))))^((Lamda_eqv2^alph_absor)-1))*(-exp((c2_planck/Lamda_eqv2)*((1/T_ap0)-(1/T_a2))))*(-c2_planck/(Lamda_eqv2*(T_ap0^2))));
    T_apn=real(T_ap0-fTn_1/derivfTn_1);
    
    if abs(T_apn-T_ap0) < TOL
        T_true=T_apn;
        %disp(['The true temperature T_true is ',num2str(T_true),' K after ',num2str(n),' iterations.']);
        %disp(['The difference between the Nth and the N-1th iteration is ',num2str(abs(T_apn-T_ap0)),' K.']);
        return
    else
        n=n+1;
        T_ap0=T_apn;
    end
    if n== N_MaxIt
        %disp(['The method failed after ',num2str(N_MaxIt),' iterations.']);
        T_true=304;
    end
end
end