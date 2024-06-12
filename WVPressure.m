function E = WVPressure(dtas)
    a1=611.21 ; a2=18.678 ; a3=234.5; a4=257.14; Ta=273.16;
    %Calculation of E saturation water vapour from Buck formula
    Tc=dtas-Ta;
    E=a1*exp((a2-Tc/a3).*Tc./(a4+Tc));
