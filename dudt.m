function Auplusb = dudt( t, u )
    global A  b    
    Auplusb = (A*u + b);
end

