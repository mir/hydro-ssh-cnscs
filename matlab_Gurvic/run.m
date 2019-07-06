m = [];
h = 0.1;

interval_begin = 0.2;
interval_end = 1;

for gamma=interval_begin:h:interval_end
    
    res = minorsGurvic(gamma);
    
    if (res < 9)
        m = [m, gamma];
        disp(gamma);
    end

end

interval_begin = min(m);
interval_end = max(m);

m = [];

for gamma=interval_begin-h:h/10:interval_begin
    
    res = minorsGurvic(gamma);
    
    if (res < 9)
        m = [m, gamma];
                disp(gamma);
    end

end

interval_begin = min(m);

m = [];

for gamma=interval_end:h/10:interval_end+h
    
    res = minorsGurvic(gamma);
    
    if (res < 9)
        m = [m, gamma];
                disp(gamma);
    end

end

interval_end = max(m);

h = h / 10;

m = [];

for gamma=interval_begin-h:h/10:interval_begin
    
    res = minorsGurvic(gamma);
    
    if (res < 9)
        m = [m, gamma];
                disp(gamma);
    end

end

interval_begin = min(m);

m = [];

for gamma=interval_end:h/10:interval_end+h
    
    res = minorsGurvic(gamma);
    
    if (res < 9)
        m = [m, gamma];
                disp(gamma);
    end

end

interval_end = max(m);

h = h / 10;


m = [];

for gamma=interval_begin-h:h/10:interval_begin
    
    res = minorsGurvic(gamma);
    
    if (res < 9)
        m = [m, gamma];
                disp(gamma);
    end

end

interval_begin = min(m);

m = [];

for gamma=interval_end:h/10:interval_end+h
    
    res = minorsGurvic(gamma);
    
    if (res < 9)
        m = [m, gamma];
                disp(gamma);
    end

end

interval_end = max(m);

h = h / 10;

disp([interval_begin, interval_end]);
