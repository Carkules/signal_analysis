function x = rectangular(A, t, T)
    x = A * sign(sin(2 * %pi * (t/T)));
endfunction

function x = triangular(A, t, T)
    x = 2*A*(abs(2*(t/T - floor(t/T + 0.5)))) - A;
endfunction

function x = sawtooth(A, t, T)
    x = 2*A*(t/T - floor(t/T + 0.5));
endfunction

function synthesis(NF, NT, x)
    //NF - liczba składników harmonicznych do przybliżania sygnału
    //NT - liczba okresów
    //x - wykorzystana funkcja fali
    
    T = 1;                  //okres
    N_per = 1000;           //liczba próbek na 1 okres
    N = N_per * NT;         // całkowita liczba próbek
    dt = (T * NT) / N;      //krok czasowy
    A = 1;                  //amplituda
    t = 0:dt:(T*NT - dt);   //oś czasu
    
    //generowanie funkcji fali z zadanymi parametrami
    xn = x(A, t, T)
    
    //obliczanie współczynników szeregu Fouriera
    a0 = sum(xn)*dt/(T*NT);
    an = zeros(1, NF);
    bn = zeros(1, NF);  

    for n = 1:NF
        an(n) = 2*sum(xn.*cos(2*%pi*n*t/T))*dt/(T*NT);
        bn(n) = 2*sum(xn.*sin(2*%pi*n*t/T))*dt/(T*NT);
    end

    //tworzenie sygnału za pomocą szeregu Fouriera
    y = zeros(1, N);
    for n = 1:NF
        y = y + an(n)*cos(2*%pi*n*t/T) +  bn(n)*sin(2*%pi*n*t/T);
    end

    plot(t,xn);         //rysowanie wykresu funkcji fali
    plot(t,y,'r');      //rysowanie sygnału wygenerowanego za pomocą szeregu Fouriera
    xlabel('czas (sek)');
    ylabel('Amplituda');
endfunction

function amplitude_modulation(NT, x)
    //NT - liczba okresów
    //x - wykorzystana funkcja fali
    
    ka = 0.5;               //amplituda modulacji
    T = 1;                  //okres
    N_per = 1000;           //liczba próbek na 1 okres
    N = N_per*NT;           //całkowita liczba próbek
    A = 1;                  //amplituda
    
    dt = (T * NT) / N;      //krok czasowy
    t = 0:dt:(T*NT - dt);   //oś czasu
    
    //generowanie funkcji fali z zadanymi parametrami
    xn = x(A, t, T);
    
    //generowanie sygnału modulującego
    fm = 0.1;
    m = sin(2 * %pi * fm * t);
    //m = fm*t
    
    //generowanie zmodulowanego sygnału
    s = A*(1 + ka*m) .* xn .* (xn ~= 0);
    plot(t, s);
    title('Modulacja amplitudy');
    xlabel('czas (sek)');
    ylabel('Amplituda');
endfunction

function frequency_modulation(NT, x)
    //NT - liczba okresów
    //x - wykorzystana funkcja fali
    
    N_per = 1000;           //liczba próbek na 1 okres
    A = 1;                  //amplituda
    T = 1;                  //okres
    N = N_per * NT;         //całkowita liczba próbek
    fc = 1/T;               //częstotliwość nośna
    kf = 3;                 //czułość częstotliwości
    dt = (T * NT) / N;      //krok czasowy
    t = 0:dt:(T*NT - dt);   //oś czasu
    
    //generowanie sygnału modulującego
    fm = 0.5;
    m = fm*(t)^2
    
    //obliczanie chwilowej fazy i przeskalowanie czasu
    phi = fc * t + kf * cumsum(m) * dt;
    t_fm = phi * T;
    
    //generowanie fali ze zmienną częstotliwością
    x_mod = x(A, t_fm, T);
    plot(t, x_mod);
    title('Modulacja częstotliwości');
    xlabel('czas (sek)');
    ylabel('Amplituda');
    a = gca();
    a.data_bounds = [min(t), min(x_mod)-0.1; max(t), max(x_mod)+0.1];
endfunction

function phase_modulation(NT, x)
    //NT - liczba okresów
    //x - wykorzystana funkcja fali
    
    N_per = 1000;           //liczba próbek na 1 okres
    N = N_per*NT;           //całkowita liczba próbek
    T = 1;                  //okres
    kp = %pi/2;             //czułość modulacji fazy
    A = 1;                  //amplituda
    dt = (T * NT) / N;      //krok czasowy
    t = 0:dt:(T*NT - dt);   //oś czasu
    
    //generowanie sygnału modulującego
    fm = 1;
    m = sin(2 * %pi * fm * t);
    //zmiana fazy
    phi = kp * m;
    x_pm = x(A, t + phi*T/(2*%pi), T); 
    
    plot(t, x_pm);
    title('Modulacja fazy');
    xlabel('czas (sek)');
    ylabel('Amplituda');
    a = gca();
    a.data_bounds = [min(t), min(x_pm)-0.1; max(t), max(x_pm)+0.1];
endfunction

phase_modulation(4, sawtooth)























