%Final Project v1.1

guess = [0.5 0.3]; %guess (radius (m), C_D)

%Optimize with initial guess and alpha = 2.5, eta = 0.0005;
answer = hooke_jeeves(@flight,guess,2.5,0.0005)

%Hooke-Jeeves translated from course textbook
function x = hooke_jeeves(f,x,alpha,epsilon)
num = 0;
gamma = 0.75;
y = f(x);
n = length(x)
signs = [-1 1];

while alpha>epsilon
    improved = 0;
    x_best = x;
    y_best = y;
    for i = 1:n
        for j = 1:length(signs)
            sign = signs(j);
            xprime = x+(sign*alpha*basis(i,n));
            yprime = f(xprime);
            if yprime<y_best
                x_best = xprime;
                y_best = yprime;
                improved = 1;
            end
        end
    end
    x = x_best;
    y = y_best;
    num = num + 1;
    if improved == 0
        alpha = alpha*gamma;
    end
end
num
end

%Basis Function translated from textbook
function k = basis(i, n) 
for k=1:n
    if i == k
        k=1;
    else
        k=0;
    end
end
end

%Flight function (objective)
function l = flight(guess) %Input radius output flight time
    %Create empty arrays for values
    rho = zeros(3600,1);
    x1 = zeros(3600,1);
    vx1 = zeros(3600,1);
    y1 = zeros(3600,1);
    ax1 = zeros(3600,1);
    temp = zeros(3600,1);
    q = zeros(3600,1);
    rho = zeros(3600,1);
    %Create history arrays
    r_hist = 0;
    Cd_hist = 0;
    l_hist = 0;

    r = guess(1);
    C_D = guess(2);
    ts = 0.005; %Thickness set to 5mm
    t = 1:1:3600;  %Create time vector of 3600 seconds (estimated maximum flight time)

    %Initial Conditions
        x1(1) = 54000; %Initial altitude of 54km
        vx1(1) = 0; %No initial vertical velocity
        vy1 = 200; %Initial horizontal velocity (m/s)
        ax1(1) = 8.87; %Gravity on Venus
        temp(1) = 270; %Atmospheric Temp at 54km
        tint(1) = 270; %Probe Initial Temp
        rho(1) = 7; %Atmospheric Density at 54km
        q(1) = 0; %Initial Heat Flux
        y(1) = 0; %Initial horizontal position
        m = 2 + (4/3)*3.14*(r^3)*4510 - (4/3)*pi*(r-ts)^3*4510; %hollow sphere + electronics (2 kg)
        area = pi*(r^2); %Area of sphere
        for i = 2:length(t)
            %Vehicle Dynamics
            %3rd order approximation of VIRA Density Data
            rho(i) = -0.000387*(x1(i-1)/1000)^3 + 0.0588*(x1(i-1)/1000)^2 - 3.23*(x1(i-1)/1000) + 64.6; 
            %Acceleration = (Weight - Drag)/Mass
            ax1(i) = ((m*8.87) - (C_D*rho(i-1)*(vx1(i-1)^2)*area)/2)/m;
            %Velocity
            vx1(i) = abs(vx1(i-1)+ax1(i-1));
            %Position (Y)
            x1(i) = x1(i-1)-vx1(i-1);
            %Position (X)
            y1(i) = y1(i-1)+vy1;
            %1st order approximation of VIRA Temperature Data
            temp(i) = -7.874*(x1(i-1)/1000)+737.8;
            %Heat Flux per second
            q(i) = 17*area*((temp(i-1)-tint(i-1))/ts);
            %Internal Temp -- Conduction 
            tint(i) = -(((ts*q(i-1)) - (17*area*temp(i-1)))/(17*area));
            if x1(i) < 0
                l = -i
                break
            end
        end
       heat = sum(q); %Sum of heat flux 
       maxtemp = max(tint); %Maximum Internal Temperature
       penalty = 0; %Initial Penalty Value
       
       if r <= 0 %If negative radius
           pen1 = ((0 - r)+1)^2
           penalty = penalty + 100*pen1;
       end
       if m <= 5 %If mass less than 5kg (from literature)
           pen2 = ((5 - m) + 1)^2
           penalty = penalty + pen2;
       end
       if maxtemp > 600 %If max temp exceeds 600K
           pen3 = ((maxtemp-600) + 1 )^2
           penalty = penalty + pen3;
       end
       if C_D <= 0.25 %If drag coefficient less than 0.25
           pen4 = ((0.25 - C_D)*100)^2
           penalty = penalty + 1000*pen4;
       end
       if m >= 100 %If mass greater than 100kg
           pen5 = ((m - 100) + 1)^2;
           penalty = penalty + 100*pen5;
       end
       l = l + 100*penalty; %Sum penalty times rho
       r_hist(end+1) = r; 
       Cd_hist(end+1) = C_D;
       l_hist(end+1) = l;

       scatter(r_hist,l_hist)
       xlim([0 1]);
       ylim([-1e7 3e7]);
       hold on
end