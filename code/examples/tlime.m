function [x_s,mu_s,extrnl,its] = tlime(a, gamma, eta, epsilon, x_s)
[n,~] = size(a);

delta = 10e-8;

x_s = x_s / norm(x_s);
r_norm = epsilon + 1;

doRQI = false;

mu_s = x_s'*a*x_s;
lambda = mu_s;

extrnl = true;
its = 0;

while r_norm > epsilon
    if doRQI
        mu = mu_s;
    else
        mu = gamma;
    end
       
    % updateShiftMu(mu, solver)
    
    temp = x_s;
    y = (a - mu*eye(n)) \ temp;
    omega = 1 / norm(y);
    
    mu_s_old = mu_s;
    
    if doRQI == false
        mu_s = gamma + (y' * x_s) * (omega * omega);
    else
        mu_s_update = (y' * x_s) * (omega * omega);
        mu_s = mu_s + mu_s_update;
    end
    
    if doRQI == false
        q_norm = omega;
    else
        temp = x_s;
        temp = temp + (mu_s - gamma)*y;
        q_norm = norm(temp) * omega;
    end
    
    if doRQI == false
        temp = x_s; temp = temp + (gamma - lambda)*y;
        r_norm = norm(temp) * omega;
    else
        temp = x_s; temp = temp - mu_s_update * y;
        r_norm = norm(temp) * omega;
    end
    
    x_s = y;
    x_s = x_s * omega;
    
    if (doRQI == false) && (q_norm < eta)
        extrnl = false;
        doRQI = true;
    end
    
    if (extrnl == false) && (doRQI == true) && (abs(mu_s - gamma) >= eta)
        doRQI = false;
    end
    
    if (extrnl == true) && (doRQI == false) && (abs(mu_s - mu_s_old) / abs(mu_s) < delta)
        doRQI = true;
    end
    
    its = its + 1;
end
end
