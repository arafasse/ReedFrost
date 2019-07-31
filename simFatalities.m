clear all; clc; close all;

N = 10000; %population size
length = 25; % time periods
T = 10; %number of simulations

actives = zeros(T, length);
susc = zeros(T, length);
immune = zeros(T, length);
fatals = zeros(T, length);

cont = 0.5;    %level of contagiousness
friendly = 10; %avg no. of people someone will contact over 1 time period
p = cont*friendly/(N-1); %prob. that two individuals selected at random 
                         %come into effective contact
                         
q = 1-p;

actives(1:end, 1) = 1; %set one person infected at time zero
immune(1:end, 1) = 0;
susc(1:end, 1) = N - actives(1) - immune(1);

f = 0.01; %mortality rate, (taken roughly from Easterday, 1980)

for k = 1:T
    for i = 2:length
        prob = 1 - q^actives(k, i-1);
        
        for j = 1:susc(k,i-1)
            random = rand(1,1);
            if random < prob
                actives(k,i) = actives(k,i) + 1;
            end;
        end;
        susc(k,i) = susc(k,i-1) - actives(k,i);
        immune(k,i) = immune(k,i-1) + actives(k,i-1); %because all of the 
        %active cases become immune after one time interval unless they
        %become a fatality

        %calculate fatalities based on f
        for j = 1:actives(k, i-1)
            %determine whether a fatality occurred with this individual
            randomFatal = rand(1, 1);    
            if randomFatal < f
                fatals(k, i) = fatals(k, i) + 1;
                immune(k, i) = immune(k, i) - 1;
                q = 1 - cont*friendly/(actives(k, i) + susc(k, i) + immune(k, i));
            end
        end
    end
end

figure (1)
hold on;
for k = 1:T
    plot(1:length, actives(k,:), 'k');
end;
title('Number of active cases over time for 10 different simulations for N=10,000, cont=2, friendly = 10');
xlabel('Time');
ylabel('Number of active cases');
hold off;

figure (2)
hold on;
for k = 1:T
    plot(1:length, susc(k,:), 'k');
end;
title('Number of susceptible cases over time for 10 different simulations for N=10,000, cont=2, friendly = 10');
xlabel('Time');
ylabel('Number of active cases');
hold off;

figure (3)
hold on;
for k = 1:T
    plot(1:length, immune(k,:), 'k');
end;
title('Number of immune cases over time for 10 different simulations for N=10,000, cont=2, friendly = 10');
xlabel('Time');
ylabel('Number of active cases');
hold off;

figure(4)
hold on;
for k = 1:T
    plot(1:length, fatals(k, :), 'k');
end
title('Number of fatalities over time for 10 different simulations for N=10,000, cont=2, friendly=10');
xlabel('Time');
ylabel('Number of fatalities');