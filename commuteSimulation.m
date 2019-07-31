% commute.m
% Author: Olivia Waring
% COS 323 Final Project

% This Matlab script is a microsimulation of a disease spreading through New 
% York City's five boroughs using the Reed-Frost model of infection. 
% Individual subjects are modeled as having a place of work, a place of 
% residence, of means of commuting, a health index (and accompanying time 
% threshold), and a location of work an residence. Subjects can occupy any
% four states: susceptible, infected, immune, and deceasd.

length = 10; % number of time periods
n = 5; % number of boroughs
numTypes = 50; % number of elements in 3-dimensional array 
plotSize = 10; % size of dots in scatter plot
numStart = 1; % number of infected individuals at t = 1

% defining boroughs
MANHATTAN = 1;
BRONX = 2;
BROOKLYN = 3;
QUEENS = 4;
STATEN = 5;

% defining states
SUSCEPTIBLE = 1;
INFECTED = 2;
IMMUNE = 3;
DECEASED = 4;

% defining types of transportation
PRIVATE = 1;
PUBLIC = 2; 

A = ipf; % call iterative potential fitting subroutine

N = sum(sum(A(:,:,1))) + sum(sum(A(:,:,2))); % total number of people
work_tots = sum(A(:,:,1)) + sum(A(:,:,2)); % number of workers in each borough
B = transpose(A(:,:,1));
B(:,:,2) = transpose(A(:,:,2));
home_tots = sum(B(:,:,1)) + sum(B(:,:,2)); % number of residents in each borough

N = N-27; % accomodate for the fact that 27 people are lost due to rounding errors

% initializing parallel arrays
home = zeros(N,1); % where the citizens live
work = zeros(N,1); % where the citizens work
mode = zeros(N,1); % how the citizens commute
health = zeros(N,1); % health of citizens 
duration = zeros(N,1); % how long the citizens have been in their current state
threshold = zeros(N,1); % duration after which health is reevaluated
location = zeros(N,4); % map coordinates (homeX, homeY, workX, workY)
dummy = 'a';
colors = repmat(dummy,1,N); % colors of scatter plot points
status = ones(N,length); % SUSCEPTIBLE, INFECTED, IMMUNE, or DECEASED

index = 1;
count = 1;

track = zeros(numTypes, 1); % keep track of index relationships
for k=1:2
    for j=1:n
        for i=1:n
            track(count) = index;
            count = count + 1;
            for p=round(1:A(i,j,k))
                home(index) = i;
                work(index) = j;
                mode(index) = k;
                health(index) = 5 + 2*randn(1); % normal distribution centered at 5, standard dev of 2
                threshold(index) = randi(10); % uniformly distributed random integer between 1 and 10
                colors(index) = 'b'; 
                [q r s t] = getLoc(i, j); % get home and work locations
                location(index,:) = [q r s t];
                index = index + 1;
            end  
        end
    end
end

% draw initial map
img=imread('mapNight.png');
min_x = 0;
max_x = 10;
min_y = 0;
max_y = 10;
figure(1);
imagesc([min_x max_x], [min_y max_y], flipdim(img,1));
hold on;

for a = 1:size(location(:,1))
    x = location(:,1); y = location(:,2);
    scatter(x(a), y(a), plotSize, colors(a), 'filled');
end
set(gca, 'ydir', 'normal');

susceptible = round(A); % number of susceptible people in each category

infected = zeros(n); % number of infected people in each category
infected(:,:,2) = infected; 

immune = zeros(n); % number of immune people in each categoty
immune(:,:,2) = immune;

deceased = zeros(n); % number of deceased people in each categoty
deceased(:,:,2) = deceased;

t = 1; % start the clock

% randomly choose infected individual
for i=1:numStart
    typhoid_Mary = randi(N); 
    status(typhoid_Mary,t) = INFECTED;
    colors(typhoid_Mary) = 'r';
end

% set initial values of infected and susceptible arrays
for i=1:N
    if (status(i, t) == INFECTED)
        infected(home(i), work(i), mode(i)) = infected(home(i), work(i), mode(i)) + 1;
        susceptible(home(i), work(i), mode(i)) = susceptible(home(i), work(i), mode(i)) - 1;
    end
end


p_home = 0.001; % probability that two people will come into effective contact in their home borough 
p_work = 0.002; % probability that two people will come into effective contact in their work borough
p_commute = 0.05; % probability that two people will come into effective contact on the subway
q_home = 1 - p_home;
q_work = 1 - p_work;
q_commute = 1 - p_commute;

for t = 2:length
    for i=1:n
        for j=1:n
            for k=1:2
                if (susceptible(i,j,k) > 0) % ensure there are still susceptible people in
                                            % a given subset of the
                                            % population
                    if (k == PUBLIC)
                        public_infected = infected(:,:,PUBLIC);
                        % Reed-Frost Model of infection 
                        % (calculate probability of contracting the disease 
                        % during the commute in question)
                        prob_commute = 1-q_commute^public_infected(i,j); 
                    else
                        prob_commute = 0;
                    end
                    sum_home = sum(sum(infected(i,:,:))); 
                    sum_work = sum(sum(infected(:,j,:))); 
                    % calculate the probability of contracting the disease in 
                    % the relevent home and work boroughs
                    prob_home = 1-q_home^sum_home;  
                    prob_work = 1-q_work^sum_work; 
  
                    % determine the number of "victims" to be infected
                    victims_home = round(prob_home * home_tots(i)); 
                    victims_work = round(prob_work * work_tots(j));
                    victims_commute = round(prob_commute * A(i,j,k));
                
                    % infect individuals who live in the same borough
                    for q = 1:victims_home
                        index = getHomeIndex(i,track,N); % determine the correct index 1 through N
                        % if susceptible, infect the chosen victim!
                        if((status(index,t) == SUSCEPTIBLE) && (susceptible(home(index), work(index), mode(index)) > 0))
                            status(index,t) = INFECTED;
                            colors(index) = 'r'; 
                            duration(index) = 0; % reset status duration
                            infected(home(index), work(index), mode(index)) = infected(home(index), work(index), mode(index)) + 1;
                            susceptible(home(index), work(index), mode(index)) = susceptible(home(index), work(index), mode(index)) - 1;
                        % if the victim is already infected and a sufficient
                        % amount of time has passed, determine future state
                        % based on health index
                        elseif((status(index,t) == INFECTED) && (duration(index) == threshold(index)))
                            random = randi(10);
                            if random <= health(index)
                                status(index,t) = IMMUNE;
                                colors(index) = 'g';
                            else
                                status(index,t) = DECEASED;
                                colors(index) = 'k';
                            end
                        end
                    end;    
                    
                    % infect individuals who work in the same borough
                    for q = 1:victims_work
                        index = getWorkIndex(j,track,N); % determine the correct index 1 through N
                        % if susceptible, infect the chosen victim!
                        if((status(index,t) == SUSCEPTIBLE)&&(susceptible(home(index), work(index), mode(index)) > 0))
                            status(index,t) = INFECTED;
                            colors(index) = 'r';
                            duration(index) = 0; % reset status duration
                            infected(home(index), work(index), mode(index)) = infected(home(index), work(index), mode(index)) + 1;
                            susceptible(home(index), work(index), mode(index)) = susceptible(home(index), work(index), mode(index)) - 1;
                         % if the victim is already infected and a sufficient
                         % amount of time has passed, determine future state
                         % based on health index
                         elseif((status(index,t) == INFECTED) && (duration(index) == threshold(index)))
                            random = randi(10);
                            if random <= health(index)
                                status(index,t) = IMMUNE;
                                colors(index) = 'g';
                            else
                                status(index,t) = DECEASED;
                                colors(index) = 'k';
                            end
                        end
                    end; 
                    
                    % infect individuals who share the same commuting path
                    for q = 1:victims_commute
                        index = getCommuteIndex(i,j,k,track,N); % determine the correct index 1 through N
                        % if susceptible, infect the chosen victim!
                        if((status(index,t) == SUSCEPTIBLE)&&(susceptible(home(index), work(index), mode(index)) > 0))
                            status(index,t) = INFECTED;
                            colors(index) = 'r';
                            duration(index) = 0; % reset status duration
                            infected(home(index), work(index), mode(index)) = infected(home(index), work(index), mode(index)) + 1;
                            susceptible(home(index), work(index), mode(index)) = susceptible(home(index), work(index), mode(index)) - 1;
                         % if the victim is already infected and a sufficient
                         % amount of time has passed, determine future state
                         % based on health index
                         elseif((status(index,t) == INFECTED) && (duration(index) == threshold(index)))
                            random = randi(10);
                            if random <= health(index)
                                status(index,t) = IMMUNE;
                                colors(index) = 'g';
                            else
                                status(index,t) = DECEASED;
                                colors(index) = 'k';
                            end
                        end
                    end;
                end;
           end;
        end;
   end;
   
   % copy current states into subsequent state vectors
   if (t ~= length)
        status(:,t+1) = status(:,t);
   end
   duration(:) = duration(:) + 1; 
   
   % update map
   min_x = 0;
   max_x = 10;
   min_y = 0;
   max_y = 10;
   figure(1);
   if (mod(t,2) == 0)
       img=imread('mapDay.png');
       imagesc([min_x max_x], [min_y max_y], flipdim(img,1));
       hold on;
       for a = 1:size(location(:,1))
           x = location(:,3); y = location(:,4);
           scatter(x(a), y(a), plotSize, colors(a), 'filled');
       end
   else
       img=imread('mapNight.png');
       imagesc([min_x max_x], [min_y max_y], flipdim(img,1));
       hold on;
       for a = 1:size(location(:,1))
           x = location(:,1); y = location(:,2);
           scatter(x(a), y(a), plotSize, colors(a), 'filled');
       end
   end
   set(gca, 'ydir', 'normal');
end

% accumulate longitudinal statistics for all of NYC
overallStats = zeros(2,length);
for t=1:length
    [susceptible infected immune deceased] = getStats(status,t);
    overallStats(1,t) = susceptible;
    overallStats(2,t) = infected;
    overallStats(3,t) = immune;
    overallStats(4,t) = deceased;
end

% accumulate longitudinal statistics for the Bronx
bronxStats = zeros(4,length);
for t=1:length
    [h i] = getHomeStats(status,home,t,BRONX);
    bronxStats(1,t) = h;
    bronxStats(2,t) = i;
    [h i] = getWorkStats(status,work,t,BRONX);
    bronxStats(3,t) = h;
    bronxStats(4,t) = i;
end

% accumulate longitudinal statistics for Brooklyn
brookStats = zeros(4,length);
for t=1:length
    [h i] = getHomeStats(status,home,t,BROOKLYN);
    brookStats(1,t) = h;
    brookStats(2,t) = i;
    [h i] = getWorkStats(status,work,t,BROOKLYN);
    brookStats(3,t) = h;
    brookStats(4,t) = i;
end

% accumulate longitudinal statistics for Manhattan
manhatStats = zeros(4,length);
for t=1:length
    [h i] = getHomeStats(status,home,t,MANHATTAN);
    manhatStats(1,t) = h;
    manhatStats(2,t) = i;
    [h i] = getWorkStats(status,work,t,MANHATTAN);
    manhatStats(3,t) = h;
    manhatStats(4,t) = i;
end

% accumulate longitudinal statistics for Queens
queensStats = zeros(4,length);
for t=1:length
    [h i] = getHomeStats(status,home,t,QUEENS);
    queensStats(1,t) = h;
    queensStats(2,t) = i;
    [h i] = getWorkStats(status,work,t,QUEENS);
    queensStats(3,t) = h;
    queensStats(4,t) = i;
end

% accumulate longitudinal statistics for Staten Island
statenStats = zeros(4,length);
for t=1:length
    [h i] = getHomeStats(status,home,t,STATEN);
    statenStats(1,t) = h;
    statenStats(2,t) = i;
    [h i] = getWorkStats(status,work,t,STATEN);
    statenStats(3,t) = h;
    statenStats(4,t) = i;
end

time = 1:length;

figure(2);
plot(time, overallStats(1,:), 'b', time, overallStats(2,:), 'r', time, overallStats(3,:), 'g', time, overallStats(4,:), 'k');
legend('susceptible', 'infected', 'immune', 'deceased');
xlabel('time periods');
ylabel('number of people');
title('Overall NYC Statistics');

figure(3);
plot(time, bronxStats(1,:), 'g', time, bronxStats(2,:), 'r', time, bronxStats(3,:), 'g--', time, bronxStats(4,:), 'r--');
legend('healthy residents', 'infected residents', 'healthy workers', 'infected workers');
xlabel('time periods');
ylabel('number of people');
title('Bronx Statistics');

figure(4);
plot(time, brookStats(1,:), 'g', time, brookStats(2,:), 'r', time, brookStats(3,:), 'g--', time, brookStats(4,:), 'r--');
legend('healthy residents', 'infected residents', 'healthy workers', 'infected workers');
xlabel('time periods');
ylabel('number of people');
title('Brooklyn Statistics');

figure(5);
plot(time, manhatStats(1,:), 'g', time, manhatStats(2,:), 'r', time, manhatStats(3,:), 'g--', time, manhatStats(4,:), 'r--');
legend('healthy residents', 'infected residents', 'healthy workers', 'infected workers');
xlabel('time periods');
ylabel('number of people');
title('Manhattan Statistics');

figure(6);
plot(time, queensStats(1,:), 'g', time, queensStats(2,:), 'r', time, queensStats(3,:), 'g--', time, queensStats(4,:), 'r--');
legend('healthy residents', 'infected residents', 'healthy workers', 'infected workers');
xlabel('time periods');
ylabel('number of people');
title('Queens Statistics');

figure(7);
plot(time, statenStats(1,:), 'g', time, statenStats(2,:), 'r', time, statenStats(3,:), 'g--', time, statenStats(4,:), 'r--');
legend('healthy residents', 'infected residents', 'healthy workers', 'infected workers');
xlabel('time periods');
ylabel('number of people');
title('Staten Island Statistics');



