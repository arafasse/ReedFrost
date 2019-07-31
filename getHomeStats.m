% getHomeStats.m
% Author: Olivia Waring
% COS 323 Final Project

function [ healthy infected ] = getHomeStats( status, home, time, N )
% This function returns the overall statistics of the residents of a given
% borough at a particular point in time.

    % define states
    SUSCEPTIBLE = 1;
    INFECTED = 2;
    IMMUNE = 3;
    DECEASED = 4;
    
    healthy = 0;
    infected = 0;
    [m n] = size(status);
    
    % populate arrays
    for i=1:m
        if (((status(i,time) == INFECTED) || (status(i,time) == DECEASED)) && (home(i) == N))
            infected = infected + 1;
        elseif (((status(i,time) == SUSCEPTIBLE) || (status(i,time) == IMMUNE)) && (home(i) == N))
            healthy = healthy + 1;
        end
    end
end

