% getStats.m
% Author: Olivia Waring
% COS 323 Final Project

function [ susceptible infected immune deceased ] = getStats( status, time )
% This function returns the overall statistics of the NYC population at a
% given point in time.
    
    % define states
    SUSCEPTIBLE = 1;
    INFECTED = 2;
    IMMUNE = 3;
    DECEASED = 4;
    
    susceptible = 0;
    infected = 0;
    immune = 0;
    deceased = 0;
    [m n] = size(status);
    
    % populate arrays
    for i=1:m
        status(i,time);
        if status(i,time) == INFECTED
            infected = infected + 1;
        elseif status(i,time) == SUSCEPTIBLE
            susceptible = susceptible + 1;
        elseif status(i,time) == IMMUNE
            immune = immune + 1;
        elseif status(i,time) == DECEASED
            deceased = deceased + 1;
        end
    end
end

