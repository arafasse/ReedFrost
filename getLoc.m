% getLoc.m
% Author: Olivia Waring
% COS 323 Final Project

function [ homeX homeY workX workY ] = getLoc( home, work )
% This function returns the x and y coordinates of random home and work
% locations in various boroughs of NYC. The different boroughs are
% represented by points on the map with a normal distribution radiating 
% outwards in all directions. 

    % defining boroughs
    MANHATTAN = 1;
    BRONX = 2;
    BROOKLYN = 3;
    QUEENS = 4;
    STATEN = 5;

    % generate the x and y coordinates of home
    if (home == BRONX)
        homeX = 7+0.5*randn(1);
        homeY = 8+0.5*randn(1);
    elseif (home == BROOKLYN)
        homeX = 5.5+0.5*randn(1);
        homeY = 3.5+0.5*randn(1);
    elseif (home == MANHATTAN)
        homeX = 5+0.25*randn(1);
        homeY = 6.5+0.75*randn(1);
    elseif (home == QUEENS)
        homeX = 7.5+0.5*randn(1);
        homeY = 5.5+0.5*randn(1);
    elseif (home == STATEN)
        homeX = 2+0.5*randn(1);
        homeY = 2+0.5*randn(1);
    end
    
    % generate the x and y coordinates of work
    if (work == BRONX)
        workX = 7+0.5*randn(1);
        workY = 8+0.5*randn(1);
    elseif (work == BROOKLYN)
        workX = 5.5+0.5*randn(1);
        workY = 3.5+0.5*randn(1);
    elseif (work == MANHATTAN)
        workX = 5+0.25*randn(1);
        workY = 6.5+0.75*randn(1);
    elseif (work == QUEENS)
        workX = 7.5+0.5*randn(1);
        workY = 5.5+0.5*randn(1);
    elseif (work == STATEN)
        workX = 2+0.5*randn(1);
        workY = 2+0.5*randn(1);
    end
end

