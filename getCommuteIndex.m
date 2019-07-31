% getCommuteIndex.m
% Author: Olivia Waring
% COS 323 Final Project

function [ num ] = getCommuteIndex( i, j, k, track, N )
% Given a work borough, a home borough, a mode of transporation, a tracking
% matrix, and the total number of people in the simulation, return the 
% index of a random person who fulfills all these criteria (using track to 
% convert between the 3D array of types and the 1D vector of people).

   index = 25*(k-1) + 5*(j-1) + i;
   lower = track(index);
   if (index ~= 50)
       upper = track(index+1)-1;
   else 
       upper = N;
   end
   num = randi([lower, upper]); % return a random integer between the upper 
                                % and lower bounds given by track
end

