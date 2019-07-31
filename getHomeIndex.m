% getHomeIndex.m
% Author: Olivia Waring
% COS 323 Final Project

function [ num ] = getHomeIndex( i, track, N )
% Given a particular borough, a tracking matrix, and the total number of 
% people in the simulation, return the index of a random person who resides 
% in the borough of interest (using track to convert between the 3D array 
% of types and the 1D vector of people).
   
   % generate random k and j values
   k = randi(2);
   j = randi(5);
   
   index = 25*(k-1) + 5*(j-1) + i;
   lower = track(index);
   if (index ~= 50)
       upper = track(index+1)-1;
   else 
       upper = N;
   end
   
   num = randi([lower upper]); % return a random integer between the upper 
                               % and lower bounds given by track
end

