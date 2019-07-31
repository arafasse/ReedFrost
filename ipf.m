% ipf.m
% Author: Olivia Waring
% COS 323 Final Project

function [ A ] = ipf()
% This function performs iterative proportional fitting, generating the 
% individual matrix elements of a 3-dimensional data set for which only
% aggregate data is available. 

    lenX = 5; % place of work
    lenY = 5; % place of residence
    lenZ = 2; % mode of transportation
    XZ = [92, 188; 258, 408; 1495, 594;
           167, 429; 23, 96];
    YZ = ones(lenY, lenZ);
    XY = ones(lenX, lenY);
    work = [280; 667; 2089; 596; 120];
    live = [1332; 2465; 1537; 2229; 443];
    live = live.*sum(work)/sum(live); % must scale it so their total sums are equal
    mode = [2038; 1717];
    
    total = sum(work); % should be the same as sum(live) and sum(mode)!
    
    % generating the XY marginals
    for q=1:5
        for i=1:lenX
            XY(i,:) = XY(i,:)*work(i)/sum(XY(i,:)); 
        end
        for j=1:lenY
            XY(:,j) = XY(:,j)*live(j)/sum(XY(:,j));
        end
    end
    
    % generating the YZ marginals
    for q=1:5 
        for j=1:lenY
            YZ(j,:) = YZ(j,:)*live(j)/sum(YZ(j,:)); 
        end
        for k=1:lenZ
            YZ(:,k) = YZ(:,k)*mode(k)/sum(YZ(:,k));
        end
    end
    
    % filling in the elements of the 3D grid 
    A = ones(lenX, lenY);
    B = A;
    for r=2:lenZ
        A(:,:,r) = B;
    end
    
    for q=1:5
        % XZ
        for i=1:lenX
            for k=1:lenZ
                A(i,:,k) = A(i,:,k)*XZ(i,k)/sum(A(i,:,k));
            end
        end
        % YZ
        for j=1:lenY
            for k=1:lenZ
                A(:,j,k) = A(:,j,k)*YZ(j,k)/sum(A(:,j,k));
            end
        end
        % XY
        for i=1:lenX
            for j=1:lenY
                A(i,j,:) = A(i,j,:)*XY(i,j)/sum(A(i,j,:));
            end
        end
    end
end
    