%SteepestAscentFunc   Find maximum of of 4D parameters grid.
%
%Latest edited: 
%   06-10-2020


function [ parameters ] = SteepestAscentFunc(paramsValues, StimLevels, NumPos, OutOfNum, PF)
    % ALGORITHM PARAMETERS 
    dx    = 0.001; 
    dy    = 0.001; 
    dz    = 0.001;
    dm    = 0.001; 
    alpha = 0.1; 

    % Initial Guess 
    x0 = paramsValues(1);
    y0 = paramsValues(2);
    z0 = paramsValues(3); 
    m0 = paramsValues(4);


    % PREFROM ALGORITHM
    tol = 1e-100; 
    g   = [inf;inf;inf];

    n = 0;
    while norm(g) > tol 
        % compute gradient 
        f1 = SingleValueLikelihoodFunc(StimLevels,NumPos, OutOfNum,[x0-dx/2, y0, z0, m0], PF); 
        f2 = SingleValueLikelihoodFunc(StimLevels,NumPos, OutOfNum,[x0+dx/2, y0, z0, m0], PF);
        gx = (f2 - f1)/dx; 

        f1 = SingleValueLikelihoodFunc(StimLevels,NumPos, OutOfNum,[x0, y0-dy/2, z0, m0], PF); 
        f2 = SingleValueLikelihoodFunc(StimLevels,NumPos, OutOfNum,[x0, y0+dy/2, z0, m0], PF);
        gy = (f2 - f1)/dy; 

        f1 = SingleValueLikelihoodFunc(StimLevels,NumPos, OutOfNum,[x0,y0, z0-dz/2, m0], PF); 
        f2 = SingleValueLikelihoodFunc(StimLevels,NumPos, OutOfNum,[x0,y0, z0+dz/2, m0], PF);
        gz = (f2 - f1)/dz; 

        f1 = SingleValueLikelihoodFunc(StimLevels,NumPos, OutOfNum,[x0,y0, z0, m0-dm/2], PF); 
        f2 = SingleValueLikelihoodFunc(StimLevels,NumPos, OutOfNum,[x0,y0, z0, m0-dm/2], PF);
        gm = (f2 - f1)/dm;

        g = [ gx ; gy ; gz ; gm];


        % update Postion of Guess
        x0 = x0 + alpha*gx; 
        y0 = y0 + alpha*gy;     
        z0 = z0 + alpha*gz; 
        m0 = m0 + alpha*gm; 

        % counter
        n = n+1; 
    end 
    
    parameters = [x0 y0 z0 m0];

end 
   