function UD = UpdateUD(UD, response)

trial = length(UD.x)
UD.response(trial) = response;

if trial == 1
    UD.initialUP   = true; 
    UD.initialDown = true; 
    UD.downCnt = 0; 
    UD.upCnt   = 0; 
    if response == 1
        UD.direction = -1;        
    else
        UD.direction = 1;
    end
end

if response == 1 % Down 
    UD.downCnt = UD.downCnt + 1; % 1 
    UD.upCnt   = 0; 
    
    if UD.direction == 1
        UD.initialUP = false; 
    end 
    UD.direction = -1; 

    if UD.downCnt == UD.down || UD.initialDown == true
        UD.downCnt = 0; 
        UD.upCnt   = 0; 
        UD.x(trial+1) = UD.x(trial) - UD.stepSizeDown; % 0.1 - 0.4 
    else
        UD.x(trial+1) = UD.x(trial); 
    end
    
else % Up
    UD.upCnt = UD.upCnt + 1;
    UD.downCnt = 0;
    
    if UD.direction == -1 
        UD.initialDown = false; 
    end 
    UD.direction = 1; 
    
    if UD.upCnt == UD.up || UD.initialUP == true
        UD.downCnt = 0; 
        UD.upCnt   = 0;
        UD.x(trial+1) = UD.x(trial) + UD.stepSizeUp;
        
    else
        UD.x(trial+1) = UD.x(trial);
    end 
end 

