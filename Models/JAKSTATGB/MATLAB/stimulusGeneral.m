function result = stimulusGeneral(timestm, epsilon, t, a)
%% General function for IFN stimulus
% t -- time
% epsilon -- 0.1
% timestm -- stimulus procedure cell array 
%              -- { {tBEGIN, tEND}, ... }
% a -- IFN stimulus level
%%
timestmsize = max(size(timestm));
itime = 1;
catchtime = false;
while (itime <= timestmsize && ~catchtime )
    timeinterval = timestm{itime};
    catchtime = (timeinterval(1) < t && timeinterval(2) >= t);
    if catchtime
        time = sort([timeinterval(1):1:timeinterval(2), timeinterval(1) + epsilon, timeinterval(2)  - epsilon]);
        GLNdata =a*[0, ones(1, max(size(time)) - 2), 0];
        result =  interp1(time,GLNdata,t,'pchip');
    end
    itime = itime + 1;
end

if ~catchtime 
    result = 0;
end

end
