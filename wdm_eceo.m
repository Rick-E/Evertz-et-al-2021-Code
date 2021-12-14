%% Weakly damped measure between EC/EO damping distributions

% EC and EO damping distribution are passed into the below
% the algorithm first detects peaks in the distribution
% then it selects the first peak and calculates the measure on that.

EC_WA = zeros(Sub,Sen);
EO_WA = zeros(Sub,Sen);

for i = 1:length(Sub_list)  
    for z = 1:length(Sen_list)
        
        % Weakly Damped Peak EC/EO
        
        [vals,locs] = findpeaks(squeeze(DEC(i,z,:))); % # of EC peak locations
        [vals1,locs1] = findpeaks(squeeze(DEO(i,z,:))); % # EO peak locations

        
        % EC area that falls below peak of EO and EC - peaks are present
        if isempty(locs) ~= 1 & isempty(locs1) ~= 1
            EC_WA(i,z) = trapz(SIG_INT(1:locs(1)) , DEC(i,z,1:locs(1)) ) + trapz(SIG_INT(1:locs1(1)) , DEC(i,z,1:locs1(1)));
            EO_WA(i,z) = trapz(SIG_INT(1:locs1(1)) , DEO(i,z,1:locs1(1)) ) + trapz(SIG_INT(1:locs(1)) , DEO(i,z,1:locs(1)));
            
        elseif isempty(locs) == 1 & isempty(locs1) == 1 % both states have no peaks
            locs = 0;
            locs1 = 0;
            EC_WA(i,z) = trapz(SIG_INT(1:locs(1)) , DEC(i,z,1:locs(1)) ) + trapz(SIG_INT(1:locs1(1)) , DEC(i,z,1:locs1(1)));
            EO_WA(i,z) = trapz(SIG_INT(1:locs1(1)) , DEO(i,z,1:locs1(1)) ) + trapz(SIG_INT(1:locs(1)) , DEO(i,z,1:locs(1)));
            
        elseif isempty(locs) == 1 | isempty(locs1) == 1 % One state has a peak and other does note
            
            if isempty(locs) == 1
                locs = 0;
            end
            if isempty(locs1) == 1
                locs1 = 0;
            end
            
            EC_WA(i,z) = trapz(SIG_INT(1:locs(1)) , DEC(i,z,1:locs(1)) ) + trapz(SIG_INT(1:locs1(1)) , DEC(i,z,1:locs1(1)));
            EO_WA(i,z) = trapz(SIG_INT(1:locs1(1)) , DEO(i,z,1:locs1(1)) ) + trapz(SIG_INT(1:locs(1)) , DEO(i,z,1:locs(1)));
            
        end

    end
end
