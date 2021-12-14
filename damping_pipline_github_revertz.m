%% Read in power spectra
% EEG file paths go here
lo = 15; % index corresponding to 7 Hz in power spectrum
up = 81; % index corresponding to 40 Hz in power spectrum
N = up-lo+1;
Sub_list = [1,2,3,4,5]; % Subject number list
Sen_list = linspace(1,64,64); % channel list
Sub = length(Sub_list); % number of subjects
Sen = length(Sen_list); % number of channels
%% Compute power spectra for subjects and channels
EEG_EC = zeros(Sub,Sen,N);
EEG_EO = zeros(Sub,Sen,N);
JSDIV = zeros(Sub,Sen);

for i = 1:length(Sub_list)

    for j = 1:length(Sen_list)
        [PEC,F] = pwelch(LEC{Sub_list(i),1}(Sen_list(j),:),320,[],320,160); % EC welch periodogram for ith subject and jth channel - 
        [PEO,~] = pwelch(LEO{Sub_list(i),1}(Sen_list(j),:),320,[],320,160); % EO welch periodogram for ith subject and jth channel - 

       
        F = F(lo:up); % clip frequency between low/up limits
        PEC = PEC(lo:up); % clip spectrum between low/up limits
        PEO = PEO(lo:up); % clip spectrum between low/up limits
        EEG_EC(i,j,:) = PEC; 
        EEG_EO(i,j,:) = PEO;
        JSDIV(i,j) = jesh(PEC,PEO,F); % Jensen-Shannon divergence calculation function (jesh.m file)
    end
end
Freq_rest = double(F);


%% Alpha band damping distribution estimation pipeline
% Can define as a function to pass each subjects spectra matrix to
% Requires alpha_peak.m, Noise_mag.m, and kernel_int.m function files to run

D0 = size(psd,1); % number of channels
D1 = size(psd,2); % length of spectrum
ITER = 100; % numer of regularisation points
damping_matrix = zeros(D0,D1,ITER); % damping distribution matrix - empty
power_matrix = zeros(D0,D1,ITER);% psd matrix - empty
K_CONSTANT = zeros(D0,ITER); % normalisation constant vector - empty
LAMBDA = logspace(-1,-5,ITER); % regularisation values vector

N_rest = D1; % Dimension of power spectrum
sig_int_rest = logspace(-1,log10(50),N_rest); % Damping interval vector
X_rest = linspace(min(freq),max(freq),N_rest); % Frequency interval
WR = ones(N_rest,N_rest); % Coefficient matrix

LR = eye(length(sig_int_rest));
ZEROR = zeros(size(LR,1),1);

A_peak = zeros(D0,1); % Peak alpha matrix - empty
BETA = zeros(D0,1); % spectral scaling matrix - empty

we = diff(sig_int_rest);
we = [we, we(end)];


for z = 1:D0
    YR = psd(z,:);
    YR = squeeze(YR);
    b_rest = YR';
    ddr = [b_rest; ZEROR]; % part of the general least squares matrix
    
    [~, betaa] = Noise_mag(Freq_rest,YR); % spectral scaling (\beta) estimation function (Noise_mag.m file)
    BETA(z,1) = betaa;

    [Peak_alpha, ~] = alpha_peak(freq,squeeze(psd(z,:))); % Alpha peak and fwhm estimation function (alpha_peak.m file)
    A_peak(z,1) = Peak_alpha;

    for ii = 1:N_rest
        for jj = 1:N_rest

            WR(ii,jj) = we(1,jj)*kernel_int(X_rest(ii),sig_int_rest(jj),Peak_alpha);
        end
        
    end

    for k = 1:ITER
        lambR = LAMBDA(k);
        gamr = lambR*LR;
        CCR = [WR; gamr]; % second part of general least squares matrix
        xx = lsqnonneg(CCR,ddr); % solve general leasts squares with non negative constraint
        damping_matrix(z,:,k) = xx/(trapz(sig_int_rest,xx)); % append normalized damping distribution
        K_CONSTANT(z,k) = 1/(trapz(sig_int_rest,xx)); % append 1/normalization constant
        power_matrix(z,:,k) = WR*xx; % append generated model fit spectra

    end
end
