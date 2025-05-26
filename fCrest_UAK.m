function [etac,etac_norm,Q,bins,params] = fCrest_UAK(etac_num,D,Tp,Hs)

%% Crest height mixture model (Al Khalili and Karmpadakis, 2025)
% Calculates the crest height distribution based on a transformation to
% second-order simulated crest heights

% Al Khalili (2025)


%% Inputs
% etac_num : array of second-order simulated crest heights
% D : sea-state directionality
% Tp : sea-state peak period
% Hs : sea-state significant wave height

% D, Tp, Hs must be stated in field scale and must corrrespond to one of
% the cases considered by (Al Khalili and Karmpadakis, 2025)

%% Outputs
% etac : array of modelled crest heights
% etac : array of modelled crest heights normalised by Hs
% Q : probability of exceedance
% bins : bins implemented in the model
% params : model parameters

%% Wave properties considered (in field scale)
Ds = [10,20];
Tps = [12,14,16];

if Tp ==12
    hs_cases =[22,44,67,89,112,134,157]./10;
elseif Tp==14
    hs_cases = [30,61,91,122,153,183]./10;
elseif Tp ==16
    hs_cases =[39,79,119,159,199]./10;
end

%% Select case of interest
iD = find(Ds==D);
iTp = find(Tps == Tp);
iHs = find(hs_cases==Hs);

if any([isempty(iD),isempty(iTp),isempty(iHs)])
    error('Sea-state not defined')
end

%% Load parameters
loadfol =  fileparts(mfilename('fullpath'));
loadfile = fullfile(loadfol,sprintf('Parameters.mat'));
load(loadfile,'bins','B','A');

%% Initialise variates
if ~issorted(etac_num,1,"ascend")
    etac_num = sort(etac_num,'ascend');
end
etac_num = reshape(etac_num,1,[]);
x = etac_num./Hs;

%% Model

% Number of model runs (for statistical sampling)
runs=100;

% Loop across runs
for run=1:runs

    % Set numerical data to model data
    X = x;

    % Split data into bins
    for j = 1:length(bins)-1

        % Get coefficients per bin:
        % Beta coefficients
        a_br = B(iD,iTp,iHs,j).a;
        b_br = B(iD,iTp,iHs,j).b;

        % Gamma coefficients
        a_amp = A(iD,iTp,iHs,j).a;
        b_amp = A(iD,iTp,iHs,j).b;

        % Probabilities
        p_amp = A(iD,iTp,iHs,j).p;
        p_br = B(iD,iTp,iHs,j).p;


        % Packaging
        params(j) = struct('p_amp',p_amp,'p_br',p_br,'a_br',a_br,'b_br',b_br,'a_amp',a_amp,'b_amp',b_amp);


        % Find data in bin
        id_bin = find(x>bins(j) & x<=bins(j+1));
        N_bin = length(id_bin);

        % Randomly select waves in bin
        indices = randperm(N_bin);

        % Modify amplified waves
        if p_amp>0
            N_amp = round(p_amp*N_bin);
            rat_amp = gamrnd(a_amp,b_amp,1,N_amp);
            ind_amp = indices(1:N_amp);
            X(id_bin(ind_amp)) =  rat_amp.*X(id_bin(ind_amp));
        else
            N_amp = 0;
        end

        % Modify breaking waves
        if p_br>0
            N_br = round(p_br*N_bin);
            rat_br = betarnd(a_br,b_br,1,N_br);
            ind_br = indices(N_amp+1:N_amp+N_br);
            X(id_bin(ind_br)) =  rat_br.*X(id_bin(ind_br));
        else
            N_br = 0;
        end

    end

    % Sort crest heights
    X_sort(run,:) = sort(X,'descend');

end

%% Average across runs
etac_norm=mean(X_sort,1);
etac = etac_norm.*Hs;

%% Set upper threshold
if any(etac_norm>bins(end))
    ind = find(etac_norm<bins(end));
    etac_norm = etac_norm(ind);
    etac = etac(ind);
end

%% Determine exceedance probability
etac_norm=sort(etac_norm,'descend');
etac=sort(etac,'descend');
Q = [1:length(etac)]./length(etac);


end
