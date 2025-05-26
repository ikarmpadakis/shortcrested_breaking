%% Start of Program
clear;
clc;
close all;


% Define parameters
gama_JONSWAP=2.5;
f_range = [26,160]./64;
th_range = [-45,45];
d = 0.5;
Fs = 128;
repeattime = 1024;
depth_case=50;
n=1;
% measurements are in lab scale: Tp_field/10, Hs_field/100: Length scale 1/100,time scale 1/10

% Add paths and input folders
restoredefaultpath
parentfol = fDrive_loc;
addpath(fullfile(parentfol,sprintf('Codes/Analysis')))
addpath(fullfile(parentfol,sprintf('Codes/Crests')))
addpath(fullfile(parentfol,sprintf('Codes/Figures')))
addpath(fullfile(parentfol,sprintf('Codes/Figures/matlab_fig_making')))
% Fig no
fig_no = {'a','b','c','d','e','f'};

% Define edges
bins = 0.2:0.2:1.6;

%% Tp and D cases
tp_cases = [12,14,16];
D_cases = [10,20];

%% Data folder
aligned_fol =  fullfile(parentfol,sprintf('Classification_tracked'));

%% Loop across sea-states
for iD = 1:length(D_cases)
    D = D_cases(iD);

    for iTp = 1:length(tp_cases)
        Tp = tp_cases(iTp);
        [gauge_int,hs_cases] = fgauge_selector(depth_case,Tp);

        for iHs = 1:length(hs_cases)
            % Ensure in field scale in m
            Hs = hs_cases(iHs)./10;
            Sps=0.01*iHs;

            %% Sea-state
            S(iD,iTp,iHs).D = D;
            S(iD,iTp,iHs).Tp = Tp;
            S(iD,iTp,iHs).Hs = Hs;


            % Amalgamate seeds
            matfiles = dir(fullfile(aligned_fol,sprintf(sprintf("Sp0.0%d_Tp%d_D%d_seed*.mat",iHs,Tp,D))));

            for i=1:length(matfiles)
                loadfile = fullfile(aligned_fol,matfiles(i).name);
                load(loadfile);

                % Indices
                ind_br = find(P2.dr<-0.1);
                ind_amp = find(P2.dr>0.1);
                ind_nbr = find(P2.dr>=-0.1 & P2.dr<=0.1);

                % Allocate ratios
                N = length(P2.r);
                T(i).rat_br = P2.r(ind_br);
                T(i).rat_nbr = P2.r(ind_nbr);
                T(i).rat_amp = P2.r(ind_amp);

                % Max spatial crest
                mg = size(P2.etac_sec,2)./2;
                for ii = 1:N
                    T(i).etac_lab(ii) = P2.etac_lab(ii,mg);
                    T(i).etac_num(ii) = P2.etac_sec(ii,mg);
                    T(i).etac_lab_norm(ii) = P2.etac_lab(ii,mg)./Hs_lab(mg);
                    T(i).etac_num_norm(ii) = P2.etac_sec(ii,mg)./Hs_sec(mg);
                    T(i).dr(ii) = P2.dr(ii);
                    T(i).r(ii) = P2.r(ii);
                end

                % Classify waves
                T(i).br_lab = T(i).etac_lab_norm(ind_br);
                T(i).nbr_lab = T(i).etac_lab_norm(ind_nbr);
                T(i).br_num = T(i).etac_num_norm(ind_br);
                T(i).nbr_num = T(i).etac_num_norm(ind_nbr);
                T(i).amp_lab = T(i).etac_lab_norm(ind_amp);
                T(i).amp_num = T(i).etac_num_norm(ind_amp);
            end

            %% Concatenate
            etac_lab = cat(2,T.etac_lab);
            etac_num = cat(2,T.etac_num);
            etac_lab_norm = cat(2,T.etac_lab_norm);
            etac_num_norm = cat(2,T.etac_num_norm);

            amp_lab = cat(2,T.amp_lab);
            br_lab = cat(2,T.br_lab);
            amp_num = cat(2,T.amp_num);
            br_num = cat(2,T.br_num);
            nbr_lab = cat(2,T.nbr_lab);
            nbr_num = cat(2,T.nbr_num);

            dr = cat(2,T.dr);
            r = cat(2,T.r);

            %% Save second-order simulations

            savefol = mkfol('Data/Numerical');
            savefile = fullfile(savefol,sprintf("D%d_Tp%d_Hs%1.1f.mat",D,Tp,Hs));
            save(savefile,'etac_num','etac_num_norm');

            %% Save lab
            savefol = mkfol('Data/Lab');
            savefile = fullfile(savefol,sprintf("D%d_Tp%d_Hs%1.1f.mat",D,Tp,Hs));
            save(savefile,'etac_lab','etac_lab_norm');

            %% Lengths
            N_amp = length(amp_lab);
            N_br = length(br_lab);
            N_nbr = length(nbr_lab);
            N_tot = N_amp+N_br+N_nbr;

            %% Probabilities
            P_amp = N_amp/N_tot;
            P_br = N_br/N_tot;
            P_nbr = N_nbr/N_tot;

            %% Ratios
            rat_amp = amp_lab./amp_num;
            rat_br = br_lab./br_num;
            rat_nbr = nbr_lab./nbr_num;

            %% Cap values
            rat_amp(rat_amp<1) = 1;
            rat_br(rat_br>1) = 1;

            %% Allocate data
            Pbr(iD,iTp,iHs) = P_br;
            Pamp(iD,iTp,iHs) = P_amp;
            Pnbr(iD,iTp,iHs) = P_nbr;

            Mbr(iD,iTp,iHs) = nanmean(nonzeros(rat_br));
            Mamp(iD,iTp,iHs) = nanmean(nonzeros(rat_amp));
            Mnbr(iD,iTp,iHs) = nanmean(nonzeros(rat_nbr));


            %% Get model parameters per interval
            for j = 1:length(bins)-1
                ind = find(etac_num_norm>bins(j) & etac_num_norm<=bins(j+1));
                id_amp = find(dr(ind)>=0.1);
                id_amp = ind(id_amp);
                id_br = find(dr(ind)<=-0.1);
                id_br = ind(id_br);
                id_nbr = find(dr(ind)<0.1 & dr(ind)>-0.1);
                id_nbr = ind(id_nbr);
                N_bin = length(id_amp) + length(id_br) + length(id_nbr);

                %% Histrograms
                if length(id_amp)>3
                    % Fit amp
                    tmp_amp = r(id_amp)';
                    pd_amp = fitdist(tmp_amp,'Gamma');
                    p_amp(j).a = pd_amp.a;
                    p_amp(j).b = pd_amp.b;
                    p_amp(j).P = length(id_amp)/N_bin;
                else                   
                    tmp_amp = NaN;
                    p_amp(j).a = NaN;
                    p_amp(j).b = NaN;
                    p_amp(j).P = NaN;
                end

                if length(id_br)>3
                    % Fit br
                    tmp_br = r(id_br)';
                    tmp_br = tmp_br(tmp_br<1);
                    pd_br = fitdist(tmp_br,'Beta');
                    p_br(j).a = pd_br.a;
                    p_br(j).b = pd_br.b;
                    p_br(j).P = length(id_br)/N_bin;
                else                   
                    tmp_br = NaN;
                    p_br(j).a = NaN;
                    p_br(j).b = NaN;
                    p_br(j).P = NaN;

                end

                %% Allocate data
                B(iD,iTp,iHs,j).a = p_br(j).a;
                B(iD,iTp,iHs,j).b = p_br(j).b;
                B(iD,iTp,iHs,j).p = p_br(j).P;
                B(iD,iTp,iHs,j).r = mean(tmp_br);

                A(iD,iTp,iHs,j).a = p_amp(j).a;
                A(iD,iTp,iHs,j).b = p_amp(j).b;
                A(iD,iTp,iHs,j).p = p_amp(j).P;
                A(iD,iTp,iHs,j).r = mean(tmp_amp);

                clear pd_amp pd_br

            end
            fprintf('iD = %d, iTp = %d, iHs = %d complete \n',iD,iTp,iHs);
        end
    end
end


%% Save coefficients
savefol = pwd;
savefile = fullfile(savefol,sprintf('Parameters.mat'));
save(savefile,'S','B','A',"Pbr",'Pamp','Pnbr','Mbr','Mnbr','Mamp','bins');



