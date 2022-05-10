%% Principal Gradient analysis script -- FLE vs TLE, lang & working-memory

% First draft written May 2018, edited May 2019, September 2019, November 2019, April/May 2020, Jul/Aug
% 2020, Jan/Feb 2021, and Jul/Sep 2021
% Lorenzo Caciagli (lorenzo.caciagli@gmail.com)

% Input and help from Boris C. Bernhardt, Casey Paquola, & Xiaosong He is
% here acknlowedged with gratitude

% Link to open-access manuscript for which data was used -- https://doi.org/10.1093/brain/awac150 

% NOTE TO USER -- this script is to be intended as a code snippet and should serve as a guide, but is by no means exhaustive
% ANOTHER NOTE TO USER -- some dependent functions are necessary to get this script to work:

% 1) Surfstat, adjusted surfstat functions and other mica functions are available @ https://github.com/MICA-MNI/micaopen 

% 2) permutation-based t-test and correlation code is available at the below links 
% https://www.mathworks.com/matlabcentral/fileexchange/29782-mult_comp_perm_t1-data-n_perm-tail-alpha_level-mu-reports-seed_state
% https://www.mathworks.com/matlabcentral/fileexchange/54585-mult_comp_perm_t2-data1-data2-n_perm-tail-alpha_level-mu-t_stat-reports-seed_state
% https://www.mathworks.com/matlabcentral/fileexchange/34920-mult_comp_perm_corr

% please contact Lorenzo Caciagli lorenzo.caciagli@gmail.com and/or Boris
% Bernhardt (boris.bernhardt@mcgill.ca) for further queries – thanks! :)

%% 1. Define working directory and set SurfStat path

%P = '/data_/mica1/04_data/18_londonGradientsBB/matlab';
P = '/home/lorenzo/Wellcome/Analysis_MNI/matlab/';

addpath([P '/surfstat'])
addpath(P)

%% 2. Load the surface data 

SP = SurfStatAvSurf({[P '/fsaverage5/lh.pial'],[P '/fsaverage5/rh.pial']})
SW = SurfStatAvSurf({[P '/fsaverage5/lh.white'],[P '/fsaverage5/rh.white']})

% generate a mid thickness surface -- help with visualization later
SM.coord = (SP.coord + SW.coord)./2; 
SM.tri   = SP.tri;

%% 3. Load Spreadsheet with subject data 

% load csv file that contains our participant ids, groups and IVs 
fid      = fopen([P '/DB_LC_comb_perf_june_2020_FN.csv']); % final database
C        = textscan(fid,'%s%s%s%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%s%n%n%s','Delimiter',',',...
                    'headerLines',1,'CollectOutput',1);
fclose(fid);

% Language below: (cell 1 column 1, then cell 1 column 2) -- i.e. groups
% the variables according to their properties (string vs n, etc)
ID          = C{1}(:,1); 
GROUP       = C{1}(:,2);
GROUP_N     = C{1}(:,3);
AGE         = C{2}(:,1);
AGE_ONSET   = C{2}(:,2);
DURATION    = C{2}(:,3);
IQ          = C{2}(:,4);
LETTER      = C{2}(:,5);
CATEGORY    = C{2}(:,6);
DIGIT       = C{2}(:,7);
N_BACK_0B   = C{2}(:,8);
N_BACK_2B   = C{2}(:,9);
ZERO_BACK   = C{2}(:,10);
ONE_BACK    = C{2}(:,11);
TWO_BACK    = C{2}(:,12);
NAM         = C{2}(:,13);
LIST_A15    = C{2}(:,14);
LIST_A6     = C{2}(:,15);
TMT_A       = C{2}(:,16);
TMT_BA      = C{2}(:,17);
SIDE        = C{2}(:,18);
SEX         = C{2}(:,19);
HANDEDNESS  = C{2}(:,20);
SZ_FREQ     = C{2}(:,21);
FBTCS       = C{2}(:,22);
TIME_SINCE  = C{2}(:,23);
DEP         = C{2}(:,24);
ANX         = C{2}(:,25);
LES         = C{3};
FCD         = C{4}(:,1);
AED         = C{4}(:,2);
TPM_ZNS     = C{5};

%% 4. Load anonymized subject data (after QC) for the verbal N-back working memory (WM) task

TASK_T_NB = 'NB_con_0003'

%spm -- paths to data
% data here consist of: surface-based beta (contrast of parameter estimate maps) 
% for the 2-0 Back verbal working memory (NB) task contrast

left_NB  = strcat ('/home/lorenzo/Wellcome/Analysis_MNI/coreg3_con_max5mm/', ID, '_', TASK_T_NB, '_bbr_lh_fsa5.mgh');
right_NB = strcat ('/home/lorenzo/Wellcome/Analysis_MNI/coreg3_con_max5mm/', ID, '_', TASK_T_NB, '_bbr_rh_fsa5.mgh');

% Load data into a matrix - one matrix per task contrast
[TL keepL]       = SurfStatReadDataTry(left_NB);
[TR keepR]      = SurfStatReadDataTry(right_NB);
keepNB = sum(T,2)~=0;

%% 5. Keep subjects with no missing imaging data after QC (i.e. exclude missing subjects) -- and mask out midline & subcortical structures

load ([P '/mask_fsa5.mat']) % midline mask in fsaverage space

%Do for NB - exclude missing subjects and mask out midline/subcortical
%areas that are not used in surface-based analyses
%(generates Tk_short, etc)

Tk = T(keepNB,:);
%size(Tk)
Tk_short= Tk(:,mask);

%Keep covariates for subjects retained in analysis (will obviously be task-specific)

Gk_NB       = GROUP(keepNB);        % group info with info re: seizure focus laterality
Gnsk_NB     = GROUP_N(keepNB);      % group ID with no info re: seizure focus laterality info
Ak_NB       = AGE(keepNB);          % age 
OK_NB       = AGE_ONSET(keepNB);    % age at epilepsy onset
DK_NB       = DURATION(keepNB);     % duration
IK_NB       = ID(keepNB);           % participant code 
Sek_NB      = SEX(keepNB);          % (binary) sex
Sik_NB      = SIDE(keepNB);         % seizure focus laterality
Han_NB      = HANDEDNESS(keepNB);   % handedness
LES_NB      = LES(keepNB);          % lesional status
FCD         = FCD(keepNB);          % diagnosed with frontal FCD         

IQ_NB           = IQ(keepNB);               % IQ (NART)
DIGk_NB         = DIGIT(keepNB);            % Digit span
ZERO_BACK_NB    = N_BACK_0B(keepNB);        % verbal NB task performance, 0 Back condition
TWO_BACK_NB     = N_BACK_2B(keepNB);        % verbal NB task performance, 2 Back condition
TMTA_NB         = TMT_A(keepNB);            % trail making, A
TMT_BA          = TMT_BA(keepNB);           % trail making, B-A
SZ_Fk           = SZ_FREQ(keepNB);          % Sz frequency
FBTCSk          = FBTCS(keepNB);            % FBTCS history (0/1)
TIMEk           = TIME_SINCE(keepNB);       % time since last seizure

%% 6. Load Gradients (previously delived)

load([P '/fsaverage5/conte69.mat'])
template        = G; 

load([P '/fsaverage5/func_G1_conte69_sjhparc.mat']) % principal gradient -- courtesy of Casey Paquola
s =     func_G1_conte69_sjhparc;

load ([P '/fsaverage5/margulies_gradient.mat'])

%load([P 'SM.mat'])
G           = SM; 

% map conte69 to fsaverage5  
n2 = size(template.coord,2); 
for i = 1:size(G.coord,2)
    i
    d = sqrt(sum((repmat(G.coord(1:3,i),1,n2) - template.coord).^2)); 
    index = find(d==min(d)); 
    index_us(1,i) = index; 
    SonG(i) = s(index); 
end


SonGS           = SurfStatSmooth(SonG, G, 5);
save([P '/fsaverage5/func_G1_fsa5_sjhparc.mat'],'SonGS');

% Mask out midline and subcortical structure using the previously loaded
% masks

SonGS_short         = SonGS(:,mask);

% Display Gradient -- to generate figures
% Use different color scales - choose parula for final display

% margulies colormap (same as Margulies paper)
f=figure; 
    BoSurfStatViewData(SonGS, G, ''); 
    colormap(margulies)
    print ('Gradient_margulies', '-dpng', '-r600');
  
% jet colormap   
f=figure; 

    BoSurfStatViewData(SonGS, G, ''); 
    colormap(jet)
    print ('Gradient_jet', '-dpng', '-r600');
 
% hot colormap
f=figure; 
    BoSurfStatViewData(SonGS, G, ''); 
    colormap(hot)
    print ('Gradient_hot', '-dpng', '-r600');
    
% parula colormap    
f=figure; 
    BoSurfStatViewData(SonG, G, '');   
    colormap([parula; .8 .8 .8]);
    %SurfStatColLim ([-0.20 0.20 ])
    print ('Gradient_parula', '-dpng', '-r600');
    
%% 7. Print individual surface-based task maps (sanity check) -- to be visually inspected

%Print all
T = Tk;  
for i = 1:size(Tk,1) 

    f=figure;
        BoSurfStatViewData(T(i,:), SW, [IK_NB{i} 'smooth']) %or unsmooth
        BoSurfStatColLim([-1 1])
        exportfigbo(f,['/home/lorenzo/Wellcome/Analysis_MNI/Task_maps_2020/NB/' 'check_tmap_smooth.' IK_NB{i} '.png'],'png',10) 
        % one can do the same for unsmoothed maps if preferred
    close(f); 
end

%%  8. Map task-related changes on Principal Gradient - NB 

TS = Tk_short; %use data ending with "_short" suffix, i.e. where midline has been masked out already
f =  figure; 
     imagesc(TS,[-6 6]); 

Data            =       TS; 
Map             =       SonGS_short; %use Gradient 1 where midline has been masked out, i.e. with "_short" suffix 
window          =       4; % keep this even as per Lowe, Paquola et al. Human Brain Mapping 2019
num_bins        =       20;

% Map task betas on principal gradient, for all participants 
%provide both smoothed and unsmoothed option

[TSBin, TSBin_unsmo] = slidingWindowEqualDistrib(Data, Map, num_bins, window); 

% Show binned, gradient-mapped task data for all participants

f=figure; 
    imagesc(TSBin,[-1 1])
    
f=figure;   
    imagesc(TSBin_unsmo, [-1 1])
    
%%  9. Plot data showing task-related changes on gradient for CTR - group averages

% Define groupings -- labels will be used further down below

clear Data Group stdev_TS SEM_TS CI95 yCI95 x y c

Group_ns= Gnsk_NB;
G1ns=strcmp(Group_ns,'Group1');
G2ns=strcmp(Group_ns,'Group2');
G3ns=strcmp(Group_ns,'Group3');
G1_idx = find(G1ns == 1); %FLE
G2_idx = find(G2ns == 1); %CTR
G3_idx = find(G3ns == 1); %TLE

Group_S= Gk_NB;
G1s=strcmp(Group_S,'Group1');
G2s=strcmp(Group_S,'Group2');
G3s=strcmp(Group_S,'Group3');
G4s=strcmp(Group_S,'Group4');
G5s=strcmp(Group_S,'Group5');
G6s=strcmp(Group_S,'Group6');

G1s_NB_idx = find(G1s == 1); %L_FLE
G2s_NB_idx = find(G2s == 1); %R_FLE
G3s_NB_idx = find(G3s == 1); %CTR_wellcomeFL
G4s_NB_idx = find(G4s ==1);  %CTR_wellcomeTL
G5s_NB_idx = find(G4s == 1); %L_TLE
G6s_NB_idx = find(G5s == 1); %R_TLE
G7s_NB_idx = vertcat(G3s_NB_idx, G4s_NB_idx); % Concatenate controls in 1 group

% Statistics - one-sample t-test in CTR, task-activity (betas) across
% gradient bins

for i=1:(size(TSBin,2))

[pval_perm, t_perm, ~, ~, ~]= mult_comp_perm_t1(TSBin(G2_idx,i),10000,0,0.05,0,1);
NB_stats_CTR_perm(1,i)= pval_perm;
NB_stats_CTR_perm(2,i)= t_perm; 

end

%FDR correction of p values across all (20) gradient bins
NB_stats_CTR_perm_FDR= mafdr(NB_stats_CTR_perm(1,:),'BHFDR', true);
NB_stats_CTR_perm(4,:)= NB_stats_CTR_perm_FDR;

cd ([P '/1_Images/Grad_con_06_2020'])
save('NB_gradient_CTR_stats.mat', 'NB_stats_CTR_perm');


% Plot gradient-stratified task betas in controls

Data    = TSBin; 
Group   = Gnsk_NB; % consider removing

% Display data in CTR - via function (in this case, DisplayData_CTR_NB - provided separately,see folder)

f = DisplayData_CTR_NB(Data, Group_ns, -0.4, 0.8); % plotting bin-wise mean and 95% CI
cd ([P '/1_Images/Grad_con_02_2022'])
print ('NB_Task_on_gradient_CTR_20', '-dpng', '-r600'); 

%% 10. Group comparisons and Data display, FLE & FLE-FCD vs controls

% First, prepare data by adjusting for age and sex via multiple regression,
% and taking residuals forward for analysis

X_all=[Ak_NB Sek_NB];

for i=1:size(TSBin,2)
    lm = fitlm(X_all,TSBin(:,i), 'CategoricalVars',[2]);
    TSBin_res(:,i)= lm.Residuals{:,1};   
end

% Z-score gradient-stratified task betas, using controls as reference group

TSBin_z         = (TSBin - repmat(mean(TSBin(G2_idx,:),1),size(TSBin,1),1)) ...
                    ./    repmat(std(TSBin(G2_idx,:),0,1),size(TSBin,1),1);
               
                
TSBin_z_res         = (TSBin_res - repmat(mean(TSBin_res(G2_idx,:),1),size(TSBin_res,1),1)) ...
                    ./    repmat(std(TSBin_res(G2_idx,:),0,1),size(TSBin_res,1),1);             
              

% Sanity check -- group means, z-scored data

for i=1:size(TSBin,2)               
m1_z_res(i)                = mean(TSBin_z_res(G1_idx,i)); %FLE
m2_z_res(i)                = mean(TSBin_z_res(G2_idx,i)); %CTR
m3_z_res(i)                = mean(TSBin_z_res(G3_idx,i)); %TLE
end

% Sanity check -- group stdev, z-scored data

for i=1:size(TSBin,2)               
std1_z_res(i)                = std(TSBin_z_res(G1_idx,i)); %FLE
std2_z_res(i)                = std(TSBin_z_res(G2_idx,i)); %CTR
std3_z_res(i)                = std(TSBin_z_res(G3_idx,i)); %TLE
end


% Statistics - one sample t test on z scores, age- and sex-adjusted data, FLE vs CTR 

for i=1:(size(TSBin_z_res,2))

[pval_perm, t_perm, ~, ~, ~]= mult_comp_perm_t1(TSBin_z_res(G1_idx,i),10000,0,0.05,0,1);
NB_FLE_stats_perm(1,i)= pval_perm;
NB_FLE_stats_perm(2,i)= t_perm; 

end

% FDR correction of p values across all (20) gradient bins
NB_FLE_stats_perm_FDR= mafdr(NB_FLE_stats_perm(1,:),'BHFDR', true);
NB_FLE_stats_perm(4,:)= NB_FLE_stats_perm_FDR

% Cohen's d - FLE vs CTR

for i=1:(size(m1_z_res,2))
d_FLE_NB(i)= m1_z_res(i)/std1_z_res(i);
end

cd ([P '/1_Images/Grad_con_06_2020'])
save('NB_gradient_FLE_stats_07_20.mat', 'NB_FLE_stats_perm', 'd_FLE_NB')

%%% Now display data

%Display gradient-stratified, Z-scored betas, FLE vs CTR along 20 gradient
%bins

f= DisplayDataZ2Groups_FLE_NB(TSBin_z_res, Group_ns, -2.5, 4)
cd ([P '/1_Images/Grad_con_06_2020'])
print ('NB_Z2_FLE_CTR_20', '-dpng', '-r600');


%%% Statistics, code example for sensitivity analysis -- FLE-FCD vs CTR

% Mean (z scores) and stdev (z scores) in FLE FCD
for i=1:(size(TSBin_z_res,2))
FLE_FCD(:,i) = TSBin_z_res(FCD==1,i);
m_z_FLE_FCD(:,i)= mean(FLE_FCD(:,i));
std_z_FLE_FCD(:,i)= std(FLE_FCD(:,i));
end

%Statistics - one sample t test on z scores, age- and sex-adjusted data, FLE-FCD vs CTR

for i=1:(size(TSBin_z_res,2))

[pval_perm, t_perm, ~, ~, ~]= mult_comp_perm_t1(FLE_FCD(:,i),10000,0,0.05,0,1);
NB_FLE_FCD_stats_perm(1,i)= pval_perm;
NB_FLE_FCD_stats_perm(2,i)= t_perm; 

end

%FDR correction of p values across all (20) gradient bins
NB_FLE_FCD_stats_perm_FDR= mafdr(NB_FLE_FCD_stats_perm(1,:),'BHFDR', true);
NB_FLE_FCD_stats_perm(4,:)= NB_FLE_FCD_stats_perm_FDR;

%Cohen's d – FLE-FCD
for i=1:(size(m1_z_res,2))
d_FLE_FCD_NB(i)= m_z_FLE_FCD(i)/std_z_FLE_FCD(i);
end

cd ([P '/1_Images/Grad_con_06_2020'])
save('NB_gradient_FLE_FCD_stats.mat', 'NB_FLE_FCD_stats_perm', 'd_FLE_FCD_NB')

%Display Z-scored FLE-FCD vs CTR on gradient

f = DisplayDataZ2Groups_FLE_FCD_NB(TSBin_z_res, FCD, Group_ns, -3.25, 3)
cd ([P '1_Images/Grad_con_06_2020/'])
print ('NB_Z2_FLE_FCD_20', '-dpng', '-r600');

%% 11. Group comparisons and Data display, TLE vs CTR

% %Statistics -- one sample t test on z scores, age- and sex-adjusted data, TLE vs CTR 

for i=1:(size(TSBin_z_res,2))

[pval_perm, t_perm, ~, ~, ~]= mult_comp_perm_t1(TSBin_z_res(G3_idx,i),10000,0,0.05,0,1);
NB_TLE_stats_perm(1,i)= pval_perm;
NB_TLE_stats_perm(2,i)= t_perm; 

end

%FDR correction of p values across all (20) gradient bins
NB_TLE_stats_perm_FDR= mafdr(NB_TLE_stats_perm(1,:),'BHFDR', true);
NB_TLE_stats_perm(4,:)= NB_TLE_stats_perm_FDR;

%Cohen's d – TLE vs CTR
for i=1:(size(m3_z_res,2))
d_TLE_NB(i)= m3_z_res(i)/std3_z_res(i);
end

cd ([P '/1_Images/Grad_con_06_2020'])
save('NB_gradient_TLE_stats.mat', 'NB_TLE_stats_perm')

%%Display gradient-stratified, Z-scored betas, TLE vs CTR along 20 gradient
%bins

f= DisplayDataZ2Groups_TLE_NB(TSBin_z_res, Group_ns, -2.5, 4);
cd ([P '1_Images/Grad_con-06_2020/'])
print ('NB_Z2_TLE_CTR_20', '-dpng', '-r600');

%% 12. Group comparisons and Data display, FLE vs TLE

% First, prepare data via linear regression to adjust for age, sex and side
% (remember that for comparison against controls, age- and sex- adjusted
% metrics are used -- side comes in for comparison of patient groups

X_all_FT = [Ak_NB Sek_NB Sik_NB];
X_all_FT(G2_idx(1:end),:)=[];

TSBin_FT=TSBin;
TSBin_FT(G2_idx(1:end),:)=[];

for i=1:size(TSBin_FT,2)
    lm_FT_NB = fitlm(X_all_FT,TSBin_FT(:,i), 'CategoricalVars',[2,3]);
    TSBin_FT_res(:,i)= lm_FT_NB.Residuals{:,1};   
end

% Statistics – comparison of FLE and TLE - pvals, two-sample t test, age-
% and sex- and side-adjusted gradient-stratified betas

for i=1:(size(TSBin_FT_res,2))

[pval_perm, t_perm, ~, ~, ~]= mult_comp_perm_t2(TSBin_FT_res(G1_idx(1:end),i),TSBin_FT_res((G1_idx(end)+1):size(TSBin_FT_res,1),i),10000,0,0.05,0,'t',1);
NB_TLE_FLE_stats_perm_cov(1,i)= pval_perm;
NB_TLE_FLE_stats_perm_cov(2,i)= t_perm; 

end

% FDR correction of p values across all (20) gradient bins

NB_TLE_FLE_stats_perm_cov_FDR= mafdr(NB_TLE_FLE_stats_perm_cov(1,:),'BHFDR', true);
NB_TLE_FLE_stats_perm_cov(4,:)= NB_TLE_FLE_stats_perm_cov_FDR;

for i=1:(size(m3_z_res,2))
    std_pooled(i) = sqrt(((std1_z_res(i)^2 + std3_z_res(i)^2))/2);
    d_FLE_TLE_NB(i)= (m1_z_res(i)-m3_z_res(i))/std_pooled(i);
end

cd ([P '1_Images/Grad_con_06_2020/'])
save('NB_gradient_FLE_TLE_stats.mat','NB_TLE_FLE_stats_perm_cov', 'd_FLE_TLE_NB');

% Display Z-scored FLE/TLE/CTR on gradient

f = DisplayDataZ3Groups_cb_TL_FL_CTR_NB(TSBin_z_res, Group_ns, -2.25, 3);
cd ([P '1_Images/Grad_con_02_2022/'])
print ('NB_Z3_FLE_TLE_20', '-dpng', '-r600');

%% 13. Functional data analysis for between-group comparisons of gradient stratified curves -- FLE vs CTR, TLE vs CTR

%%% Permutation-based analysis conducted as in Bassett et al., Neuroimage
%%% 2012 (https://www.sciencedirect.com/science/article/abs/pii/S1053811911011633)
%%% Analyses conducted here are adjusted for age and sex

% First - define groups (obtained from merging another script)

Group_ns= Gnsk_NB;
G1ns=strcmp(Group_ns,'Group1');
G2ns=strcmp(Group_ns,'Group2');
G3ns=strcmp(Group_ns,'Group3');
G1_idx = find(G1ns == 1); %FLE
G2_idx = find(G2ns == 1); %CTR
G3_idx = find(G3ns == 1); %TLE

%Create group ID for permutation - after converting cell array, cutting out
%'Group', and converting to vector
Gf_NB=string(Gnsk_NB);
Gf_NB_new= erase(Gf_NB, "Group");
Gf_NB_newch= convertStringsToChars(Gf_NB_new);
Gf_NB_newch= cell2mat(Gf_NB_newch);
Gfa_NB= str2num(Gf_NB_newch);

%Get separate groups -- FLE/CTR, TLE/CTR, TLE/FLE (will be necessary to
%swap group allocation in the context of permutation-based method
Gfa_NB_FLE_CTR= Gfa_NB(Gfa_NB~=3);
Gfa_NB_TLE_CTR= Gfa_NB(Gfa_NB~=1);
Gfa_NB_FLE_TLE= Gfa_NB(Gfa_NB~=2);

% Regress out age and sex -- take residuals for further analysis

X_all=[Ak_NB Sek_NB];

for i=1:size(TSBin,2)
    lm = fitlm(X_all,TSBin(:,i), 'CategoricalVars',[2]);
    TSBin_res_ALL(:,i)= lm.Residuals{:,1};   
end

% Get separate NB arrays for TLE, FLE, CTR and combine
TSBin_res_CTR = TSBin_res_ALL(G2_idx,:);
TSBin_res_FLE = TSBin_res_ALL(G1_idx,:);
TSBin_res_TLE = TSBin_res_ALL(G3_idx,:);

TSBin_res_FLE_CTR = [TSBin_res_FLE; TSBin_res_CTR];
TSBin_res_TLE_CTR = [TSBin_res_CTR; TSBin_res_TLE];

%%% Now -- get to the stats, part 1! Compute Area between Groups according to Bassett et al. 2012

% Get real AbC - FLE vs CTR

for i=1:size(TSBin_res_ALL,2)
    Int_NB_res_FLE_CTR(i)=(abs(mean(TSBin_res_ALL(G2_idx,i))-mean(TSBin_res_ALL(G1_idx,i))));
end

AbC_NB_res_FLE_CTR= sum(Int_NB_res_FLE_CTR);

% Get real AbC - TLE vs CTR

for i=1:size(TSBin_res_ALL,2)
    Int_NB_res_TLE_CTR(i)=(abs(mean(TSBin_res_ALL(G2_idx,i))-mean(TSBin_res_ALL(G3_idx,i))));
end

AbC_NB_res_TLE_CTR= sum(Int_NB_res_TLE_CTR);

%%% Stats, part 2 -- Get permuted AbC values, 10000 permutations, and then
%%% compute permuted statistic for comparison of real vs permuted AbC

% Get permuted AbC - FLE vs CTR

for k=1:10000
    
    idx_rd = randperm(length(Gfa_NB_FLE_CTR));
    idx_rd_CTR = find(Gfa_NB_FLE_CTR(idx_rd)== 2);
    idx_rd_FLE = find(Gfa_NB_FLE_CTR(idx_rd)== 1);
    
    for i=1:size(TSBin_res_FLE_CTR,2)
    Int_NB_res_FLE_CTR_perm(i,k)=(abs(mean(TSBin_res_FLE_CTR(idx_rd_CTR,i))-mean(TSBin_res_FLE_CTR(idx_rd_FLE,i))));
    end
end


for k=1:10000
AbC_NB_res_FLE_CTR_perm(k)= sum(Int_NB_res_FLE_CTR_perm(1:i,k));
end 

% Compute p value for real AbC vs permuted - FLE vs CTR

Treal = AbC_NB_res_FLE_CTR;
Trand = AbC_NB_res_FLE_CTR_perm(1,:);
AbC_P_NB_res_FLE_CTR = sum(abs(Trand)>abs(Treal))/length(Trand);


% Get permuted AbC - TLE vs CTR


 for k=1:10000
    idx_rd = randperm(length(Gfa_NB_TLE_CTR));
    idx_rd_CTR = find(Gfa_NB_TLE_CTR(idx_rd)== 2);
    idx_rd_TLE = find(Gfa_NB_TLE_CTR(idx_rd)== 3);

    for i=1:size(TSBin_res_TLE_CTR,2)
 
   
    Int_NB_res_TLE_CTR_perm(i,k)=(abs(mean(TSBin_res_TLE_CTR(idx_rd_CTR,i))-mean(TSBin_res_TLE_CTR(idx_rd_TLE,i))));
    end 
end


for k=1:10000
AbC_NB_res_TLE_CTR_perm(k)= sum(Int_NB_res_TLE_CTR_perm(:,k));
end 

% Compute p value for real AbC vs permuted - FLE vs CTR

Treal = AbC_NB_res_TLE_CTR;
Trand = AbC_NB_res_TLE_CTR_perm(1,:);
AbC_P_NB_res_TLE_CTR = sum(abs(Trand)>abs(Treal))/length(Trand);

%% 14. Functional data analysis for between-group comparisons of gradient stratified curves -- FLE vs TLE


% For FLE and TLE - Regress out age, sex, and side -- take residuals for further analysis

X_all_FT = [Ak_NB Sek_NB Sik_NB];
X_all_FT(G2_idx(1:end),:)=[];
TSBin_FT=TSBin;
TSBin_FT(G2_idx(1:end),:)=[];

for i=1:size(TSBin_FT,2)
    lm_FT = fitlm(X_all_FT,TSBin_FT(:,i), 'CategoricalVars',[2,3]);
    TSBin_res_FT(:,i)= lm_FT.Residuals{:,1};   
end


% Get real AbC - FLE vs TLE (here contrast to be read as FLE - TLE)

for i=1:size(TSBin_res_FT,2)
    Int_NB_res_FLE_TLE(i)=(abs(mean(TSBin_FT_res(G1_idx(1:end),i))-mean(TSBin_FT_res((G1_idx(end)+1):size(TSBin_FT_res,1),i))));
end

AbC_NB_res_FLE_TLE= sum(Int_NB_res_FLE_TLE);


% Get permuted  AbC - FLE vs TLE (here contrast to be read as FLE - TLE)

for k=1:10000
    idx_rd = randperm(length(Gfa_NB_FLE_TLE));
    idx_rd_FLE = find(Gfa_NB_FLE_TLE(idx_rd)== 1);
    idx_rd_TLE = find(Gfa_NB_FLE_TLE(idx_rd)== 3);
    
    for i=1:size(TSBin_res_FT,2)
    Int_NB_res_FLE_TLE_perm(i,k)=(abs(mean(TSBin_res_FT(idx_rd_FLE,i))-mean(TSBin_res_FT(idx_rd_TLE,i))));
    end 
end


for k=1:10000
AbC_NB_res_FLE_TLE_perm(k)= sum(Int_NB_res_FLE_TLE_perm(:,k));
end 

%Calculate p value for real vs permuted AbC - FLE vs TLE

Treal = AbC_NB_res_FLE_TLE;
Trand = AbC_NB_res_FLE_TLE_perm(1,:);
AbC_P_NB_res_FLE_TLE= sum(abs(Trand)>abs(Treal))/length(Trand);

save ('NB_gradient_FDA_2021.mat', 'AbC_NB_res_FLE_CTR', 'AbC_NB_res_FLE_TLE', 'AbC_NB_res_TLE_CTR', ...
      'AbC_P_NB_res_FLE_CTR','AbC_P_NB_res_FLE_TLE', 'AbC_P_NB_res_TLE_CTR',...
      'AbC_NB_res_TLE_CTR_perm', 'AbC_NB_res_FLE_CTR_perm', 'AbC_NB_res_FLE_TLE_perm');
 
%% 15. Correlations of fMRI measures across systems with task/cognitive test performance

% Prepare data subset for permutation-based correlation analyses
% (also includes group allocation as a sanity check)

% NB fMRI data and Digit Span (DS) scores

NB_DS = horzcat(TSBin_res, DIGk_NB, Gfa_NB);
NB_DS(any(isnan(NB_DS),2),:)= [];

% Get correlation values across gradient bins

R_comb_NB_DS= zeros(2,(size(NB_DS,2)-2));

for i=1:size(TSBin_res,2)
    
[pval_perm, R_perm, ~, ~, ~]= mult_comp_perm_corr_orig(NB_DS(:,i),NB_DS(:,size(NB_DS,2)-1),10000,0,0.05,'linear',1); % change 'linear' to 'rank' if you want Spearman's rank correlations
R_comb_NB_DS(1,i)= R_perm;
R_comb_NB_DS(2,i)= pval_perm;

end

% FDR correction of p values across gradient bins

R_comb_NB_DS_FDR= mafdr(R_comb_NB_DS(2,:),'BHFDR', true);
R_comb_NB_DS(3,:)= R_comb_NB_DS_FDR;

% Plot correlational data - create figure

x       = 1:20;
y       = R_comb_NB_DS(1,:); % Correlation values (R) per bin, for all 20 bins
c       = 1:20;

%%%% Display data, DS & NB (gradient-stratified)

f = figure ('Position', [10 10 600 450]);
H.mainLine = color_line(x, y, c); 
colormap(parula)
H.mainLine.LineWidth = 8.0;
shading interp
xlim([1 20])
ax = gca;
ax.FontName = 'Gill Sans';
ax.FontSize = 18;
ax.XLim = [1 20]; % reflects number of bins
ax.YLim = [-0.1 0.4];
ax.YTick = -0.1:0.1:0.4;
ax.XTick = 2:2:20;

hold on

plot(1:20, zeros(1,20),'k','linewidth',1.10)

hold on 

% Plot asterisks to highlight statistically significant associations (both
% corrected and uncorrected for mc)

star1 = find(R_comb_NB_DS(3,:)<0.05);
star2 = find(R_comb_NB_DS(2,:)<0.05);

uncorr_only = setxor(star1, star2);

for k=1:numel(star1(1,:))
        
        if isempty(star1(k))
           plot(nan,nan)
        else
           plot(star1(k),0.35, '*', 'Color', [0 0 0], 'MarkerSize', 11, 'MarkerEdgeColor',[0 0 0], 'LineWidth',2);
        end
    hold on
    end
    
    
   for k=1:numel(uncorr_only(1,:))
        if isempty(uncorr_only(k))
           plot(nan,nan)
        else
           plot(uncorr_only(k), 0.35, 'v', 'Color', [0 0 0], 'MarkerSize', 11, 'MarkerEdgeColor',[0 0 0], 'LineWidth',2);
        end
    hold on
   end
     
cd ([P '1_Images/Grad_con_02_2022/'])
print ('Correl_NB_grad_DS_20', '-dpng', '-r600');   

%Scatterplot for correlation of DS scores and gradient-stratified NB fMRI activity in bin of
%maximum R value (here, 15) -- remember variable ordering of prepared data subset as above

f=figure('Position', [10 10 1300 900]);
f.Color=[1 1 1];

hd =plot(NB_DS(:,21),NB_DS(:,15),'.','MarkerSize',85,'Color',[0.75 0 0.93]); % PLOT EVERYONE;
hold on

set(gca, 'FontSize',54,'FontName','Gill Sans');
hd=gca; 
hd.XLim=([0 20]); 
hd.XTick=0:5:30; 
hd.YLim=[-1 1];
hd.YTick=-1:0.5:2;

hd1 = lsline(hd); % this command will generate the linear fit. 
hd1.LineWidth=10;
hd1.LineStyle='-.';
hd1.Color=[0.5 0.5 0.5];
hd1.AlignVertexCenters='off';

G1_l = find(NB_DS(:,22)==1);
G2_l = find(NB_DS(:,22)==2);
G3_l = find(NB_DS(:,22)==3);

hdc=plot(NB_DS(G2_l,21),NB_DS(G2_l,15),'.','MarkerSize',85,'Color',[0.238 0.304 0.434]); % Group2 - CTR
hdf=plot(NB_DS(G1_l,21),NB_DS(G1_l,15),'.','MarkerSize',85, 'Color',[0.758 0 0]); % Group1 - FLE
hdt=plot(NB_DS(G3_l,21),NB_DS(G3_l,15),'.','MarkerSize',85,'Color',[0.949 0.508 0.316]); % Group3 - TLE

cd ([P '1_Images/Grad_con_06_2020/'])
print ('Scatter_NB_DS_bin15_11_2020', '-dpng', '-r600');

close all
clear f hd hd1 hdc hdf hdt

%%%% Correlation of NB fMRI data with NB (2 Back) task performance –permutation-based

% Prepare data subset (slightly different from the above - group allocation
% not included this time, not really needed this time as there is no associated scatterplot)

NB_NB = horzcat(TSBin_res, TWO_BACK_NB);
NB_NB(any(isnan(NB_NB),2),:)= [];

% Compute permutation-based correlations

R_comb_NB_NB= zeros(2,(size(NB_NB,2)-1));

for i=1:size(TSBin_res,2)
[pval_perm, R_perm, ~, ~, ~]= mult_comp_perm_corr_orig(NB_NB(:,i),NB_NB(:,size(NB_NB,2)),10000,0,0.05,'linear',1);
R_comb_NB_NB(1,i)= R_perm;
R_comb_NB_NB(2,i)= pval_perm;
end

%FDR correction of p values across gradient bins

R_comb_NB_NB_FDR= mafdr(R_comb_NB_NB(2,:),'BHFDR', true);
R_comb_NB_NB(3,:)= R_comb_NB_NB_FDR;

% Plot correlational data - create figure

x       = 1:20;
y       = R_comb_NB_NB(1,:); 
c       = 1:20;

f = figure ('Position', [10 10 600 450]);
H.mainLine = color_line(x, y, c); 
colormap(margulies)
H.mainLine.LineWidth = 8.0;
shading interp
xlim([1 20])
ax = gca;
ax.FontName = 'Gill Sans';
ax.FontSize = 18;
ax.XLim = [1 20];
ax.YLim = [-0.1 0.4];
ax.YTick = -0.1:0.1:0.4;
ax.XTick = 2:2:20;

hold on

plot(1:20, zeros(1,20),'k','linewidth',1.10)

hold on 

% Plot asterisks to highlight statistically significant associations (both
% corrected and uncorrected for mc)

star1 = find(R_comb_NB_NB(3,:)<0.05);
star2 = find(R_comb_NB_NB(2,:)<0.05);

uncorr_only = setxor(star1, star2);


for k=1:numel(star1(1,:))
        
        if isempty(star1(k))
           plot(nan,nan)
        else
           plot(star1(k),0.35, '*', 'Color', [0 0 0], 'MarkerSize', 10);
        end
    hold on
end    
    

    
   for k=1:numel(uncorr_only(1,:))
        if isempty(uncorr_only(k))
           plot(nan,nan)
        else
           plot(uncorr_only(k), 0.35, 'v', 'Color', [0 0 0], 'MarkerSize', 8);
        end
    hold on
   end
     
cd ([P '1_Images/Grad_con_06_2020/'])
print ('Correl_NB_grad_NB_20', '-dpng', '-r600');   

save('NB_correl_perf_grad.mat', 'R_comb_NB_DS', 'R_comb_NB_NB');

%% Correlations of gradient-stratified task fMRI activity with clinical variables

% Here -- example with age at onset and FLE (same to be done for duration of epilepsy, seizure
% frequency, FBTCS, time since last seizure)

% Prepare data subset - FLE (includes regressing out effects of sex and side)

NB_Ao_FL= horzcat(TSBin(G1_idx,:), OK_NB(G1_idx), Sik_NB(G1_idx), Sek_NB(G1_idx));
NB_Ao_FL(any(isnan(NB_Ao_FL),2),:)= [];

for i=1:size(NB_Ao_FL,2)-3
    lm_NB_Ao = fitlm(NB_Ao_FL(:,22:23),NB_Ao_FL(:,i), 'CategoricalVars',[1,2]);
    NB_Ao_FL_res(:,i)= lm_NB_Ao.Residuals{:,1};   
end

% Compute correlations

R_comb_NB_Ao_FL= zeros(2,(size(NB_Ao_FL,2)-3));

for i=1:size(TSBin_res,2)
[pval_perm, R_perm, ~, ~, ~]= mult_comp_perm_corr_orig(NB_Ao_FL_res(:,i),NB_Ao_FL(:,size(NB_Ao_FL,2)-2),10000,0,0.05,'linear',1) % change 'linear' to 'rank' if you want Spearman's rank correlations
R_comb_NB_Ao_FL(1,i)= R_perm
R_comb_NB_Ao_FL(2,i)= pval_perm
end

%FDR correction of p values across gradient bins

R_comb_NB_Ao_FL_FDR= mafdr(R_comb_NB_Ao_FL(2,:),'BHFDR', true);
R_comb_NB_Ao_FL(3,:)= R_comb_NB_Ao_FL_FDR

% Display correlation data - age at onset and NB fMRI, FLE

x       = 1:20;
y       = R_comb_NB_Ao_FL(1,:); %N-Back on gradient,for controls only to start with -- so 54:75, in this case
c       = 1:20;


f = figure ('Position', [10 10 600 450]);

H.mainLine = color_line(x, y, c); 
colormap(parula)
H.mainLine.LineWidth = 8.0;
shading interp
xlim([1 20])
ax = gca;
ax.FontName = 'Gill Sans';
ax.FontSize = 18;
ax.XLim = [1 20];
ax.YLim = [-0.1 0.5];
ax.XTick = 2:2:20;

hold on

plot(1:20, zeros(1,20),'k','linewidth',1.10)

hold on 

% Plot asterisks to highlight statistically significant associations (both
% corrected and uncorrected for mc)

star1 = find(R_comb_NB_Ao_FL(3,:)<0.05);
star2 = find(R_comb_NB_Ao_FL(2,:)<0.05);

uncorr_only = setxor(star1, star2);


for k=1:numel(star1(1,:))
        
        if isempty(star1(k))
           plot(nan,nan)
        else
           plot(star1(k),0.45, '*', 'Color', [0 0 0], 'MarkerSize', 11, 'MarkerEdgeColor',[0 0 0], 'LineWidth',2);
        end
    hold on
    end
    
    
   for k=1:numel(uncorr_only(1,:))
        if isempty(uncorr_only(k))
           plot(nan,nan)
        else
           plot(uncorr_only(k), 0.45, 'v', 'Color', [0 0 0], 'MarkerSize', 11, 'MarkerEdgeColor',[0 0 0], 'LineWidth',2);
        end
    hold on
   end


cd ([P '1_Images/Grad_con_02_2022/'])
print ('Correl_NB_Ao_FLE_20_sex_side_cov', '-dpng', '-r600');

% Scatter plot - age at onset and gradient-stratified fMRI activity, 16th
% gradient bin

f=figure('Position', [10 10 1300 900]);
f.Color=[1 1 1];

hd = plot(NB_Ao_FL(:,21),NB_Ao_FL_res(:,16),'.','MarkerSize',85,'Color',[0.758 0 0]); % PLOT EVERYONE;
hold on

set(gca, 'FontSize',54,'FontName','Gill Sans');
hd=gca; 
hd.XLim=([0 41]); 
hd.XTick=0:10:50; 
hd.YLim=[-1 1];
hd.YTick=-1:0.5:2;

hd1 = lsline(hd); % generate the linear fit
hd1.LineWidth=10;
hd1.LineStyle='-.'
hd1.Color=[0.5 0.5 0.5];
hd1.AlignVertexCenters='off';

cd ([P '1_Images/Grad_con_06_2020/'])
print ('Scatter_NB_FLE_Ao_11_2020_16thbin', '-dpng', '-r600');


% Plot FLE

x       = 1:20;
y       = R_comb_NB_Ao_FL(1,:);
c       = 1:20;

f = figure ('Position', [10 10 600 450]);
H.mainLine = color_line(x, y, c); 
colormap(margulies)
H.mainLine.LineWidth = 8.0;
shading interp
xlim([1 20])
ax = gca;
ax.FontName = 'Gill Sans';
ax.FontSize = 18;
ax.XLim = [1 20];
ax.YTick = -0.1:0.1:0.4;
ax.XTick = 2:2:20;

hold on

plot(1:20, zeros(1,20),'k','linewidth',1.10)

hold on 

% Plot asterisks to highlight statistically significant associations (both
% corrected and uncorrected for mc)

star1 = find(R_comb_NB_Ao_FL(3,:)<0.05)
star2 = find(R_comb_NB_Ao_FL(2,:)<0.05)

uncorr_only = setxor(star1, star2)


for k=1:numel(star1(1,:))
        
        if isempty(star1(k))
           plot(nan,nan)
        else
           plot(star1(k),0.375, '*', 'Color', [0 0 0], 'MarkerSize', 10);
        end
    hold on
    end
    
    
   for k=1:numel(uncorr_only(1,:))
        if isempty(uncorr_only(k))
           plot(nan,nan)
        else
           plot(uncorr_only(k), 0.375, 'v', 'Color', [0 0 0], 'MarkerSize', 8);
        end
    hold on
   end
   
cd ([P '1_Images/Grad_con_06_2020/'])
print ('Correl_NB_grad_Ao_FL', '-dpng', '-r600');

save('R_NB_Ao_grad_sex_side_cov.mat','R_comb_NB_Ao_FL')

%%%% End of code snippet %%%%