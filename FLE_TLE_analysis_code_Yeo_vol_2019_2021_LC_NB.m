%% Yeo systems analysis script -- FLE vs TLE, lang & working-memory

% First written 19/11/2019, edited December 2019, April/May 2020, Jul/Aug
% 2020, Jan/Feb 2021, and Jul/Sep 2021
% Lorenzo Caciagli (lorenzo.caciagli@gmail.com)

% Input and help from Boris C. Bernhardt, Casey Paquola, & Xiaosong He is
% here acknlowedged with gratitude

% Link to open-access manuscript https://doi.org/10.1093/brain/awac150 

% NOTE TO USER -- this script is to be intended as a code snippet and should serve as a guide, but is by no means exhaustive
% ANOTHER NOTE TO USER -- some dependent functions are necessary to get this script to work:

% 1) mica-spider and other mica functions are available @ https://github.com/MICA-MNI/micaopen 

% 2) permutation-based t-test and correlation code is available at the below links 
% https://www.mathworks.com/matlabcentral/fileexchange/29782-mult_comp_perm_t1-data-n_perm-tail-alpha_level-mu-reports-seed_state
% https://www.mathworks.com/matlabcentral/fileexchange/54585-mult_comp_perm_t2-data1-data2-n_perm-tail-alpha_level-mu-t_stat-reports-seed_state
% https://www.mathworks.com/matlabcentral/fileexchange/34920-mult_comp_perm_corr

% please contact Lorenzo Caciagli lorenzo.caciagli@gmail.com and/or Boris
% Bernhardt (boris.bernhardt@mcgill.ca) for further queries – thanks! :)

%% 1. Define working directory and set relevant paths

P = '/home/lorenzo/Wellcome/Analysis_MNI/matlab/';

addpath([P '/surfstat'])
addpath(P)

%% 2. Load anonymized database

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

%% 3. Load anonymized subject list (after QC) for the verbal N-back working memory (WM) task

TASK_T_NB = 'NB_con_0003';

%spm -- paths to QC info

left_NB  = strcat ('/home/lorenzo/Wellcome/Analysis_MNI/coreg3_con_max5mm/', ID, '_', TASK_T_NB, '_bbr_lh_fsa5.mgh');
right_NB = strcat ('/home/lorenzo/Wellcome/Analysis_MNI/coreg3_con_max5mm/', ID, '_', TASK_T_NB, '_bbr_rh_fsa5.mgh');

%% 4. Keep subjects with no missing imaging data after QC (i.e. exclude missing subjects)

Gk_NB       = GROUP(keepNB);                % group info with info re: seizure focus laterality
Gnsk_NB     = GROUP_N(keepNB);              % group ID with no info re: seizure focus laterality info
Ak_NB       = AGE(keepNB);                  % age
OK_NB       = AGE_ONSET(keepNB);            % age at epilepsy onset
DK_NB       = DURATION(keepNB);             % duration
IK_NB       = ID(keepNB);                   % participant code
Sek_NB      = SEX(keepNB);                  % (binary) sex
Sik_NB      = SIDE(keepNB);                 % seizure focus laterality
Han_NB      = HANDEDNESS(keepNB);           % handedness
LES_NB      = LES(keepNB);                  % lesional status
FCD         = FCD(keepNB);                  % diagnosed with frontal FCD

IQ_NB           = IQ(keepNB);               % IQ (NART)
DIGk_NB         = DIGIT(keepNB);            % Digit span
ZERO_BACK_NB    = N_BACK_0B(keepNB);        % verbal NB task performance, 0 Back condition
TWO_BACK_NB     = N_BACK_2B(keepNB);        % verbal NB task performance, 2 Back condition
TMTA_NB         = TMT_A(keepNB);            % trail making, A
TMT_BA          = TMT_BA(keepNB);           % trail making, B-A
SZ_Fk           = SZ_FREQ(keepNB);          % Sz frequency
FBTCSk          = FBTCS(keepNB);            % FBTCS history (0/1)
TIMEk           = TIME_SINCE(keepNB);       % time since last seizure

%% 5. Load Task-related betas (contrasts of parameter estimates) in volume space 
% Concatenated as: LFLE, RFLE, CTR, CTR_TL, LTLE, RTLE

% task-related betas (2-0 Back) extracted using the Schaefer Atlas
% (https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal) 
% MNI 200 parcel partition with parcel-wise Yeo/Krienen 7 systems IDs

Tk  = load('/home/lorenzo/Wellcome/Analysis_MNI/Con_maps_volume/Yeo/NB/Extracted_con_NB_Schaefer200.txt');

% Use the below when repeating analyses with the Schaefer Atlas (https://pubmed.ncbi.nlm.nih.gov/28981612/) 200
% parcels (MNI) for Yeo/Krienen 17 systems

% Tk  = load('/home/lorenzo/Wellcome/Analysis_MNI/Con_maps_volume/Yeo/NB/Extracted_con_NB_Schaefer200_17.txt');

%% 6. Map task-related changes on Yeo networks (using mica-spyder for spider/radar plots) -- NB

%Define groupings -- labels will be used further down below

Group_ns= Gnsk_NB;
G1ns=strcmp(Group_ns,'Group1');
G2ns=strcmp(Group_ns,'Group2');
G3ns=strcmp(Group_ns,'Group3');
G1_NB_idx = find(G1ns == 1);        % FLE
G2_NB_idx = find(G2ns == 1);        % CTR
G3_NB_idx = find(G3ns == 1);        % TLE

Group_S= Gk_NB;
G1s=strcmp(Group_S,'Group1');
G2s=strcmp(Group_S,'Group2');
G3s=strcmp(Group_S,'Group3');
G4s=strcmp(Group_S,'Group4');
G5s=strcmp(Group_S,'Group5');
G6s=strcmp(Group_S,'Group6');

G1s_NB_idx = find(G1s == 1);                    % L_FLE
G2s_NB_idx = find(G2s == 1);                    % R_FLE
G3s_NB_idx = find(G3s == 1);                    % CTR_wellcomeFL
G4s_NB_idx = find(G4s ==1);                     % CTR_wellcomeTL
G5s_NB_idx = find(G4s == 1);                    % L_TLE
G6s_NB_idx = find(G5s == 1);                    % R_TLE
G7s_NB_idx = vertcat(G3s_NB_idx, G4s_NB_idx);   % Concatenate controls in 1 group


% Derive betas across Schaefer 200 parcels, 7 Yeo/Krienen systems - for
% all participants

% Order to be displayed in spider/radar plots which (remember) does not coincide with systems labeling in the atlas:
% 1 dorsal attn (DAN), 2 frontoparietal control (FPN), 3
% default-mode (DMN), 4 Limbic (Lim), 5 Visual (Vis), 6 somatomotor (SM), ventral attention/salience (SAL)

for i=1:size(Tk,1)

Tk_yeo(i,1) = mean(Tk(i,[31:43,135:147]));  %DAN
Tk_yeo(i,2) = mean(Tk(i,[61:73,165:181]));  %FPN
Tk_yeo(i,3) = mean(Tk(i,[74:100,182:200])); %DMN
Tk_yeo(i,4) = mean(Tk(i,[55:60,159:164]));  %Lim
Tk_yeo(i,5) = mean(Tk(i,[1:14,101:115]));   %Vis
Tk_yeo(i,6) = mean(Tk(i,[15:30,116:134]));  %SM
Tk_yeo(i,7) = mean(Tk(i,[44:54,148:158]));  %SAL

end

% Derive betas across Schaefer 200 parcels, 17 Yeo/Krienen systems -
% for all participants
% Order as above = 1 DAN A/B, FPN A/B/C, DMN A/B/C, Limb A/B, Vis A/B, SM
% A/B, SAL A/B

% for i=1:size(Tk,1);
%  
% Tk_yeo(i,1) = mean(Tk(i,[29:34,131:136]));    %DAN A     --5
% Tk_yeo(i,2) = mean(Tk(i,[35:39,137:141]));    %DAN B     --6
% Tk_yeo(i,3) = mean(Tk(i,[57:66,165:170]));    %FPN A     --11
% Tk_yeo(i,4) = mean(Tk(i,[67:71,171:180]));    %FPN B     --12
% Tk_yeo(i,5) = mean(Tk(i,[72:74,181:183]));    %FPN C     --13
% Tk_yeo(i,6) = mean(Tk(i,[75:82,184:189]));    %DMN A     --14
% Tk_yeo(i,7) = mean(Tk(i,[83:95,190:193]));    %DMN B     --15
% Tk_yeo(i,8) = mean(Tk(i,[96:98,194:196]));    %DMN C     --16
% Tk_yeo(i,9) = mean(Tk(i,[99:100,197:200]));   %TP        --17
% Tk_yeo(i,10) = mean(Tk(i,[53:56,161:164]));   %Limbic A  --10
% Tk_yeo(i,11) = mean(Tk(i,[51:52,157:160]));   %Limbic B  --9
% Tk_yeo(i,12) = mean(Tk(i,[1:6,101:106]));     %VisA      --1
% Tk_yeo(i,13) = mean(Tk(i,[7:12,107:112]));    %VisB      --2
% Tk_yeo(i,14) = mean(Tk(i,[13:20,113:123]));   %SM A      --3
% Tk_yeo(i,15) = mean(Tk(i,[21:28,124:130]));   %SM B      --4
% Tk_yeo(i,16) = mean(Tk(i,[40:46,142:150]));   %Sal A     --7
% Tk_yeo(i,17) = mean(Tk(i,[47:50,151:156]));   %Sal B     --8
% 
% end

% Derive group-wise mean betas across systems 

for i=1:size(Tk_yeo,2)
m1(i)   = mean(Tk_yeo(G1_NB_idx,i)); %FLE
m2(i)   = mean(Tk_yeo(G2_NB_idx,i)); %CTR
m3(i)   = mean(Tk_yeo(G3_NB_idx,i)); %TLE
end


% Derive betas (residuals) across sytems - after adjustment for age and
% sex via multiple regression

X_all=[Ak_NB Sek_NB];

for i=1:size(Tk_yeo,2)
    lm = fitlm(X_all,Tk_yeo(:,i), 'CategoricalVars',[2]);
    Tk_yeo_res(:,i)= lm.Residuals{:,1};   
end


% Derive group-wise, age- and sex-adjusted betas 

for i=1:size(Tk_yeo,2)
m1_res(i)   = mean(Tk_yeo_res(G1_NB_idx,i)); %FLE
m2_res(i)   = mean(Tk_yeo_res(G2_NB_idx,i)); %CTR
m3_res(i)   = mean(Tk_yeo_res(G3_NB_idx,i)); %TLE
end

for i=1:size(Tk_yeo,2)
m1a_res(i)      = mean(Tk_yeo_res(G1s_NB_idx,i));   %L_FLE
m1b_res(i)      = mean(Tk_yeo_res(G2s_NB_idx,i));   %R_FLE
m3a_res(i)      = mean(Tk_yeo_res(G5s_NB_idx,i));   %L_TLE
m3b_res(i)      = mean(Tk_yeo_res(G6s_NB_idx,i));   %R_TLE
end


% Derive z-scored betas (raw)
Tk_yeo_z         = (Tk_yeo - repmat(mean(Tk_yeo(G2_NB_idx,:),1),size(Tk_yeo,1),1)) ...
                    ./    repmat(std(Tk_yeo(G2_NB_idx,:),0,1),size(Tk_yeo,1),1);
    
%mean  (sanity check)              
for i=1:size(Tk_yeo,2)              
m1_z(i)                = mean(Tk_yeo_z(G1_NB_idx,i)); %FLE
m2_z(i)                = mean(Tk_yeo_z(G2_NB_idx,i)); %CTR
m3_z(i)                = mean(Tk_yeo_z(G3_NB_idx,i)); %TLE
end

%stdev (canity check)
for i=1:size(Tk_yeo,2)
std1_z(i)              = std(Tk_yeo_z(G1_NB_idx,i)); %FLE
std2_z(i)              = std(Tk_yeo_z(G2_NB_idx,i)); %CTR
std3_z(i)              = std(Tk_yeo_z(G3_NB_idx,i)); %TLE
end 

% Derive z-scored betas (on residuals)

Tk_yeo_z_res         = (Tk_yeo_res - repmat(mean(Tk_yeo_res(G2_NB_idx,:),1),size(Tk_yeo_res,1),1)) ...
                    ./    repmat(std(Tk_yeo_res(G2_NB_idx,:),0,1),size(Tk_yeo_res,1),1);
                
                
%mean (sanity check)                
for i=1:size(Tk_yeo_z_res,2)               
m1_z_res(i)                = mean(Tk_yeo_z_res(G1_NB_idx,i)); %FLE
m2_z_res(i)                = mean(Tk_yeo_z_res(G2_NB_idx,i)); %CTR
m3_z_res(i)                = mean(Tk_yeo_z_res(G3_NB_idx,i)); %TLE
end

%stdev (sanity check)

for i=1:size(Tk_yeo,2)
std1_z_res(i)              = std(Tk_yeo_z_res(G1_NB_idx,i)); %FLE
std2_z_res(i)              = std(Tk_yeo_z_res(G2_NB_idx,i)); %CTR
std3_z_res(i)              = std(Tk_yeo_z_res(G3_NB_idx,i)); %TLE
end 

%mean (sanity check)
for i=1:size(Tk_yeo_z_res,2)
m1a_z_res(i)               = mean(Tk_yeo_z_res(G1s_NB_idx,i)); %L_FLE
m1b_z_res(i)               = mean(Tk_yeo_z_res(G2s_NB_idx,i)); %R_FLE
m3a_z_res(i)               = mean(Tk_yeo_z_res(G5s_NB_idx,i)); %L_TLE
m3b_z_res(i)               = mean(Tk_yeo_z_res(G6s_NB_idx,i)); %R_TLE
end

%stdeev (sanity check)
for i=1:size(Tk_yeo_z_res,2)
std1a_z_res(i)              = std(Tk_yeo_z_res(G1s_NB_idx,i)); %L_FLE
std1b_z_res(i)              = std(Tk_yeo_z_res(G2s_NB_idx,i)); %R_FLE
std3a_z_res(i)              = std(Tk_yeo_z_res(G5s_NB_idx,i)); %L_TLE
std3b_z_res(i)              = std(Tk_yeo_z_res(G6s_NB_idx,i)); %R_TLE
end

%% 7. Display data and carry out stats -- NB, CTR

%%%% DISPLAY IN CTR ONLY -- NB
        
network.names_17      = {'DAN_A','DAN_B','FP_A','FP_B', 'FP_C','DMN_A', 'DMN_B', 'DMN_C','TP','Lim_A','Lim_B','Vis_A','Vis_B','SM_A','SM_B','SAL_A','SAL_B'};
network.names         = {'DAN','FP','DMN','LIM','VIS','SM','VAN'};

f = figure;

% Plot data in CTR

[ax,H]              = mica_spider([m2(1:end)'], [], repmat([-0.15 0.3],7,1), ...
                    network.names(1:end),[0 0 0],gca) ;
                    set(gcf,'color','w');
                    set(H.labels,'FontSize',22,'FontName', 'Gill Sans');keepR
                    set(H.tickMarks,'FontSize',19,'FontName', 'Gill Sans');
                    H.tickMarks(1).String = '-0.15';
                    
                                     
cd ([P '/1_Images/Yeo_06_2020']) 
print('Spider_NB_CTR', '-dpng', '-r600'); % Option for Yeo/Krienen 17 systems
%print('Spider_NB_CTR_17', '-dpng', '-r600'); 

% Stats -- one sample t-tests in controls

for i=1:(size(Tk_yeo,2))

[pval_perm, t_perm, ~, ~, ~]= mult_comp_perm_t1(Tk_yeo(G2_NB_idx,i),10000,0,0.05,0,1);
NB_stats_CTR_perm(1,i)= pval_perm;
NB_stats_CTR_perm(2,i)= t_perm; 

end

%FDR correction of p values across all systems

NB_stats_CTR_perm_FDR= mafdr(NB_stats_CTR_perm(1,:),'BHFDR', true);
NB_stats_CTR_perm(4,:)= NB_stats_CTR_perm_FDR;

save('NB_CTR_stats.mat', 'NB_stats_CTR_perm');
%save('NB_CTR_stats_17.mat', 'NB_stats_CTR_perm'); use the below when analyzing 17 systems



% Stats-- sanity check, display in all subjects (FLE, CTR, TLE - not used
% for display

[ax,H]              = mica_spider([mean(Tk_yeo(:,1:end))'], [], repmat([-0.15 0.3],7,1), ...
                    network.names(1:end),[0 0 0],gca) ;
                    set(gcf,'color','w');
                    set(H.labels,'FontSize',22,'FontName', 'Gill Sans');keepR
                    set(H.tickMarks,'FontSize',19,'FontName', 'Gill Sans');
                    H.tickMarks(1).String = '-0.15';
                    
                                     
cd ([P '/1_Images/Yeo_06_2020']) 
print('Spider_NB_CTR', '-dpng', '-r600');
%print('Spider_NB_CTR_17', '-dpng', '-r600'); use the below when analyzing 17 systems
 
for i=1:(size(Tk_yeo,2))

[pval_perm, t_perm, ~, ~, ~]= mult_comp_perm_t1(Tk_yeo(:,i),10000,0,0.05,0,1);
NB_stats_all_perm(1,i)= pval_perm;
NB_stats_all_perm(2,i)= t_perm; 

end

%FDR correction of p values across all systems
NB_stats_all_perm_FDR= mafdr(NB_stats_all_perm(1,:),'BHFDR', true);
NB_stats_all_perm(4,:)= NB_stats_all_perm_FDR;



%% 8. Group comparisons -- NB, FLE vs CTR
         
%%%% Spider data and Display adjusted z scores in 2 groups (FLE vs CTR)
       
network.names       = {'DAN','FP','DMN','LIM','VIS','SM','VAN'}; 
f = figure;
[ax,H]              = mica_spider_z(m1_z_res(1:end)', [], repmat([-1 0.5],7,1), ...
                    network.names(1:end),[0.7 0 0],gca);
                    set(gcf,'color','w');
                    % Change labels size & style
                    set(H.labels,'FontSize',22,'FontName', 'Gill Sans')
                    % Change tickMark size & style
                    set(H.tickMarks,'FontSize',19,'FontName', 'Gill Sans')
           
cd ([P '/1_Images/Yeo_06_2020'])          
print('Spider_NB_Z2_FLE_CTR', '-dpng', '-r600');
%print('Spider_NB17_Z2_FLE_CTR', '-dpng', '-r600'); %option for Yeo/Krienen
%17 systems
               

%Statistics -- one sample t test on z scores, age- and sex-adjusted data, FLE vs CTR     

for i=1:(size(Tk_yeo_z_res,2))

[pval_perm, t_perm, ~, ~, ~]= mult_comp_perm_t1(Tk_yeo_z_res(G1_NB_idx,i),10000,0,0.05,0,1);
NB_FLE_stats_perm(1,i)= pval_perm;
NB_FLE_stats_perm(2,i)= t_perm; 

end

%FDR correction of p values across all systems
NB_FLE_stats_perm_FDR= mafdr(NB_FLE_stats_perm(1,:),'BHFDR', true);
NB_FLE_stats_perm(4,:)= NB_FLE_stats_perm_FDR;

%Cohen's D in FLE

for i=1:(size(m1_z_res,2))
d_FLE_NB(i)= m1_z_res(i)/std1_z_res(i);
end

%%%% Statistics, code example for sensitivity analysis -- FLE-FCD vs CTR

% mean (z scores) in FLE FCD
for i=1:(size(Tk_yeo_z_res,2))
FLE_FCD(:,i) = Tk_yeo_z_res(FCD==1,i);
m_z_FLE_FCD(:,i)= mean(FLE_FCD(:,i));
end

% std in FLE FCD
for i=1:size(FLE_FCD,2)
std_FCD(:,i)    = std(FLE_FCD(:,i));
end

%  one sample t test on z scores, age- and sex-adjusted data, FLE-FCD vs CTR 

for i=1:(size(Tk_yeo_z_res,2))

[pval_perm, t_perm, ~, ~, ~]= mult_comp_perm_t1(FLE_FCD(:,i),10000,0,0.05,0,1);
NB_FLE_FCD_stats_perm(1,i)= pval_perm;
NB_FLE_FCD_stats_perm(2,i)= t_perm; 

end

%Cohen's d in FLE-FCD
for i=1:(size(m_z_FLE_FCD,2))
d_FLE_FCD_NB(i)= m_z_FLE_FCD(i)/std_FCD(i);
end

%FDR correction of p values across all systems
NB_FLE_FCD_stats_perm_FDR= mafdr(NB_FLE_FCD_stats_perm(1,:),'BHFDR', true);
NB_FLE_FCD_stats_perm(4,:)= NB_FLE_FCD_stats_perm_FDR;


% spider data - FLE-FCD

network.names       = {'DAN','FP','DMN','LIM','VIS','SM','VAN'}; 
f = figure;
[ax,H]              = mica_spider_z(m_z_FLE_FCD(1:end)', [], repmat([-1 0.5],7,1), ...
                    network.names(1:end),[0.7 0 0],gca);
                    set(gcf,'color','w');
                    % Change labels size & style
                    set(H.labels,'FontSize',22,'FontName', 'Gill Sans')
                    % Change tickMark size & style
                    set(H.tickMarks,'FontSize',19,'FontName', 'Gill Sans')
                    
cd ([P '/1_Images/Yeo_06_2020'])          
print('Spider_NB_Z2_FLE_FCD_CTR', '-dpng', '-r600');
save('NB_FCD.mat', 'm_z_FLE_FCD', 'std_FCD', 'NB_FLE_FCD_stats_perm');

%% 9. Group Comparisons, TLE vs CTR                    

%%%% Spider data and Display Z2 groups (TLE vs CTR)

f = figure;
[ax,H]              = mica_spider_z(m3_z_res(1:end)', [], repmat([-1 0.5],7,1), ...
                    network.names(1:end),[1 0.6 0.1; 0 0 0],gca ) ;
                    set(gcf,'color','w');
                    % Change labels size & style
                    set(H.labels,'FontSize',22,'FontName', 'Gill Sans')
                    % Change tickMark size & style
                    set(H.tickMarks,'FontSize',19,'FontName', 'Gill Sans')

cd ([P '/1_Images/Yeo_06_2020'])          
print('Spider_NB_Z2_TLE_CTR', '-dpng', '-r600');
%print('Spider_NB17_Z2_TLE_CTR', '-dpng', '-r600'); %option for Yeo/Krienen
%17 systems


% Statistics -- TLE vs CTR, one-sample t tests on z scores, age- and
% sex-adjusted data

for i=1:(size(Tk_yeo_z_res,2))

[pval_perm, t_perm, ~, ~, ~]= mult_comp_perm_t1(Tk_yeo_z_res(G3_NB_idx,i),10000,0,0.05,0,1);
NB_TLE_stats_perm(1,i)= pval_perm;
NB_TLE_stats_perm(2,i)= t_perm; 

end

% FDR correction of p values across all systems
NB_TLE_stats_perm_FDR= mafdr(NB_TLE_stats_perm(1,:),'BHFDR', true);
NB_TLE_stats_perm(4,:)= NB_TLE_stats_perm_FDR;

% Cohen's d in FLE

for i=1:(size(m3_z_res,2))
d_TLE_NB(i)= m3_z_res(i)/std3_z_res(i);
end

%% 10. Group Comparisons, FLE vs TLE

%%%% First, Display z scores for 3 groups - Z3, these will be the paper's figures (FLE vs TLE vs CTR)

f = figure;
[ax,H]              = mica_spider_z_17([m3_z_res(1:end)', m1_z_res(1:end)'], [], repmat([-1 0.5],17,1), ...
                    network.names_17(1:end),[1 0.6 0.1; 0.7 0 0],gca ) ;
                    set(gcf,'color','w');
                    % Change labels size & style
                    set(H.labels,'FontSize',22,'FontName', 'Gill Sans')
                    % Change tickMark size & style
                    set(H.tickMarks,'FontSize',19,'FontName', 'Gill Sans')
                                        
cd ([P '/1_Images/Yeo_06_2020'])    
print('Spider_NB_Z3_FLE_TLE_CTR', '-dpng', '-r600');
%print('Spider_NB_Z3_FLE_TLE_CTR_17', '-dpng', '-r600'); %Option for
%Yeo/Krienen 17 systems


%Comparison of FLE and TLE - pvals, two-sample t test, age- and sex-adjusted measures, not
%adjusted for side

for i=1:(size(Tk_yeo_res,2))

[pval_perm, t_perm, ~, ~, ~]= mult_comp_perm_t2(Tk_yeo_res(G1_NB_idx,i),Tk_yeo_res(G3_NB_idx,i),10000,0,0.05,0,'t',1);
NB_TLE_FLE_stats_perm(1,i)= pval_perm;
NB_TLE_FLE_stats_perm(2,i)= t_perm; 

end

%FDR correction of p values across all systems (not adjusted for side)
NB_TLE_FLE_stats_perm_FDR= mafdr(NB_TLE_FLE_stats_perm(1,:),'BHFDR', true);
NB_TLE_FLE_stats_perm(4,:)= NB_TLE_FLE_stats_perm_FDR;


% Linear regression in FLE & TLE – for comparison of FLE/TLE for age/sex/side

X_all_FT = [Ak_NB Sek_NB Sik_NB];
X_all_FT(G2_NB_idx(1:end),:)=[];
Tk_yeo_FT=Tk_yeo;
Tk_yeo_FT(G2_NB_idx(1:end),:)=[];

for i=1:size(Tk_yeo_FT,2);
    lm_FT_NB = fitlm(X_all_FT,Tk_yeo_FT(:,i), 'CategoricalVars',[2,3]);
    Tk_yeo_FT_res(:,i)= lm_FT_NB.Residuals{:,1};   
end


%Stats -- two-sample t test, p vals, comparison of FLE/TLE adjusted for
%age, sex and side

for i=1:(size(Tk_yeo_FT_res,2))

[pval_perm, t_perm, ~, ~, ~]= mult_comp_perm_t2(Tk_yeo_FT_res(G1_NB_idx(1:end),i),Tk_yeo_FT_res((G1_NB_idx(end)+1):size(Tk_yeo_FT_res,1),i),10000,0,0.05,0,'t',1);
NB_TLE_FLE_stats_perm_cov(1,i)= pval_perm;
NB_TLE_FLE_stats_perm_cov(2,i)= t_perm; 

end

%FDR correction of p values across all systems (adjusted for side)
NB_TLE_FLE_stats_perm_cov_FDR= mafdr(NB_TLE_FLE_stats_perm_cov(1,:),'BHFDR', true);
NB_TLE_FLE_stats_perm_cov(4,:)= NB_TLE_FLE_stats_perm_cov_FDR;

%Cohen's d FLE/TLE

for i=1:(size(m3_z_res,2))
    std_pooled(i) = sqrt(((std1_z_res(i)^2 + std3_z_res(i)^2))/2);
    d_FLE_TLE_NB(i)= (m1_z_res(i)-m3_z_res(i))/std_pooled(i);
end


save('NB_stats_17.mat', 'NB_FLE_stats_perm', 'NB_TLE_stats_perm', 'NB_TLE_FLE_stats_perm',...
     'NB_TLE_FLE_stats_perm_cov', 'd_FLE_NB', 'd_TLE_NB', 'd_FLE_TLE_NB');

%% 11. Correlations of fMRI measures across systems with task/cognitive test performance

%Correlations of NB fMRI data with Digit Span (DS) scores - permutation-based

%The below is to create separate group labels - after converting cell array, cutting
%'Group', and converting to vector, and prepar data subset for correlation
% analyses

Gf_NB=string(Gnsk_NB);
Gf_NB_new= erase(Gf_NB, "Group");
Gf_NB_newch= convertStringsToChars(Gf_NB_new);
Gf_NB_newch= cell2mat(Gf_NB_newch);
Gfa_NB= str2num(Gf_NB_newch);

NB_DS = horzcat(Tk_yeo_res, DIGk_NB,Gfa_NB);

NB_DS(any(isnan(NB_DS),2),:)= [];

% Get correlation values – across systems

R_comb_NB= zeros(2,(size(NB_DS,2)-2));

for i=1:size(Tk_yeo_res,2)
[pval_perm, R_perm, ~, ~, ~]= mult_comp_perm_corr_orig(NB_DS(:,i),NB_DS(:,8),10000,0,0.05,'linear',1); %change 'linear' to 'rank' if you want Spearman's rank correlations
R_comb_NB(1,i)= R_perm;
R_comb_NB(2,i)= pval_perm;
end

%FDR correction of p values across systems
R_comb_NB_FDR= mafdr(R_comb_NB(2,:),'BHFDR', true);
R_comb_NB(3,:)= R_comb_NB_FDR;

% Display correlation data - DS & NB, using spider plot

network.names       = {'DAN','FP','DMN','LIM','VIS','SM','VAN'}; 
f = figure;

[ax,H] = mica_spider_corr([R_comb_NB(1,(1:end))'], [], repmat([-0.2 0.4],7,1), ...
                    network.names(1:end),[0 0 0],gca);
                    set(gcf,'color','w');
                    set(H.labels,'FontSize',22,'FontName', 'Gill Sans');
                    set(H.tickMarks,'FontSize',19, 'FontName', 'Gill Sans');
                    H.tickMarks(1).String = '-0.2';
                    H.tickMarks(2).String = '0';
                    
cd ([P '/1_Images/Yeo_06_2020'])          
print('Spider_NB_correl_DS', '-dpng', '-r600');


% Scatter plot for correlation of NB fMRI and DS data -- frontoparietal
% (FP) control network

f=figure('Position', [10 10 1300 900]);
f.Color=[1 1 1];

hd =plot(NB_DS(:,8),NB_DS(:,2),'.','MarkerSize',85,'Color',[0.75 0 0.93]); % PLOT EVERYONE;
hold on

set(gca, 'FontSize',54,'FontName','Gill Sans');
hd=gca; 
hd.XLim=([0 20]); 
hd.XTick=0:5:30; 
hd.YLim=[-1 0.75];
hd.YTick=-1:0.5:2;

hd1 = lsline(hd);% generate the linear fit. 
hd1.LineWidth=10;
hd1.LineStyle='-.';
hd1.Color=[0.5 0.5 0.5];
hd1.AlignVertexCenters='off';

% Color scheme to distinguish groups in plot
% Colors scheme conversions -- FLE =red -- 253,0,0, gets 0.699 0.258 0.199; TLE=orange 243, 130, 81, gets 0.820 0.558 0.199; CTR green-ish 61 108 11 gets 0.238 0.304 0.434

G1_l = find(NB_DS(:,9)==1);
G2_l = find(NB_DS(:,9)==2);
G3_l = find(NB_DS(:,9)==3);

hdc=plot(NB_DS(G2_l,8),NB_DS(G2_l,2),'.','MarkerSize',85,'Color',[0.238 0.304 0.434]); % Group2 - CTR
hdf=plot(NB_DS(G1_l,8),NB_DS(G1_l,2),'.','MarkerSize',85, 'Color',[0.758 0 0]); % Group1 - FLE
hdt=plot(NB_DS(G3_l,8),NB_DS(G3_l,2),'.','MarkerSize',85,'Color',[0.949 0.508 0.316]); % Group3 - TLE

print ([P '/1_Images/Yeo_06_2020/' 'Scatter_NB_DS_FP_11_2020'], '-dpng', '-r600');

close all
clear f hd hd1 hdc hdf hdt


%%%% Correlation of NB fMRI data with NB (2 Back) task performance –permutation-based

% Prepare data subset

NB_NB = horzcat(Tk_yeo_res, TWO_BACK_NB);
NB_NB(any(isnan(NB_NB),2),:)= [];

% Compute permutation-based correlations

R_comb_NB_NB= zeros(2,(size(NB_NB,2)-1));
for i=1:size(Tk_yeo_res,2)
[pval_perm, R_perm, ~, ~, ~]= mult_comp_perm_corr_orig(NB_NB(:,i),NB_NB(:,8),10000,0,0.05,'rank',1);
R_comb_NB_NB(1,i)= R_perm;
R_comb_NB_NB(2,i)= pval_perm;
end

%FDR correction of p values across sytems

R_comb_NB_NB_FDR= mafdr(R_comb_NB_NB(2,:),'BHFDR', true);
R_comb_NB_NB(3,:)= R_comb_NB_NB_FDR;


% Display correlational data - NB fMRI and NB (2 Back) task performance

network.names       = {'DAN','FP','DMN','LIM','VIS','SM','VAN'}; 
f = figure;

[ax,H] = mica_spider_corr([R_comb_NB_NB(1,(1:end))'], [], repmat([-0.2 0.4],7,1), ...
                    network.names(1:end),[0.7 0 0],gca);
                    set(gcf,'color','w');
                    set(H.labels,'FontSize',22,'FontName', 'Gill Sans');
                    set(H.tickMarks,'FontSize',19, 'FontName', 'Gill Sans');
                    H.tickMarks(1).String = '-0.2';
                    H.tickMarks(2).String = '0';
                    
cd ([P '/1_Images/Yeo_06_2020'])          
print('Spider_NB_correl_DS', '-dpng', '-r600');

save('NB_correl_perf.mat','R_comb_NB','R_comb_NB_NB'); 

%% 12. Correlations of fMRI measures (systems-wise) with clinical variables

% Here -- example with age at onset and FLE (same to be done for duration of epilepsy, seizure
% frequency, FBTCS, time since last seizure)

% Prepare data subset - FLE (includes regressing out effects of sex and side)

NB_Ao_FL= horzcat(Tk_yeo(G1_NB_idx,:), OK_NB(G1_NB_idx), Sik_NB(G1_NB_idx), Sek_NB(G1_NB_idx));
NB_Ao_FL(any(isnan(NB_Ao_FL),2),:)= [];


for i=1:size(NB_Ao_FL,2)-3
    lm_NB_Ao = fitlm(NB_Ao_FL(:,9:10),NB_Ao_FL(:,i), 'CategoricalVars',[1,2]);
    NB_Ao_FL_res(:,i)= lm_NB_Ao.Residuals{:,1};   
end

% Compute correlations

R_comb_NB_Ao_FL= zeros(2,(size(NB_Ao_FL,2)-3));

for i=1:size(Tk_yeo_res,2)
[pval_perm, R_perm, ~, ~, ~]= mult_comp_perm_corr_orig(NB_Ao_FL_res(:,i),NB_Ao_FL(:,size(NB_Ao_FL,2)-2),10000,0,0.05,'linear',1); %change 'linear' to 'rank' if you want Spearman's rank correlations
R_comb_NB_Ao_FL(1,i)= R_perm;
R_comb_NB_Ao_FL(2,i)= pval_perm;
end

%FDR correction of p values across sytems

R_comb_NB_Ao_FL_FDR= mafdr(R_comb_NB_Ao_FL(2,:),'BHFDR', true);
R_comb_NB_Ao_FL(3,:)= R_comb_NB_Ao_FL_FDR;


% Display correlation data - age at onset and NB fMRI, FLE

network.names       = {'DAN','FP','DMN','LIM','VIS','SM','VAN'}; 
f = figure;

[ax,H] = mica_spider_corr([R_comb_NB_Ao_FL(1,(1:end))'], [], repmat([-0.2 0.4],7,1), ...
                    network.names(1:end),[0.758 0 0],gca);
                    set(gcf,'color','w');
                    set(H.labels,'FontSize',22,'FontName', 'Gill Sans');
                    set(H.tickMarks,'FontSize',19, 'FontName', 'Gill Sans');
                    H.tickMarks(1).String = '-0.2';
                    H.tickMarks(2).String = '0';
                    
cd ([P '/1_Images/Yeo_06_2020'])          
print('Spider_NB_correl_ao_FL', '-dpng', '-r600');

%Scatter plot - age at onset and fMRI activity across frontoparietal
%control network

f=figure('Position', [10 10 1300 900]);
f.Color=[1 1 1];

hd =plot(NB_Ao_FL(:,8),NB_Ao_FL_res(:,2),'.','MarkerSize',85,'Color',[0.758 0 0]); % PLOT EVERYONE;
hold on

set(gca, 'FontSize',54,'FontName','Gill Sans');
hd=gca; 
hd.XLim=([0 41]); 
hd.XTick=0:10:50; 
hd.YLim=[-1 1];
hd.YTick=-1:0.5:2;

hd1 = lsline(hd); % this will generate the linear fit
hd1.LineWidth=10;
hd1.LineStyle='-.';
hd1.Color=[0.5 0.5 0.5];
hd1.AlignVertexCenters='off';

cd ([P '/1_Images/Yeo_06_2020'])          
print ('Scatter_NB_FLE_Ao_11_2020_FP', '-dpng', '-r600');

save('R_NB_Ao.mat','R_comb_NB_Ao_FL');

%%%% End of code snippet %%%%