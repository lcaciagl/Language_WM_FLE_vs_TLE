function f = FLE_TLE_analysis_code_DisplayGradientData_CTR_NB(Data, Group_ns, minD, maxD)
    % Data  - num_sub * num_bins 
    % Group - 1,2,3 num_sub 
    
if ~iscell(Group_ns)==1
Group_ns = var2fac(Group_ns,'Group_ns');
end
Group_ns;
G1ns=strcmp(Group_ns,'Group1'); %FLE
G2ns=strcmp(Group_ns,'Group2'); %CTR
G3ns=strcmp(Group_ns,'Group3'); %TLE
        
n =      size(Data,2);

stdev_CTR   = std(Data(G2ns,:));
SEM_TS      = std(Data(G2ns,:)) / sqrt(length(Data(G2ns,:)));
CI95        = tinv([0.025 0.975], (length(Data(G2ns,:))-1));
yCI95       = bsxfun(@times, SEM_TS, CI95(:)); 

%% Plot CTR data

P = '/home/lorenzo/Wellcome//Analysis_MNI/matlab/';
load ([P '/fsaverage5/margulies_gradient.mat'])

f = figure; 
% declare figure (prev 600 450)
%('Position', [10 10 800 600])
hold on;

yG2ns= mean(Data(G2ns,:));
xG2ns=(1:size(yG2ns,2));
cG2ns=(1:size(yG2ns,2));

H.mainLine = color_line(xG2ns, yG2ns, cG2ns); 
colormap(parula)
H.mainLine.LineWidth = 4.0;
shading interp

uE = yG2ns + yCI95(2,:); % if one wants to plot SEM (or std) -- this should be done uE = y + SEM_TS; 
lE = yG2ns + yCI95(1,:); % if one wants to plot SEM (or std) -- this should be done lE = y - SEM_TS;

% uE = yG2ns + SEM_TS; % if one wants to plot SEM (or std) -- this should be done uE = y + SEM_TS; 
% lE = yG2ns - SEM_TS; % if one wants to plot SEM (or std) -- this should be done lE = y - SEM_TS;

yP = [lE,fliplr(uE)];
xP = [xG2ns,fliplr(xG2ns)];
cP = [cG2ns,fliplr(cG2ns)];

%Shaded area
H.patch = patch(xP,yP,cP);
H.patch.FaceAlpha = 0.2; %Modify to change the transparency
H.patch.EdgeColor = 'none';

%Edges (comment if you don't want them)
H.edge(1) = color_line(xG2ns,lE,cG2ns);
H.edge(2) = color_line(xG2ns,uE,cG2ns);
H.edge(1).LineWidth = 0.05; %Modify to change the tickness of the edge
H.edge(2).LineWidth = 0.05; %Modify to change the tickness of the edge

%%%% Color Bar Addition 

H.cb = colorbar;
H.cb.LineWidth= 0.1;
H.cb.Ticks = '';
H.cb.Location = 'south';
H.cb.AxisLocation = 'out';
H.cb.Position = [0.1311 0.1098 0.7737 0.0431];

uistack(H.mainLine,'top')

%%%% Axes settings
ax = gca;
ax.FontName = 'Gill Sans';
ax.FontSize = 16;
ax.XLim = [1 20]; 
ax.YLim = ([minD maxD]);
ax.XTick = 2:2:20; 

hold on
     
%% Statistics to add stars             

load ([P '/1_Images/Grad_con_06_2020/NB_gradient_CTR_stats.mat'])

star1 = find(NB_stats_CTR_perm(4,:)<0.05);
star2 = find(NB_stats_CTR_perm(1,:)<0.05);

uncorr_only = setxor(star1, star2);


for k=1:numel(star1(1,:))
        
        if isempty(star1(k))
           plot(nan,nan)
        else
           plot(star1(k),0.7, '*', 'Color', [0 0 0], 'MarkerSize', 11, 'MarkerEdgeColor',[0 0 0], 'LineWidth',2);
        end
    hold on
    end
    
    
   for k=1:numel(uncorr_only(1,:))
        if isempty(uncorr_only(k))
           plot(nan,nan)
        else
           plot(uncorr_only(k), 0.7, 'v', 'Color', [0 0 0], 'MarkerSize', 9, 'MarkerEdgeColor',[0 0 0], 'LineWidth',2);
        end
    hold on
   end
        
