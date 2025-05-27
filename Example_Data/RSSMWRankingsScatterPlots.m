clear; 
clc
close all

fig=figure;
set(fig, 'Position', [1,126,1314,867]);
set(gcf, 'Color', 'w'); % Set background color to white

for jj=1:3 %change for the number of grains that are in the sample
% 
% Sample 1 - Grain 1 in powerpoint and grain 3 in slack green grain
% if jj==2
% v = [-0.406052; -0.053809; 0.129460];
% winners = [56,55]; %56 ranked 1st by MW and 2nd by RSS 55 ranked 10th by MW and 11th by RSS 
% end
% if jj==1
% % Sample 1 - Grain 2 in powerpoint, grain 1 in slack blue grain
% v = [0.341501; -0.178892; 0.059482];
% winners = 57; %ranked 10th by MW 11th by RSS
% end

% % Sample 2 - Grain 1 in powerpoint and slack green grain
    % if jj==1
    % v = [-0.076822; 0.389908; 0.090307];
    % winners = 18; %ranked 9th by MW and 11th by RSS
    % end
    % if jj==2
    % % Sample 2 - Grain 2 in powerpoint and slack yellow grain
    % v = [-0.033486; -0.149606; -0.204647];
    % winners = [54, 56]; % 54 ranked 9th in MW and 10th by RSS, 56 is ranked 10th by MW and 12th by RSS
    % end

% % S3   
    if jj==1
% Sample 3 - Grain 1 in powerpoint and slack light purple
v = [0.294993; 0.180003; 0.044470];
winners = 31; %ranked 5th by work and 6th by RSS
    end
    if jj==2
% Sample 3 - Grain 2 in powerpoint and slack pink grain
v = [0.059719; -0.102081; 0.121441];
winners = [58,60 ]; %58 is ranked 3rd by MW and 8th by RSS 60 is 4th by MW and 9th by RSS
    end
    if jj==3
% Sample 3 - Grain 3 in powerpoint and slack small purple grain
v = [-0.298081; 0.166393; 0.012648];
winners = 19;%ranked 10th by MW 11th by RSS
    end


R_crystal = Rodrot(norm(v), unit(v));
cd /Users/janicemoya/Library/CloudStorage/OneDrive-Personal/Research/Matlab/Example_Data/
CTM = readmatrix("Shield_CTM_CuAlNi_Results.csv");

HPVNum = CTM(:,1); 
b = CTM(:,4:6);
m = CTM(:,7:9);
type=CTM(:,10);
R_S2C = R_crystal;

F_sample = cell(1,96);
F_crystal = cell(1,96);    
schmid = cell(1,96);
schmid2 = cell(1,96);
transform_strain = cell(1,96);
work_transform=zeros(1,96);
tau_r=zeros(1,96);
sigma = [0 0 0; 0 0 0; 0 0 6/0.0005^2*10^(-6)]; % MPa % 6 N for Sample 2 and 3 but 10N for Sample 1
sigma_prime = R_S2C*sigma*R_S2C';

for ii = 1 : 96

    btmp = b(ii,:)';
    mtmp = m(ii,:)';
    btmpunit = (btmp) / (norm(b(ii,:)));

    schmid{ii} = btmp * transpose(mtmp);
    F_crystal{ii} = eye(3) + schmid{ii};
    F_sample{ii} = R_S2C' * F_crystal{ii} * R_S2C;

    transform_strain{ii} = 1/2 * (F_sample{ii}' * F_sample{ii} - eye(3));
    work_transform(ii) = trace(sigma * transform_strain{ii}');
    % work_transform(ii) = transform_strain{ii}(3,3);

    schmid2{ii} = mtmp * transpose(btmpunit);
    tau_r(ii) = trace(schmid2{ii} * sigma_prime');

end


% CPFEM results
% kk=jj+1;
% work = readmatrix('CPFEM_SimpleT.xlsx', 'Sheet', kk);
% work_transform=work(1:96,1);
% tau_r=readmatrix('CPFEM_02.xlsx', 'Sheet', kk);
% HPVNum=work(:,2);

[workValues, ind] = sort(work_transform, 'descend');
HPVNum_Ranked_MW1 = HPVNum(ind);
% 
[tau_rValues, ind2] = sort(tau_r, 'descend');
HPVNum_Ranked_RSS1 = HPVNum(ind2);

top10_MW1 = HPVNum_Ranked_MW1(1:10); 

table(top10_MW1,(1:10)')

rankingWin1 = find(ind==winners(1));
rankingWin1RSS = find(ind2==winners(1));
% F_sampleWin1 = F_sample{winners(1)};
% transform_strainWin1 = transform_strain{winners(1)};
% m_winner=m(winners(1),:);
% b_winner=b(winners(1),:);

if length(winners)==2
    rankingWin2 = find(ind==winners(2));
    rankingWin2RSS = find(ind2==winners(2));
    % F_sampleWin2 = F_sample{winners(2)};
    % transform_strainWin2 = transform_strain{winners(2)};
    % m_winner2=m(winners(2),:);
    % b_winner2=b(winners(2),:);

end

%%
hold on;
if jj==1
    h1=scatter(1:96, workValues, 80, 'r', 'filled', 'DisplayName', 'Grain 2'); % MW rankings
    h2=scatter(rankingWin1, workValues(rankingWin1), 350, 'r', 'filled', 'LineWidth',3,'MarkerEdgeColor','k','DisplayName', 'Grain 2 HPV'); % Highlight winner
end

if jj==2
    h4=scatter(1:96, workValues, 80, 'b', 'filled','DisplayName', 'Grain 1'); % MW rankings
    h5=scatter(rankingWin1, workValues(rankingWin1),  350, 'b', 'filled', 'LineWidth',3,'MarkerEdgeColor','k', 'DisplayName', 'Grain 1 Primary HPV'); % Highlight winner

    if exist('rankingWin2', 'var')
        h6=scatter(rankingWin2,workValues(rankingWin2) , 350, 'b', 'filled', 'LineWidth',3,'MarkerEdgeColor','k', 'DisplayName', 'Grain 1 Crossing HPV'); 
    end


end

if jj==3
    h7=scatter(1:96, workValues-5, 80, 'k', 'filled', 'DisplayName', 'Grain 3'); % MW rankings
    h8=scatter(rankingWin1, workValues(rankingWin1)-5, 350, 'k', 'filled','LineWidth',3,'MarkerEdgeColor','k', 'DisplayName', 'Grain 3 HPV'); % Highlight winner
end

xlabel('Maximum Work Ranking', 'FontSize', 30);
ylabel('Maximum Work (J)', 'FontSize', 30);
set(gca, 'FontSize', 28, 'Color', 'w'); % Set axis background color to white
xlim([0, 100])
xticks([0, 20,40,60,80,100])
ylim([-60, 60])


if jj==3
legend([h4,h1,h7], 'FontSize', 24);
end

if jj==2
legend([h4,h1], 'FontSize', 24);
end


%% RSS Scatter plot rankings
% if jj==1
%     h1=scatter(1:96, tau_rValues, 80, 'r', 'filled', 'DisplayName', 'Grain 2'); % MW rankings
%     h2=scatter(rankingWin1RSS, tau_rValues(rankingWin1RSS), 350, 'r', 'filled','LineWidth',3,'MarkerEdgeColor','k', 'DisplayName', 'Grain 2 HPV'); % Highlight winner
% 
% end
% 
% if jj==2
%     h4=scatter(1:96, tau_rValues, 80, 'b', 'filled','DisplayName', 'Grain 1'); % MW rankings
%     h5=scatter(rankingWin1RSS, tau_rValues(rankingWin1RSS), 350, 'b', 'filled','LineWidth',3,'MarkerEdgeColor','k', 'DisplayName', 'Grain 1 Primary HPV'); % Highlight winner
% 
%     if exist('rankingWin2', 'var')
%         h6=scatter(rankingWin2RSS,tau_rValues(rankingWin2RSS) , 350, 'b', 'filled','LineWidth',3,'MarkerEdgeColor','k', 'DisplayName', 'Grain 1 Crossing HPV'); 
%     end
% 
% 
% 
% end
% 
% if jj==3
%     h7=scatter(1:96, tau_rValues-0.5, 80, 'k', 'filled', 'DisplayName', 'Grain 3'); % MW rankings
%     h8=scatter(rankingWin1RSS, tau_rValues(rankingWin1RSS)-0.5, 350, 'k', 'filled','LineWidth',3,'MarkerEdgeColor','k', 'DisplayName', 'Grain 3 HPV'); % Highlight winner
% 
% 
% end
% % 
% xlabel('RSS Ranking', 'FontSize', 30);
% ylabel('RSS (MPa)', 'FontSize', 30);
% % title('Sample 1 RSS Rankings', 'FontSize', 36, 'FontWeight', 'bold');
% set(gca, 'FontSize', 28, 'Color', 'w'); % Set axis background color to white
% xlim([0, 100])
% xticks([0, 20,40,60,80,100])
% if jj==3
% legend([h4,h1,h7], 'FontSize', 24);
% end
% 
% if jj==2
% legend([h4,h1], 'FontSize', 24);
% end
% ylim([-22, 20])
% ylim([-14, 13.5])
% ylim([-1.5, 1.5])
% ylim([-2.5, 2.5])
end

