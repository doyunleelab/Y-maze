%%
% Supplementary code for
% Adverse effects of acute social isolation on spatial working memory: A quantitative analysis of spontaneous alternation behaviors on a Y-maze
% by Joowon Kim and Doyun Lee

%% Load data
load('C:\Users\supplementary_data.mat');

%% Calculate normalized positions of the head and the tailbase along arms A, B and C
head_abc = zeros(size(t,1),3);
for i_placein = 1:3
  placein = [];
  if (i_placein == 1)
    p0 = center;
    p1 = A;
  elseif (i_placein == 2)
    p0 = center;
    p1 = B;
  elseif (i_placein == 3)
      p0 = center;
      p1 = C;
  end

  for i_time = 1:size(head,1)

    q = head(i_time,1:2);
    a = [-q(1)*(p1(1)-p0(1)) - q(2)*(p1(2)-p0(2)); ...
        -p0(2)*(p1(1)-p0(1)) + p0(1)*(p1(2)-p0(2))];
    b = [p1(1) - p0(1), p1(2) - p0(2);...
        p0(2) - p1(2), p1(1) - p0(1)];
    ProjPoint = -(b\a)';

    p10 = p1 - p0;
    if (p10(1) > 0 & p10(2) > 0)
        if (ProjPoint(1) > p0(1) & ProjPoint(2) > p0(2))
            placein = [placein; sqrt(sum(   (ProjPoint - p0).^2   )) / sqrt(sum(   (p1 - p0).^2   ))];
        else
            placein = [placein; 0];
        end
    end
    if (p10(1) > 0 & p10(2) <= 0)
        if (ProjPoint(1) > p0(1) & ProjPoint(2) <= p0(2))
            placein = [placein; sqrt(sum(   (ProjPoint - p0).^2   )) / sqrt(sum(   (p1 - p0).^2   ))];
        else
            placein = [placein; 0];
        end
    end
    if (p10(1) <= 0 & p10(2) > 0)
        if (ProjPoint(1) <= p0(1) & ProjPoint(2) > p0(2))
            placein = [placein; sqrt(sum(   (ProjPoint - p0).^2   )) / sqrt(sum(   (p1 - p0).^2   ))];
        else
            placein = [placein; 0];
        end
    end
    if (p10(1) <= 0 & p10(2) <= 0)
        if (ProjPoint(1) <= p0(1) & ProjPoint(2) <= p0(2))
            placein = [placein; sqrt(sum(   (ProjPoint - p0).^2   )) / sqrt(sum(   (p1 - p0).^2   ))];
        else
            placein = [placein; 0];
        end
    end

  end %i_time

  head_abc(:, i_placein) = placein;
end %i_placein


tailbase_abc = zeros(size(t,1),3);
for i_placein = 1:3
  placein = [];
  if (i_placein == 1)
        p0 = center;
        p1 = A;
  end
  if (i_placein == 2)
        p0 = center;
        p1 = B;
  end
  if (i_placein == 3)
        p0 = center;
        p1 = C;
  end

  for i_time = 1:size(tailbase,1)

    q = tailbase(i_time,1:2);
    a = [-q(1)*(p1(1)-p0(1)) - q(2)*(p1(2)-p0(2)); ...
        -p0(2)*(p1(1)-p0(1)) + p0(1)*(p1(2)-p0(2))];
    b = [p1(1) - p0(1), p1(2) - p0(2);...
        p0(2) - p1(2), p1(1) - p0(1)];
    ProjPoint = -(b\a)';

    p10 = p1 - p0;

    if (p10(1) > 0 & p10(2) > 0)
      if (ProjPoint(1) > p0(1) & ProjPoint(2) > p0(2))
          placein = [placein; sqrt(sum(   (ProjPoint - p0).^2   )) / sqrt(sum(   (p1 - p0).^2   ))];
      else
          placein = [placein; 0];
      end
    end
    if (p10(1) > 0 & p10(2) <= 0)
        if (ProjPoint(1) > p0(1) & ProjPoint(2) <= p0(2))
            placein = [placein; sqrt(sum(   (ProjPoint - p0).^2   )) / sqrt(sum(   (p1 - p0).^2   ))];
        else
            placein = [placein; 0];
        end
    end
    if (p10(1) <= 0 & p10(2) > 0)
        if (ProjPoint(1) <= p0(1) & ProjPoint(2) > p0(2))
            placein = [placein; sqrt(sum(   (ProjPoint - p0).^2   )) / sqrt(sum(   (p1 - p0).^2   ))];
        else
            placein = [placein; 0];
        end
    end
    if (p10(1) <= 0 & p10(2) <= 0)
        if (ProjPoint(1) <= p0(1) & ProjPoint(2) <= p0(2))
            placein = [placein; sqrt(sum(   (ProjPoint - p0).^2   )) / sqrt(sum(   (p1 - p0).^2   ))];
        else
            placein = [placein; 0];
        end
    end

  end %i_time

  tailbase_abc(:, i_placein) = placein;
end %i_placein


%% Calculate correct rate along varying arm thresholds
k_range = [25:70];
N_frame = size(head_abc, 1);
arm_where = zeros(N_frame, 2, length(k_range));
behav_rate_k = zeros(1,length(k_range));

k_id = 0;
for k = k_range
  k_id = k_id + 1;
  center_thres = 0.25;
  arm_thres = 0.01 * k;

  % Position on the Y maze
  for fr = 1:N_frame
    for arm_idx = 1:3
      if fr == 1
        if head_abc(fr,arm_idx) > center_thres && tailbase_abc(fr,arm_idx) > center_thres
          arm_where(fr,1,k_id) = arm_idx;
          if head_abc(fr,arm_idx) > arm_thres && tailbase_abc(fr,arm_idx) > arm_thres
            arm_where(fr,2,k_id) = 1; % exceed arm thres
          end
        end
      elseif fr > 1
        if head_abc(fr,arm_idx) > center_thres && tailbase_abc(fr,arm_idx) > center_thres
          arm_where(fr,1,k_id) = arm_idx;
          if head_abc(fr,arm_idx) > arm_thres && tailbase_abc(fr,arm_idx) > arm_thres
            arm_where(fr,2,k_id) = 1; % exceed arm thres
          elseif (head_abc(fr,arm_idx) > arm_thres) ~= (tailbase_abc(fr,arm_idx) > arm_thres)
            if arm_where(fr-1,2,k_id) == 1 % if previously beyond arm thres
              arm_where(fr,2,k_id) = 1;
            end
          end
        elseif (head_abc(fr,arm_idx) > center_thres) ~= (tailbase_abc(fr,arm_idx) > center_thres)
          if arm_where(fr-1,1,k_id) == arm_idx
            arm_where(fr,1,k_id) = arm_idx;
          end
        end
      end
    end
  end

  % Generate trial_index [col1=what arm, col2=t_enter, col3=t_exit]
  diff_idx = diff(arm_where(:,1,k_id));
  t_enter_idx = find(diff_idx > 0) + 1;
  t_exit_idx = find(diff_idx < 0) + 1;

  if t_enter_idx(1) > t_exit_idx(1)
    t_enter_idx = [1; t_enter_idx];
  end

  if t_enter_idx(end) > t_exit_idx(end)
    t_exit_idx = [t_exit_idx; size(arm_where, 1)];
  end

  trial_index = [];
  for j = 1:size(t_enter_idx,1)
    if any(arm_where(t_enter_idx(j):t_exit_idx(j),2,k_id))
      trial_index = [trial_index; arm_where(t_enter_idx(j),1,k_id), t(t_enter_idx(j)), t(t_exit_idx(j))];
    end
  end
  trial_index_k{k_id} = trial_index;

  % Behavioral performance on the Y maze
  arm_seq = trial_index(:,1);
  behav_correct = zeros(length(arm_seq)-2, 1);
  for v = 1:length(arm_seq)-2
    if arm_seq(v+2) ~= arm_seq(v+1)
      if arm_seq(v+2) ~= arm_seq(v)
        if arm_seq(v+1) ~= arm_seq(v)
          behav_correct(v) = 1; % correct; mice made a triad
        else
          behav_correct(v) = 0; % neutral; mice had two choices
        end
      else
        behav_correct(v) = -1; % incorrect; mice visited the 2nd last arm again.
      end
    else
      behav_correct(v) = -2; % incorrect; mice visited the last arm again.
    end
  end

  N_triad = length(find(behav_correct == 1));
  behav_rate = length(find(behav_correct == 1))/length(find(behav_correct));

behav_rate_k(k_id) = behav_rate;

end

%% Plot correct rate along varying arm thresholds
figure('pos',[100,100,500,400]);
object = behav_rate_k;
xx = 0.25:0.01:0.7;
plot(xx,mean(object,1), 'b');
ylim([0.45, 0.75]);
box off
ylabel 'Correct rate'
xlabel 'Arm threshold'
