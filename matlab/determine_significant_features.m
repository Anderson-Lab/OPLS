function [sig_inxs,not_sig_inxs,significant,p_permuted,pvalues] = determine_significant_features(X,Y,orig_model,num_permutations,talpha,inside_num_permutations,inside_talpha)
[num_samples,num_variables] = size(X);
p_permuted = cell(num_variables,1);
% if ~exist('variables')
%     variables = 1:num_variables;
% end

%% Determine the significant features for a model
P_original = orig_model.p; % Grab the original P
% for each variable, permute the labels N times and recalculate P
max_num_permutations = factorial(num_samples);
fprintf('Maximum number of permutations is %d\n',max_num_permutations);
% fprintf('The number of permutations is %d\n',N);
% fprintf('The test alpha is %f\n',talpha);

%% Inner permutation test
variables = 1:num_variables;
N = min([inside_num_permutations,max_num_permutations]);
inner_P_permuted = run(X,Y,orig_model,N,variables);

% Determine the outright not significant bins
significant = ones(1,num_variables)*NaN;
sig_inxs = [];
not_sig_inxs = [];
pvalues = ones(1,num_variables)*NaN;
for v = 1:num_variables
    sorted = sort(inner_P_permuted(v,:),'descend');
    inside_ix = max([1,round(N*inside_talpha/2)]); % Two tailed
    inside_thres1 = sorted(inside_ix);
    inside_thres2 = sorted(end-inside_ix+1);
    % Needed for p-value
    ix = max([1,round(N*talpha/2)]); % Two tailed
    thres1 = sorted(ix);
    thres2 = sorted(end-ix+1);    
    if inside_thres1 >= P_original(v) && P_original(v) >= inside_thres2 % Not close to the boundary
        p_permuted{v} = inner_P_permuted(v,:);
        not_sig_inxs(end+1) = v;
        significant(v) = false;
        pL = length(find(sorted <= min(P_original(v),-P_original(v))))/N;
        pR = length(find(sorted >= max(P_original(v),-P_original(v))))/N;
        pvalues(v) = pL + pR;
    end
end

%% Border permutation test
variables = find(isnan(significant));
N = min([num_permutations,max_num_permutations]);
border_P_permuted = run(X,Y,orig_model,N,variables);

% Determine the outright (not) significant bins
sig_inxs = [];
not_sig_inxs = [];
for v = variables
    sorted = sort(border_P_permuted(v,:),'descend');
    ix = max([1,round(N*talpha/2)]); % Two tailed
    thres1 = sorted(ix);
    thres2 = sorted(end-ix+1);    
    if P_original(v) >= thres1
        p_permuted{v} = border_P_permuted(v,:);
        significant(v) = true;
        sig_inxs(end+1) = v;
    elseif P_original(v) <= thres2
        p_permuted{v} = border_P_permuted(v,:);
        significant(v) = true;
        sig_inxs(end+1) = v;
    else
        p_permuted{v} = border_P_permuted(v,:);
        not_sig_inxs(end+1) = v;
        significant(v) = false;
    end
    pL = length(find(sorted <= min(P_original(v),-P_original(v))))/N;
    pR = length(find(sorted >= max(P_original(v),-P_original(v))))/N;
    pvalues(v) = pL + pR;
end

% p_permuted = P_permuted;

function P_permuted = run(X,Y,orig_model,N,variables)
[num_samples,num_variables] = size(X);
P_original = orig_model.p; % Grab the original P

num_OPLS_fact = orig_model.num_OPLS_fact;
P_permuted = NaN*ones(N,num_variables);
for v = variables
    v_P_permuted = NaN*ones(N,1);
    parfor n = 1:N
        X_permuted = X;
        inxs = randperm(num_samples);
        X_permuted(:,v) = X(inxs,v);
        [model,stats] = opls(X_permuted,Y,num_OPLS_fact);
        v_P_permuted(n) = model.p(v);
        % Make sure the direction of the vector is the same
        model.p(v) = 0; % Remove from calculation
        P_test = P_original;
        P_test(v) = 0;
        err1 = sum((P_test - model.p).^2);
        err2 = sum((P_test - (-1*model.p)).^2);
        if err2 < err1
            v_P_permuted(n) = -1 * v_P_permuted(n);
        end
    end
    P_permuted(:,v) = v_P_permuted;
    %fprintf('Finished %d\n',v);
end
P_permuted = P_permuted';
