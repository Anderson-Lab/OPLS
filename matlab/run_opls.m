function [model,stats] = run_opls(X,Y,num_permutations,CV)

%%%%%
% start of cross validation
%%%%%
CV_array = {};
if length(CV) > 1 % Random bootstrap with replacement
    num_times = CV(1);
    min_num_samples = CV(2);
    max_num_samples = CV(3);
    num_samples_in_each_test_set = round(min_num_samples + (max_num_samples-min_num_samples).*rand(num_times,1));
    for i = 1:length(num_samples_in_each_test_set)
        num_samples_in_test_set = num_samples_in_each_test_set(i);
        CV_array{i} = round(1 + (length(Y)-1).*rand(num_samples_in_test_set,1));
    end
else
    fold = CV;
    if fold == -1
       fold = length(Y); 
    end
    
    j=1;
    while(1)
        for i=1:fold
            fold_array(j) = i;
            j=j+1;
            if j > length(Y)
                break;
            end
        end
        if j > length(Y)
            break;
        end
    end
    CV_array = {};
    for i = 1:fold
        CV_array{i} = find(fold_array == i);
    end
end


num_OPLS_fact = 0;
[Q2,Q2s,press,accuracy,AUC] = opls_CV(X,Y,num_OPLS_fact,CV_array);
while num_OPLS_fact < length(X(1,:))
    [next_Q2,next_Q2s,next_press,next_accuracy,next_AUC] = opls_CV(X,Y,num_OPLS_fact+1,CV_array);
    R = next_press/press;
    if R >= 1        
        break;
    else
        Q2 = next_Q2;
        Q2s = next_Q2s;
        AUC = next_AUC;
        press = next_press;
        accuracy = next_accuracy;
        num_OPLS_fact = num_OPLS_fact + 1;
    end
end

[model,stats] = opls(X,Y,num_OPLS_fact);
stats.Q2 = Q2;
stats.Q2s = Q2s;
stats.AUC = AUC;
stats.accuracy = accuracy;

[stats.permutation_Q2s,stats.permutation_AUCs] = permutation_test(X,Y,model.num_OPLS_fact,num_permutations,CV_array);
stats.Q2_pvalue = max([1/num_permutations,length(find(stats.permutation_Q2s >= Q2))/num_permutations]);
stats.AUC_pvalue = max([1/num_permutations,length(find(stats.permutation_AUCs >= AUC))/num_permutations]);
