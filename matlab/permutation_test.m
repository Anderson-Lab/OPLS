function [permutation_Q2s,permutation_AUCs] = permutation_test(X,Y,num_OPLS_fact,num_permutations,CV_array)
%permutation
permutation_Q2s = zeros(1,num_permutations);
permutation_AUCs = zeros(1,num_permutations);
for i=1:num_permutations
    %permute labels
    Y_rand_idx = randperm(length(Y));
    Y_rand = Y(Y_rand_idx);

    %run OPLS on permuted data    
    [permutation_Q2s(i),Q2s,press,accuracy,permutation_AUCs(i)] = opls_CV(X,Y_rand,num_OPLS_fact,CV_array);
end
