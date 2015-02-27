function [t,t_ortho,Y_pred] = apply_opls_model(X,Y,model,new_X)
t_ortho = [];

Xres = bsxfun(@minus,new_X,mean(X));
for iter = 1:model.num_OPLS_fact
    %run OSC filter on Xres
    t_ortho(:,iter) = Xres*model.w_ortho(:,iter) / (model.w_ortho(:,iter)'*model.w_ortho(:,iter));
    Xres = Xres - t_ortho(:,iter)*model.p_ortho(:,iter)';
end

t = Xres*model.w/(model.w'*model.w);

% Original space
Wstar=model.w*inv(model.p'*model.w);
B_pls=Wstar*diag(model.b_l)*model.c';
Xres = bsxfun(@minus,new_X,mean(X));
z=Xres;
% filter out OPLS components
for filter=1:model.num_OPLS_fact
  z = (z - (z*model.w_ortho(:,filter)/(model.w_ortho(:,filter)'*model.w_ortho(:,filter)))*model.p_ortho(:,filter)');
end
Y_pred = z*B_pls + mean(Y);

