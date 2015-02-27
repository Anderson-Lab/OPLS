function [model,stats] = opls(X,Y,num_OPLS_fact)
%%%%%%%%%%%%
% OPLS on full data
%%%%%%%%%%%%
w_ortho = [];
t_ortho = [];
p_ortho = [];

Xres = bsxfun(@minus,X, mean(X));
Yres = (Y)-mean(Y);
SS_Y=sum(sum(Yres.^2));
SS_X=sum(sum(Xres.^2));

for iter=1:num_OPLS_fact
    %find PLS component
    w = (Yres'*Xres / (Yres'*Yres))';
    w = w / norm(w);
    t = Xres*w / (w'*w);
    p = (t'*Xres / (t'*t))';

    %run OSC filter on Xres
    w_ortho(:,iter) = p - (w'*p / (w'*w)) * w;
    w_ortho(:,iter) = w_ortho(:,iter) / norm(w_ortho(:,iter));
    t_ortho(:,iter) = Xres*w_ortho(:,iter) / (w_ortho(:,iter)'*w_ortho(:,iter));
    p_ortho(:,iter) = (t_ortho(:,iter)'*Xres / (t_ortho(:,iter)'*t_ortho(:,iter)))';
    Xres = Xres - t_ortho(:,iter)*p_ortho(:,iter)';
end;

%%%%%%%%%%
% PLS on full data
%%%%%%%%%%
%find PLS component
w = (Yres'*Xres / (Yres'*Yres))';
w = w / norm(w);
t = Xres*w / (w'*w);
c = (t'*Yres / (t'*t))';
u = Yres*c / (c'*c);
p = (t'*Xres / (t'*t))';
% b coef
b_l=((t'*t)^(-1))*(u'*t);

%save model params
model = {};
model.b=b_l;
model.b_l=b_l;
model.c=c;
model.p=p;
model.w=w;
model.t=t;
model.u = u;
model.t_ortho = t_ortho;
model.p_ortho = p_ortho;
model.w_ortho = w_ortho;
model.num_OPLS_fact = num_OPLS_fact;
% Original space
Wstar=w*inv(p'*w);
B_pls=Wstar*diag(b_l)*c';
m=mean(X);
Xres = bsxfun(@minus,X, m);
z=Xres;
% filter out OPLS components
for filter=1:num_OPLS_fact
  z = (z - (z*w_ortho(:,filter)/(w_ortho(:,filter)'*w_ortho(:,filter)))*p_ortho(:,filter)');
end
%predict
model.Y_pred = z*B_pls + mean(Y);

stats = {};
stats.R2_X=(model.t'*model.t)*(model.p'*model.p)./SS_X;
stats.R2_Y=(model.t'*model.t)*(model.b.^2)*(model.c'*model.c)./SS_Y;
