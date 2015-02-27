function [Q2,Q2s,press,accuracy,AUC] = opls_CV(X,Y,num_OPLS_fact,CV_array)

press = 0;
y_sum = 0;
errors = 0;
Q2s = NaN*ones(1,length(CV_array));
Ys_predicted = [];
Ys_actual = [];
for CV_count=1:length(CV_array)
    %%%%
    % set up CV data
    %%%%        
    mask = ones(size(Y)); % Start with all
    mask(CV_array{CV_count}) = 0;
    inxs = find(mask == 1);
    Xtemp = X(inxs,:);
    Ytemp = Y(inxs,:);
    %zscore this CV data set
    m=mean(Xtemp);
    m_Y=mean(Ytemp);
    Xres = bsxfun(@minus,Xtemp, m);
    Yres = Ytemp - m_Y;
    
    [model,stats] = opls(Xres,Yres,num_OPLS_fact);
    
    %%%%%%
    % calc partial press
    %%%%%%
    X_leftOut = X(CV_array{CV_count},:);
    X_leftOut = bsxfun(@minus,X_leftOut, m);
    Y_leftOut = Y(CV_array{CV_count},:)-m_Y;
    temp_press = 0;
    temp_y_sum = 0;
    for cpp = 1:length(Y_leftOut)
        Wstar=model.w*inv(model.p'*model.w);
        B_pls=Wstar*diag(model.b_l)*model.c';
        z=(X_leftOut(cpp,:));
        % filter out OPLS components
        for filter=1:num_OPLS_fact
            z = (z - (z*model.w_ortho(:,filter)/(model.w_ortho(:,filter)'*model.w_ortho(:,filter)))*model.p_ortho(:,filter)');
        end
        %predict
        Y_pred = z*B_pls;
        Ys_predicted(end+1) = Y_pred;
        Ys_actual(end+1) = Y_leftOut(cpp);
        temp_press = temp_press + (Y_pred - Y_leftOut(cpp))^2;
        temp_y_sum = temp_y_sum + (Y_leftOut(cpp))^2;
        correct_Y = Y_pred - (Y_leftOut(cpp));
        for k=1:length(Y)
            if (abs(Y_pred - (Y(k))) < abs(correct_Y))
                errors = errors+1;
                break;
            end
        end
    end
    Q2s(CV_count) = 1 - temp_press/temp_y_sum;
    press = press + temp_press;
    y_sum = y_sum + temp_y_sum;
end

Q2 = 1 - press/y_sum;
accuracy = (length(Y)-errors) / length(Y);

[X,Y,T,AUC] = perfcurve(Ys_actual,Ys_predicted,max(Ys_actual));