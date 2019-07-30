function [validity, invert] = train_model(data_filename)
    global epsilon_auc;
    global epsilon_mu;

    invert = 0;
    
    ret_vals = py.ml_svm.train_wrapper(data_filename);  % TODO 11/1 train the model

    auc=ret_vals(1);   % The area under ROC curve for the test data.
    auc = auc{1};
    mu=ret_vals(2);
    mu = mu{1};
    fprintf('Returned AUC for the training data: %.4f.\n', auc);
    fprintf('Returned misclassification error for the training data: %.4f (cutoff is %.4f).\n', mu, epsilon_mu);
    
%     if isnan(auc)
%         validity = 0;
%     else
    
%         if (auc < 0.5)
% Note if AUC < 0.5, it normally means that the test data is not from the
% same distribution as the training data.
%             auc = 1 - auc;
%             invert = 1;
%         end

        if (mu <= epsilon_mu) 
            validity = 1;
        else
            validity = 0;
        end 
%     end
end 