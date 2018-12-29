% ----------------------------------------------------------------------- %
%                           R O C    C U R V E                            %
% ----------------------------------------------------------------------- %
% Function "roc_curve" calculates the Receiver Operating Characteristic   %
% curve, which represents the 1-specificity and sensitivity, of two classes
% of data, called class_1 and class_2.                                    %
%                                                                         %
%   Input parameters                                                      %
%       * class_1:  Column vector that represents the data of the first   %
%                   class.                                                %
%       * class_2:  Column vector that represents the data of the second  %
%                   class.                                                %
%       * dispp:    (Optional) If dispp is 1, the ROC Curve will be disp- %
%                   ayed inside the active figure. If dispp is 0, no figure
%                   will be displayed.                                    %
%       * dispt:    (Optional) If dispt is 1, the optimum threshold para- %
%                   meters obtained will be displayed on the MATLAB log.  %
%                   Otherwise, if dispt is 0, no parameters will be disp- %
%                   ayed there.                                           %
%                                                                         %
%   Output variables                                                      %
%       * ROC_data: Struct that contains all the curve parameters.        %
%           - param:    Struct that contains the cuantitative parameters  %
%                       of the obtained curve, which are:                 %
%               + Threshold:Optimum threshold calculated in order to maxi-%
%                           mice the sensitivity and specificity values,  %
%                           which is colocated in the nearest point to    %
%                           (0,1).                                        %
%               + AROC:     Area under ROC curve.                         %
%               + Accuracy: Accuracy.                                     %
%               + Sensi:    Sensitivity (i.e., recall, hit rate, or true  %
%                           positive rate).                               %
%               + Speci:    Specificity (i.e., selectivity, or true       %
%                           negative rate).                               %
%               + PPV:      Positive predicted value (i.e., precision).   %
%               + NPV:      Negative predicted value.                     %
%               + FNR:      False negative rate (i.e., miss rate).        %
%               + FPR:      False positive rate (i.e., fall-out).         %
%               + FDR:      False discovery rate.                         %
%               + FOR:      False omission rate.                          %
%               + F1_score: F1 score (harmonic mean of precision and      %
%                           sensitivity).                                 %
%               + MCC:      Matthews correlation coefficient.             %
%               + BM:       Bookmaker informedness.                       %
%               + MK:       Markedness.                                   %
%           - curve:    Matrix that contains the specificity and specifi- %
%                       city of each threshold point in columns.          %
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       class_1 = 0.5*randn(100,1);                                       %
%       class_2 = 0.5+0.5*randn(100,1);                                   %
%       roc_curve(class_1, class_2);                                      %
% ----------------------------------------------------------------------- %
%   Version log:                                                          %
%       - 1.0: Original script (02/09/2015).                              %
%       - 2.0: Different sizes of data_1 and data_2 are allowed as long as%
%              one column is filled with NaNs (26/06/2018).               %
%       - 2.1: Classes are now indicated separately (10/07/2018).         %
%       - 3.0: More parameters are indicated as an output (13/12/2018).   %
%       - 3.1: All parameters are now displayed in CMD (13/12/2018).      %
% ----------------------------------------------------------------------- %
%   Author:  Víctor Martínez Cagigal                                      %
%   Date:    13/12/2018                                                   %
%   E-mail:  vicmarcag (dot) gmail (dot) com                              %
%                                                                         %
%   Sources: Fawcett (2006), Powers (2011), and Ting (2011)               %
% ----------------------------------------------------------------------- %
function ROC_data = roc_curve(class_1, class_2, dispp, dispt)

    % Setting default parameters and detecting errors
    if(nargin<4), dispt = 1;    end
    if(nargin<3), dispp = 1;    end
    if(nargin<2), error('Params "class_1" or "class_2" are not indicated.'); end
    class_1 = class_1(:);
    class_2 = class_2(:);
    
    % Calculating the threshold values between the data points
    s_data = unique(sort([class_1; class_2]));          % Sorted data points
    s_data(isnan(s_data)) = [];                 % Delete NaN values
    d_data = diff(s_data);                      % Difference between consecutive points
    if(isempty(d_data)), error('Both class data are the same!'); end
    d_data(length(d_data)+1,1) = d_data(length(d_data));% Last point
    thres(1,1) = s_data(1) - d_data(1);                 % First point
    thres(2:length(s_data)+1,1) = s_data + d_data./2;   % Threshold values
        
    % Calculating the sensibility and specificity of each threshold
    curve = zeros(size(thres,1),2);
    distance = zeros(size(thres,1),1);
    for id_t = 1:1:length(thres)
        TP = length(find(class_2 >= thres(id_t)));    % True positives
        FP = length(find(class_1 >= thres(id_t)));    % False positives
        FN = length(find(class_2 < thres(id_t)));     % False negatives
        TN = length(find(class_1 < thres(id_t)));     % True negatives
        
        curve(id_t,1) = TP/(TP + FN);   % Sensitivity
        curve(id_t,2) = TN/(TN + FP);	% Specificity
        
        % Distance between each point and the optimum point (0,1)
        distance(id_t)= sqrt((1-curve(id_t,1))^2+(curve(id_t,2)-1)^2);
    end
    
    % Optimum threshold and parameters
    [~, opt] = min(distance);
    TP = length(find(class_2 >= thres(opt)));    % No. true positives
    FP = length(find(class_1 >= thres(opt)));    % No. false positives 
    FN = length(find(class_2 < thres(opt)));     % No. false negatives                                 
    TN = length(find(class_1 < thres(opt)));     % No. true negatives       
    
    % Output parameters
    param.Threshold = thres(opt);               % Optimum threshold position
    param.Sensi = curve(opt,1);                 % Sensitivity
    param.Speci = curve(opt,2);                 % Specificity
    param.AROC  = abs(trapz(1-curve(:,2),curve(:,1))); % Area under curve
    param.Accuracy = (TP+TN)/(TP+TN+FP+FN);     % Aaccuracy
    param.PPV   = TP/(TP+FP);                   % Positive predictive value
    param.NPV   = TN/(TN+FN);                   % Negative predictive value
    param.FNR   = FN/(FN+TP);                   % False negative rate
    param.FPR   = FP/(FP+TN);                   % False positive rate
    param.FDR   = FP/(FP+TP);                   % False discovery rate
    param.FOR   = FN/(FN+TN);                   % False omission rate
    param.F1_score = 2*TP/(2*TP+FP+FN);         % F1 score
    param.MCC   = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));  % Matthews correlation coefficient
    param.BM    = param.Sensi+param.Speci-1;    % Informedness
    param.MK    = param.PPV+param.NPV-1;        % Markedness
    
    param.TP = TP;    % No. true positives
    param.FP = FP;    % No. false positives 
    param.FN = FN;    % No. false negatives                                 
    param.TN = TN;    % No. true negatives  
    
    % Plotting if required
    if(dispp == 1)
        fill_color = [11/255, 208/255, 217/255];
        fill([1-curve(:,2); 1], [curve(:,1); 0], fill_color,'FaceAlpha',0.5);
        hold on; plot(1-curve(:,2), curve(:,1), '-b', 'LineWidth', 2);
        hold on; plot(1-curve(opt,2), curve(opt,1), 'or', 'MarkerSize', 10);
        hold on; plot(1-curve(opt,2), curve(opt,1), 'xr', 'MarkerSize', 12);
        hold off; axis square; grid on; xlabel('1 - specificity'); ylabel('sensibility');
        title(['AROC = ' num2str(param.AROC)]);
    end
    
    % AROC warning
    if param.AROC < 0.5
        warning('Since AROC is less than 0.5, you should swap the classes: roc_curve(class_2,class_1).');
    end
    
    % Log screen parameters if required
    if(dispt == 1)
        fprintf('\n ROC CURVE PARAMETERS\n');
        fprintf(' ------------------------------\n');
        fprintf('  - Distance:     %.4f\n', distance(opt));
        fprintf('  - Threshold:    %.4f\n', param.Threshold);
        fprintf('  - Sensitivity:  %.4f\n', param.Sensi);
        fprintf('  - Specificity:  %.4f\n', param.Speci);
        fprintf('  - AROC:         %.4f\n', param.AROC);
        fprintf('  - Accuracy:     %.4f\n', param.Accuracy);
        fprintf('  - PPV:          %.4f\n', param.PPV);
        fprintf('  - NPV:          %.4f\n', param.NPV);
        fprintf('  - FNR:          %.4f\n', param.FNR);
        fprintf('  - FPR:          %.4f\n', param.FPR);
        fprintf('  - FDR:          %.4f\n', param.FDR);
        fprintf('  - FOR:          %.4f\n', param.FOR);
        fprintf('  - F1 score:     %.4f\n', param.F1_score);
        fprintf('  - MCC:          %.4f\n', param.MCC);
        fprintf('  - BM:           %.4f\n', param.BM);
        fprintf('  - MK:           %.4f\n', param.MK);
        fprintf(' \n');
    end
    
    % Assinging parameters and curve data
    ROC_data.param = param;
    ROC_data.curve = curve;
end