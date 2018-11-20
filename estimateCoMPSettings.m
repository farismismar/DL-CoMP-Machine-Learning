function [validity, inverted, model] = estimateCoMPSettings(simulation_data, isStatic)
    validity = false;
    inverted = false;
    model = [];
    
    global epsilon;
    global model_choice;
    
    if ~isStatic
        CQI = [];         
        RSRP = [];
        BLER = [];
        TBSINR = [];
         for i=1:length(simulation_data) 
            CQI = [CQI; simulation_data.TBCQI];           
            RSRP = [RSRP; simulation_data.RX_Power_TB];
            BLER = [BLER; simulation_data.BLER];
            TBSINR = [TBSINR; simulation_data.TBSINR];
         end

         % Vectorize and
         % Prepare data and clean it up.
         BLER = BLER(:);
         CQI = CQI(:);       
         RSRP = 10*log10(RSRP(:));
         TBSINR = TBSINR(:);
         
         % Find rows with infinity
        inclusion = union(find(RSRP > -Inf), ...
            find(TBSINR ~= 0));  % Get rid of -inf dBm and 0 SINR which is invalid.

        RSRP = RSRP(inclusion);
        BLER = BLER(inclusion);
        CQI = CQI(inclusion);
        TBSINR = TBSINR(inclusion);
         
        if (model_choice == 'svm')
            [err,model] = SVM(TBSINR,RSRP,BLER);
            
            if (err > epsilon)
                validity = false;
            else
                validity = true;
            end
                
        end
        
        if (model_choice == 'dnn')
            % Create a structure
            data.RSRP = RSRP;
            data.BLER = BLER;
            data.CQI = CQI;
            data.TBSINR = TBSINR;

            filename = 'measurements.mat';
            save(filename, '-struct', 'data')
            ret_vals = py.dnn.train_wrapper(filename);  % TODO 11/1 train the model

            auc=ret_vals(3);   % The area under ROC curve based on a split test data from the data collected.

            auc_roc = auc{1};
            %fprintf('Test AUC ROC =  %.3f.\n', double(auc_roc));
            
            if (auc_roc > 0)
                % Inverted decisions.
                if (auc_roc < 0.5) 
                    auc_roc = 1 - auc_roc;
                    inverted = true;
                    fprintf('Warning: Test AUC ROC < 0.5 is complemented and the decision is inverted.\n');
                end

                if (auc_roc < epsilon) 
                    validity = false;
                else
                    validity = true;
                end
            else
                fprintf('Warning: Model is not valid due to reported error.\n');
                validity = false;
            end
        end
    end 