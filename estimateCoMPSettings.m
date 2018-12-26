function [validity, inverted, model] = estimateCoMPSettings(simulation_data, isStatic, is_dnn)
    validity = false;
    inverted = false;
    model = [];
    
    global epsilon;
        
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
         
        if (~is_dnn)
            [err,model] = SVM(TBSINR,RSRP,BLER,false);
            
            if (err > epsilon)
                validity = false;
            else
                validity = true;
            end
                
        end
        
        if (is_dnn)
            % Create a structure
            data.RSRP = RSRP;
            data.BLER = BLER;
            data.CQI = CQI;
            data.TBSINR = TBSINR;

            [err,model] = SVM(TBSINR,RSRP,BLER,true);
            
            if (err > epsilon)
                validity = false;
            else
                validity = true;
            end
                
        end
    end 