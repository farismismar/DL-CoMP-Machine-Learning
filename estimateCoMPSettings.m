function [validity,model] = estimateCoMPSettings(simulation_data, isStatic)
    model = [];
    validity = 0;
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
        inclusion = find(RSRP > -Inf);  % Get rid of -inf dBm

        RSRP = RSRP(inclusion);
        BLER = BLER(inclusion);
        CQI = CQI(inclusion);
        TBSINR = TBSINR(inclusion);
           
        [err,model] = SVM(TBSINR,RSRP,BLER);
        if (err > epsilon) 
            validity = 0;            
        else
            validity = 1;
        end     
    end 