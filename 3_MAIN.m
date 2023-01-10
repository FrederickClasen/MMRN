%% FLUX SIMULATIONS PRESENTED IN CLASEN ET AL              %%%%%%%%%%%%%
% COBRA and RAVEN has to be added to the path to run the simulations

clear;
setRavenSolver('mosek');          % SEE MOSEK WEBSITE FOR MORE DETAIL
% RER AND OBJECTIVES
RERRxns = {'MMRNR10362','MMRNR10338'}; % O2 uptake and CO2 production    
objective = {'MMRN_Biomass'};          % OBJECTIVE FUNCTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RUN FLUX SIMULATIONS FOR GENERIC AND CONSTRAINT-BASED GSMMs %% 

% load diet information
conditions = {'nonDEN_Liver_CD','nonDEN_Liver_WD','DEN_Liver_CD','DEN_AdjLiver_WD','DEN_Tumour_WD'};
CD = importModel('model/MMRNHep/CD/MMRNHep-CD.xml',false);
WD = importModel('model/MMRNHep/WD/MMRNHep-WD.xml',false);
models = {CD,WD};
fluxVector = zeros(length(WD.rxns),(length(conditions)*2)+2); 
diets = {'CD','WD'};

counter = 0;
for i=1:length(diets)
    counter = counter + 1;
    model = models{i};
    % SET RER AND OBJECTIVE
    model = setParam(model,'lb',RERRxns,[17,13]);
    model = setParam(model,'ub',RERRxns,[1000,1000]);
    model = setParam(model,'obj',objective,1);
    fba = solveLP(model,1); %RUN FBA
    fluxVector(:,counter) = fba.x;
    
    for j=1:length(conditions)
        counter = counter + 1;
        fc = 'data/Eflux/' + string(conditions(j)) + '.csv';
        constraints = importdata(fc);
        csModel = addConstraints(model,constraints);  % EFLUX CONSTRAINT
        
        if contains(conditions(j),'nonDEN')
            csModel = setParam(csModel,'lb',RERRxns,[15,11]); % RER CONSTRAINT
            csModel = setParam(csModel,'ub',RERRxns,[1000,1000]);
        else
            csModel = setParam(csModel,'lb',RERRxns,[17,13]); % RER CONSTRAINT
            csModel = setParam(csModel,'ub',RERRxns,[1000,1000]);
        end
        
        csModel = setParam(csModel,'obj',objective,1);
        fba = solveLP(csModel,1);
        fluxVector(:,counter) = fba.x;
        
        modelName = string(diets{i}) + '_' + string(conditions{j});
        
    end
    
end





