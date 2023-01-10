% TESTS FOR ATP HYDROLYSIS IN CLOSED CS-GSMMs
% ALL MODELS RETURN INFEASIBLE SOLUTION

clear;
setRavenSolver('mosek');
RERRxns = {'MMRNR10362','MMRNR10338'}; % O2 uptake and CO2 production    
objective = {'MMRNR06619'};            % Objective (ATP Hydrolysis)

% RUN FLUX SIMULATIONS FOR CONSTRAINT-BASED GSMMs %% 

% load diet models
conditions = {'nonDEN_Liver_CD','nonDEN_Liver_WD','DEN_Liver_CD','DEN_AdjLiver_WD','DEN_Tumour_WD'};
CD = importExcelModel('model/MMRNHep/CD/MMRNHep-CD.xlsx',false);
WD = importExcelModel('model/MMRNHep/WD/MMRNHep-WD.xlsx',false);
models = {CD,WD};
diets = {'CD','WD'};

for i=1:length(diets)
    model = models{i};
    % set objective as atp hydrolysis
    model = setParam(model,'obj',objective,1);
    % block exhanges
    [selExc, selUpt] = findExcRxns(model);
    exch = union(model.rxns(selExc),model.rxns(selUpt));
    model = setParam(model,'eq',exch,0);
    
    for j=1:length(conditions)
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
        
        fba = solveLP(csModel) %fba will return an infeasible solution
        
    end
    
end

