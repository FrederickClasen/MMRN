function outModel = addConstraints(model,constraints)

constraints.textdata(1,:) = [];
constraints.textdata(:,2) = [];
constraints.textdata(:,2) = [];
constraints.rxn = constraints.textdata;
constraints.data = 2.5*(constraints.data); 
constraints.lb = constraints.data(:,[1]);
constraints.lb(constraints.lb < -1000) = -1000;
constraints.ub = constraints.data(:,[2]);
constraints.ub(constraints.ub > 1000) = 1000;

outModel = setParam(model,'ub',constraints.rxn,constraints.ub);
outModel = setParam(outModel,'lb',constraints.rxn,constraints.lb);
