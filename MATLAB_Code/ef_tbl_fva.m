function [ tbl ] = ef_tbl_fva( pct,model,tbl,insert )
%ef_tbl_fva performs flux variablility analysis with inputs of a required
%percentage of the objective function, cobra model, and table to add
%results of analysis

tbl = struct(tbl);
[minF,maxF] = fluxVariability(model,pct);

fva_req = (double(minF > 1*10^-9) + double(maxF < -1*10^-9)) > 0;
fva_on = (double(abs(maxF) > (1*10^-9)) + double(abs(minF) > 1*10^-9)) > 0;

if(insert ~= 0)
    start = length(tbl.fva_lb)+1;
    fin = length(tbl.fva_lb)+length(minF);
    tbl.fva_lb(start:fin) = minF;
    tbl.fva_ub(start:fin) = maxF;
    tbl.fva_req(start:fin) = fva_req;
    tbl.fva_on(start:fin) = fva_on;
    tbl.fva_pct(start:fin) = zeros(length(minF),1) + pct;
    tbl.rxn_id(start:fin) = tbl.rxn_id(1:length(fva_req));
    tbl.rxn_name(start:fin) = tbl.rxn_name(1:length(fva_req));
    tbl.lb(start:fin) = tbl.lb(1:length(fva_req));
    tbl.ub(start:fin) = tbl.ub(1:length(fva_req));
    tbl.rxn_formula(start:fin) = tbl.rxn_formula(1:length(fva_req));
else
    tbl.fva_lb = minF;
    tbl.fva_ub = maxF;
    tbl.fva_req = fva_req;
    tbl.fva_on = fva_on;
    tbl.fva_pct = zeros(length(minF),1) + pct;

end
    

end

