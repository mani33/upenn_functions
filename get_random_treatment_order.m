function varargout = get_random_treatment_order(nMicePerCond)
% get the treatment order in which saline and sch or any other drug is
% given during field recording experiment.
% Mani Subramaniyan: 2017-12-14

ind = randperm(nMicePerCond*2);
cond = [zeros(1,nMicePerCond),ones(1,nMicePerCond)];
treatOrder = cond(ind);
if nargout
    varargout{1} = treatOrder;
end
disp(treatOrder)

