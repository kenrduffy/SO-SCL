function [L_APP, L_E, chat, chat_list, p_incorrect, p_notinList] = SOSCL(llr, frz, dCons, L)
[~, chat_list, PM, p_NL] = SOSCL_mex(llr, frz, dCons, L, 0);
p_NL = max(eps, exp(-p_NL));
%% bitSO
assert( size(chat_list,2) == length(PM), "List size error !")
% channel observation
pp1 = 1./(1+exp(llr)); 
pp0 = 1 - pp1;
% avoid numerical problem
pp1 = max(pp1, eps); pp1 = min(pp1, 1-eps);
pp0 = max(pp0, eps); pp0 = min(pp0, 1-eps);
%
p = exp(-PM);
p1 = sum(chat_list.*p,2);
p0 = sum((1-chat_list).*p,2);

p0 = p0 + p_NL * pp0;
p1 = p1 + p_NL * pp1;

L_APP = log(p0) - log(p1);
L_E = L_APP - llr;
%% blkSO
i_ml = find(PM == min(PM));
chat = chat_list(:, i_ml);
pC = p_NL + sum(p);
p_incorrect = 1 - (p(i_ml) / pC);
p_notinList = 1 - (sum(p) / pC);
end