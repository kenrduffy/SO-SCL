% All code is subject to license:
% GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf

% Simulation assesses the use of of blockwise SO from SO-SCL to manage
% undetected error rate.

% SO-SCL
% P. Yuan, K. R. Duffy & M. Médard. "Near-optimal generalized decoding of 
% Polar-like codes.", IEEE ISIT, 2024. 
% P. Yuan, K. R. Duffy & M. Médard. "Soft-output successive cancellation 
% list decoding", IEEE Transactions on Information Theory, 71 (2), 
% 1007–1017, 2025.




clear;
%% Monte-Carlo parameters
EbN0dB = 2;
maxIt = 10^5;
minlogp = 5;
step = 0.5;
r = 10.^-(0:step:minlogp);
%% Code parameters
load("code/dRM_64_42.mat")
%% Decoder parameters
L = 1;
%% 
R = k / n;
EsN0dB = EbN0dB + 10 * log10(2*k/n);
scal = sqrt(10^(EsN0dB / 10));
%% Loop over SNRs
n_sample = zeros(1, length(r));
n_negative = zeros(1, length(r));
n_p = zeros(1, length(r));
for ntx = 1:maxIt
    msg = randsrc(k,1,[0 1]);
    u = preencode_dpolar(msg, frz, dCons);
    c = polarTrans(u, 0);

    x = (1 - 2 * c) * scal;
    y = x + randn([n, 1]);
    llr = 2 * scal * y;

    % confidence
    [~, ~, ~, chat_list, ~, p_notInList] = SOSCL(llr, frz, dCons, L);

    ind = ceil(-(log10(p_notInList) / step));
    ind = max(1, ind);
    ind = min(ind, length(r));
    n_p(ind) = n_p(ind) + p_notInList;
    n_sample(ind) = n_sample(ind) + 1;

    %
    flag = 0;
    for l = 1:L
        if isequal(chat_list(:, l), c)
            flag = 1;
            break
        end
    end
    %
    if flag == 0
        n_negative(ind) = n_negative(ind) + 1;
    end
    %
    if mod(ntx, maxIt/10) == 0
        disp(['@ ', num2str(ntx), ' / ', num2str(maxIt)])
    end
end

figure(1)
clf
loglog([10^(-minlogp - 1), 1], [10^(-minlogp - 1), 1], 'k','LineWidth',2)
hold on
grid on
p = n_negative ./ n_sample;
pp = n_p ./ n_sample;
loglog(pp, p, 'rx-','LineWidth',2)
grid 'on'
xlabel('Predicted')
ylabel('Conditional')
set(gca,'FontSize',16)
