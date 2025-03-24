% All code is subject to license:
% GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf

% Simulation assesses the calibration of blockwise SO from SO-SCL.

% SO-SCL
% P. Yuan, K. R. Duffy & M. Médard. "Near-optimal generalized decoding of 
% Polar-like codes.", IEEE ISIT, 2024. 
% P. Yuan, K. R. Duffy & M. Médard. "Soft-output successive cancellation 
% list decoding", IEEE Transactions on Information Theory, 71 (2), 
% 1007–1017, 2025.


clear;
%% Code parameters
load("code/dRM_64_42.mat")
%% Decoder parameters
L           = 4;
p_e         = 0.1;
%% Monte-Carlo parameters
EbN0dB      = 1:0.5:6;
NoErrors    = 50/p_e;
maxIt       = 10^9;
minIt       = 10^4;
%% Code and channel
EsN0dB      = EbN0dB + 10 * log10(2*k/n);
numSNR      = length(EsN0dB);

%% Loop over SNRs
BLER        = ones(1, numSNR);
UER         = ones(1, numSNR);
ER          = ones(1, numSNR);
for sp = 1:numSNR
    BlockError              = 0;
    undetectedBlockError    = 0;
    Erasure                 = 0;
    ntx                     = 0;
    scal = sqrt(10^(EsN0dB(sp) / 10));

    while ((BlockError < NoErrors && ntx < maxIt) || ntx < minIt)
        ntx = ntx + 1;
        msg = randsrc(k,1,[0 1]);
        u = preencode_dpolar(msg, frz, dCons);
        c = polarTrans(u, 0);
        x = (1 - 2 * c) * scal;
        y = x + randn([n, 1]);
        llr = 2 * scal * y;

        [~, ~, chat, ~, p_incorrect, ~] = SOSCL(llr, frz, dCons, L);

        if p_incorrect > p_e
            flag = 0;
        else
            flag = 1;
        end

        if flag == 0
            BlockError = BlockError + 1;
            Erasure = Erasure + 1;
        elseif (~isequal(c, chat))
            BlockError = BlockError + 1;
            undetectedBlockError = undetectedBlockError + 1;
        end
    end
    BLER(sp)    = BlockError / ntx;
    UER(sp)   = undetectedBlockError / ntx;
    ER(sp)      = Erasure / ntx;
    disp(['---[' num2str(n) ',' num2str(k) ']---'])
    disp([num2str(EbN0dB(sp)), 'dB: BLER             = ', num2str(BLER(sp))]);
    disp([num2str(EbN0dB(sp)), 'dB: UER              = ', num2str(UER(sp))]);
    disp([num2str(EbN0dB(sp)), 'dB: ER               = ', num2str(ER(sp))]);
    save(['./results/uer-soscl-' num2str(n) '-' num2str(k) '-' num2str(L) '-' num2str(p_e) '.mat'], 'EbN0dB', 'BLER', 'UER', 'ER');

end
