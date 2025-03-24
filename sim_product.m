% All code is subject to license:
% GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf

% Simulation peforms turbo product code decoding with SO-SCL as the SISO
% component decoder.

% SO-SCL
% P. Yuan, K. R. Duffy & M. Médard. "Near-optimal generalized decoding of 
% Polar-like codes.", IEEE ISIT, 2024. 
% P. Yuan, K. R. Duffy & M. Médard. "Soft-output successive cancellation 
% list decoding", IEEE Transactions on Information Theory, 71 (2), 
% 1007–1017, 2025.


clear;
%% Monte-Carlo parameters
EbN0dB      = 2:0.25:3.5;
NoErrors    = 20;
maxIt       = 10^6;
minIt       = 10^2;
%% Code parameters
load("code/RM_64_57.mat")
%% Decoder parameters
L           = 4;  % Maximum list size
Imax        = 20; % maximum number of iterations
alpha       = 0.5 * ones(1, 50); % Extrinsic LLR scaling
%% Code and channel
R           = (k / n)^2;
EsN0dB      = EbN0dB + 10 * log10(2*R);
numSNR      = length(EsN0dB);
%% Loop over SNRs
BLER        = zeros(1, numSNR);
BER         = zeros(1, numSNR);
Iavg        = zeros(1, numSNR);
for sp = 1:numSNR
    BlockError  = 0;
    BitError    = 0;
    n_iter      = 0;
    ntx         = 0;
    sigma = 1 / sqrt(10^(EsN0dB(sp) / 10));
    while ((BlockError < NoErrors && ntx < maxIt) || ntx < minIt)
        ntx = ntx + 1;
        c = zeros(n, n); % To store codeword.
        %% Enc
        u = randsrc(k, k, [0, 1]);
        for row = 1:k
            c(row, :) = mod(u(row, :)*G, 2);
        end
        for col = 1:n
            c(:, col) = mod(c(1:k, col)'*G, 2);
        end
        %% binary input AWGN channel
        x = 1 - 2 * c;
        y = x + sigma * randn([n, n]);
        L_channel = 2 * y / (sigma^2);
        %% Decoding
        L_APP = zeros(size(L_channel));
        L_E = zeros(size(L_channel));
        L_A = zeros(size(L_channel));
        c_HD = 0.5 * (1 - sign(L_channel));
        for i = 1:Imax
            %% row
            n_iter = n_iter + 0.5;
            L_A = alpha(2*i-1) * L_E;
            input = L_channel + L_A;
            for row = 1:n
                [L_APP(row, :), L_E(row, :),] = SOSCL(input(row, :)', frz, dCons, L);
            end
            % early-termination
            c_HD = 0.5 * (1 - sign(L_APP));
            s1 = mod(c_HD*H', 2);
            s2 = mod(c_HD'*H', 2);
            if sum(s1(:)) == 0 && sum(s2(:)) == 0
                break
            end
            %% column
            n_iter = n_iter + 0.5;
            L_A = alpha(2*i) * L_E;
            input = L_channel + L_A;
            for col = 1:n
                [L_APP(:, col), L_E(:, col)] = SOSCL(input(:, col), frz, dCons, L);
            end
            %% early-termination
            c_HD = 0.5 * (1 - sign(L_APP));
            s1 = mod(c_HD*H', 2);
            s2 = mod(c_HD'*H', 2);
            if sum(s1(:)) == 0 && sum(s2(:)) == 0
                break
            end
        end
        %% error collection
        if (~isequal(c, c_HD))
            BlockError  = BlockError + 1;
            BitError    = BitError + sum(c(:) ~= c_HD(:));
        end
    end
    disp(['---[' num2str(n) ',' num2str(k) ']^2---Eb/N0 dB ', num2str(EbN0dB(sp)), ' dB:---'])
    BLER(sp)    = BlockError / ntx;
    BER(sp)     = BitError / (ntx * n * n);
    Iavg(sp)    = n_iter / ntx;
    disp([' BLER             = ', num2str(BLER(sp))]);
    disp([' BER              = ', num2str(BER(sp))]);
    disp([' Iavg             = ', num2str(Iavg(sp))]);

    save(['./results/prod-soscl-' num2str(n) '-' num2str(k) '-' num2str(L) '.mat'], 'EbN0dB', 'BLER', 'BER', 'G', 'H', 'dCons', 'frz');
end
