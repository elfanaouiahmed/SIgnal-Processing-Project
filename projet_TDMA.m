%% ============================================================
%  TP MF-TDMA — Démodulation AM + décodage (avec Toolbox Signal)
% ============================================================

clear; clc;

%% 1) Chargement du signal reçu MF-TDMA
data = load('Signal_MF_TDMA.mat');
signal = data.signal(:);     % on force un vecteur colonne (plus robuste)

%% 2) Paramètres du sujet
Fe = 200e3;                  % fréquence d'échantillonnage [Hz]
Te = 1/Fe;                   % période d'échantillonnage [s]

F1 = 15e3;                   % porteuse allouée à l'utilisateur 1 [Hz]
F2 = 45e3;                   % porteuse allouée à l'utilisateur 2 [Hz]

Ts = 50e-6;                  % durée d'un bit [s]
Ns = round(Ts*Fe);           % nb d'échantillons par bit (normalement 10)

Tslot = 24e-3;               % durée d'un slot [s]
Nslot = round(Tslot*Fe);     % nb d'échantillons par slot (normalement 4800)
t = (0:Nslot-1)' * Te;       % vecteur temps sur un slot

%% 3) Extraction des slots (d'après l'énoncé)
% Utilisateur 1 : slot n°2, porteuse F1
% Utilisateur 2 : slot n°5, porteuse F2
slotU1 = signal(Nslot+1 : 2*Nslot);
slotU2 = signal(4*Nslot+1 : 5*Nslot);

%% 4) Conception d'un filtre passe-bas (sinc fenêtrée Hamming)
Fc = 10e3;                   % fréquence de coupure [Hz]
fc_tilde = Fc/Fe;            % fréquence normalisée (par rapport à Fe)

M = 201;                     % longueur impaire du filtre
n = double(-(M-1)/2 : (M-1)/2);

% réponse impulsionnelle idéale d'un passe-bas (sinc)
h = (2*fc_tilde) * sinc(2*fc_tilde*n);

% fenêtrage Hamming pour réduire les oscillations (Gibbs)
w = hamming(M).';
h_pb = h .* w;

% normalisation : gain DC = 1 (évite de changer l'amplitude moyenne)
h_pb = h_pb / sum(h_pb);

%% 5) Retour en bande de base (démodulation AM cohérente) + filtrage
% On multiplie par cos(2πF t), puis on passe-bas pour enlever le terme à 2F.
bbU1 = slotU1 .* cos(2*pi*F1*t);
bbU2 = slotU2 .* cos(2*pi*F2*t);

% filtfilt : filtrage aller-retour => pas de retard de phase
messageU1 = filtfilt(h_pb, 1, bbU1);
messageU2 = filtfilt(h_pb, 1, bbU2);

%% 6) Décodage sans bruit
txtU1 = Demod_BdB(messageU1, Ns);
txtU2 = Demod_BdB(messageU2, Ns);

fprintf('Sans bruit; Message utilisateur 1 : %s\n', txtU1);
fprintf('Sans bruit; Message utilisateur 2 : %s\n', txtU2);

%% 7) Ajout de bruit AWGN + décodage (plusieurs SNR)
SNR_db_list = [20, 10, 0, -5];

% puissance moyenne du signal reçu
P_signal = mean(abs(signal).^2);

for SNR_db = SNR_db_list
    % puissance du bruit correspondant au SNR demandé
    P_bruit = P_signal / (10^(SNR_db/10));

    % bruit gaussien blanc (même taille que signal)
    bruit = sqrt(P_bruit) * randn(size(signal));
    signal_bruite = signal + bruit;

    % on ré-extrait les mêmes slots
    slotU1_b = signal_bruite(Nslot+1 : 2*Nslot);
    slotU2_b = signal_bruite(4*Nslot+1 : 5*Nslot);

    % retour bande de base + filtrage
    bbU1_b = slotU1_b .* cos(2*pi*F1*t);
    bbU2_b = slotU2_b .* cos(2*pi*F2*t);

    messageU1_b = filtfilt(h_pb, 1, bbU1_b);
    messageU2_b = filtfilt(h_pb, 1, bbU2_b);

    % décodage
    txtU1_b = Demod_BdB(messageU1_b, Ns);
    txtU2_b = Demod_BdB(messageU2_b, Ns);

    fprintf('\nSNR : %d dB\n', SNR_db);
    fprintf('Avec bruit; Message utilisateur 1 : %s\n', txtU1_b);
    fprintf('Avec bruit; Message utilisateur 2 : %s\n', txtU2_b);
end