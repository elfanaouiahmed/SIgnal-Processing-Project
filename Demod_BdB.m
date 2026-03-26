%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Démodulation bande de base avec les hypothèses suivantes :
% - mapping binaire à moyenne nulle
% - mise en forme rectangulaire de durée Ts=NsTe, Te étant la période
% d'échantillonnage et Ns le nombre d'échantillon par symbole
% - symboles embrouillés de manière connue
% 
% Et reconstruction du texte à partir de la suite de bits retrouvés
%
% Paramètres d'entrée :
% - signal à démoduler
% - Ns
% Paramètre de sortie : le texte retrouvé
% 
% Nathalie Thomas, Janvier 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function texte=Demod_BdB(message, Ns)

%Filtrage adapté 
signal_sortie_FA=filter(ones(1,Ns),1,message);

%échantillonnage à Ns
signal_sortie_FA_ech=signal_sortie_FA(Ns:Ns:end);

%Désembrouillage des symboles
mat_desembrouillage=reshape(signal_sortie_FA_ech,20,24);
signal_sortie_FA_ech=reshape(mat_desembrouillage.',1,480);

%décisions sur les symboles
symboles_retrouves=sign(signal_sortie_FA_ech);

%Demapping
bits_retrouves=(symboles_retrouves+1)/2;

%Reconstruction du message utilisateur1 à partir de la suite de bits retrouvée
texte=bin2str(bits_retrouves);