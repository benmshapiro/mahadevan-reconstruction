% Description: Makes a synthetic dataset for testing metabolic rearrangment in a 
% community of binned species
% Author: Benjamin Shapiro, University of Oregon
% Date: 2015-01-23
% Modified: 2015-01-30
% Comments: Use in conjunction with MahadevanSimpleModel.m

mets = { %.mets, .metNames = .metFormulas = .mets
    'Glcxt'; %1
    'Ac'   ; %2
    'O2'   ; %3
    'X'  ;}; %4
 
rxns = { % .rxns, .rxnNames,                        .rev,   .lb,    .ub
    'v1' '39.43 Ac + 35 O2 to X'                    0       0       1000 ;
    'v2' '9.46 Glcxt + 12.92 O2 to X'               0       0       1000 ;
    'v3' '9.84 Glcxt + 12.73 O2 to 1.24 Ac + X'     0       0       1000 ;
    'v4' '19.23 Glcxt to 12.12 Ac + X'              0       0       1000 ;};

S = zeros(4,4);
S(2,1) = -39.43; S(3,1) = -35;  S(4,1) = 1; %v1
S(1,2) = -9.46; S(3,2) = -12.92; S(4,2) = 1; %v2
S(1,3) = -9.84; S(3,3) = -12.73; S(2,3) = 1.24; S(4,3) = 1; %v3
S(1,4) = -19.23; S(2,4) = 12.12; S(4,4) = 1; %v4

%Objective function is the X sink (v5)
msm.mets = mets(:,1);
msm.metNames = mets(:,1);
msm.metFormulas = mets(:,1);
r1 = 1:4;
msm.rxns = {rxns{r1,1}}';
msm.rxnNames = {rxns{r1,2}}';
msm.rev = [rxns{r1,3}]';
msm.lb = [rxns{r1,4}]';
msm.ub = [rxns{r1,5}]';
msm.S = sparse([S(:,r1)]);
msm.description = 'Mahadevan Simple Model';