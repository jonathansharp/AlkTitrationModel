function [Result,Result_Headers] = TITRATE(SAL,TMP,DIC,PH,TP,TSI,NH4,...
                                       H2S,MSW,CA,TORG1,PKORG1,TORG2,...
				                       PKORG2,TORG3,PKORG3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INFORMATION
%
%  This function will provide measured total alkalinity values of a
%  solution with the specified salinity, temperature, dissolved inorganic
%  carbon, pH, and total concentrations of inorganic chemical species,
%  along with organic chemical species with the specified total
%  concentrations and acid dissociation constants. Five different methods
%  of titration and data analysis are used in this function to provide five
%  distinct measured alkalinity values.
%
%  ***SYNTAX:
%   [Result,Result_Headers] = 
%   TITRATE(SAL,TMP,DIC,PH,TP,TSI,NH4,H2S,MSW,CA,TORG1,PKORG1,TORG2,PKORG2,TORG3,PKORG3)
%
%  ***SYNTAX EXAMPLE:
%   [Result,Result_Headers] =
%   TITRATE(35,25,2000,8.1,0.2,8.0,0,0,200,0.2,5,4.5,5,9.0,0,0.0)
%
%  ***INPUT:
%    SAL:   Salnity
%    TMP:    Titration temperature (degrees Celcius)
%    DIC:    Total dissolved inorganic carbon (umol/kgSW)
%    PH:     Total scale pH
%    TP:     Total dissolved phosphate (umol/kgSW)
%    TSI:    Total dissolved silicate (umol/kgSW)
%    NH4:    Total dissolved ammonia (umol/kgSW)
%    H2S:    Total dissolved hydrogen sulfide (umol/kgSW)
%    MSW:    Mass of titrated seawater sample
%    CA:     Acid concentration (mol/kg)
%    TORGx:  Total organic concentration x (umol/kgSW)
%    PKORGx: Organic acid dissociation constant x
%
%  ***OUTPUT:
%    Result = array containing the following values (one row per sample):
%    COL  VALUE
%    1    Inorganic Alkalinity defined by Dickson (1981)
%    2    Organic Alkalinity (organic bases treated using pK = 4.5 as the
%         cutoff for the zero level of protons)
%    3    Open Cell Gran Alkalinity
%    4    Closed Cell Nonlinear Fit Alkalinity
%    5    Open Cell Nonlinear Fit Alkalinity
%    6    Difference Derivative Alkalinity
%    7    Single-Step Alkalinity
%         * Values 3-7 result from a given treatment of the simulated
%           titration data
%    Result_Headers = cell array containing descriptions of 'Result'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUT CONDITIONING

% Determine lengths of input vectors
veclengths = ...
[ length(SAL)   length(TMP)    length(DIC)    length(PH)        length(TP)...
  length(TSI)   length(NH4)    length(H2S)    length(MSW)       length(CA)...
  length(TORG1) length(PKORG1) length(TORG2)  length(PKORG2)    length(TORG3)...
  length(PKORG3)];

if length(unique(veclengths))>2
	disp(' '); disp('*** INPUT ERROR: Input vectors must all be of same length, or of length 1. ***'); disp(' '); return
end

% Make column vectors of all input vectors
SAL    = SAL    (:);
TMP    = TMP    (:);
DIC    = DIC    (:);
PH     = PH     (:);
TP     = TP     (:);
TSI    = TSI    (:);
NH4    = NH4    (:);
H2S    = H2S    (:);
MSW    = MSW    (:);
CA     = CA     (:);
TORG1  = TORG1  (:);
PKORG1 = PKORG1 (:);
TORG2  = TORG2  (:);
PKORG2 = PKORG2 (:);
TORG3  = TORG3  (:);
PKORG3 = PKORG3 (:);

% Find the longest column vector:
ntps   = max(veclengths);
Result = zeros(ntps,7);

% Populate column vectors
SAL(1:ntps,1)    = SAL(:);
TMP(1:ntps,1)    = TMP(:);
DIC(1:ntps,1)    = DIC(:);
PH(1:ntps,1)     = PH(:);
TP(1:ntps,1)     = TP(:);
TSI(1:ntps,1)    = TSI(:);
NH4(1:ntps,1)    = NH4(:);
H2S(1:ntps,1)    = H2S(:);
MSW(1:ntps,1)    = MSW(:);
CA(1:ntps,1)     = CA(:);
TORG1(1:ntps,1)  = TORG1(:);
PKORG1(1:ntps,1) = PKORG1(:);
TORG2(1:ntps,1)  = TORG2(:);
PKORG2(1:ntps,1) = PKORG2(:);
TORG3(1:ntps,1)  = TORG3(:);
PKORG3(1:ntps,1) = PKORG3(:);

%% DEFINE RATIOS OF CONSERVATIVE CONSTITUENTS

% BORON from Lee et al. (2010)
  TB  = (0.0004326.*(SAL./35)).*10.^6;
% SULFUR from Morris and Riley (1966)
  TS  = ((0.14./96.062).*(SAL./1.80655)).*10.^6;
% FLOURIDE from Riley (1965)
  TF  = ((0.000067./18.998).*(SAL./1.80655)).*10.^6;

%% DEFINE EQUILIBRIUM CONSTANTS (ON THE TOTAL SCALE UNLESS SPECIFIED)

% Boric acid dissociation constant from Dickson (1990); S=5:45, T=0:45
  KB = exp((-8966.9-2890.53.*sqrt(SAL)-77.942.*SAL+1.728.*sqrt(SAL).*SAL-...
       0.0996.*SAL.^2)./(TMP+273.15)+148.0248+137.1942.*sqrt(SAL)+1.62142.*...
       SAL+(-24.4344-25.085.*sqrt(SAL)-0.2474.*SAL).*log(TMP+273.15)+...
       0.053105.*sqrt(SAL).*(TMP+273.15));
% Sulfuric acid dissociation constant from Dickson (1990)
  KS = exp(-4276.1./(TMP+273.15)+141.328-23.093.*log(TMP+273.15)+(-13856./...
       (TMP+273.15)+324.57-47.986.*log(TMP+273.15)).*sqrt((19.924.*SAL./...
       (1000-1.005.*SAL)))+(35474./(TMP+273.15)-771.54+114.723.*log(TMP+...
       273.15)).*(19.924.*SAL./(1000-1.005.*SAL))+(-2698./(TMP+273.15)).*...
       sqrt((19.924.*SAL./(1000-1.005.*SAL))).*(19.924.*SAL./(1000-1.005.*...
       SAL))+(1776./(TMP+273.15)).*(19.924.*SAL./(1000-1.005.*SAL)).^2).*...
       (1-0.001005.*SAL); % FREE SCALE
% Define factor (Z) to convert from free scale to total scale
  Z = (1+(TS.*10.^-6)./KS);
% Hydrogen flouride dissociation constant from Dickson and Riley (1979)
  KF = exp(1590.2./(TMP+273.15)-12.641+1.525.*(19.924.*SAL./(1000-...
       1.005.*SAL)).^0.5).*(1-0.001005.*SAL); % FREE SCALE
% Define factor (Y) to convert from seawater scale to total scale
  Y = Z./(1+(TS.*10.^-6)./KS+(TF.*10.^-6)./KF);
% Define KF on total scale
  KF = KF.*Z;

for w = 1:ntps;
% Carbonic acid dissociation constants
if  SAL(w) >= 20
% Carbonic acid dissociation constants from Leuker et al. (2000)
% S=19:43 T=2:35
  KC1(w) = 10.^-(3633.86./(TMP(w)+273.15)-61.2172+9.6777.*log(TMP(w)+273.15)-...
           0.011555.*SAL(w)+0.0001152.*SAL(w).^2)';     
  KC2(w) = 10.^-(471.78./(TMP(w)+273.15)+25.929-3.16967.*log(TMP(w)+273.15)-...
           0.01781.*SAL(w)+0.0001122.*SAL(w).^2)';
elseif   SAL(w) > 0 && SAL(w) < 20
% Carbonic acid dissociation constants from Waters et al. (2014), which fixed
% inconsistencies with Millero (2010) identified by Orr et al. (2015)
% S=1:50 T=0:50
  pK10(w) = -126.34048 + 6320.813./(TMP(w)+273.15) + 19.568224.*log(TMP(w)+273.15);
    A1(w) = 13.409160.*SAL(w).^0.5 + 0.031646.*SAL(w) - 5.1895e-5.*SAL(w).^2;
    B1(w) = -531.3642.*SAL(w).^0.5 - 5.713.*SAL(w);
    C1(w) = -2.0669166.*SAL(w).^0.5;
   pK1(w) = pK10(w) + A1(w) + B1(w)./(TMP(w)+273.15) + C1(w).*log(TMP(w)+273.15);
   KC1(w) = (10.^-pK1(w)).*Y(w);
  pK20(w) =  -90.18333 + 5143.692./(TMP(w)+273.15) + 14.613358.*log(TMP(w)+273.15);
    A2(w) = 21.225890.*SAL(w).^0.5 + 0.12450870.*SAL(w) - 3.7243e-4.*SAL(w).^2;
    B2(w) = -779.3444.*SAL(w).^0.5 - 19.91739.*SAL(w);
    C2(w) = -3.3534679.*SAL(w).^0.5;
   pK2(w) = pK20(w) + A2(w) + B2(w)./(TMP(w)+273.15) + C2(w).*log(TMP(w)+273.15);
   KC2(w) = (10.^-pK2(w)).*Y(w);
elseif   SAL(w) == 0
% Carbonic acid dissociation constants from Millero (1979)
% S=0 T=0:50
  lnK1(w) = 290.9097 - 14554.21./(TMP(w)+273.15) - 45.0575.*log(TMP(w)+273.15);
   KC1(w) = exp(lnK1(w));
  lnK2(w) = 207.6548 - 11843.79./(TMP(w)+273.15) - 33.6485.*log(TMP(w)+273.15);
   KC2(w) = exp(lnK2(w));
end

end

% Dissociaton of water from Millero (1995)
  KW = (exp(148.9802-13847.26./(TMP+273.15)-23.6521.*log(TMP+273.15)+...
       (-5.977+118.67./(TMP+273.15)+1.0495.*log(TMP+273.15)).*(SAL.^0.5)-...
       0.01615*SAL)).*Y;
% First phosphoric acid dissociation constant from Yao and Millero (1995)
  KP1 = (exp(-4576.752./(TMP+273.15)+115.54-18.453.*log(TMP+273.15)+...
        (-106.736./(TMP+273.15)+0.69171).*(SAL.^0.5)+(-0.65643./...
        (TMP+273.15)-0.01844).*SAL)).*Y;
% Second phosphoric acid dissociation constant from Yao and Millero (1995)
  KP2 = (exp(-8814.715./(TMP+273.15)+172.1033-27.927.*log(TMP+273.15)+...
        (-160.34./(TMP+273.15)+1.3566).*(SAL.^0.5)+(0.37335./(TMP+273.15)-...
        0.05778).*SAL)).*Y;
% Third phosphoric acid dissociation constant from Yao and Millero (1995)
  KP3 = (exp(-3070.75./(TMP+273.15)-18.126+(17.27039./(TMP+273.15)+...
        2.81197).*SAL.^0.5+(-44.99486./(TMP+273.15)-0.09984).*SAL)).*Y;
% Silicic acid dissociation constant from Yao and Millero (1995)
  KSI = (exp(-8904.2./(TMP+273.15)+117.4-19.334.*log(TMP+273.15)+(-458.79./...
        (TMP+273.15)+3.5913).*sqrt(19.924.*SAL./(1000-1.005.*SAL))+(188.74./...
        (TMP+273.15)-1.5998).*(19.924.*SAL./(1000-1.005.*SAL))+(-12.1652./...
        (TMP+273.15)+0.07871).*(19.924.*SAL./(1000-1.005*SAL)).^2).*(1-...
        0.001005.*SAL)).*Y;
% Ammonium dissociation constant from Yao and Millero (1995)
  KNH4 = (exp(-6285.33./(TMP+273.15)+0.0001635.*(TMP+273.15)-0.25444+...
         (0.46532-123.7184./(TMP+273.15)).*SAL.^0.5+(-0.01992+3.17556./(TMP+...
         273.15)).*SAL)).*Y;
% First hydrogen sulfide dissociation constant from Millero et al. (1988)
  KS1 = (exp(225.838-13275.3./(TMP+273.15)-34.6435.*log(TMP+273.15)+0.3449.*...
        SAL.^0.5-0.0274.*SAL)).*Y;
% Defined organic acid dissociation constant
  KORG1 = 10.^-PKORG1; KORG2 = 10.^-PKORG2; KORG3 = 10.^-PKORG3;
  
%% IDEAL GAS CONSTANT (kCal/(mol*K)) AND FARADAY CONSTANT (kCal/mol)

   Rg = 8.3144621; F = 96.4853365;
  
%% START LOOP FOR EACH LINE OF INPUT DATA

for x = 1:ntps

% Set organic proton acceptors and donors to zero
ORGa1  = 0; ORGd1  = 0; ORGap1  = 0; ORGdp1  = 0;
ORGa2  = 0; ORGd2  = 0; ORGap2  = 0; ORGdp2  = 0;
ORGa3  = 0; ORGd3  = 0; ORGap3  = 0; ORGdp3  = 0;

%% CALCULATE TA FROM pH AND TOTAL CARBON
HTOT   = 10.^(-PH(x));
HFr    = HTOT./Z(x);
HCO3   = (DIC(x).*10.^-6)./(1+HTOT./KC1(x)+KC2(x)./HTOT);
CO3    = (DIC(x).*10.^-6)./(1+HTOT.^2./(KC1(x).*KC2(x))+HTOT./KC2(x));
BOH4   = (TB(x).*10.^-6)./(1+HTOT./KB(x));
OH     = (KW(x)./HTOT);
PO4    = (TP(x).*10.^-6)./(1+HTOT./KP3(x)+(HTOT.^2)./(KP2(x).*KP3(x))+...
         (HTOT.^3)./(KP1(x).*KP2(x).*KP3(x)));
HPO4   = (TP(x).*10.^-6)./(1+HTOT./KP2(x)+(HTOT.^2)./(KP1(x).*KP2(x))+KP3(x)./HTOT);
SiO4H3 = (TSI(x).*10.^-6)./(1+HTOT./KSI(x));
NH3    = (NH4(x).*10.^-6)./(1+HTOT./KNH4(x));
HS     = (H2S(x).*10.^-6)./(1+HTOT./KS1(x));
HSO4   = (TS(x).*10.^-6)./(1+(KS(x)./HFr));
HF     = (TF(x).*10.^-6)./(1+(KF(x)./HTOT));
H3PO4  = (TP(x).*10.^-6)./(1+KP1(x)./HTOT+(KP1(x).*KP2(x))./(HTOT.^2)+...
         (KP1(x).*KP2(x).*KP3(x))./(HTOT.^3));
if        PKORG1(x) >= 4.5; ORGa1 = (TORG1(x).*10.^-6)./(1+HTOT./KORG1(x));
elseif    PKORG1(x) < 4.5;  ORGd1 = (TORG1(x).*10.^-6)./(1+KORG1(x)./HTOT); end
if        PKORG2(x) >= 4.5; ORGa2 = (TORG2(x).*10.^-6)./(1+HTOT./KORG2(x));
elseif    PKORG2(x) < 4.5;  ORGd2 = (TORG2(x).*10.^-6)./(1+KORG2(x)./HTOT); end
if        PKORG3(x) >= 4.5; ORGa3 = (TORG3(x).*10.^-6)./(1+HTOT./KORG3(x));
elseif    PKORG3(x) < 4.5;  ORGd3 = (TORG3(x).*10.^-6)./(1+KORG3(x)./HTOT); end
TA       = (10.^6).*(HCO3+2.*CO3+BOH4+OH+2.*PO4+HPO4+SiO4H3+NH3+HS+ORGa1+...
                     ORGa2+ORGa3-HFr-HSO4-HF-H3PO4-ORGd1-ORGd2-ORGd3);
TA_inorg = (10.^6).*(HCO3+2.*CO3+BOH4+OH+2.*PO4+HPO4+SiO4H3+NH3+HS-...
                     HFr-HSO4-HF-H3PO4);
TA_org   = (10.^6).*(ORGa1+ORGa2+ORGa3-ORGd1-ORGd2-ORGd3);

%% PREPARE ACIDIMETRIC TITRATION

% Declare initial acid addition number, pH tolerance
  I = (0.0035.*MSW)./CA; R = I./1000; S = 1; pHTol = 0.0000001;
  DATA = zeros(S,4); DATAp = zeros(S,4);

%% PERFORM STEP-WISE ACIDIMETRIC TITRATION

for MA = 0:R(x):I(x);

% Use either initial pH Guess of 8 or pH from last titration step
if exist('pHtit','var'); pH=pHtit; else pH=PH(x); end
% Seawater mass and acid mass in kg
kgSW = MSW(x).*(10.^-3); kgA = MA.*(10.^-3); kgTOT = kgSW+kgA;

%% CALCULATE pH DURING EACH TITRATION STEP
% Declare initial delta pH
deltapH = pHTol+1;
% Iterate to determine step-wise pH until deltapH is below tolerance limit
while any(abs(deltapH) > pHTol)
    % Concentrations of all species relevant to TA
    HTOT      = (10.^(-pH));
    HFr       = HTOT./Z(x);
    HCO3      = (DIC(x).*10.^-6)./(1+HTOT./KC1(x)+KC2(x)./HTOT);
    CO3       = (DIC(x).*10.^-6)./(1+HTOT.^2./(KC1(x).*KC2(x))+HTOT./KC2(x));
    BOH4      = (TB(x).*10.^-6)./(1+HTOT./KB(x));
    OH        = (KW(x)./HTOT);
    PO4       = (TP(x).*10.^-6)./(1+HTOT./KP3(x)+(HTOT.^2)./(KP2(x).*KP3(x))+...
                (HTOT.^3)./(KP1(x).*KP2(x).*KP3(x)));
    HPO4      = (TP(x).*10.^-6)./(1+HTOT./KP2(x)+(HTOT.^2)./(KP1(x).*KP2(x))+...
                KP3(x)./HTOT);
    SiO4H3    = ((TSI(x).*10.^-6))./(1+HTOT./KSI(x));
    NH3       = (NH4(x).*10.^-6)./(1+HTOT./KNH4(x));
    HS        = (H2S(x).*10.^-6)./(1+HTOT./KS1(x));
    HSO4      = (TS(x).*10.^-6)./(1+(KS(x)./HFr));
    HF        = (TF(x).*10.^-6)./(1+(KF(x)./HTOT));
    H3PO4     = (TP(x).*10.^-6)./(1+KP1(x)./HTOT+(KP1(x).*KP2(x))./(HTOT.^2)+...
                (KP1(x).*KP2(x).*KP3(x))./(HTOT.^3));
    if     PKORG1(x)  >= 4.5; ORGa1  = (TORG1(x).*10.^-6)./(1+HTOT./KORG1(x));
    elseif PKORG1(x)   < 4.5; ORGd1  = (TORG1(x).*10.^-6)./(1+KORG1(x)./HTOT); end
    if     PKORG2(x)  >= 4.5; ORGa2  = (TORG2(x).*10.^-6)./(1+HTOT./KORG2(x));
    elseif PKORG2(x)   < 4.5; ORGd2  = (TORG2(x).*10.^-6)./(1+KORG2(x)./HTOT); end
    if     PKORG3(x)  >= 4.5; ORGa3  = (TORG3(x).*10.^-6)./(1+HTOT./KORG3(x));
    elseif PKORG3(x)   < 4.5; ORGd3  = (TORG3(x).*10.^-6)./(1+KORG3(x)./HTOT); end
    
    % Difference between initial TA and calculated TA at each titration step
    Residual  = (TA.*10.^-6).*kgSW-HCO3.*kgSW-2.*CO3.*kgSW-BOH4.*kgSW-...
                OH.*kgTOT-2.*PO4.*kgSW-HPO4.*kgSW-SiO4H3.*kgSW-NH3.*kgSW-...
                HS.*kgSW+HFr.*kgTOT+HSO4.*kgSW+HF.*kgSW+H3PO4.*kgSW-...
                CA(x).*kgA-ORGa1.*kgSW+ORGd1.*kgSW-ORGa2.*kgSW+ORGd2.*kgSW-...
		ORGa3.*kgSW+ORGd3.*kgSW;
    
    % find Slope dTA/dpH
    Slope     = log(10).*(((DIC(x).*10.^-6).*kgSW.*KC1(x).*HTOT.*(HTOT.^2+...
                KC1(x).*KC2(x)+4.*HTOT.*KC2(x))./((HTOT.^2+KC1(x).*HTOT+...
                KC1(x).*KC2(x)).^2))+... % CAlk
                (BOH4.*HTOT.*kgSW./(HTOT+KB(x)))+... % BAlk
                ((TP(x).*10.^-6).*HTOT.*KP1(x).*kgSW.*(HTOT.^4+4.*KP2(x).*...
                HTOT.^3+(9.*KP2(x).*KP3(x)+KP1(x).*KP2(x)).*HTOT.^2+4.*...
                KP1(x).*KP2(x).*KP3(x).*HTOT+KP1(x).*KP2(x).^2.*KP3(x))./...
                (HTOT.^3+KP1(x).*HTOT.^2+KP1(x).*KP2(x).*HTOT+KP1(x).*...
                KP2(x).*KP3(x)).^2)+... % PAlk
                (SiO4H3.*HTOT.*kgSW./(HTOT+KSI(x)))+... % SiAlk
                ORGa1.*HTOT.*kgSW./(HTOT+KORG1(x))+...
                ORGd1.*HTOT.*kgSW./(HTOT+KORG1(x))+...
		ORGa2.*HTOT.*kgSW./(HTOT+KORG2(x))+...
                ORGd2.*HTOT.*kgSW./(HTOT+KORG2(x))+...
		ORGa3.*HTOT.*kgSW./(HTOT+KORG3(x))+...
                ORGd3.*HTOT.*kgSW./(HTOT+KORG3(x))+... % OrgAlk
                NH3.*HTOT.*kgSW./(HTOT+KNH4(x))+... % NH3Alk
                HS.*HTOT.*kgSW./(HTOT+KS1(x))+... % HSAlk
                HSO4.*HFr.*kgSW./(HFr+KS(x))+... % HSO4Alk
                HF.*HTOT.*kgSW./(HTOT+KF(x))+... % HFAlk
                OH.*kgTOT+HFr.*kgTOT+CA(x).*kgA); % OHAlk, HAlk, & Acid
    
    % Newton's method
    deltapH = Residual./Slope;
    % To keep the jump from being too big
    while any(abs(deltapH) > 1)
    FF = abs(deltapH)>1; deltapH(FF)=deltapH(FF)./2;
    end
    % Determine pH
    pH = pH + deltapH;
end;

% Log titration pH and estimate e.m.f.
pHtit = pH;
Etit  = 400+(Rg.*(TMP(x)+273.15)./F).*log(10.^-pHtit);

% Populate data vector with added acid, pH values, and species concentrations
DATA(S,1)=MA; DATA(S,2)=pHtit; DATA(S,3)=HTOT; DATA(S,4)=Etit;

%% CALCULATE pH DURING EACH TITRATION STEP AFTER PURGING OF CO2
% Declare initial delta pH
deltapH = pHTol+1;
% Iterate to determine stepwise purged pH until deltapH is below tolerance limit
while any(abs(deltapH) > pHTol)  
    % Concentrations of all species relevant to TA
    % (carbonate species have been bubbled out)
    HTOTp      = (10.^(-pH));
    HFrp       = HTOTp./Z(x);
    BOH4p      = (TB(x).*10.^-6)./(1+HTOTp./KB(x));
    OHp        = (KW(x)./HTOTp);
    PO4p       = (TP(x).*10.^-6)./(1+HTOTp./KP3(x)+(HTOTp.^2)./(KP2(x).*KP3(x))+...
                 (HTOTp.^3)./(KP1(x).*KP2(x).*KP3(x)));
    HPO4p      = (TP(x).*10.^-6)./(1+HTOTp./KP2(x)+(HTOTp.^2)./(KP1(x).*KP2(x))+...
                  KP3(x)./HTOTp);
    SiO4H3p    = ((TSI(x).*10.^-6))./(1+HTOTp./KSI(x));
    NH3p       = (NH4(x).*10.^-6)./(1+HTOTp./KNH4(x));
    HSp        = (H2S(x).*10.^-6)./(1+HTOTp./KS1(x));
    HSO4p      = (TS(x).*10.^-6)./(1+(KS(x)./HFrp));
    HFp        = (TF(x).*10.^-6)./(1+(KF(x)./HTOTp));
    H3PO4p     = (TP(x).*10.^-6)./(1+KP1(x)./HTOTp+(KP1(x).*KP2(x))./(HTOTp.^2)+...
                 (KP1(x).*KP2(x).*KP3(x))./(HTOTp.^3));
    if     PKORG1(x)  >= 4.5; ORGap1  = (TORG1(x).*10.^-6)./(1+HTOTp./KORG1(x));
    elseif PKORG1(x)   < 4.5; ORGdp1  = (TORG1(x).*10.^-6)./(1+KORG1(x)./HTOTp); end
    if     PKORG2(x)  >= 4.5; ORGap2  = (TORG2(x).*10.^-6)./(1+HTOTp./KORG2(x));
    elseif PKORG2(x)   < 4.5; ORGdp2  = (TORG2(x).*10.^-6)./(1+KORG2(x)./HTOTp); end
    if     PKORG3(x)  >= 4.5; ORGap3  = (TORG3(x).*10.^-6)./(1+HTOTp./KORG3(x));
    elseif PKORG3(x)   < 4.5; ORGdp3  = (TORG3(x).*10.^-6)./(1+KORG3(x)./HTOTp); end
    
    % Difference between initial TA and calculated TA at each titration step
    Residual  = (TA.*10.^-6).*kgSW-BOH4p.*kgSW-OHp.*kgTOT-2.*PO4p.*kgSW-...
                HPO4p.*kgSW-SiO4H3p.*kgSW-NH3p.*kgSW-HSp.*kgSW+HFrp.*kgTOT+...
                HSO4p.*kgSW+HFp.*kgSW+H3PO4p.*kgSW-CA(x).*kgA-ORGap1.*kgSW+...
                ORGdp1.*kgSW-ORGap2.*kgSW+ORGdp2.*kgSW-ORGap3.*kgSW+ORGdp3.*kgSW;

    % find Slope dTA/dpH
    Slope     = log(10).*(BOH4p.*HTOTp.*kgSW./(HTOTp+KB(x))+... % BAlk
                ((TP(x).*10.^-6).*HTOTp.*KP1(x).*kgSW.*(HTOTp.^4+4.*KP2(x).*...
                HTOTp.^3+(9.*KP2(x).*KP3(x)+KP1(x).*KP2(x)).*HTOTp.^2+4.*...
                KP1(x).*KP2(x).*KP3(x).*HTOTp+KP1(x).*KP2(x).^2.*KP3(x))./...
                (HTOTp.^3+KP1(x).*HTOTp.^2+KP1(x).*KP2(x).*HTOTp+KP1(x).*...
                KP2(x).*KP3(x)).^2)+... % PAlk
                SiO4H3p.*HTOTp.*kgSW./(HTOTp+KSI(x))+... % SiAlk
                ORGap1.*HTOTp.*kgSW./(HTOTp+KORG1(x))+...
                ORGdp1.*HTOTp.*kgSW./(HTOTp+KORG1(x))+...
		ORGap2.*HTOTp.*kgSW./(HTOTp+KORG2(x))+...
                ORGdp2.*HTOTp.*kgSW./(HTOTp+KORG2(x))+...
		ORGap3.*HTOTp.*kgSW./(HTOTp+KORG3(x))+...
                ORGdp3.*HTOTp.*kgSW./(HTOTp+KORG3(x))+... % OrgAlk
                NH3p.*HTOTp.*kgSW./(HTOTp+KNH4(x))+... % NH3Alk
                HSp.*HTOTp.*kgSW./(HTOTp+KS1(x))+... % HSAlk
                HSO4p.*HFrp.*kgSW./(HFrp+KS(x))+... % HSO4Alk
                HFp.*HTOTp.*kgSW./(HTOTp+KF(x))+... % HFAlk
                OHp.*kgTOT+HFrp.*kgTOT+CA(x).*kgA); % OHAlk, HAlk, & Acid
 
    % Newton's method
    deltapH = Residual./Slope;
    % To keep the jump from being too big
    while any(abs(deltapH) > 1)
    FF=abs(deltapH)>1; deltapH(FF)=deltapH(FF)./2;
    end
    % Determine total scale pH
    pH = pH + deltapH;
end
% Log purged pH and estimate purged pH
pHp = pH;
Ep  = 400+(Rg.*(TMP(x)+273.15)./F).*log(10.^-pHp);

% Populate vector with added acid, pH values, and species concentrations
DATAp(S,1)=MA; DATAp(S,2)=pHp; DATAp(S,3)=HTOTp; DATAp(S,4)=Ep;

S=S+1;

end

%% EXTRACT ACID MASS, pH, H+, AND E FROM 'DATA' FOR CALCULATIONS
MA  =  DATA(:,1);
kgA =  DATA(:,1).*10.^-3;
pH  =  DATA(:,2);
pHp =  DATAp(:,2);
H   =  DATA(:,3);
Hp  =  DATAp(:,3);
E   =  DATA(:,4);
Ep  =  DATAp(:,4);

%% CALCULATE ALKALINITY USING THE SECOND GRAN FUNCTION (OPEN-CELL)
% Index to a pH Range of 3 to 3.5
GRidx = pHp > 3.0 & pHp < 3.5;
Ep_GR  = Ep(GRidx); kgA_GR = kgA(GRidx);
% Obtain initial estimates of F2, TA, E0
F_est_GR     = (kgA_GR+kgSW).*exp(Ep_GR./((Rg.*(TMP(x)+273.15))./F));
A0_est_GR    = fitlm(F_est_GR,kgA_GR);
A0_est_GR    = (A0_est_GR.Coefficients.Estimate(1).*CA(x))./kgSW;
Ep0_est_GR   = Ep_GR - (((Rg.*(TMP(x)+273.15))./F).*...
               log((kgA_GR.*CA(x)-kgSW.*A0_est_GR)./(kgA_GR+kgSW)));
Ep0_est_GR   = mean(Ep0_est_GR);
% Set initial values for E0 and TA
Ep0      = Ep0_est_GR;
TA       = A0_est_GR;
deltaTA  = 1;
% Iteratively determine TA
while any(abs(deltaTA) > 0.00001)
   HTOTg  = exp((Ep_GR - Ep0)./(Rg.*(TMP(x)+273.15)./F));
   HFrg   = HTOTg./Z(x);
   HSO4g  = ((TS(x).*10.^-6)./(1+KS(x)./HFrg));
   HFg    = ((TF(x).*10.^-6)./(1+KF(x)./HTOTg));
   H3PO4g = ((TP(x).*10.^-6)./(1+KP1(x)./HTOTg+(KP1(x).*KP2(x))./(HTOTg.^2)+...
            (KP1(x).*KP2(x).*KP3(x))./(HTOTg.^3)));
   F2     = (kgSW+kgA_GR).*HFrg+kgSW.*(HSO4g+HFg+H3PO4g);
   
   TA0     = TA;
   fit2    = regstats(F2,kgA_GR,'linear','beta');
   TA      = (-fit2.beta(1)./fit2.beta(2)).*(CA(x)./kgSW);
   deltaTA = TA - TA0;
   
   deltaE0 = Ep0 - (Ep_GR-((Rg.*(TMP(x)+273.15))./F).*...
             log((kgA_GR.*CA(x)-TA.*kgSW-...
             kgSW.*(HSO4g+HFg+H3PO4g))./...
             (kgA_GR+kgSW))); deltaE0 = mean(deltaE0);
   Ep0     = Ep0-deltaE0;
end

TA_Gran = TA.*10.^6;

%% CALCULATE ALKALINITY (AND TOTAL CARBON) USING NON-LINEAR CURVE FITTING (CLOSED-CELL)
% Index to low-pH range for initial estimates
CLidx = pH > 3 & pH < 3.5;
E_CL  = E(CLidx); kgA_CL = kgA(CLidx);
% Obtain initial estimates for alkalinity and E0
F_est_CL    = (kgA_CL+kgSW).*exp(E_CL./((Rg.*(TMP(x)+273.15))./F));
A0_est_CL   = fitlm(F_est_CL,kgA_CL);
A0_est_CL   = (A0_est_CL.Coefficients.Estimate(1).*CA(x))./kgSW;
E0_est_CL   = E_CL - (((Rg.*(TMP(x)+273.15))./F).*...
              log((kgA_CL.*CA(x)-kgSW.*A0_est_CL)./(kgA_CL+kgSW)));
E0_est_CL   = mean(E0_est_CL);
% Obtain initial estimate for H
HTOT_est_CL = exp((E-E0_est_CL)./((Rg.*(TMP(x)+273.15))./F));
f           = 1;
% Fit full titration curve
nlin0 = [A0_est_CL 0.002 f];
Eq = @(nlin,w)(nlin(1,1)-... % TA
   nlin(1,2).*((KC1(x).*nlin(1,3).*HTOT_est_CL+2.*KC1(x).*KC2(x))./...
   ((nlin(1,3).*HTOT_est_CL).^2+KC1(x).*nlin(1,3).*HTOT_est_CL+KC1(x).*KC2(x)))-... % CAlk
   (TB(x).*10.^-6)./(1+(nlin(1,3).*HTOT_est_CL)./KB(x))-... % BAlk
   (TSI(x).*10.^-6)./(1+(nlin(1,3).*HTOT_est_CL)./KSI(x))-... % SiAlk
   TP(x).*10.^-6.*((KP1(x).*KP2(x).*nlin(1,3).*HTOT_est_CL+...
   2.*KP1(x).*KP2(x).*KP3(x)-(nlin(1,3).*HTOT_est_CL).^3)./...
   ((nlin(1,3).*HTOT_est_CL).^3+KP1(x).*(nlin(1,3).*HTOT_est_CL).^2+...
   KP1(x).*KP2(x).*(nlin(1,3).*HTOT_est_CL)+KP1(x).*KP2(x).*KP3(x)))+... % PAlk
   (TS(x).*10.^-6)./(1+KS(x)./((nlin(1,3).*HTOT_est_CL)./Z(x)))+... % HSO4
   (TF(x).*10.^-6)./(1+KF(x)./(nlin(1,3).*HTOT_est_CL))+... % HF
   ((kgSW+kgA)./kgSW).*(((nlin(1,3).*HTOT_est_CL)./Z(x))-... % HFr
   (KW(x)./(nlin(1,3).*HTOT_est_CL)))-... % OH
   (kgA./kgSW).*CA(x));
[nlin1]  = nlinfit(HTOT_est_CL,zeros(size(HTOT_est_CL,1),1),Eq,nlin0);

TA_nlin1  = nlin1(1,1).*10.^6;

%% CALCULATE ALKALINITY USING NON-LINEAR CURVE FITTING ACROSS GRAN RANGE (OPEN-CELL)
% Index to low-pH range
OPidx = pHp > 3.0 & pHp < 3.5;
Ep_OP    = Ep(OPidx); kgA_OP = kgA(OPidx);
% Obtain initial estimates for alkalinity and E0
F_est_OP    = (kgA_OP+kgSW).*exp(Ep_OP./((Rg.*(TMP(x)+273.15))./F));
A0_est_OP   = fitlm(F_est_OP,kgA_OP);
A0_est_OP   = (A0_est_OP.Coefficients.Estimate(1).*CA(x))./kgSW;
Ep0_est_OP  = Ep_OP - (((Rg.*(TMP(x)+273.15))./F).*...
              log((kgA_OP.*CA(x)-kgSW.*A0_est_OP)./(kgA_OP+kgSW)));
Ep0_est_OP  = mean(Ep0_est_OP);
% Obtain initial estimate for H
HTOT_est_OP = exp((Ep_OP-Ep0_est_OP)./((Rg.*(TMP(x)+273.15))./F));
f           = 1;
% Fit open titration curve
nlin0 = [A0_est_OP f];
Eq2 = @(nlin,w)(nlin(1,1)+... % TA
       (TS(x).*10.^-6)./(1+(KS(x)./((nlin(1,2).*HTOT_est_OP)./Z(x))))+... % HSO4
       ((TF(x).*10.^-6)./(1+(KF(x)./(nlin(1,2).*HTOT_est_OP))))+... % HF
       ((kgSW+kgA_OP)./kgSW).*(((nlin(1,2).*HTOT_est_OP)./Z(x)))-... % HFr
       (kgA_OP./kgSW).*CA(x));
[nlin2] = nlinfit(HTOT_est_OP,zeros(size(HTOT_est_OP,1),1),Eq2,nlin0);

TA_nlin2 = nlin2(1,1).*10.^6;

%% CALCULATE ALKALINITY USING DIFFERENCE DERIVATIVE
Diffs   = -diff(E)./diff(MA);
kgA0    = kgA(2:end);
sp      = spline(kgA0,Diffs);
deriv   = fnder(sp);
zer2    = fnzeros(deriv);
TA_Diff = (10.^6).*((zer2(end,end).*CA(x))./kgSW);

%% CALCULATE ALKALINITY USING SINGLE-STEP METHOD
% Index to pH endpoint which will be used to minimize excess acid
spH       = spline(pH,kgA);
kgAi      = ppval(spH,4.2);
sHp       = spline(pH,Hp);
Hpi       = ppval(sHp,4.2);
% Calculate TA
TA_YaoByrne = (10.^6).*((CA(x).*kgAi-...                            % HCl
              ((TS(x).*10.^-6)./(1+(KS(x)./(Hpi./Z(x))))).*kgSW-... % HSO4
              ((TF(x).*10.^-6)./(1+(KF(x)./Hpi))).*kgSW-...         % HF
              (Hpi./Z(x)).*(kgAi+kgSW))...                    % HFr
               ./kgSW);

%% LOG THE RESULTS
Result(x,:) = [TA_inorg,TA_org,TA_Gran,TA_nlin1,TA_nlin2,TA_Diff,TA_YaoByrne];

end

Result_Headers = {'Inorg. TA Input','Org. TA Input','Open-Cell Gran TA','Closed-Cell Nonlinear Fit TA',...
                  'Open-Cell Nonlinear Fit TA','Difference Derivative TA','Single-Step TA'};
