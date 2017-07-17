Here are the files which were generated with the GUI. This exmple is for the oneCSTR model.
The inputs were as follows:
name:
  oneCSTR
path:
  /home/users/pme/Downloads/MapleModules/
number of states:
  6
State names:
  cA1
  cB1
  T1 
  cAF
  cBF
  TF
number delaytimes:
  2
expression for delay:
  5*350/(Fr)
  5*350/(Fr+F01)
number of parameter:
  8
parameter:
  Fr  	50
  beta	45
  T01	400
  F01	60
  Tc1	480
  V1	1000
  TcF	480
  VF	1000
eqautions:
  (F01/3600)/(V1/1000)*(1-cA1)+(Fr/3600)/(V1/1000)*(alphaA*cAFtau1/(alphaA*cAFtau1+alphaB*cBFtau1+alphaC*(1-cAFtau1-cBFtau1))-cA1)-kA*exp(-EA/(R*T1))*cA1
  (F01/3600)/(V1/1000)*(0-cB1)+(Fr/3600)/(V1/1000)*(alphaB*cBFtau1/(alphaA*cAFtau1+alphaB*cBFtau1+alphaC*(1-cAFtau1-cBFtau1))-cB1)+kA*exp(-EA/(R*T1))*cA1-kB*exp(-EB/(R*T1))*cB1
  (F01/3600)/(V1/1000)*(T01-T1)+(Fr/3600)/(V1/1000)*(TFtau1-T1)+BA*kA*exp(-EA/(R*T1))*cA1+BB*kB*exp(-EB/(R*T1))*cB1-(beta/10)/((V1/1000)^(1/3))*(T1-Tc1)
  ((F01/3600)+(Fr/3600))/(VF/1000)*(cA1tau2-cAF)-((Fr/3600)+(Fp/3600))/(VF/1000)*(alphaA*cAF/(alphaA*cAF+alphaB*cBF+alphaC*(1-cAF-cBF))-cAF)
  ((F01/3600)+(Fr/3600))/(VF/1000)*(cB1tau2-cBF)-((Fr/3600)+(Fp/3600))/(VF/1000)*(alphaB*cBF/(alphaA*cAF+alphaB*cBF+alphaC*(1-cAF-cBF))-cBF)
  ((F01/3600)+(Fr/3600))/(VF/1000)*(T1tau2-TF)-(beta/10)/((VF/1000)^(1/3))*(TF-TcF)
additional explicit AEs:
  p1=3.5, 
  p2=1, 
  p3=0.5, 
  p4=50000, 
  p5=60000, 
  p6=4.2, 
  p7=-60000, 
  p8=-70000, 
  p9=14285.7143, 
  p10=16666.6667, 
  p11=8.314, 
  p12=2770, 
  p13=2500, 
  p14=0, 
  alphaA=p1, 
  alphaB=p2, 
  alphaC=p3, 
  EA=p4, 
  EB=p5, 
  Cp=p6, 
  dHA=p7, 
  dHB=p8, 
  BA=p9, 
  BB=p10, 
  R=p11, 
  kA=p12, 
  kB=p13, 
  Fp=p14
Manifold:
  fold