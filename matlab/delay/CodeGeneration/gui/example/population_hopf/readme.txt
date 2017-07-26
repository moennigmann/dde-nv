Here are the files which were  generated with the GUI. This exmple is for the popualtion model.
The inputs were as follows:
name:
  population_model_2
path:
  /home/users/pme/Downloads/MapleModules/
number of states:
  2
State names:
  xJuv
  xMat
number delaytimes:
  1
expression for delay:
  1.5-0.5*exp(-xMat)
number of parameter:
  2
parameter:
  alpha	5
  gamma	2
eqautions:
  alpha*xMat-gamma*xJuv-alpha*exp(-gamma*(tau[1]))*xMattau1
  alpha*exp(-gamma*(tau[1]))*xMattau1-beta*xMat^2
additional explicit AEs:
  p1=1
  beta=p1
Manifold:
  hopf