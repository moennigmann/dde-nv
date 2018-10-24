close all
clear
clc


open('crit_Manifold_3L_fixed_K_clear.fig')
box on
xlabel('p_1')
ylabel('p_2')
zlabel('p_3')
view([-35, 10])
drawnow

matlab2tikz('3LasymManis.tex')