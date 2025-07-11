clear all;
clc;
%% part i
x = [0 0.02 0.04 0.06 0.08 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1.0];
y = [0 0.134 0.23 0.304 0.365 0.418 0.579 0.665 0.729 0.779 0.825 0.87 0.915 0.958 0.979 1];
func = @(variables,xdata) variables(1).*xdata./(1+variables(2).*xdata+variables(3).*xdata.^2);
initial_guess = [1,1,1];
variables = lsqcurvefit(func,initial_guess,x,y);
figure(1);
plot(x,y,'r-.');
hold on;
plot(x,x,'g');
grid on;
plot(x,func(variables,x),'b--');
legend('given data','y=x line','fitted data',Location='best');
func_eqm = @(x) variables(1).*x./(1+variables(2).*x+variables(3).*x.^2);
%% part iii
syms D W
F = 500;
S = 50;
eqns = [F == D+W+S;0.45*F == 0.97*D+0.02*W+0.7*S];
D_and_W = solve(eqns,[D W]);
%% part iv-viii
figure(2);
p1 = plot(x,x,'cyan','linewidth',1.1);
hold on;
grid on;
x_eqm = linspace(0,1,100);
func1 = @(xdata) variables(1).*xdata./(1+variables(2).*xdata+variables(3).*xdata.^2);
p2 = plot(x_eqm,func1(x_eqm),'g','linewidth',1.1);
xd = 0.97;
xw = 0.02;
xf = 0.45;
xs = 0.70;
q = 0.8;
scatter(xd,xd,'black','o','filled');
scatter(xw,xw,'black','o','filled');
scatter(xs,xs,'black','o','filled');
scatter(xf,xf,'black','o','filled');
func2 = @(x) xf+(x-xf).*(q/(q-1));
x_feed = [xf,xf-0.1];
y_feed = func2(x_feed);
p3 = plot(x_feed,y_feed);
x_side = [xs,xs];
y_side = [xs,xs+0.25];
p7 = plot(x_side,y_side);
pinch_point = func(variables,x(12));
scatter(x(12),pinch_point,'black','o','filled');
pp_slope = (pinch_point-xd)/(xs-xd);
R_min = pp_slope/(1-pp_slope);
R_actual = 2.5*R_min;
slope_actual = R_actual/(1+R_actual);
start_point = xd/(1+R_actual);
%% section-1
L = vpa(D_and_W.D*R_actual);
V = vpa((R_actual+1)*D_and_W.D);
%% section-2
L_sec2 = L-S;
V_sec2 = V;
slope_sideline = L_sec2/V_sec2;
intersection_X1 = xs;
intersection_Y1 = xd+slope_actual*(intersection_X1-xd);
op_line_RS_x = [intersection_X1 xd];
op_line_RS_y = [intersection_Y1 xd];
p5 = plot(op_line_RS_x,op_line_RS_y,'magenta');
dottedX = [0 intersection_X1];
func_RS = @(x) R_actual*x/(R_actual+1)+xd/(R_actual+1);
dottedY = [start_point func_RS(dottedX(2))];
plot(dottedX,dottedY,'m--');
func_1 = @(x) (x-intersection_X1)*slope_sideline+intersection_Y1;
op_sideline_x = [xs-0.31 xs];
op_sideline_y = [func_1(op_sideline_x(1)) func_1(op_sideline_x(2))];
p4 = plot(op_sideline_x,op_sideline_y,'black');
%% section-3
L_bar = 0.8*F+L_sec2;
V_bar = V_sec2+0.2*F;
syms x2
eqn = (x2-intersection_X1)*slope_sideline+intersection_Y1 == xf+(x2-xf).*(q/(q-1));
intersection_X2 = solve(eqn,x2);
intersection_Y2 = func2(intersection_X2);
op_line_SS_x = [xw intersection_X2];
op_line_SS_y = [xw intersection_Y2];
p6 = plot(op_line_SS_x,op_line_SS_y,'red');
func3 = @(x) (L_bar/vpa(L_bar-D_and_W.W)).*x - vpa(D_and_W.W/(L_bar-D_and_W.W)).*xw;
dotted2x = [intersection_X2, xw+0.55];
dotted2y = func3(dotted2x);
plot(dotted2x,dotted2y,'red--');
dotted3x = [0, intersection_X2];
dotted3y = func_1(dotted3x);
plot(dotted3x,dotted3y,'black--');
slope_SS = (op_line_SS_y(2)-op_line_SS_y(1))/(op_line_SS_x(2)-op_line_SS_x(1));
%% part ix (McCabe-Thiele Construction)
% till intersection_X1
nt = -1;
Y = xd;
X = xd;
func3 = @(x) (L_bar/vpa(L_bar-D_and_W.W)).*x - vpa(D_and_W.W/(L_bar-D_and_W.W)).*xw;
while X>=intersection_X1
    Xold = X;
    Yold = Y;
    func4 = @(X) Yold-((variables(1)*X)/(1+variables(2)*X+variables(3)*X^2));
    X = fsolve(func4,0.9);
    Y = Yold + (X-Xold)*slope_actual;
    nt = nt+1;
    line([Xold X],[Yold Yold],'color','#A0F','linewidth',1.1);
    line([X X],[Yold Y],'color','#A0F','linewidth',1.1);
end
% till intersection_X2
while X>=intersection_X2
    Xold = X;
    Yold = Y;
    syms X
    func4 = ((variables(1)*X)/(1+variables(2)*X+variables(3)*X^2)) == Y;
    X = vpasolve(func4,X,[0 0.7]);
    Y = Yold + (X-Xold)*slope_sideline;
    nt = nt+1;
    line([Xold X],[Yold Yold],'color','#A0F','linewidth',1.1);
    line([X X],[Yold Y],'color','#A0F','linewidth',1.1);
end
line([X X],[Yold func3(X)],'color','#A0F','linewidth',1.1);
Y = func3(X);
% till xw
while X>=xw
    Xold = X;
    Yold = Y;
    syms X;
    func4 = ((variables(1)*X)/(1+variables(2)*X+variables(3)*X^2)) == Y;
    X = (vpasolve(func4,X,[0 0.26]));

    Y = Yold + (X-Xold)*slope_SS;
    nt = nt+1;
    line([Xold X],[Yold Yold],'color','#A0F','linewidth',1.1);
    if X>=0 && Y>=xw
        line([X X],[Yold Y],'color','#A0F','linewidth',1.1);
    else
        line([X X],[Yold 0.005],'color','#A0F','linewidth',1.1);
    end
end
%% calculations for ideal number of trays
nt = round((nt+(Yold-xw)/(Yold-Y)),4)
%% by observation:
% number of ideal stages = 9.8147 trays
% feed introduced in 7th tray
% sideline removed from 4th tray
%% plotting codes
scatter(intersection_X1,intersection_Y1,'black','o','filled');
scatter(intersection_X2,intersection_Y2,'black','o','filled');
line([xd xd],[0 xd],'linestyle','-.','color','black');
line([xw xw],[0 xw],'linestyle','-.','color','black');
line([xf xf],[0 xf],'linestyle','-.','color','black');
line([xs xs],[0 xs],'linestyle','-.','color','black');
line([intersection_X2 intersection_X2],[0 intersection_Y2],'linestyle','-.','color','black');
text(xw+0.015,0.02,'boiler fraction',Rotation=0);
text(xf+0.015,0.02,'feed fraction',Rotation=90);
text(xd-0.045,0.02,'distillate fraction',Rotation=90);
text(xs+0.035,0.02,'sideline removal',Rotation=90);
text(intersection_X2-0.025,0.02,'feed introduction',Rotation=90);
legend([p1 p2 p3 p4 p5 p6 p7],{'y=x line','fitted data','feed line','section 2','op. line RS','op. line SS','sideline'},Location='northwest');

