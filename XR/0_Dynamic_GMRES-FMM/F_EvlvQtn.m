function [QtnEps1,QtnEps2,QtnEps3,QtnEta, Eps1rateold,Eps2rateold,...
    Eps3rateold,Etarateold] = F_EvlvQtn(dt,QtnEps1,QtnEps2,QtnEps3,QtnEta,...
    Eps1rateold,Eps2rateold,Eps3rateold,Etarateold,rot_pf)
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
    Eps1rate = 0.5 * (QtnEta.*rot_pf(:,1)-QtnEps3.*rot_pf(:,2)+...
        QtnEps2.*rot_pf(:,3));
    Eps2rate = 0.5 * (QtnEps3.*rot_pf(:,1)+QtnEta.*rot_pf(:,2)-...
        QtnEps1.*rot_pf(:,3));
    Eps3rate = 0.5 * (-QtnEps2.*rot_pf(:,1)+QtnEps1.*rot_pf(:,2)+...
        QtnEta.*rot_pf(:,3));
    Etarate = 0.5 * (-QtnEps1.*rot_pf(:,1)-QtnEps2.*rot_pf(:,2)-...
        QtnEps3.*rot_pf(:,3));
    
    deps1 = 1.5*Eps1rate-0.5*Eps1rateold;
    deps2 = 1.5*Eps2rate-0.5*Eps2rateold;
    deps3 = 1.5*Eps3rate-0.5*Eps3rateold;
    deta = 1.5*Etarate-0.5*Etarateold;
    
    QtnEps1 = QtnEps1 + deps1*dt;
    QtnEps2 = QtnEps2 + deps2*dt;
    QtnEps3 = QtnEps3 + deps3*dt;
    QtnEta = QtnEta + deta*dt;
end
