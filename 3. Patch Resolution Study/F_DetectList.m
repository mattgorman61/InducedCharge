function [DetectList] = F_DetectList(npars,rpar,a,b,c)
%F_COLLIDELIST �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    Distance = pdist(rpar); % Interparticle distance
    DistMatrx = squareform(Distance); % Distance matrix
    
    AxisLong = max(a,max(b,c)); % Use longest axis as criterion
    DetectList = (DistMatrx <= (3*AxisLong)); % Label close pair
    
    DetectList = DetectList - eye(npars); % Remove particle itself
end

