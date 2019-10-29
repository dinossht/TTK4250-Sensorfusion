function qprod = quatProd(ql, qr)
% Equation (10.34)

    if numel(ql) == 3 % assume pure quat
        ql = [0; ql];
    end
    
    if numel(qr) == 3 % assume pure quat
        qr = [0; qr];
    end
    
    epsl = ql(2:4);
    qprod = (ql(1)*eye(4) + [0,-epsl';epsl,crossProdMat(epsl)]) * qr;
end