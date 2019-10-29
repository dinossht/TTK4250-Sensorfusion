function eul = quat2eul(q)
% Equation (10.38) & (10.39)
    qSquared = q.^2;

    phi = atan2(2*(q(4)*q(3)+q(1)*q(2)), qSquared(1)-qSquared(2)-qSquared(3)+qSquared(4));
    theta = asin(2*(q(1)*q(3)-q(2)*q(4)));
    psi = atan2(2*(q(2)*q(3)+q(1)*q(4)), qSquared(1)+qSquared(2)-qSquared(3)-qSquared(4));

    eul = [phi; theta; psi];
end