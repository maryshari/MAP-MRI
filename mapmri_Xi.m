function Xi = mapmri_Xi(J,L,M,u0,q,bvecs)

mag_q = sqrt(q(1,:).^2 + q(2,:).^2 + q(3,:).^2);
mag_q = mag_q';

% mag_q1 = mag_q;
% mag_q1(mag_q1 == 0) = 1;
% bvecs = [q(1,:)'./mag_q1 q(2,:)'./mag_q1 q(3,:)'./mag_q1];

Xi = sqrt(4*pi)*1i^(-L)*(2*pi^2*u0^2*mag_q.^2).^(L/2).*exp(-2*pi^2*u0^2*mag_q.^2).*...
                    laguerreL(J-1,(L+1/2),4*pi^2*u0^2*mag_q.^2).*spharm(L,M,bvecs');

end



