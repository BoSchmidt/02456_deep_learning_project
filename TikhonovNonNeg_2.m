function [xalpha, alpha] = TikhonovNonNeg(transfer_matrix,G,vpara,vperp,L)
ax=[-10.5 10.5 -0.5 20.5];
n = size(transfer_matrix,2);
L0 = eye(n);
[L1vpara,L1vperp] = gradient_v_space_matrix(vpara,vperp,'custom');
scaling_factor = 1/max(max(transfer_matrix));
transfer_matrix = transfer_matrix*scaling_factor;

switch L
    case 0
        H = L0'*L0;
        % alpha = logspace(-3,3,30)';
        % alpha = linspace(0.0035, 0.35, 10);
        % alpha = 0.05;
        % alpha = linspace(0.001, 0.1, 50);
        % alpha = logspace(-4,-1,16)'
        % alpha = vect;
        alpha = vect(itime);
        norm_operator = 1;
    case 1
        H = L1vperp'*L1vperp + L1vpara'*L1vpara;
        alpha=([1e-3 3e-3 1e-2 3e-2 1e-1 3e-1 1 3 10]').^2;
    otherwise
        error('L must be 0 or 1 (0th or 1st order penalty function)')
end

xalpha = zeros(n,length(alpha));

for i = 1:length(alpha)
   i
   switch L
       case 0
           GalphaL0=zeros(size(L0,2),1);
           WalphaL=double([transfer_matrix; sqrt(alpha(i))*L0]);
           GalphaL=double([G; GalphaL0]);
           xalpha(:,i) = (transfer_matrix'*transfer_matrix + alpha(i)*H)\transfer_matrix'*G;
       case 1 
           GalphaL1=zeros(2*size(L1vpara,2),1);
           WalphaL=double([transfer_matrix; sqrt(alpha(i))*L1vpara; sqrt(alpha(i))*L1vperp]);
           GalphaL=double([G; GalphaL1]);
           xalpha(:,i) = lsqnonneg(WalphaL,GalphaL);
   end
end  

xi = zeros(length(alpha),1);
xalpha_norm = zeros(length(alpha),1);
xi_grad_test = zeros(length(alpha),1);
rho_grad_test = zeros(length(alpha),1);
f= zeros(length(alpha),1);
rho= zeros(length(alpha),1);

for i=1:length(alpha)
 f(i) = norm(transfer_matrix*xalpha(:,i) - G);
    if L == 0
        xalpha_norm(i) = norm(xalpha(:,i));
        xi(i) = xalpha_norm(i).^2;
    elseif L ==1
        xalpha_norm(i) = 1/sqrt(2)*norm(L1E*xalpha(:,i)) + 1/sqrt(2)*norm(L1p*xalpha(:,i));
        xi(i) = norm(L1E*xalpha(:,i))^2 + norm(L1p*xalpha(:,i))^2;
    end
end

rho = f.^2;
xi_grad_test = zeros(length(alpha),1);
rho_grad_test = zeros(length(alpha),1);

for i = 2:length(alpha)-1
    xi_grad_test(i) = (xi(i+1)-xi(i-1))/(sqrt(alpha(i+1))-sqrt(alpha(i-1)));
    rho_grad_test(i) = (rho(i+1)-rho(i-1))/(sqrt(alpha(i+1))-sqrt(alpha(i-1)));
end

curvature = -2*xi.*rho./xi_grad_test.*(alpha.*xi_grad_test.*rho + 2*sqrt(alpha).*xi.*rho + alpha.^2.*xi.*xi_grad_test)./(alpha.^2.*xi.^2 + rho.^2).^(3/2);
curvature_MY_OWN = (alpha.*rho.^2.*xi + 2*sqrt(alpha).*rho.^2.*xi.^2./xi_grad_test + alpha.^2.*rho.*xi.^2)./(alpha.^2.*xi.^2 + rho.^2).^(3/2);

L_curve_index = find(max(curvature) == curvature);

if 0
figure(100);clf;hold on;
subplot(1,2,1)
loglog(f.^2,xalpha_norm.^2,'s')
hold all
loglog(f(L_curve_index).^2,xalpha_norm(L_curve_index).^2,'s','linewidth',2)
subplot(1,2,2)
semilogx(alpha,curvature,'s')
end

xalpha = xalpha*scaling_factor;