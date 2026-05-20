function gp = gbern_upper2lower(f, g, solseries, fUHP, pUHP)
    
    B = solseries.B;
    mu = solseries.mu;    
    
    % Note if we desire f with real(f) > 0, we need interpolated values in
    % real(-f) < 0 and imag(-f) > 0
    if real(f) > 0
        p = interp2(real(fUHP.left), imag(fUHP.left), pUHP.left, -real(f), -imag(f), 'spline');
    elseif real(f) <= 0
        p = interp2(real(fUHP.right), imag(fUHP.right), pUHP.right, -real(f), -imag(f), 'spline');
    end
    
%     minusf = -f;
%     if real(minusf) < 0
%         [~,indphi] = min(abs(real(fUHP.left(1,:)) - real(minusf)));
%         p = interp1(imag(fUHP.left(:,1)), pUHP.left(:,indphi), imag(minusf));        
%     elseif real(minusf) >= 0
%         [~,indphi] = min(abs(real(fUHP.right(1,:)) - real(minusf)));
%         p = interp1(imag(fUHP.right(:,1)), pUHP.right(:,indphi), imag(minusf));
%     end
    
    gp = mu./(2*p.*(mu*B + 1i*pi*g)) - p;  
%     disp(['f = ', num2str(f), ', gp = ', num2str(gp)]);
end