function r = erfcxf_implementation(x)

    a = max(x, single(0) - x);
    
    m = a - single(2);
    p = a + single(2);
    r = single(1)./p;
    q = m.*r;
    
    t = fma(q+single(1), -single(2), a);
    e = fma(q, -a, t);
    q = fma(r, e, q);

    p =           single( 5.92470169e-5);
    p = fma(p, q, single( 1.61224554e-4));
    p = fma(p, q, single(-3.46481771e-4));
    p = fma(p, q, single(-1.39681227e-3));
    p = fma(p, q, single( 1.20588380e-3));
    p = fma(p, q, single( 8.69014394e-3));
    p = fma(p, q, single(-8.01387429e-3));
    p = fma(p, q, single(-5.42122945e-2));
    p = fma(p, q, single( 1.64048523e-1));
    p = fma(p, q, single(-1.66031078e-1));
    p = fma(p, q, single(-9.27637145e-2));
    p = fma(p, q, single( 2.76978403e-1));

    d = a + single(0.5);
    r = single(1) ./ d;
    r = r .* single(0.5);
    q = fma(p, r, r); 
    t = q + q;
    e = (p - q) + fma(t, -a, single(1));
    r = fma(e, r, q);

    r(a > realmax('single')) = single(0);
        
    k = x < single(0);
    if any(k)
        s = x(k).*x(k);
        d = fma(x(k), x(k), -s);
        e = exp(s);
        R = e - r(k);
        R = fma(e, d + d, R);
        R = R + e;
        R(e > realmax('single')) = e(e > realmax('single'));
        r(k) = R;
    end
end

% function f = fma(a,b,c)
%     f = a.*b + c;
% end