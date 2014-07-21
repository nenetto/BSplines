function y = oMom3(x)

a = abs(x);
y(a<1) = 0.5*a(a<1).^3 -a(a<1).^2 + a(a<1)/14 + 13/21;
y(and((a>=1),(a<2))) = -(1/6)*a(and((a>=1),(a<2))).^3 + a(and((a>=1),(a<2))).^2 -(85/42)*a(and((a>=1),(a<2))) + 29/21;
y(a>=2) = 0;

end