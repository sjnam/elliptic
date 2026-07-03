% ecfig.mp -- 타원 곡선 y^2 = x^3 - 2x + 1 위의 현(덧셈)과 접선(두 배) 그림.
% mpost ecfig.mp 하면 ecfig.1(그림 1), ecfig.2(그림 2)가 나온다.

u := 30; % 단위 길이(bp)

vardef yy(expr x) = sqrt(max(0, x*x*x - 2*x + 1)) enddef;

% [xa,xb] 구간에서 곡선의 위쪽 절반을 잇는 경로
vardef arcpath(expr xa, xb, n) =
  (xa*u, yy(xa)*u)
  for i=1 upto n: .. ((xa+(xb-xa)*i/n)*u, yy(xa+(xb-xa)*i/n)*u) endfor
enddef;

def drawaxes =
  drawarrow (-2.1u, 0)--(2.55u, 0) withcolor .55white;
  drawarrow (0, -2.15u)--(0, 2.35u) withcolor .55white;
enddef;

def drawcurve =
  begingroup
    save ovalu, branchu;
    path ovalu, branchu;
    ovalu := arcpath(-1.61803, 0.61803, 48);
    branchu := arcpath(1, 2.25, 40);
    pickup pencircle scaled 1bp;
    draw ovalu .. reverse(ovalu reflectedabout ((0,0),(1,0))) .. cycle;
    draw (reverse(branchu reflectedabout ((0,0),(1,0)))) .. branchu;
    pickup defaultpen;
  endgroup
enddef;

def spot(expr z) = fill fullcircle scaled 3.6bp shifted z enddef;

beginfig(1) % 덧셈: P + Q
  drawaxes; drawcurve;
  pair P, Q, R, S;
  P := (-1, 1.41421)*u;     % P
  Q := (0, 1)*u;            % Q
  R := (1.17157, 0.51472)*u;   % 현의 셋째 교점 R
  S := (1.17157, -0.51472)*u;  % P+Q = R의 반사
  % 현: y = 1 - 0.41421x 를 양쪽으로 늘인 것
  draw (-1.65u, 1.68345u)--(1.65u, 0.31655u) withcolor .3white;
  draw R--S dashed evenly withcolor .3white;
  spot(P); spot(Q); spot(R); spot(S);
  label.ulft(btex $P$ etex, P);
  label.ulft(btex $Q$ etex, Q);
  label.urt(btex $R$ etex, R);
  label.lrt(btex $P{+}Q$ etex, S);
endfig;

beginfig(2) % 두 배: 2P
  drawaxes; drawcurve;
  pair P, R, S;
  P := (1.2, 0.57271)*u;       % P
  R := (1.70274, 1.59093)*u;   % 접선의 둘째 교점 R
  S := (1.70274, -1.59093)*u;  % 2P = R의 반사
  % P에서의 접선: 기울기 2.02551
  draw (0.85u, -0.13622u)--(1.9u, 1.99057u) withcolor .3white;
  draw R--S dashed evenly withcolor .3white;
  spot(P); spot(R); spot(S);
  label.lrt(btex $P$ etex, P);
  label.ulft(btex $R$ etex, R);
  label.lrt(btex $2P$ etex, S);
endfig;

end
