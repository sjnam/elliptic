@* 나눗셈 다항식.
Schoof 알고리즘으로 건너가기 전에 다리 하나를 놓아야 한다. $E$의
$n$-등분점($n$-torsion point), 곧 $nP=\cal O$인 점들을 통째로 붙드는
다항식이다. {\it 나눗셈 다항식}(division polynomial) $\psi_n\in F_p[x]$은
정확히 그 일을 한다: $\cal O$ 아닌 점 $P=(x,y)$에 대해
$$nP=\cal O\iff\psi_n(x)=0.$$
``점 $P$를 $n$으로 나누는 문제''에서 나온 이름이다. 홀수 $n$에서 $\psi_n$의
차수는 $(n^2-1)/2$인데, 이는 $n$-등분점이 $n^2-1$개($\cal O$ 제외)이고
$x$좌표를 $\pm y$가 나눠 가진다는 사실과 정확히 맞아떨어진다. $n^2$이라는
숫자에 잠깐 눈길을 주자---실수 위 곡선이라면 $n$등분점이 고리 하나에
$n$개뿐일 텐데, 복소수(그리고 유한체)의 타원 곡선은 도넛 모양이라 고리가
두 방향으로 있어 $n\times n$개가 된다. Schoof 알고리즘 비용의 대부분이 이
$n^2$에서 나온다.

고맙게도 $\psi_n$은 점화식으로 술술 나온다. $f=x^3+Ax+B$라 할 때
$$\psi_{2m+1}=\psi_{m+2}\psi_m^3-\psi_{m-1}\psi_{m+1}^3\cdot(16f^2)^{\pm1},
\qquad
\psi_{2m}={\psi_m\over\psi_2}\,(\psi_{m+2}\psi_{m-1}^2-\psi_{m-2}\psi_{m+1}^2)$$
꼴이다(원래 $\psi_n$은 $n$이 짝수일 때 $y$의 홀수 차수 항을 가지는데,
$y^2=f$로 짝수 차수만 남기고 $y$ 한 개를 약속으로 떼어 둔 것이 위 식의
$16f^2$ 곱셈·나눗셈으로 나타난다). 셈은 |divPoly|가 하고, 재귀가 같은 항을
거듭 찾으므로 캐시를 둔다.

점 세기 엔진의 곡선은 |uint64| 체 위에 산다. |big.Int| 층의 |Curve|와 별개로
작은 곡선 타입을 하나 두고, 캐시도 여기에 딸려 보낸다---Schoof가 소수마다
고루틴을 띄울 때 곡선을 하나씩 나눠 가지면 잠금 없이 병렬이 된다.
@<나눗셈 다항식@>=
type fpCurve struct {
	f    *Fp
	A, B uint64
	dp   map[int64]FpPoly // 나눗셈 다항식 캐시
}

func (c *fpCurve) poly() FpPoly {
	return c.f.Poly(c.B, c.A, 0, 1) // $x^3+Ax+B$
}

@ 밑돌 다섯 장은 손으로 놓는다. $\psi_0=0$, $\psi_1=1$이고
$$\psi_2=4f,\qquad \psi_3=3x^4+6Ax^2+12Bx-A^2,$$
$$\psi_4=(8x^6+40Ax^4+160Bx^3-40A^2x^2-32ABx-8A^3-64B^2)\cdot f$$
($\psi_2$와 $\psi_4$는 본디 $2y$, $4y(\cdots)$인데 위 약속대로 $y\mapsto y^2=f$로
바꿔 둔 것이다).
@<나눗셈 다항식@>=
func (c *fpCurve) divPoly(n int64) FpPoly {
	f := c.f
	if c.dp == nil {
		c.dp = make(map[int64]FpPoly)
	}
	if d, ok := c.dp[n]; ok {
		return d
	}
	cache := func(dp FpPoly) FpPoly { c.dp[n] = dp; return dp }
	A, B := c.A, c.B
	switch n {
	case 0:
		return cache(f.Poly(0))
	case 1:
		return cache(f.Poly(1))
	case 2:
		return cache(c.poly().MulScalar(4))
	case 3:
		@<$\psi_3$를 만들어 돌려준다@>@;
	case 4:
		@<$\psi_4$를 만들어 돌려준다@>@;
	}
	@<점화식으로 $\psi_n$을 만들어 돌려준다@>@;
}

@ @<$\psi_3$를 만들어 돌려준다@>=
a2 := f.mul(A, A)
return cache(f.Poly(
	f.neg(a2),                // $-A^2$
	f.mul(f.fromInt(12), B),  // $12B$
	f.mul(f.fromInt(6), A),   // $6A$
	0,
	f.fromInt(3),             // $3x^4$
))

@ @<$\psi_4$를 만들어 돌려준다@>=
a2 := f.mul(A, A)
a3 := f.mul(a2, A)
b2 := f.mul(B, B)
ab := f.mul(A, B)
c0 := f.sub(f.neg(f.mul(f.fromInt(8), a3)), f.mul(f.fromInt(64), b2))
return cache(f.Poly(
	c0,                             // $-8A^3-64B^2$
	f.neg(f.mul(f.fromInt(32), ab)), // $-32AB$
	f.neg(f.mul(f.fromInt(40), a2)), // $-40A^2$
	f.mul(f.fromInt(160), B),        // $160B$
	f.mul(f.fromInt(40), A),         // $40A$
	0,
	f.fromInt(8),                    // $8x^6$
).Mul(c.poly()))

@ 점화식 부분. $m=\lfloor n/2\rfloor$ 언저리의 다항식 다섯을 모아 조립한다.
홀수 $n=2m+1$이면 $\psi_{m+2}\psi_m^3-\psi_{m-1}\psi_{m+1}^3$인데, 두 항 중
짝수 첨자 쪽이 $(4f)$ 인수를 세제곱으로 세 개 들고 있으니 $16f^2$로 나눠
균형을 맞춘다(어느 쪽인지는 $m$의 홀짝이 정한다). 짝수 $n=2m$이면 전체를
$\psi_2$로 나눈다. 나눗셈은 언제나 나누어떨어진다---점화식이 보증하는
정체성이라, 나머지는 버려도 좋은 것이 아니라 애초에 $0$이다.
@<점화식으로 $\psi_n$을 만들어 돌려준다@>=
m := n / 2
p2m := c.divPoly(m - 2)
p1m := c.divPoly(m - 1)
pm := c.divPoly(m)
pm1 := c.divPoly(m + 1)
pm2 := c.divPoly(m + 2)
var dp FpPoly
if n&1 == 1 {
	den := c.poly().Mul(c.poly()).MulScalar(16)
	t1 := pm2.Mul(pm.Mul(pm).Mul(pm))
	t2 := p1m.Mul(pm1.Mul(pm1).Mul(pm1))
	if m&1 == 0 {
		t1, _ = t1.DivMod(den)
	} else {
		t2, _ = t2.DivMod(den)
	}
	dp = t1.Sub(t2)
} else {
	dp = pm.Mul(pm2.Mul(p1m.Mul(p1m)).Sub(p2m.Mul(pm1.Mul(pm1))))
	dp, _ = dp.DivMod(c.divPoly(2))
}
return cache(dp)

@ 차수 공식으로 조립을 검산한다. $n$이 홀수면 $\deg\psi_n=(n^2-1)/2$다. 짝수는
표준 나눗셈 다항식이 인수 $2y$를 품는데, 이 구현은 그 $y$를 $y^2=f$로 갈아
$\psi_2=4f$처럼 두므로 $x$-차수가 $(n^2-4)/2$에 $\deg f=3$이 더 붙어
$(n^2+2)/2$가 된다($n=2$면 $3$, $n=4$면 $9$).
@(elliptic_test.go@>=
func TestDivPolyDeg(t *testing.T) {
	c := &fpCurve{f: NewFp(10007), A: 3, B: 8}
	for n := int64(3); n <= 30; n++ {
		want := int((n*n - 1) / 2)
		if n%2 == 0 {
			want = int((n*n + 2) / 2)
		}
		if got := c.divPoly(n).Deg(); got != want {
			t.Fatalf("deg psi_%d = %d, want %d", n, got, want)
		}
	}
}
