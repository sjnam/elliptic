@* Schoof 알고리즘.
곡선 $E/F_p$ 위의 점은 몇 개일까? 암호학자에게 이 물음은 한가한 호기심이
아니다 — 군의 위수 $\#E(F_p)$에 큰 소인수가 없으면 이산 로그가
Pohlig--Hellman에게 낱낱이 쪼개지므로(다음 장), 곡선을 고르려면 먼저 세어야
한다. $x$마다 $x^3+Ax+B$가 제곱잉여인지 보면 $O(p)$에 셀 수 있지만, $p$가
256비트면 우주가 먼저 식는다.

1933년 Hasse가 과녁을 크게 좁혀 두었다: $\#E(F_p)=p+1-t$로 쓰면
$$\vert t\vert\le2\sqrt p.$$
$x$마다 해가 평균 하나쯤(제곱잉여일 확률 절반에 해가 둘)이라는 직관의
엄밀한 형태로, 오차항 $t$를 {\it 프로베니우스의 자취}(trace of Frobenius)라
부른다. 그래도 후보가 $4\sqrt p$개 — 여전히 지수적으로 많다.

1985년 René Schoof가 처음으로 다항 시간 알고리즘을 내놓았다. 발상은
CRT식 분할이다: {\it $t$ 전체는 몰라도 $t\bmod\ell$은 작은 무대에서 알아낼
수 있다.} 작은 소수 $\ell$을 곱이 $4\sqrt p$를 넘을 때까지 모아 각각에서
$t\bmod\ell$을 구하면, 중국인의 나머지 정리가 $t$를 하나로 꿰어 준다.
$\ell$은 $O(\log p)$개면 족하다. 뒷이야기 하나: 이 논문은 처음에 ``우아하나
비실용적''이라는 평을 들었다. 큰 $\ell$에서 $\deg\psi_\ell=(\ell^2-1)/2$짜리
다항식을 주무르는 비용이 만만치 않았던 탓인데, Elkies와 Atkin이 $\psi_\ell$을
차수 $(\ell-1)/2$의 인수로 갈아 끼우는 개량(SEA 알고리즘)을 내놓아 실용의
문턱을 넘었다. 이 문서는 원조 Schoof를 구현한다 — 학습에는 원형이 낫다.

@ 작은 무대의 주인공은 프로베니우스 사상 $\pi(x,y)=(x^p,y^p)$다.
$F_p$의 원소는 $a^p=a$로 꿈쩍 않지만(페르마), 확대체에 사는 등분점들은
$\pi$가 뒤섞는다. 핵심 정리는 $\pi$가 $E$ 위에서 특성 방정식
$$\pi^2-t\,\pi+p=0$$
을 만족한다는 것이다. 이를 $\ell$-등분점에 국한하면, $\bar t=t\bmod\ell$과
$\bar p=p\bmod\ell$에 대해 $\pi^2(P)+\bar p\,P=\bar t\,\pi(P)$가 된다.
좌변을 계산해 놓고 $\bar t=0,1,\ldots,\ell-1$을 하나씩 우변에 넣어 맞는
것을 찾으면 $t\bmod\ell$이 나온다.

``$\ell$-등분점에 국한''을 코드로 말하면 이렇다: 점의 좌표를 구체적 수가
아니라 잉여환
$$R=F_p[x]/(\psi_\ell(x))$$
의 원소로 다룬다. $\psi_\ell$로 나눈 나머지만 기억하는 세계에서는 ``모든
$\ell$-등분점에서 동시에 성립하는 등식''만 참이 되므로, $R$에서의 등식
검사가 곧 원하는 국한이다. 등분점을 하나도 명시적으로 찾지 않고 전부를
한꺼번에 다루는 것 — 이것이 Schoof의 마술이다.
@<Schoof 알고리즘@>=
type qring struct {
	h FpPoly // 법: $\psi_\ell$ 또는 그 인수
	f *Fp
}

func (qr *qring) reduce(p FpPoly) FpPoly { return p.Mod(qr.h) }

@ $R$ 위의 ``점''은 사상(endomorphism)이다. $x$좌표는 $R$의 원소 $a(x)$,
$y$좌표는 $b(x)\cdot y$ 꼴로 두되 $y$는 장부에만 있는 유령으로 다룬다 —
$y^2$이 나올 때마다 곡선 방정식으로 $f(x)=x^3+Ax+B$로 바꿔치기하면 $y$는
식에서 영영 지수 $1$을 넘지 않으므로, $b$만 저장하면 된다. 항등 사상은
$(x,\,1\cdot y)$, 프로베니우스는 $(x^p,\,f^{(p-1)/2}\cdot y)$다
($y^p=y\,(y^2)^{(p-1)/2}=y\,f^{(p-1)/2}$).

$\cal O$는 \GO/의 |nil|로 나타낸다. 하늘의 점에 이만큼 어울리는 값도 없다.
@<Schoof 알고리즘@>=
type endo struct {
	qr   *qring
	x, y FpPoly // 점 $(x,\ y\cdot y_{\rm 유령})$
}

func newEnd(qr *qring, x, y FpPoly) *endo {
	return &endo{qr, qr.reduce(x), qr.reduce(y)}
}

func endoEq(pe, qe *endo) bool {
	return pe.x.equal(qe.x) && pe.y.equal(qe.y)
}

func endoNeg(pe *endo) *endo {
	if pe == nil {
		return nil
	}
	return newEnd(pe.qr, pe.x, pe.y.Neg())
}

@ 사상들의 덧셈은 곡선 장의 현-접선 공식 그대로다 — 수 대신 다항식이
들어갈 뿐이다. 다만 기울기의 분모를 ``나누기''가 문제다. $R$는 체가 아니라
환이라 역원이 없을 수도 있다. 역원 계산이 실패하면 어떻게 하나?

여기서 이 알고리즘의 가장 아름다운 반전이 나온다. |ModInverse|가 실패했다는
것은 분모와 $\psi_\ell$의 최대공약수가 자명하지 않다는 뜻 — 곧 방금
$\psi_\ell$의 {\it 진약수를 공짜로 주웠다}는 뜻이다. 그 인수로 법을 갈아
끼우면 무대는 좁아지고 셈은 빨라진다. 실패가 아니라 횡재다. 오류 값에
주운 인수를 실어 보내고, 바깥 고리가 법을 좁혀 처음부터 다시 돌게 한다.
@<Schoof 알고리즘@>=
type zeroDivError struct{ factor FpPoly }

func (e *zeroDivError) Error() string { return "0으로 나눗셈: 법의 인수 발견" }

@ 서로 다른 두 사상의 덧셈. $x$좌표가 같은데 $y$까지 같으면 두 배로 넘기고,
$y$가 다르면(반드시 서로 반수다) $\cal O$다. 기울기는
$m=(b_2-b_1)/(a_2-a_1)\cdot y$ 꼴인데, $x_3$에 쓰이는 $m^2$에서 $y^2=f$가
튀어나오므로 셋째 교점 공식이 $x_3=f\,m^2-a_1-a_2$로 변형된다($m$은 $y$를
뗀 몫만 저장). 인자 $A$와 $fx$는 여기서는 안 쓰이지만 |endoDouble|과
서명을 맞춘다.
@<Schoof 알고리즘@>=
func endoAdd(pe, qe *endo, A uint64, fx FpPoly) (*endo, error) {
	if pe == nil {
		return qe, nil
	}
	if qe == nil {
		return pe, nil
	}
	qr := pe.qr
	a1, b1 := pe.x, pe.y
	a2, b2 := qe.x, qe.y
	if a1.equal(a2) {
		if b1.equal(b2) {
			return endoDouble(pe, A, fx)
		}
		return nil, nil
	}
	adif := a2.Sub(a1)
	inv, ok := adif.ModInverse(qr.h)
	if !ok {
		return nil, &zeroDivError{adif}
	}
	m := qr.reduce(b2.Sub(b1).Mul(inv))
	m2 := qr.reduce(m.Mul(m))
	a3 := qr.reduce(fx.Mul(m2)).Sub(a1.Add(a2))
	b3 := qr.reduce(m.Mul(a1.Sub(a3))).Sub(b1)
	return newEnd(qr, a3, b3), nil
}

@ 두 배. 접선 기울기 $(3a_1^2+A)/(2b_1y)$의 분모에서도 $y$를 $f$로
바꿔치기한다.
@<Schoof 알고리즘@>=
func endoDouble(pe *endo, A uint64, fx FpPoly) (*endo, error) {
	if pe == nil {
		return nil, nil
	}
	qr := pe.qr
	a1, b1 := pe.x, pe.y
	num := qr.reduce(a1.Mul(a1)).MulScalar(3).addConst(A)
	den := qr.reduce(b1.Mul(fx)).MulScalar(2)
	inv, ok := den.ModInverse(qr.h)
	if !ok {
		return nil, &zeroDivError{den}
	}
	m := qr.reduce(num.Mul(inv))
	a3 := qr.reduce(fx.Mul(m.Mul(m))).Sub(a1.MulScalar(2))
	b3 := qr.reduce(m.Mul(a1.Sub(a3))).Sub(b1)
	return newEnd(qr, a3, b3), nil
}

@ 스칼라 곱은 곡선 장과 같은 두 배-덧셈이다. 여기서 $n$은 $p\bmod\ell$이라
|uint64|로 족하다.
@<Schoof 알고리즘@>=
func endoScalarMul(pe *endo, n, A uint64, fx FpPoly) (*endo, error) {
	if n == 0 {
		return nil, nil
	}
	re := newEnd(pe.qr, pe.x, pe.y)
	started := false
	for i := 63; i >= 0; i-- {
		var err error
		if started {
			if re, err = endoDouble(re, A, fx); err != nil {
				return nil, err
			}
		}
		if (n>>uint(i))&1 == 1 {
			if !started {
				started = true // |re|가 이미 $1\cdot pe$다
			} else if re, err = endoAdd(re, pe, A, fx); err != nil {
				return nil, err
			}
		}
	}
	return re, nil
}

@ $\ell=2$는 따로 논다. $t\equiv\#E\equiv p+1-t\pmod2$에서 $t$가 짝수인
것은 $E$에 위수 $2$의 점이 있을 때, 곧 $f(x)=x^3+Ax+B$가 $F_p$에 근을 가질
때다. 삼차식이 근을 가지는지는 기약성 검사로 판별한다: $x^p-x$는 $F_p$의
원소 전부를 근으로 가지는 다항식이므로, $\gcd(x^p-x,f)=1$이면 근이 없어
$t\equiv1$, 아니면 $t\equiv0\pmod2$다.
@<Schoof 알고리즘@>=
func irreducible(qr *qring, q uint64) bool {
	x := qr.f.Poly(0, 1)
	xq := x.PowMod(q, qr.h).Sub(x)
	return xq.GCD(qr.h).equal(qr.f.Poly(1))
}

@ $t\bmod\ell$을 구하는 본체다. 잉여환을 차리고, 프로베니우스와 그 제곱을
만들고, 특성 방정식의 우변 후보를 차례로 대 본다. 0으로 나눗셈이 인수를
물어다 주면 법을 좁혀 재시도한다.
@<Schoof 알고리즘@>=
var errNoCharPoly = errors.New("프로베니우스가 특성 방정식을 만족하지 않는다")

func traceMod(f *Fp, A, B uint64, ell int64) (int64, error) {
	c := &fpCurve{f: f, A: A, B: B}
	fx := c.poly()
	q := f.p
	if ell == 2 {
		if irreducible(&qring{fx, f}, q) {
			return 1, nil
		}
		return 0, nil
	}
	qr := &qring{c.divPoly(ell).Monic(), f}
	qModEll := q % uint64(ell)
	qHalf := q / 2
	x := f.Poly(0, 1)
	var err error
	for {
		@<오류를 살핀다: 인수를 주웠으면 법을 좁히고, 가망 없으면 접는다@>@;
		@<프로베니우스 $\pi$와 $\pi^2$을 만든다@>@;
		@<특성 방정식에 $\bar t$ 후보를 대 본다@>@;
	}
}

@ 첫 순회에서는 |err|가 |nil|이라 그냥 지나간다. 이후로는 0으로 나눗셈이
남긴 인수로 $h\gets\gcd(h,\hbox{인수})$ 하고 다시 돈다.
@<오류를 살핀다: 인수를 주웠으면 법을 좁히고, 가망 없으면 접는다@>=
var zd *zeroDivError
switch {
case errors.As(err, &zd):
	qr.h = qr.h.GCD(zd.factor)
case errors.Is(err, errNoCharPoly):
	return 0, err
}
err = nil

@ $\pi=(x^q,\,f^{(q-1)/2}y)$이고, $\pi^2$은 지수를 $q^2$으로 키우는 대신
이미 셈한 $x^q$과 $y$부분 $f^{q/2}$를 다시 $q$제곱, $q+1$제곱해 얻는다 —
$(x^q)^q=x^{q^2}$이고 $f^{\lfloor q/2\rfloor(q+1)}=f^{(q^2-1)/2}$이다.
$q^2$은 |uint64|를 넘칠 수 있으니 이 우회가 멋이 아니라 필수다.
@<프로베니우스 $\pi$와 $\pi^2$을 만든다@>=
xq := x.PowMod(q, qr.h)
yq := fx.PowMod(qHalf, qr.h)
pi := newEnd(qr, xq, yq)
pi2 := newEnd(qr, xq.PowMod(q, qr.h), yq.PowMod(q+1, qr.h))

@ 좌변 $S=\pi^2+\bar q\,{\rm id}$를 만들고 $S=\bar t\,\pi$가 되는 $\bar t$를
찾는다. $S=\cal O$면 $\bar t=0$, $S=\pm\pi$면 $\bar t=\pm1$이고, 나머지는
$\pi$를 거듭 더해 가며 맞춘다. 사상 연산이 |err|를 내면 |continue|로 바깥
고리에 돌아가 법을 좁히고 재시도한다. 후보가 다 떨어지면 이론상 있을 수
없는 일이니 |errNoCharPoly|로 접는다.
@<특성 방정식에 $\bar t$ 후보를 대 본다@>=
id := newEnd(qr, x, f.Poly(1))
var Q, S *endo
if Q, err = endoScalarMul(id, qModEll, A, fx); err != nil {
	continue
}
if S, err = endoAdd(pi2, Q, A, fx); err != nil {
	continue
}
if S == nil {
	return 0, nil
}
if endoEq(S, pi) {
	return 1, nil
}
if endoEq(endoNeg(S), pi) {
	return -1, nil
}
P := newEnd(qr, pi.x, pi.y)
for t := int64(2); t < ell-1; t++ {
	if P, err = endoAdd(P, pi, A, fx); err != nil {
		break
	}
	if endoEq(P, S) {
		return t, nil
	}
}
if err == nil {
	err = errNoCharPoly
}

@ 남은 것은 지휘부다. 그 전에 연장 두 개 — 다음 소수와 중국인의 나머지
정리. |NextPrime|은 Miller--Rabin 판정(|ProbablyPrime|)으로 홀수를 걸러
올라간다.
@<Schoof 알고리즘@>=
func NextPrime(n *big.Int) *big.Int {
	if n.Cmp(big.NewInt(2)) < 0 {
		return big.NewInt(2)
	}
	if n.Cmp(big.NewInt(2)) == 0 {
		return big.NewInt(3)
	}
	p := new(big.Int).Set(n)
	if p.Bit(0) == 0 {
		p.Add(p, big.NewInt(1))
		if p.ProbablyPrime(20) {
			return p
		}
	}
	for {
		p.Add(p, big.NewInt(2))
		if p.ProbablyPrime(20) {
			return p
		}
	}
}

@ 손자(孫子)의 《손자산경》에 ``셋씩 세면 둘 남고 다섯씩 세면 셋 남고
일곱씩 세면 둘 남는 물건''을 묻는 문제가 있다 — 답은 23이고, 서양이 이를
{\it 중국인의 나머지 정리}라 부르는 연유다. 서로소인 법 $n_i$들에 대한
잉여 $a_i$들에서 $x\bmod\prod n_i$를 복원한다. $q_i=\prod n/n_i$의
역원을 확장 유클리드(|GCD|)로 얻는 표준 공식이다.
@<Schoof 알고리즘@>=
func CRT(a, n []*big.Int) *big.Int {
	if a == nil || n == nil {
		return nil
	}
	prod := big.NewInt(1)
	for _, x := range n {
		prod.Mul(prod, x)
	}
	var c, q, s, z big.Int
	for i, x := range n {
		q.Div(prod, x)
		z.GCD(nil, &s, x, &q)
		if z.Int64() != 1 {
			return nil
		}
		c.Add(&c, s.Mul(a[i], s.Mul(&s, &q)))
	}
	return c.Mod(&c, prod)
}

@ 지휘부 |Schoof|는 |Curve|의 메서드다. 점 세기 엔진이 |uint64| 체라
$p<2^{61}$만 받는다($q+1$제곱 같은 지수가 |uint64| 안에 머물도록 여유를
두었다). 60비트 소수의 곡선이면 초 단위로 세어 준다 — |math/big|으로 짰던
옛 구현이 며칠 밤을 새우던 크기다.
@<Schoof 알고리즘@>=
func (c *Curve) Schoof() (*big.Int, error) {
	if c.P == nil || c.P.Sign() <= 0 || c.P.BitLen() > 61 {
		return nil, errors.New("elliptic: Schoof는 0 < p < 2^61 이 필요하다")
	}
	p := c.P.Uint64()
	f := NewFp(p)
	A := reduceBig(c.A, p)
	B := reduceBig(c.B, p)
	@<곱이 하세 구간을 덮을 때까지 소수 $\ell$을 모은다@>@;
	@<소수 $\ell$마다 고루틴 하나가 $t\bmod\ell$을 구한다@>@;
	@<CRT로 $t$를 복원하고 $\#E=p+1-t$를 낸다@>@;
}

@ 하세 구간의 길이가 $4\sqrt p$이므로 $\prod\ell$이 그보다 커야 $t$가
유일하게 잡힌다. |big.Int|의 |Sqrt|는 내림이니 1을 더해 안전쪽으로 잡는다.
@<곱이 하세 구간을 덮을 때까지 소수 $\ell$을 모은다@>=
fsq := new(big.Int).Sqrt(c.P)
fsq.Add(fsq, big.NewInt(1)).Mul(fsq, big.NewInt(4))
var ells []*big.Int
M := big.NewInt(1)
for l := big.NewInt(2); M.Cmp(fsq) <= 0; l = NextPrime(l) {
	ells = append(ells, new(big.Int).Set(l))
	M.Mul(M, l)
}

@ $\ell$끼리는 완전히 독립이니 병렬은 공짜다. |traceMod|가 곡선과 캐시를
제 것으로 만들므로 공유하는 것은 불변인 체 |f|뿐 — 잠금이 필요 없다.
@<소수 $\ell$마다 고루틴 하나가 $t\bmod\ell$을 구한다@>=
traces := make([]*big.Int, len(ells))
errs := make([]error, len(ells))
var wg sync.WaitGroup
for i := range ells {
	wg.Add(1)
	go func() {
		defer wg.Done()
		t, err := traceMod(f, A, B, ells[i].Int64())
		if err != nil {
			errs[i] = err
			return
		}
		traces[i] = big.NewInt(t)
	}()
}
wg.Wait()
for _, err := range errs {
	if err != nil {
		return nil, err
	}
}

@ CRT는 $[0,M)$의 대표를 주지만 하세의 $t$는 음수일 수 있으니, $M/2$를
넘으면 $M$을 빼서 대칭 구간으로 옮긴다.
@<CRT로 $t$를 복원하고 $\#E=p+1-t$를 낸다@>=
t := CRT(traces, ells)
if t.Cmp(new(big.Int).Rsh(M, 1)) >= 0 {
	t.Sub(t, M)
}
res := new(big.Int).Sub(c.P, t)
return res.Add(res, big.NewInt(1)), nil

@ 음수일 수 있는 |big.Int| 계수를 $[0,p)$로 데려오는 잔심부름.
@<Schoof 알고리즘@>=
func reduceBig(v *big.Int, p uint64) uint64 {
	P := new(big.Int).SetUint64(p)
	return new(big.Int).Mod(v, P).Uint64()
}

@ 시험은 삼심제다. 1심: 작은 소수들에서 무식한 전수 세기와 대조한다.
전수 세기는 제곱잉여 표를 만들어 $x$마다 $f(x)$가 표에 있는지 본다 —
$O(p)$ 시간, $O(p)$ 메모리의 정직한 셈이다.
@(elliptic_test.go@>=
func naiveCount(p, a, b uint64) *big.Int {
	f := NewFp(p)
	isQR := make([]bool, p)
	for x := uint64(0); x < p; x++ {
		isQR[f.mul(x, x)] = true
	}
	n := uint64(1) // 무한원점
	for x := uint64(0); x < p; x++ {
		r := f.add(f.mul(f.mul(x, x), x), f.add(f.mul(a, x), b))
		if r == 0 {
			n++
		} else if isQR[r] {
			n += 2
		}
	}
	return new(big.Int).SetUint64(n)
}

@ @(elliptic_test.go@>=
func TestSchoofSmall(t *testing.T) {
	for _, p := range []uint64{7, 11, 101, 1009, 10007} {
		f := NewFp(p)
		for range 3 {
			a, b := rng.Uint64()%p, rng.Uint64()%p
			disc := f.add(f.mul(4, f.mul(f.mul(a, a), a)), f.mul(27, f.mul(b, b)))
			if disc == 0 {
				continue // 특이 곡선은 타원 곡선이 아니다
			}
			c := &Curve{
				P: new(big.Int).SetUint64(p),
				A: new(big.Int).SetUint64(a),
				B: new(big.Int).SetUint64(b),
			}
			got, err := c.Schoof()
			if err != nil {
				t.Fatal(err)
			}
			if want := naiveCount(p, a, b); got.Cmp(want) != 0 {
				t.Fatalf("p=%d a=%d b=%d: Schoof=%v, 전수=%v", p, a, b, got, want)
			}
		}
	}
}

@ 2심: 전수 세기가 아직 가능한 중간 크기 $p=1000003$.
@(elliptic_test.go@>=
func TestSchoofMedium(t *testing.T) {
	const p = 1000003
	c := &Curve{
		P: big.NewInt(p),
		A: big.NewInt(31337),
		B: big.NewInt(271828),
	}
	got, err := c.Schoof()
	if err != nil {
		t.Fatal(err)
	}
	if want := naiveCount(p, 31337, 271828); got.Cmp(want) != 0 {
		t.Fatalf("Schoof=%v, 전수=%v", got, want)
	}
}

@ 3심: 전수 세기가 무리인 메르센 소수 $p=2^{31}-1$. 이제 대조할 답이
없으니 라그랑주 정리를 판사로 모신다 — 위수 $n$이 맞다면 곡선 위의 어떤
점이든 $nP=\cal O$여야 한다. $p\equiv3\pmod4$라 제곱근이
$r^{(p+1)/4}$로 바로 나오니 점 찾기도 쉽다.
@(elliptic_test.go@>=
func TestSchoofOrder(t *testing.T) {
	const p = 2147483647 // $2^{31}-1$
	c := &Curve{P: big.NewInt(p), A: big.NewInt(7), B: big.NewInt(11)}
	n, err := c.Schoof()
	if err != nil {
		t.Fatal(err)
	}
	f := NewFp(p)
	var px, py uint64
	for x := uint64(1); ; x++ {
		r := f.add(f.mul(f.mul(x, x), x), f.add(f.mul(7, x), 11))
		if y := f.pow(r, (p+1)/4); f.mul(y, y) == r {
			px, py = x, y
			break
		}
	}
	X, Y := c.ScalarMult(new(big.Int).SetUint64(px), new(big.Int).SetUint64(py), n)
	if !isInf(X, Y) {
		t.Fatalf("#E=%v인데 nP != O", n)
	}
}
