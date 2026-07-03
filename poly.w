@* 다항식.
Schoof 알고리즘은 곡선 위의 점이 아니라 다항식 위에서 논다. 프로베니우스
사상이니 $\ell$-등분점이니 하는 등장인물이 모두 $F_p[x]$의 원소로 분장하고
나오기 때문에, 먼저 다항식 산술을 갖추어야 한다. 이 장의 목표는 사칙과
합동식 산술, 그리고 — 이 프로그램의 속도를 결정짓는 — 빠른 곱셈이다.

계수는 낮은 차수부터 |[]uint64| 한 판에 담는다. 첨자가 곧 차수라서
$3x^3+2x+1$은 |{1, 2, 0, 3}|이다. 필기 순서와 반대라 처음에는 어색하지만,
$i$차 항이 |c[i]|에 있다는 규칙은 곱셈에서 |c[i+j]|처럼 첨자 셈이 그대로
지수 법칙이 되는 즐거움을 준다. 연산은 언제나 새 다항식을 돌려주고 피연산자를
건드리지 않으며, 최고차 계수가 $0$이 아니도록(상수항만은 예외로 남겨)
다듬어진 상태를 유지한다. 그래서 |Deg|는 길이에서 1을 빼면 끝이고, 영
다항식의 차수는 편의상 $0$으로 친다.
@<다항식@>=
type FpPoly struct {
	f *Fp
	c []uint64
}

func (f *Fp) Poly(coeffs ...uint64) FpPoly {
	c := make([]uint64, len(coeffs))
	for i, v := range coeffs {
		c[i] = f.reduce(v)
	}
	if len(c) == 0 {
		c = []uint64{0}
	}
	return FpPoly{f, c}.trim()
}

func (p FpPoly) Deg() int     { return len(p.c) - 1 }
func (p FpPoly) isZero() bool { return len(p.c) == 1 && p.c[0] == 0 }

@ 덧셈·뺄셈 뒤에는 최고차 계수가 상쇄되어 $0$이 될 수 있으니, |trim|이 웃자란
$0$들을 쳐낸다. 상수항까지 치지 않는 것은 영 다항식도 길이 $1$은 가져야 하기
때문이다.
@<다항식@>=
func (p FpPoly) trim() FpPoly {
	n := len(p.c)
	for n > 1 && p.c[n-1] == 0 {
		n--
	}
	return FpPoly{p.f, p.c[:n]}
}

func (p FpPoly) clone() FpPoly {
	return FpPoly{p.f, append([]uint64(nil), p.c...)}
}

func (p FpPoly) equal(q FpPoly) bool {
	if len(p.c) != len(q.c) {
		return false
	}
	for i := range p.c {
		if p.c[i] != q.c[i] {
			return false
		}
	}
	return true
}

@ 덧셈은 짧은 쪽을 긴 쪽에 포갠다.
@<다항식@>=
func (p FpPoly) Add(q FpPoly) FpPoly {
	f := p.f
	a, b := p.c, q.c
	if len(a) < len(b) {
		a, b = b, a
	}
	r := make([]uint64, len(a))
	for i := range b {
		r[i] = f.add(a[i], b[i])
	}
	copy(r[len(b):], a[len(b):])
	return FpPoly{f, r}.trim()
}

@ 뺄셈은 피감수가 짧을 수도 있어 자리마다 있는 만큼만 꺼내 뺀다.
@<다항식@>=
func (p FpPoly) Sub(q FpPoly) FpPoly {
	f := p.f
	n := max(len(p.c), len(q.c))
	r := make([]uint64, n)
	for i := range n {
		var a, b uint64
		if i < len(p.c) {
			a = p.c[i]
		}
		if i < len(q.c) {
			b = q.c[i]
		}
		r[i] = f.sub(a, b)
	}
	return FpPoly{f, r}.trim()
}

@ 부호 반전, 상수배, 상수항 더하기는 계수별 손질이다. |addConst|는 Schoof
장에서 $3x^2+A$ 같은 식을 만들 때 쓴다.
@<다항식@>=
func (p FpPoly) Neg() FpPoly {
	f := p.f
	r := make([]uint64, len(p.c))
	for i, v := range p.c {
		r[i] = f.neg(v)
	}
	return FpPoly{f, r}
}

func (p FpPoly) MulScalar(a uint64) FpPoly {
	f := p.f
	a = f.reduce(a)
	r := make([]uint64, len(p.c))
	for i, v := range p.c {
		r[i] = f.mul(v, a)
	}
	return FpPoly{f, r}.trim()
}

func (p FpPoly) addConst(a uint64) FpPoly {
	r := p.clone()
	r.c[0] = p.f.add(r.c[0], p.f.reduce(a))
	return r
}

@ 곱셈은 두 벌을 마련한다. 교과서식 $O(n^2)$ 곱셈과, 다음 절부터 펼쳐질
수론적 변환(NTT)의 $O(n\log n)$ 곱셈이다. 작은 다항식에는 교과서가 이긴다 —
변환의 상차림 비용이 본 요리보다 비싸기 때문이다. 문턱은 결과 길이 128로
잡았다. 반대쪽 끝의 안전판도 하나 둔다. 뒤에 보겠지만 우리 변환은 길이
$2^{23}$까지만 지원하므로, 그보다 긴(이 프로그램에서는 나올 일 없는) 곱은
교과서로 되돌린다.
@<다항식@>=
const mulThreshold = 128 // 이보다 짧은 곱은 교과서식이 더 빠르다

func (p FpPoly) Mul(q FpPoly) FpPoly {
	n := len(p.c) + len(q.c) - 1
	if n < mulThreshold || n > 1<<23 {
		return p.mulSchool(q)
	}
	return FpPoly{p.f, p.mulNTT(q, n)}.trim()
}

@ 교과서식 곱셈. $0$인 계수를 건너뛰는 것은 나눗셈 다항식처럼 듬성듬성한
피연산자에서 쏠쏠하다.
@<다항식@>=
func (p FpPoly) mulSchool(q FpPoly) FpPoly {
	f := p.f
	r := make([]uint64, len(p.c)+len(q.c)-1)
	for i, ai := range p.c {
		if ai == 0 {
			continue
		}
		for j, bj := range q.c {
			if bj == 0 {
				continue
			}
			r[i+j] = f.add(r[i+j], f.mul(ai, bj))
		}
	}
	return FpPoly{f, r}.trim()
}

@* 수론적 변환.
빠른 곱셈의 사연은 1805년으로 거슬러 오른다. 가우스는 소행성 궤도를 맞추다가
삼각급수 보간을 빨리 하는 법을 노트에 적어 두었는데, 그게 바로 고속 푸리에
변환(FFT)이다 — 푸리에가 열 방정식 논문을 내기 두 해 {\it 전}이다. 노트는
사후 전집에 라틴어로 묻혀 있다가, 1965년 Cooley와 Tukey가 같은 것을 재발견해
계산의 역사를 바꾼 뒤에야 재조명되었다. Tukey의 동기가 재미있다: 소련 핵실험을
지진계로 감시하려면 긴 신호의 스펙트럼을 그때그때 계산해야 했던 것이다.
냉전이 알고리즘 하나를 캐냈다.

다항식 곱셈이 여기 걸리는 이유는 곱셈이 곧 계수열의 합성곱이고, 변환이
합성곱을 점별 곱으로 바꿔 주기 때문이다. 차수 $n$ 다항식을 $2n$개 지점에서
값매김하고(변환), 값끼리 곱한 뒤(점별 곱 $O(n)$), 도로 보간하면(역변환) 곱이
나온다. 값매김·보간을 $1$의 $2^k$제곱근들에서 하면 분할정복으로 $O(n\log n)$ —
이것이 FFT 곱셈이다.

다만 복소수 FFT를 그대로 쓰면 부동소수점 오차가 계수에 스며든다. 우리 계수는
$F_p$의 원소라 반올림을 한 톨도 용납할 수 없다. 다행히 FFT에 필요한 것은
복소수 자체가 아니라 ``$1$의 원시 $2^k$제곱근이 있는 환''이라는 무대 장치뿐이다.
$m=c\cdot2^s+1$ 꼴 소수의 $F_m$이 정확히 그런 무대다: 곱셈군의 위수가
$m-1=c\cdot2^s$이니 원시근 $g$의 거듭제곱 $g^{(m-1)/2^k}$이 $1$의 원시
$2^k$제곱근이 된다($k\le s$). 이렇게 유한체에서 하는 FFT를 수론적
변환(number-theoretic transform, NTT)이라 부른다. 오차는 원천 봉쇄, 셈은
전부 정수다.

@ 그런데 무대가 하나로는 모자란다. 진짜 계수 — 법으로 줄이기 전의 합성곱
값 — 는 최대 $n\,(p-1)^2$까지 자라는데, $p$가 $2^{63}$ 언저리면 이는
$2^{149}$쯤 되어 어떤 64비트 법 $m$으로도 담을 수 없다. 처방은 오래된 것이다:
서로소인 법 세 개로 각각 셈하고 중국인의 나머지 정리(CRT)로 합치면, 세 법의
곱 $m_1m_2m_3\approx2^{152.5}$ 미만의 값이 유일하게 복원된다.
$n\le2^{23}$에서 $n\,(p-1)^2<2^{23}\cdot2^{126}=2^{149}$이므로 여유 있게
들어간다 — 앞서 |Mul|이 $2^{23}$에서 손을 뗀 이유가 이것이다.

세 소수는 경시 프로그래밍판에서 오래 굴러 검증된 명물들로 골랐다. 각각의
원시근도 함께 적는다. $m-1$이 $2$의 큰 거듭제곱을 인수로 가진다는 것, 즉
지원하는 변환 길이가 각각 $2^{23}$, $2^{56}$, $2^{57}$이라는 것이 자격
요건이다(셋의 최솟값이 전체의 한계가 된다).
@<다항식@>=
const (
	ntt1, ntt1g = 998244353, 3           // $119\cdot2^{23}+1$
	ntt2, ntt2g = 1945555039024054273, 5 // $27\cdot2^{56}+1$
	ntt3, ntt3g = 4179340454199820289, 3 // $29\cdot2^{57}+1$
)

@ 변환 본체는 제자리(in-place) 반복문 꼴의 Cooley--Tukey다. 먼저 첨자를
비트 뒤집기 순서로 재배열하고, 길이 $2,4,8,\ldots,n$의 나비(butterfly)
연산을 겹겹이 쌓는다. 나비 하나는 $(u,v)\mapsto(u+\omega^jv,\ u-\omega^jv)$ —
길이 $2$의 변환이다. 역변환은 회전 인자 $\omega$를 역원으로 갈고 끝에
$n^{-1}$을 곱하면 된다는 것이 변환의 잘 알려진 대칭이다.
@<다항식@>=
func ntt(a []uint64, m, g uint64, invert bool) {
	n := len(a)
	@<첨자를 비트 뒤집기 순서로 재배열한다@>@;
	for length := 2; length <= n; length <<= 1 {
		w := powmod(g, (m-1)/uint64(length), m)
		if invert {
			w = powmod(w, m-2, m)
		}
		half := length >> 1
		for i := 0; i < n; i += length {
			wn := uint64(1)
			for j := i; j < i+half; j++ {
				u, v := a[j], mulmod(a[j+half], wn, m)
				a[j] = addmod(u, v, m)
				a[j+half] = submod(u, v, m)
				wn = mulmod(wn, w, m)
			}
		}
	}
	if invert {
		ninv := powmod(uint64(n), m-2, m)
		for i := range a {
			a[i] = mulmod(a[i], ninv, m)
		}
	}
}

@ 재배열은 첨자 $i$의 비트를 뒤집은 자리 $j$와 맞바꾸는 것인데, $j$를 매번
새로 뒤집지 않고 ``이진수 덧셈을 최상위 비트부터 하는'' 증가 트릭으로
이어간다. 뒤집힌 세계에서 $1$을 더하는 셈이다.
@<첨자를 비트 뒤집기 순서로 재배열한다@>=
for i, j := 1, 0; i < n; i++ {
	bit := n >> 1
	for ; j&bit != 0; bit >>= 1 {
		j ^= bit
	}
	j |= bit
	if i < j {
		a[i], a[j] = a[j], a[i]
	}
}

@ 법 하나에서의 합성곱. 두 계수열을 길이 $n$으로 늘려 변환하고, 점별로
곱하고, 역변환한다. 들어오는 계수는 $p$로 줄어 있지 국소 법 $m$으로 줄어
있지는 않으니 먼저 $m$으로 줄인다($m_1$은 $p$보다 작을 수 있다).
@<다항식@>=
func nttConv(a, b []uint64, n int, m, g uint64) []uint64 {
	fa := make([]uint64, n)
	for i, v := range a {
		fa[i] = v % m
	}
	fb := make([]uint64, n)
	for i, v := range b {
		fb[i] = v % m
	}
	ntt(fa, m, g, false)
	ntt(fb, m, g, false)
	for i := range fa {
		fa[i] = mulmod(fa[i], fb[i], m)
	}
	ntt(fa, m, g, true)
	return fa
}

@ NTT 곱셈의 지휘부. 결과 길이 이상의 $2$의 거듭제곱 $n$을 잡고, 세 법에서
따로 합성곱을 구한 뒤 계수마다 CRT로 복원한다.
@<다항식@>=
func (p FpPoly) mulNTT(q FpPoly, rlen int) []uint64 {
	n := 1
	for n < rlen {
		n <<= 1
	}
	r1 := nttConv(p.c, q.c, n, ntt1, ntt1g)
	r2 := nttConv(p.c, q.c, n, ntt2, ntt2g)
	r3 := nttConv(p.c, q.c, n, ntt3, ntt3g)
	@<Garner의 방법으로 세 나머지에서 계수를 복원한다@>@;
	return r
}

@ 복원은 Garner의 방법을 쓴다. 참값 $x<m_1m_2m_3$를 혼합기수
$x=a_1+t_2\,m_1+t_3\,m_1m_2$ 꼴로 풀면
$$t_2=(a_2-a_1)\,m_1^{-1}\bmod m_2,\qquad
t_3=\bigl(a_3-(a_1+t_2\,m_1)\bigr)\,(m_1m_2)^{-1}\bmod m_3$$
이고, 이 전개를 $p$로 줄이며 그대로 읽으면 원하는 계수 $x\bmod p$다.
128비트 큰 수를 실제로 만들 필요가 전혀 없다는 것이 이 방법의 미덕이다 —
모든 셈이 |mulmod| 몇 번으로 끝난다. $m_1<m_2<m_3$이라 중간값들이 항상
해당 법 미만임도 눈여겨보라.
@<Garner의 방법으로 세 나머지에서 계수를 복원한다@>=
inv12 := powmod(ntt1, ntt2-2, ntt2)         // $m_1^{-1}\bmod m_2$
m12 := mulmod(ntt1, ntt2%ntt3, ntt3)        // $m_1m_2\bmod m_3$
inv123 := powmod(m12, ntt3-2, ntt3)         // $(m_1m_2)^{-1}\bmod m_3$
pm := p.f.p
m1p := ntt1 % pm
m12p := mulmod(m1p, ntt2%pm, pm)
r := make([]uint64, rlen)
for i := range r {
	t2 := mulmod(submod(r2[i], r1[i], ntt2), inv12, ntt2)
	x12 := addmod(r1[i], mulmod(t2, ntt1, ntt3), ntt3)
	t3 := mulmod(submod(r3[i], x12, ntt3), inv123, ntt3)
	v := addmod(r1[i]%pm, mulmod(t2%pm, m1p, pm), pm)
	r[i] = addmod(v, mulmod(t3%pm, m12p, pm), pm)
}

@* 다항식 나눗셈과 그 친구들.
나눗셈은 초등학교 세로셈 그대로다: 남은 것의 최고차 항을 제수의 최고차
항으로 나눠 몫의 한 자리를 얻고, 그만큼 덜어 낸다. 제수가 모닉(최고차 계수
$1$)이면 자리마다 하던 역원 곱이 통째로 사라지므로 따로 우대한다 — 나눗셈
다항식 계산이 바로 이 우대의 단골손님이다.
@<다항식@>=
func (p FpPoly) DivMod(q FpPoly) (quo, rem FpPoly) {
	f := p.f
	if len(p.c) < len(q.c) {
		return f.Poly(0), p.clone()
	}
	qd := q.Deg()
	monic := q.c[qd] == 1
	var qInv uint64
	if !monic {
		qInv = f.inv(q.c[qd])
	}
	quoC := make([]uint64, len(p.c)-len(q.c)+1)
	remC := append([]uint64(nil), p.c...)
	@<몫의 자리를 하나씩 정하며 나머지를 깎는다@>@;
	return FpPoly{f, quoC}.trim(), FpPoly{f, remC}.trim()
}

@ 안쪽 고리. 나머지의 차수가 제수 아래로 떨어지면 끝난다.
@<몫의 자리를 하나씩 정하며 나머지를 깎는다@>=
for {
	td := len(remC) - 1
	rd := td - qd
	if rd < 0 || (len(remC) == 1 && remC[0] == 0) {
		break
	}
	coef := remC[td]
	if !monic {
		coef = f.mul(coef, qInv)
	}
	quoC[rd] = coef
	for i, qi := range q.c {
		if qi == 0 {
			continue
		}
		remC[rd+i] = f.sub(remC[rd+i], f.mul(qi, coef))
	}
	n := len(remC)
	for n > 1 && remC[n-1] == 0 {
		n--
	}
	remC = remC[:n]
}

@ 나머지만 필요할 때의 준말과, 모닉으로 만드는 손질.
@<다항식@>=
func (p FpPoly) Mod(q FpPoly) FpPoly {
	_, r := p.DivMod(q)
	return r
}

func (p FpPoly) Monic() FpPoly {
	d := p.Deg()
	if p.c[d] == 1 {
		return p.clone()
	}
	return p.MulScalar(p.f.inv(p.c[d]))
}

@ 합동식 거듭제곱 $p^e\bmod m$. Schoof 알고리즘의 심장 박동이다 —
프로베니우스 $x\mapsto x^q$를 계산한다는 것이 곧 이 함수로 $e=q$를 먹이는
일이고, 그때마다 NTT 곱셈 예순 번쯤이 연달아 뛴다.
@<다항식@>=
func (p FpPoly) PowMod(e uint64, m FpPoly) FpPoly {
	r := p.f.Poly(1)
	if e == 0 {
		return r
	}
	base := p.Mod(m)
	for e > 0 {
		if e&1 == 1 {
			r = r.Mul(base).Mod(m)
		}
		e >>= 1
		if e > 0 {
			base = base.Mul(base).Mod(m)
		}
	}
	return r
}

@ 최대공약수는 유클리드 호제법이다. 기원전 300년의 알고리즘이 21세기의
암호 코드에 그대로 앉아 있다 — 이보다 오래 현역인 코드는 없다. 답은 상수배
차이가 남으니 모닉으로 다듬어 대표를 정한다.
@<다항식@>=
func (p FpPoly) GCD(q FpPoly) FpPoly {
	a, b := p, q
	for !b.isZero() {
		a, b = b, a.Mod(b)
	}
	return a.Monic()
}

@ 법 $h$에 대한 역원은 확장 유클리드로 구한다. 베주 계수 $t$를 나란히
끌고 가면 마지막에 $tp\equiv r\pmod h$가 남는데, $r$가 $0$ 아닌 상수면
$p^{-1}=t/r$이고, $r$의 차수가 남아 있으면 $p$와 $h$가 서로소가 아니라
역원이 없다. 실패를 |ok|로 알리는 까닭은 Schoof 장에서 밝혀진다 — 놀랍게도
이 ``실패''가 거기서는 횡재다.
@<다항식@>=
func (p FpPoly) ModInverse(h FpPoly) (FpPoly, bool) {
	f := p.f
	r, newR := h, p.Mod(h)
	t, newT := f.Poly(0), f.Poly(1)
	for !newR.isZero() {
		quo, _ := r.DivMod(newR)
		r, newR = newR, r.Sub(quo.Mul(newR))
		t, newT = newT, t.Sub(quo.Mul(newT))
	}
	if r.Deg() > 0 {
		return f.Poly(0), false
	}
	return t.MulScalar(f.inv(r.c[0])), true
}

@ 시험대. 무작위 다항식을 만드는 손이 먼저 필요하다. 최고차 계수는 $0$이
되지 않게 한다.
@(elliptic_test.go@>=
func randPoly(f *Fp, n int) FpPoly {
	c := make([]uint64, n)
	for i := range c {
		c[i] = rng.Uint64() % f.p
	}
	if c[n-1] == 0 {
		c[n-1] = 1
	}
	return FpPoly{f, c}
}

@ NTT 소수 세 개가 정말 명물인지부터 재확인한다. 소수성은 |big.Int|의
Miller--Rabin에게 묻고, 원시근은 ``$g^{(m-1)/q}\ne1$이 $m-1$의 모든 소인수
$q$에 대해 성립''이라는 판정으로 확인한다.
@(elliptic_test.go@>=
func TestNTTModuli(t *testing.T) {
	cases := []struct {
		m, g uint64
		qs   []uint64 // $m-1$의 소인수들
	}{
		{ntt1, ntt1g, []uint64{2, 7, 17}},
		{ntt2, ntt2g, []uint64{2, 3}},
		{ntt3, ntt3g, []uint64{2, 29}},
	}
	for _, c := range cases {
		if !new(big.Int).SetUint64(c.m).ProbablyPrime(32) {
			t.Fatalf("%d은 소수가 아니다", c.m)
		}
		for _, q := range c.qs {
			if (c.m-1)%q != 0 {
				t.Fatalf("%d은 %d-1의 인수가 아니다", q, c.m)
			}
			if powmod(c.g, (c.m-1)/q, c.m) == 1 {
				t.Fatalf("%d은 %d의 원시근이 아니다", c.g, c.m)
			}
		}
	}
}

@ 곱셈 두 벌이 같은 답을 내는지 — 작은 체와 큰 체에서 각각 — 대조한다.
길이는 일부러 문턱 너머로 잡아 |Mul|이 NTT 길로 가게 한다.
@(elliptic_test.go@>=
func TestPolyMul(t *testing.T) {
	for _, p := range []uint64{97, p61} {
		f := NewFp(p)
		a, b := randPoly(f, 230), randPoly(f, 179)
		if !a.Mul(b).equal(a.mulSchool(b)) {
			t.Fatalf("p=%d: NTT 곱셈과 교과서 곱셈이 다르다", p)
		}
	}
}

@ 나눗셈은 $a=qb+r$와 $\deg r<\deg b$를, 역원은 $aa^{-1}\equiv1\pmod h$를
확인한다.
@(elliptic_test.go@>=
func TestPolyDiv(t *testing.T) {
	f := NewFp(10007)
	for range 50 {
		a := randPoly(f, 1+int(rng.Uint64()%40))
		b := randPoly(f, 1+int(rng.Uint64()%20))
		q, r := a.DivMod(b)
		if !r.isZero() && r.Deg() >= b.Deg() {
			t.Fatal("나머지의 차수가 제수보다 크거나 같다")
		}
		if !q.Mul(b).Add(r).equal(a) {
			t.Fatal("a != q*b + r")
		}
	}
}

func TestPolyInverse(t *testing.T) {
	f := NewFp(101)
	h := f.Poly(1, 1, 0, 0, 1) // $x^4+x+1$
	for range 50 {
		a := randPoly(f, 1+int(rng.Uint64()%4))
		inv, ok := a.ModInverse(h)
		if !ok {
			continue
		}
		if !a.Mul(inv).Mod(h).equal(f.Poly(1)) {
			t.Fatal("a * a^{-1} != 1 (mod h)")
		}
	}
}
