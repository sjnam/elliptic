@* 이산 로그 문제.
암호의 판을 뒤집어 볼 차례다. 밑점 $P$와 그 배수 $H=kP$가 주어졌을 때
$k$를 찾는 것이 타원 곡선 이산 로그 문제(ECDLP)다. ``로그''라는 이름은
곱셈군 버전 $g^k=h$에서 $k=\log_g h$인 데서 왔다---우리 군은 덧셈으로
쓰니 사실 나눗셈이지만, 이름은 이미 굳었다. 앞 장에서 본 비대칭을
기억하자: $k$에서 $kP$로는 오백 걸음, 거꾸로는 (알려진 최선이) $\sqrt N$
걸음. 이 낙차가 타원 곡선 암호의 밑천이고, 이 장은 그 ``알려진 최선''들을
직접 만들어 밑천의 크기를 손으로 재 본다. 모두 c의 밑점 위수 |c.N|이
아니라 임의의 위수 $n$을 받는 안짝 함수로 만드는데, 마지막 손님
Pohlig--Hellman이 부분군들에서 재활용하기 위해서다.

@ 첫 손님은 Shanks의 아기걸음-거인걸음(baby-step giant-step)이다. 주인공
Daniel Shanks는 학위 없이 해군 연구소에서 일하다 마흔넷에 박사가 된
늦깎이로, 원주율을 십만 자리까지 처음 계산한 사람이기도 하다. 1971년의 이
알고리즘은 시간-공간 맞바꿈의 원형이다. $m=\lceil\sqrt n\rceil$로 잡고
$k=a+mb$ ($1\le a\le m$, $0\le b\le m$)로 쪼개면
$$H=(a+mb)P\iff H-b\,(mP)=aP.$$
좌변을 모르니 우변을 몽땅 준비한다: 아기걸음으로 $aP$들을 사전에 담아
두고($m$개), 거인걸음으로 $H,\ H-mP,\ H-2mP,\ldots$를 성큼성큼 밟으며
사전에 있나 물으면 된다. $O(\sqrt n)$ 시간에 $O(\sqrt n)$ 공간---사전
열쇠는 좌표를 문자열로 이어 붙여 만든다.
@<이산 로그 문제@>=
func ptKey(x, y *big.Int) string { return x.Text(62) + "," + y.Text(62) }

func (c *Curve) Shank(px, py, hx, hy *big.Int) *big.Int {
	return c.shank(px, py, hx, hy, c.N)
}

func (c *Curve) shank(px, py, hx, hy, n *big.Int) *big.Int {
	m := new(big.Int).Sqrt(n)
	m.Add(m, big.NewInt(1))
	@<아기걸음: $aP$를 사전에 담는다@>@;
	@<거인걸음: $H-b(mP)$가 사전에 있는지 밟아 본다@>@;
	return nil
}

@ @<아기걸음: $aP$를 사전에 담는다@>=
baby := make(map[string]*big.Int)
rx, ry := new(big.Int), new(big.Int)
for a := big.NewInt(1); a.Cmp(m) <= 0; a.Add(a, big.NewInt(1)) {
	rx, ry = c.Add(rx, ry, px, py)
	baby[ptKey(rx, ry)] = new(big.Int).Set(a)
}

@ 걸음마다 $-mP$를 더하면 $H-b(mP)$가 차례로 나온다. 맞았다면 $k=a+mb$.
@<거인걸음: $H-b(mP)$가 사전에 있는지 밟아 본다@>=
sx, sy := c.Neg(px, py)
sx, sy = c.ScalarMult(sx, sy, m)
rx, ry = new(big.Int).Set(hx), new(big.Int).Set(hy)
for b := new(big.Int); b.Cmp(m) <= 0; b.Add(b, big.NewInt(1)) {
	if a, ok := baby[ptKey(rx, ry)]; ok {
		return new(big.Int).Add(a, new(big.Int).Mul(m, b))
	}
	rx, ry = c.Add(rx, ry, sx, sy)
}

@ 둘째 손님, Pollard의 $\rho$. 사전이 $\sqrt n$개나 되는 것이 아까웠던
John Pollard가 1978년에 공간을 상수로 줄였다. 발상은 생일 역설이다:
군 안에서 ``충분히 무작위한'' 걸음을 걸으면 $O(\sqrt n)$걸음 안에 같은
점을 두 번 밟을 공산이 크다. 걸음의 궤적은 꼬리 달린 고리---그리스 문자
$\rho$ 모양이라 이름이 $\rho$다. 걸음마다 현재 점을 $aP+bH$ 꼴로 장부에
적어 두면, 충돌 $a_1P+b_1H=a_2P+b_2H$에서
$$k\equiv(a_1-a_2)(b_2-b_1)^{-1}\pmod n$$
이 나온다.

걸음 함수는 $x$좌표를 3으로 나눈 나머지로 갈림길을 정한다: $P$를 더하거나,
두 배 하거나, $H$를 더한다. 어느 갈래든 장부($a$, $b$)를 같이 갱신한다.
이 함수가 진짜 무작위는 아니지만 무작위 흉내로는 충분하다는 것이 경험칙이다.
@<이산 로그 문제@>=
func (c *Curve) PollardRho(px, py, hx, hy *big.Int) *big.Int {
	return c.pollardRho(px, py, hx, hy, c.N)
}

func (c *Curve) pollardRho(px, py, hx, hy, n *big.Int) *big.Int {
	step := func(x, y, a, b *big.Int) (*big.Int, *big.Int, *big.Int, *big.Int) {
		switch new(big.Int).Mod(x, big.NewInt(3)).Int64() {
		case 0:
			x, y = c.Add(x, y, px, py)
			a = new(big.Int).Add(a, big.NewInt(1))
			a.Mod(a, n)
		case 1:
			x, y = c.Double(x, y)
			a = new(big.Int).Lsh(a, 1)
			a.Mod(a, n)
			b = new(big.Int).Lsh(b, 1)
			b.Mod(b, n)
		default:
			x, y = c.Add(x, y, hx, hy)
			b = new(big.Int).Add(b, big.NewInt(1))
			b.Mod(b, n)
		}
		return x, y, a, b
	}
	@<무작위 출발점에서 토끼와 거북을 달리게 한다@>@;
	return nil
}

@ 충돌 감지는 Floyd의 토끼와 거북이다: 거북은 한 걸음, 토끼는 두 걸음.
고리에 들어서면 토끼가 거북을 반드시 따라잡는다---같은 점을 몽땅 저장하는
대신 둘만 기억하는 것이 공간 $O(1)$의 비결이다. 궤적이 불운하면(이를테면
$b_1=b_2$거나 분모가 $n$과 서로소가 아니면) 새 출발점으로 다시 달린다.
기대 걸음 수는 $\sqrt{\pi n/2}$쯤이니 여유를 두어 $8\sqrt n$걸음까지 달리고
접는다.
@<무작위 출발점에서 토끼와 거북을 달리게 한다@>=
limit := new(big.Int).Sqrt(n)
limit.Mul(limit, big.NewInt(8)).Add(limit, big.NewInt(16))
li := int64(1) << 62
if limit.IsInt64() {
	li = limit.Int64()
}
for range 64 {
	@<무작위 $a_1P+b_1H$에서 출발한다@>@;
	x2, y2, a2, b2 := x1, y1, a1, b1
	for j := int64(0); j < li; j++ {
		x1, y1, a1, b1 = step(x1, y1, a1, b1)
		x2, y2, a2, b2 = step(step(x2, y2, a2, b2))
		if x1.Cmp(x2) == 0 && y1.Cmp(y2) == 0 {
			@<충돌 장부에서 $k$를 풀어 본다; 성공이면 답한다@>@;
			break
		}
	}
}

@ @<무작위 $a_1P+b_1H$에서 출발한다@>=
a1, _ := rand.Int(rand.Reader, n)
b1, _ := rand.Int(rand.Reader, n)
vx, vy := c.ScalarMult(px, py, a1)
wx, wy := c.ScalarMult(hx, hy, b1)
x1, y1 := c.Add(vx, vy, wx, wy)

@ 분모 $b_2-b_1$의 역원이 없으면($n$이 합성수면 생긴다) 이 궤적은 버린다.
풀린 $k$는 $kP=H$로 검산하고 내보낸다.
@<충돌 장부에서 $k$를 풀어 본다; 성공이면 답한다@>=
den := new(big.Int).Sub(b2, b1)
den.Mod(den, n)
if den.Sign() == 0 || den.ModInverse(den, n) == nil {
	break
}
k := new(big.Int).Sub(a1, a2)
k.Mul(k, den).Mod(k, n)
tx, ty := c.ScalarMult(px, py, k)
if tx.Cmp(hx) == 0 && ty.Cmp(hy) == 0 {
	return k
}

@ 셋째 손님 앞에 연장이 하나 필요하다: 정수의 소인수분해다. 여기서도
Pollard의 $\rho$가 나온다---같은 생일 역설을 $x\mapsto x^2+1\bmod n$의
궤적에 적용하면 $n$의 소인수 $q$를 $O(\sqrt q)$에 뽑아내는 인수분해법이
된다(궤적을 $q$로 줄여 보면 더 일찍 고리에 드는데, 그 순간이
$\gcd(x_i-x_j,n)>1$로 드러난다). 한 사람이 같은 곡조로 두 문제를 푼 셈이다.
뽑아낸 인수가 소수라는 보장은 없으니---$\rho$는 $12=2\cdot6$처럼 아무
약수나 물어 온다---재귀로 양쪽을 마저 쪼갠다. 소수를 만나면 목록에 얹는다.
@<이산 로그 문제@>=
func factorize(n *big.Int) []*big.Int {
	var factors []*big.Int
	var split func(m *big.Int)
	split = func(m *big.Int) {
		if m.Cmp(big.NewInt(1)) == 0 {
			return
		}
		if m.ProbablyPrime(20) {
			factors = append(factors, m)
			return
		}
		d := rhoFactor(m)
		split(d)
		split(new(big.Int).Div(m, d))
	}
	split(new(big.Int).Set(n))
	return factors
}

@ 인수 하나를 찾는 $\rho$ 본체. 거북(|xs|)을 박아 두고 토끼(|x|)를 달리게
한 뒤 주기적으로 갱신하는 Floyd 변형이다. 걸음 함수 $x\mapsto x^2+c$의
상수 $c$를 잘못 골라 궤적이 통째로 한 고리에 갇히면($5^k$ 같은 소수
거듭제곱의 단골 사고다) 자명한 인수 $n$만 나오는데, 그럴 땐 $c$를 갈아
다시 달린다. 짝수는 $\rho$가 서툴러 $2$를 먼저 벗겨 넘긴다. 합성수에는
반드시 $\sqrt n$ 이하의 소인수가 있으니 언젠가 걸린다.
@<이산 로그 문제@>=
func rhoFactor(n *big.Int) *big.Int {
	if n.Bit(0) == 0 {
		return big.NewInt(2)
	}
	one := big.NewInt(1)
	for c := int64(1); ; c++ {
		xs := big.NewInt(2)
		x := big.NewInt(2)
		factor := big.NewInt(1)
		cycle := uint64(2)
		for factor.Cmp(one) == 0 {
			for j := uint64(0); j < cycle && factor.Cmp(one) == 0; j++ {
				x.Mul(x, x).Add(x, big.NewInt(c)).Mod(x, n)
				factor.GCD(nil, nil, new(big.Int).Sub(x, xs), n)
			}
			cycle <<= 1
			xs.Set(x)
		}
		if factor.Cmp(n) != 0 {
			return factor // 자명하지 않은 인수를 잡았다
		}
	}
}

@ 셋째 손님은 Pohlig--Hellman이다. 1978년, Diffie--Hellman 논문의 잉크가
마르기도 전에 나온 이 공격의 전언은 서늘하다: {\it 이산 로그의 어려움은
군 위수 $n$이 아니라 $n$의 가장 큰 소인수가 정한다.} $n=\prod q_i$
($q_i$는 서로소인 소수 거듭제곱)로 쪼개면, $t=n/q_i$ 배 한 점들
$P'=tP$, $H'=tH$는 위수 $q_i$의 부분군에 살고 거기서 $k\bmod q_i$가
풀린다. 조각들은 CRT가 꿰매 준다. 각 조각은 Shanks나 $\rho$로 $\sqrt{q_i}$
걸음이니, $n$이 매끄러우면(smooth---작은 소인수뿐이면) 전체가 와르르
무너진다. 실무에서 곡선의 위수를 소수로 고르는 이유가 정확히 이것이다 —
secp256k1의 $N$도 소수다.
@<이산 로그 문제@>=
func (c *Curve) PohligHellman(px, py, hx, hy *big.Int) *big.Int {
	N := c.N
	@<위수를 소수 거듭제곱들로 쪼갠다@>@;
	@<부분군마다 작은 이산 로그를 푼다@>@;
	return CRT(dLogs, qs)
}

@ |factorize|가 준 소인수 목록을 정렬해 같은 소수끼리 거듭제곱으로 뭉친다.
@<위수를 소수 거듭제곱들로 쪼갠다@>=
factors := factorize(new(big.Int).Set(N))
sort.Slice(factors, func(i, j int) bool {
	return factors[i].Cmp(factors[j]) < 0
})
var qs []*big.Int
for i, j := 0, 0; i < len(factors); i = j {
	q := new(big.Int).Set(factors[i])
	for j = i + 1; j < len(factors) && factors[j].Cmp(factors[i]) == 0; j++ {
		q.Mul(q, factors[i])
	}
	qs = append(qs, q)
}

@ 조각이 작으면 사전을 쓰는 Shanks가, 크면 공간 걱정 없는 $\rho$가 낫다.
경계는 $2^{32}$로 잡았다---사전 항목 $2^{16}$개쯤은 껌값이다.
@<부분군마다 작은 이산 로그를 푼다@>=
var dLogs []*big.Int
for _, q := range qs {
	t := new(big.Int).Div(N, q)
	gx, gy := c.ScalarMult(px, py, t)
	bx, by := c.ScalarMult(hx, hy, t)
	var k *big.Int
	if q.BitLen() <= 32 {
		k = c.shank(gx, gy, bx, by, q)
	} else {
		k = c.pollardRho(gx, gy, bx, by, q)
	}
	if k == nil {
		return nil
	}
	dLogs = append(dLogs, k)
}

@ 장난감 곡선($N=19$)에서 Shanks와 $\rho$가 숨긴 $k$를 도로 찾는지 본다.
Shanks의 답 $a+mb$는 $N$을 넘을 수 있으니 $\bmod N$으로 견준다.
@(elliptic_test.go@>=
func TestShankAndRho(t *testing.T) {
	c := toyCurve()
	k := big.NewInt(13)
	hx, hy := c.ScalarBaseMult(k)
	for name, dlp := range map[string]func(a, b, x, y *big.Int) *big.Int{
		"Shank": c.Shank, "PollardRho": c.PollardRho,
	} {
		got := dlp(c.Gx, c.Gy, hx, hy)
		if got == nil || new(big.Int).Mod(got, c.N).Cmp(k) != 0 {
			t.Fatalf("%s: k=13을 %v로 잘못 찾았다", name, got)
		}
	}
}

@ Pohlig--Hellman은 위수가 매끄러운 제물이 필요하니 $F_{1019}$ 위에서
직접 사냥한다: $x$를 훑으며 곡선 위의 점을 찾고($1019\equiv3\pmod4$라
제곱근이 거듭제곱 한 방이다), 그 위수를 전수 덧셈으로 재서 합성수 위수인
점을 밑점으로 삼는다. 그런 다음 무작위 $k$를 숨기고 도로 찾게 한다.
@(elliptic_test.go@>=
func TestPohligHellman(t *testing.T) {
	const p = 1019
	c := &Curve{P: big.NewInt(p), A: big.NewInt(5), B: big.NewInt(7), BitSize: 10}
	f := NewFp(p)
	for x := uint64(1); x < p; x++ {
		r := f.add(f.mul(f.mul(x, x), x), f.add(f.mul(5, x), 7))
		y := f.pow(r, (p+1)/4)
		if f.mul(y, y) != r {
			continue
		}
		@<이 점의 위수를 재고, 합성수면 밑점으로 채택한다@>@;
	}
	if c.N == nil {
		t.Fatal("합성수 위수의 점을 찾지 못했다")
	}
	k := new(big.Int).Mod(big.NewInt(int64(rng.Uint64())), c.N)
	hx, hy := c.ScalarBaseMult(k)
	got := c.PohligHellman(c.Gx, c.Gy, hx, hy)
	if got == nil || got.Cmp(k) != 0 {
		t.Fatalf("N=%v에서 k=%v를 %v로 잘못 찾았다", c.N, k, got)
	}
}

@ 위수 재기는 $P,2P,3P,\ldots$를 $\cal O$가 나올 때까지 미련하게 더한다.
$p$가 천 남짓이라 하세 상한으로도 몇천 걸음이다.
@<이 점의 위수를 재고, 합성수면 밑점으로 채택한다@>=
gx := new(big.Int).SetUint64(x)
gy := new(big.Int).SetUint64(y)
d := int64(1)
for tx, ty := new(big.Int).Set(gx), new(big.Int).Set(gy); !isInf(tx, ty); d++ {
	tx, ty = c.Add(tx, ty, gx, gy)
}
if d > 4 && !big.NewInt(d).ProbablyPrime(20) {
	c.Gx, c.Gy, c.N = gx, gy, big.NewInt(d)
	break
}
