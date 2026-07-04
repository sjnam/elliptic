@* 타원 곡선.
드디어 주인공이다. 체 $K$ 위의 타원 곡선은
$$E:\ y^2=x^3+Ax+B,\qquad 4A^3+27B^2\ne0$$
꼴 방정식의 해집합에 무한원점 $\cal O$를 하나 보탠 것이다(표수가 $2$나 $3$이
아닐 때는 어떤 삼차 곡선도 이 바이어슈트라스 표준형으로 옮길 수 있다).
조건 $4A^3+27B^2\ne0$은 판별식이 $0$이 아니라는 말로, 곡선이 스스로를
찌르거나(첨점) 겹치지(교차점) 않고 매끈하다는 보증이다.

이름부터 짚고 가자. 타원 곡선은 타원이 {\it 아니다}. 이름의 사연은 이렇다.
타원의 둘레를 구하려던 18세기 수학자들은 $\int dx/\sqrt{3차식}$ 꼴의 적분
앞에서 멈춰 섰다 — 초등함수로는 안 풀리는 이 부류가 {\it 타원 적분}이라는
이름을 얻었다. 아벨과 야코비가 그 역함수({\it 타원 함수})를 연구했고, 타원
함수가 매개화하는 곡선이 {\it 타원 곡선}으로 불리게 된 것이다. 말하자면
증조부의 직업이 성씨가 된 격이다. 혈통은 유서 깊다: 디오판토스의
《산술》에도 사실상의 타원 곡선 문제가 나오고, 페르마가 무한강하법을 벼린
곳도 여기다. 20세기에는 모델(Mordell)이 유리점들이 유한생성군을 이룸을
보였고, 와일스가 페르마의 마지막 정리를 무너뜨린 결정타도 타원 곡선의
모듈러성이었다. 그리고 1985년, Koblitz와 Miller가 서로 모르게 동시에
``이 군을 암호에 쓰자''고 제안하면서 이 곡선은 순수 수학의 귀족에서 만인의
호주머니(스마트폰의 TLS, 비트코인 지갑) 속 일꾼이 되었다.

@ 암호가 탐낸 것은 곡선 위의 점들이 이루는 {\it 군}(group)이다. 군 법칙은
자와 컴퍼스로 그릴 수 있을 만큼 기하적이다: 두 점 $P$, $Q$를 지나는 직선은
곡선과 반드시 세 번째 점 $R$에서 만나고(삼차 방정식이니 근이 셋), $P+Q$는
그 $R$를 $x$축에 대해 반사한 점이다.

\medskip
\centerline{\pdfpic{ecfig-1.pdf}}
\smallskip
\centerline{그림 1: $y^2=x^3-2x+1$ 위의 덧셈. 현이 곡선과 만나는 셋째 점
$R$를 반사하면 $P+Q$다.}
\medskip

``셋째 교점을 그냥 합이라 하지, 왜 굳이 반사까지?''라는 물음이 당연히
나온다. 반사가 없으면 결합법칙이 깨진다. 반사 규칙 아래에서는 한 직선 위의
세 점이 언제나 $P+Q+R=\cal O$를 만족하는데, 이 대칭적인 진술이 군 공리를
모두 끌어낸다. 항등원은 무한원점 $\cal O$ — 수직선이 곡선과 만나는 ``하늘
끝의 점''이다. $P=(x,y)$의 역원은 수직으로 마주 보는 $-P=(x,-y)$이고,
결합법칙은... 직접 좌표로 확인하려 들면 종이 몇 장이 순식간에 사라진다.
(제대로 하려면 사영 기하의 베주 정리나 복소 해석의 바이어슈트라스 $\wp$가
필요하다. 여기서는 시험 코드가 무작위 점 삼십 벌로 확인해 주는 것으로
만족하자 — 증명은 아니지만 잠은 잘 온다.)

이 프로그램은 점을 아핀 좌표 $(x,y)$ 그대로 다루고 무한원점은 $(0,0)$으로
표기한다. $B\ne0$인 곡선에서 $(0,0)$은 곡선 위에 없으므로 안심하고 쓸 수
있는 자리다(진지한 라이브러리는 나눗셈을 아끼려 사영 좌표를 쓰지만, 학습이
목적인 우리는 그림과 똑같이 생긴 공식이 더 값지다). 매개변수는 256비트
소수까지 감당해야 하니 이 층은 |big.Int|다.
@<타원 곡선@>=
type Curve struct {
	P       *big.Int // 바탕 체 $F_p$의 소수 $p$
	A, B    *big.Int // 곡선 방정식 $y^2=x^3+Ax+B$의 계수
	Gx, Gy  *big.Int // 밑점(생성원) $G$
	N       *big.Int // $G$의 위수
	H       *big.Int // 여인수(cofactor)
	BitSize int      // $p$의 비트 수
	Name    string   // 곡선의 통칭
}

@ 곡선 판정. 우변 $x^3+Ax+B$는 여러 곳에서 쓰이니 |rhs|로 떼어 둔다.
관례대로 무한원점 $(0,0)$에는 |false|를 답한다 — ``곡선 위의 점이냐''는
물음에 ``점이긴 한데 하늘에 있다''고 답할 수는 없으니.
@<타원 곡선@>=
func (c *Curve) rhs(x *big.Int) *big.Int {
	r := new(big.Int).Mul(x, x)
	r.Mul(r, x)
	r.Add(r, new(big.Int).Mul(c.A, x))
	r.Add(r, c.B)
	return r.Mod(r, c.P)
}

func (c *Curve) IsOnCurve(x, y *big.Int) bool {
	if x.Sign() < 0 || x.Cmp(c.P) >= 0 || y.Sign() < 0 || y.Cmp(c.P) >= 0 {
		return false
	}
	y2 := new(big.Int).Mul(y, y)
	y2.Mod(y2, c.P)
	return c.rhs(x).Cmp(y2) == 0
}

@ 무한원점 판별과 역원. $-P$는 $x$축 반사이니 $y$의 부호만 갈면 된다.
@<타원 곡선@>=
func isInf(x, y *big.Int) bool { return x.Sign() == 0 && y.Sign() == 0 }

func (c *Curve) Neg(x, y *big.Int) (*big.Int, *big.Int) {
	ny := new(big.Int).Neg(y)
	ny.Mod(ny, c.P)
	return new(big.Int).Set(x), ny
}

@ 덧셈. 그림 1을 좌표로 옮기면 된다. $P=(x_1,y_1)$과 $Q=(x_2,y_2)$를 잇는
현의 기울기가 $m$일 때, 직선 $y=m(x-x_1)+y_1$을 곡선 방정식에 넣으면
$x^3-m^2x^2+\cdots=0$이 되고, 세 근의 합이 $m^2$이라는 비에타 공식에서
$$x_3=m^2-x_1-x_2,\qquad y_3=m\,(x_1-x_3)-y_1$$
이 나온다($y_3$의 꼴에 반사가 이미 반영되어 있다). 특수한 경우가 셋이다:
한쪽이 $\cal O$면 다른 쪽이 답이고, $x_1=x_2$인데 $y_1+y_2\equiv0$이면 현이
수직이라 $\cal O$이며, 두 점이 같으면 현 대신 접선을 써야 하니 |Double|에
넘긴다.
@<타원 곡선@>=
func (c *Curve) Add(x1, y1, x2, y2 *big.Int) (*big.Int, *big.Int) {
	if isInf(x1, y1) {
		return new(big.Int).Set(x2), new(big.Int).Set(y2)
	}
	if isInf(x2, y2) {
		return new(big.Int).Set(x1), new(big.Int).Set(y1)
	}
	if x1.Cmp(x2) == 0 {
		s := new(big.Int).Add(y1, y2)
		if s.Mod(s, c.P).Sign() == 0 {
			return new(big.Int), new(big.Int)
		}
		return c.Double(x1, y1)
	}
	d := new(big.Int).Sub(x2, x1)
	d.Mod(d, c.P)
	d.ModInverse(d, c.P)
	m := new(big.Int).Sub(y2, y1)
	m.Mul(m, d).Mod(m, c.P)
	return c.chord(m, x1, y1, x2)
}

@ 기울기 $m$에서 셋째 교점을 구해 반사하는 마무리는 덧셈과 두 배에서
똑같으니 함수로 나눠 가진다.
@<타원 곡선@>=
func (c *Curve) chord(m, x1, y1, x2 *big.Int) (*big.Int, *big.Int) {
	x3 := new(big.Int).Mul(m, m)
	x3.Sub(x3, x1).Sub(x3, x2).Mod(x3, c.P)
	y3 := new(big.Int).Sub(x1, x3)
	y3.Mul(y3, m).Sub(y3, y1).Mod(y3, c.P)
	return x3, y3
}

@ 두 배. $Q$가 $P$로 다가가는 극한에서 현은 접선이 된다. 곡선 방정식을
음함수 미분하면 $2y\,y'=3x^2+A$, 곧 접선의 기울기는 $m=(3x^2+A)/2y$다.

\medskip
\centerline{\pdfpic{ecfig-2.pdf}}
\smallskip
\centerline{그림 2: 두 배. $P$에서의 접선이 곡선과 다시 만나는 점 $R$를
반사하면 $2P$다.}
\medskip

$y=0$이면 접선이 수직이라 $2P=\cal O$다 — 그런 $P$는 위수 $2$의 점이다.
$(0,0)$ 관례 덕에 이 검사가 무한원점까지 한꺼번에 처리한다.
@<타원 곡선@>=
func (c *Curve) Double(x, y *big.Int) (*big.Int, *big.Int) {
	if y.Sign() == 0 {
		return new(big.Int), new(big.Int)
	}
	d := new(big.Int).Lsh(y, 1)
	d.Mod(d, c.P)
	d.ModInverse(d, c.P)
	m := new(big.Int).Mul(x, x)
	m.Mul(m, big.NewInt(3)).Add(m, c.A)
	m.Mul(m, d).Mod(m, c.P)
	return c.chord(m, x, y, x)
}

@ 스칼라 곱 $kP$가 암호의 심장이다. $P$를 $k$번 더하면 되지만 $k$가
$2^{256}$쯤 되면 우주의 나이로도 모자라다. 두 배 셈이 있으니 이진법이
구원한다: $k$의 비트를 최상위부터 읽으며 ``두 배, (비트가 $1$이면) 더하기''를
반복하면 $\log_2k$걸음이다. $2^{256}$이 오백 걸음 남짓으로 준다 — 지수
함수를 상대로 거둔 로그의 승리이고, 뒤에 볼 이산 로그 문제의 어려움과
정확히 짝을 이루는 비대칭이다: 곱하기는 오백 걸음, 되돌리기는 $2^{128}$걸음.
$k$의 부호는 무시하고 절댓값을 쓴다.
@<타원 곡선@>=
func (c *Curve) ScalarMult(px, py, k *big.Int) (*big.Int, *big.Int) {
	x, y := new(big.Int), new(big.Int)
	kk := new(big.Int).Abs(k)
	for i := kk.BitLen() - 1; i >= 0; i-- {
		x, y = c.Double(x, y)
		if kk.Bit(i) == 1 {
			x, y = c.Add(x, y, px, py)
		}
	}
	return x, y
}

@ 밑점 $G$에 대한 준말과, ECDSA 검증이 쓰는 이중 스칼라 곱 $mG+nQ$.
@<타원 곡선@>=
func (c *Curve) ScalarBaseMult(k *big.Int) (*big.Int, *big.Int) {
	return c.ScalarMult(c.Gx, c.Gy, k)
}

func (c *Curve) CombinedMult(qx, qy, m, n *big.Int) (*big.Int, *big.Int) {
	x1, y1 := c.ScalarBaseMult(m)
	x2, y2 := c.ScalarMult(qx, qy, n)
	return c.Add(x1, y1, x2, y2)
}

@ 시험용 장난감 곡선은 암호 교과서의 고전인 $F_{17}$ 위의
$y^2=x^3+2x+2$다. 점 $G=(5,1)$의 위수는 $19$ — 군 전체가 열아홉 점짜리
순환군이라 손으로도 다 그려 볼 수 있는 크기다.
@(elliptic_test.go@>=
func toyCurve() *Curve {
	return &Curve{
		P:  big.NewInt(17),
		A:  big.NewInt(2),
		B:  big.NewInt(2),
		Gx: big.NewInt(5),
		Gy: big.NewInt(1),
		N:  big.NewInt(19),
		H:  big.NewInt(1), BitSize: 5, Name: "toy17",
	}
}

@ 군다운지 심문한다: $19G=\cal O$, $18G=-G$, 그리고 무작위 점들로
교환·결합법칙.
@(elliptic_test.go@>=
func TestCurveGroup(t *testing.T) {
	c := toyCurve()
	if !c.IsOnCurve(c.Gx, c.Gy) {
		t.Fatal("G가 곡선 위에 없다")
	}
	x, y := c.ScalarBaseMult(c.N)
	if !isInf(x, y) {
		t.Fatalf("19G != O (%v, %v)", x, y)
	}
	x, y = c.ScalarBaseMult(big.NewInt(18))
	nx, ny := c.Neg(c.Gx, c.Gy)
	if x.Cmp(nx) != 0 || y.Cmp(ny) != 0 {
		t.Fatal("18G != -G")
	}
	@<무작위 배수들로 교환법칙과 결합법칙을 확인한다@>@;
}

@ $P=iG$, $Q=jG$, $R=kG$를 뽑아 $P+Q=Q+P$와 $(P+Q)+R=P+(Q+R)$를 재 본다.
@<무작위 배수들로 교환법칙과 결합법칙을 확인한다@>=
for range 30 {
	pt := func() (*big.Int, *big.Int) {
		return c.ScalarBaseMult(big.NewInt(int64(rng.Uint64() % 19)))
	}
	px, py := pt()
	qx, qy := pt()
	rx, ry := pt()
	x1, y1 := c.Add(px, py, qx, qy)
	x2, y2 := c.Add(qx, qy, px, py)
	if x1.Cmp(x2) != 0 || y1.Cmp(y2) != 0 {
		t.Fatal("P+Q != Q+P")
	}
	x1, y1 = c.Add(x1, y1, rx, ry)
	x2, y2 = c.Add(qx, qy, rx, ry)
	x2, y2 = c.Add(px, py, x2, y2)
	if x1.Cmp(x2) != 0 || y1.Cmp(y2) != 0 {
		t.Fatal("(P+Q)+R != P+(Q+R)")
	}
}
