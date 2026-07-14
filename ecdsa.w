@* ECDSA.
마지막 장은 앞의 모든 것을 조립해 실제로 돈이 오가는 물건을 만든다.
타원 곡선 디지털 서명 알고리즘(ECDSA)이다. 서명의 요구는 셋이다: 개인 키
$d$를 가진 사람만 서명을 만들 수 있고, 공개 키 $Q=dG$만 아는 누구든 검증할
수 있고, 서명에서 $d$가 새어 나가면 안 된다. 이산 로그가 어려운 한 셋째가
지켜진다는 것이 앞 장의 교훈이다.

먼저 열쇠부터. 개인 키는 $[1,N-1]$의 무작위 정수, 공개 키는 그 배수다.
암호학적 난수원을 밖에서 주입받는다---시험에서는 |crypto/rand|를 꽂는다.
@<ECDSA@>=
func (c *Curve) GenerateKey(rnd io.Reader) (priv, x, y *big.Int, err error) {
	nMinus1 := new(big.Int).Sub(c.N, big.NewInt(1))
	if priv, err = rand.Int(rnd, nMinus1); err != nil {
		return
	}
	priv.Add(priv, big.NewInt(1))
	x, y = c.ScalarBaseMult(priv)
	return
}

@ 서명 대상은 메시지의 해시다. FIPS 186-4의 지시대로 해시의 왼쪽
비트들을 위수의 비트 길이에 맞춰 자른다.
@<ECDSA@>=
func (c *Curve) hashToInt(hash []byte) *big.Int {
	orderBits := c.N.BitLen()
	orderBytes := (orderBits + 7) / 8
	if len(hash) > orderBytes {
		hash = hash[:orderBytes]
	}
	ret := new(big.Int).SetBytes(hash)
	if excess := len(hash)*8 - orderBits; excess > 0 {
		ret.Rsh(ret, uint(excess))
	}
	return ret
}

func FermatInverse(k, n *big.Int) *big.Int {
	return new(big.Int).Exp(k, new(big.Int).Sub(n, big.NewInt(2)), n)
}

@ 서명. 일회용 비밀 $k$를 뽑아 $R=kG$의 $x$좌표를 $r$로 삼고,
$$s=k^{-1}(z+rd)\bmod N$$
을 짝지으면 $(r,s)$가 서명이다. $r$나 $s$가 $0$이면 다시 뽑는다.

여기서 잔소리를 늘어놓을 대목이 있다. $k$는 서명마다 새로, 예측 불가능하게
뽑아야 한다. 같은 $k$를 두 번 쓰면 두 서명에서 $k$가, 이어서 개인 키
$d=(sk-z)r^{-1}$가 초등 대수로 풀려 나온다. 농담 같지만 2010년 소니
플레이스테이션 3가 정확히 이 사고를 쳤다---펌웨어 서명에 $k$를 상수로
박아 둔 것이다. 해커 그룹 fail0verflow는 연말 발표회에서 칠판 두 줄로
소니의 마스터 키를 유도해 보였다. 백억 달러짜리 보안이 중학교 연립방정식에
무너진 순간이니, 아래 |GenerateKey| 호출 한 줄에는 그 수업료가 들어 있는
셈이다.
@<ECDSA@>=
func (c *Curve) Sign(priv *big.Int, hash []byte) (r, s *big.Int) {
	N := c.N
	z := c.hashToInt(hash)
	for {
		k, x, _, err := c.GenerateKey(rand.Reader)
		if err != nil {
			continue
		}
		r = x.Mod(x, N)
		if r.Sign() == 0 {
			continue
		}
		s = new(big.Int).Mul(r, priv)
		s.Add(s, z)
		s.Mul(s, FermatInverse(k, N))
		s.Mod(s, N)
		if s.Sign() != 0 {
			return
		}
	}
}

@ 검증은 서명식을 뒤집은 것이다. $u_1=zs^{-1}$, $u_2=rs^{-1}$로 놓고
$u_1G+u_2Q$를 계산하면, 서명이 진짜일 때 이 점이 $kG$와 같아지므로 그
$x$좌표가 $r$와 일치해야 한다. 대수 한 줄 검산:
$u_1+u_2d=(z+rd)s^{-1}=k\pmod N$.
@<ECDSA@>=
func (c *Curve) Verify(qx, qy *big.Int, hash []byte, r, s *big.Int) bool {
	N := c.N
	if r.Sign() <= 0 || s.Sign() <= 0 || r.Cmp(N) >= 0 || s.Cmp(N) >= 0 {
		return false
	}
	sInv := FermatInverse(s, N)
	u1 := c.hashToInt(hash)
	u1.Mul(u1, sInv).Mod(u1, N)
	u2 := new(big.Int).Mul(r, sInv)
	u2.Mod(u2, N)
	x, y := c.CombinedMult(qx, qy, u1, u2)
	if isInf(x, y) {
		return false
	}
	return x.Mod(x, N).Cmp(r) == 0
}

@ 무대에 세울 곡선은 secp256k1이다. NIST 표준 곡선들(P-256 등)의 계수가
출처 불명의 상수를 품은 것과 달리 이 곡선은 $A=0$, $B=7$---``소매에 아무것도
숨기지 않은'' 담백한 선택이라, NSA가 표준 난수 생성기 Dual\_EC\_DRBG에
뒷문을 심었다는 사실이 드러난 뒤로 부쩍 사랑받았다. 사토시 나카모토가
비트코인에 채택하면서 지금은 지구에서 가장 부지런히 일하는 곡선이 되었다 —
이 순간에도 초당 수천 번씩 서명하고 검증한다. $p=2^{256}-2^{32}-977$이고
위수 $N$은 소수다(Pohlig--Hellman 장의 교훈 그대로).
@<ECDSA@>=
var S256 = &Curve{
	Name:    "secp256k1",
	P:       bigFromHex("fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f"),
	A:       big.NewInt(0),
	B:       big.NewInt(7),
	Gx:      bigFromHex("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798"),
	Gy:      bigFromHex("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8"),
	N:       bigFromHex("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141"),
	H:       big.NewInt(1),
	BitSize: 256,
}

func bigFromHex(s string) *big.Int {
	b, ok := new(big.Int).SetString(s, 16)
	if !ok {
		panic("elliptic: 잘못된 16진수 상수")
	}
	return b
}

@ secp256k1이 제대로 박혔는지($G$가 곡선 위에, $NG=\cal O$), 그리고
서명-검증 왕복이 도는지 본다. 위조 검사로는 다른 메시지의 해시를 들이민다.
@(elliptic_test.go@>=
func TestS256(t *testing.T) {
	c := S256
	if !c.IsOnCurve(c.Gx, c.Gy) {
		t.Fatal("G가 곡선 위에 없다")
	}
	if x, y := c.ScalarBaseMult(c.N); !isInf(x, y) {
		t.Fatal("NG != O")
	}
}

func TestECDSA(t *testing.T) {
	priv, qx, qy, err := S256.GenerateKey(rand.Reader)
	if err != nil {
		t.Fatal(err)
	}
	hash := sha256.Sum256([]byte("타원 곡선과 그 암호"))
	r, s := S256.Sign(priv, hash[:])
	if !S256.Verify(qx, qy, hash[:], r, s) {
		t.Fatal("참 서명이 검증에 떨어졌다")
	}
	forged := sha256.Sum256([]byte("위조된 메시지"))
	if S256.Verify(qx, qy, forged[:], r, s) {
		t.Fatal("위조 서명이 검증을 통과했다")
	}
}
