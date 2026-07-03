@* 유한체 $F_p$.
타원 곡선 암호의 무대는 실수 평면이 아니라 유한체다. 소수 $p$에 대해
$\{0,1,\ldots,p-1\}$에 ``$p$로 나눈 나머지'' 셈을 얹으면 덧셈·뺄셈·곱셈은
물론 $0$ 아닌 수의 나눗셈까지 되는 체 $F_p$가 된다. 시계가 좋은 비유다:
12시에서 세 시간을 더 가면 3시이듯, $F_{13}$에서 $12+3=2$다. 다만 시계와
달리 곱셈의 역원까지 있다는 점이 요긴하다 — $F_{13}$에서 $5\times8=40=1$이니
$1/5=8$이다.

이 장은 $p<2^{63}$인 소수에 대한 $F_p$를 |uint64| 하나로 구현한다. 원소는
$[0,p)$의 대표 잉여로 저장하고, 힙 할당이 전혀 없다. 굳이 이렇게 하는 데는
쓰라린 사연이 있다. 전에 같은 것을 |math/big|으로 짰더니 Schoof 알고리즘이
소수 하나 세는 데 하세월이었다. |big.Int|는 아무리 작은 수라도 힙에 눕고,
덧셈 한 번에도 포인터를 좇는다. 반면 |uint64|는 레지스터에서 태어나
레지스터에서 죽는다. 수백만 번 반복되는 안쪽 고리에서 이 차이는 수십 배로
벌어진다. 물론 256비트 소수를 쓰는 ECDSA 층은 |big.Int|로 갈 수밖에 없지만
(그쪽은 연산 횟수가 몇백 번뿐이라 아무래도 좋다), 다항식을 수십만 번 곱하는
점 세기 엔진은 이 작고 빠른 체 위에 세운다.

$p<2^{63}$이라는 제한은 게으름이 아니라 설계다. 두 원소가 $2^{63}$ 미만이면
합이 |uint64|를 넘치지 않아 덧셈에 128비트 중간값이 필요 없다.
@<유한체@>=
type Fp struct {
	p uint64
}

func NewFp(p uint64) *Fp {
	if p < 2 || p>>63 != 0 {
		panic("elliptic: Fp의 법은 2 <= p < 2^63 이어야 한다")
	}
	return &Fp{p: p}
}

@ 낱개 연산 넷은 손가락 셈이다. 덧셈은 넘치면 $p$를 한 번 빼고, 뺄셈은
모자라면 $p$를 한 번 꾼다. |fromInt|는 음수일 수도 있는 정수를 $[0,p)$로
데려온다 — \GO/의 |%|는 피제수의 부호를 따르므로 음수에는 $p$를 한 번 더
얹어야 한다.
@<유한체@>=
func (f *Fp) reduce(a uint64) uint64 { return a % f.p }

func (f *Fp) fromInt(v int64) uint64 {
	r := v % int64(f.p)
	if r < 0 {
		r += int64(f.p)
	}
	return uint64(r)
}

func (f *Fp) add(a, b uint64) uint64 { return addmod(a, b, f.p) }
func (f *Fp) sub(a, b uint64) uint64 { return submod(a, b, f.p) }

func (f *Fp) neg(a uint64) uint64 {
	if a == 0 {
		return 0
	}
	return f.p - a
}

@ 진짜 일꾼은 법이 매개변수인 다음 네 함수다. |Fp|의 메서드로 만들지 않고
따로 둔 것은 뒤에 나올 NTT가 자기만의 법 세 개를 들고 이 함수들을 다시 찾아올
것이기 때문이다.

곱셈이 볼거리다. $a,b<p<2^{64}$의 곱은 128비트까지 자라는데, \GO/에는
128비트 정수가 없다. 대신 |math/bits|가 두 마디를 내주는 곱셈 |Mul64|와 두
마디를 받는 나눗셈 |Div64|를 갖추고 있어, 곱하고 나머지를 취하는 일이 정확히
두 명령이다. |Div64|는 상위 마디가 법 이상이면 몫이 넘친다고 패닉하지만,
$a,b<m$이면 $ab<m\cdot2^{64}$라 상위 마디가 늘 $m$ 미만이니 안전하다.
@<유한체@>=
func mulmod(a, b, m uint64) uint64 {
	hi, lo := bits.Mul64(a, b)
	_, r := bits.Div64(hi, lo, m)
	return r
}

func addmod(a, b, m uint64) uint64 {
	s := a + b // $a,b<m<2^{63}$이니 넘치지 않는다
	if s >= m {
		s -= m
	}
	return s
}

func submod(a, b, m uint64) uint64 {
	if a >= b {
		return a - b
	}
	return m - b + a
}

@ 거듭제곱은 지수의 이진 표현을 따라 제곱하며 오르는 정석 그대로다. 지수가
$e=e_0+2e_1+4e_2+\cdots$일 때 $a^e=\prod_{e_i=1}a^{2^i}$이므로, 비트를 하나씩
밀며 밑을 제곱해 간다. $e<2^{64}$라도 예순네 걸음이면 끝난다.
@<유한체@>=
func powmod(a, e, m uint64) uint64 {
	r := uint64(1) % m
	a %= m
	for ; e > 0; e >>= 1 {
		if e&1 == 1 {
			r = mulmod(r, a, m)
		}
		a = mulmod(a, a, m)
	}
	return r
}

@ 나눗셈은 페르마의 작은 정리에 기댄다. 소수 $p$와 $p$의 배수가 아닌 $a$에
대해 $a^{p-1}\equiv1\pmod p$이므로 $a^{-1}=a^{p-2}$다. 1640년에 페르마가
친구 프레니클에게 보낸 편지에서 ``증명을 보내 주고 싶지만 너무 길어질까
두렵다''며 증명 없이 알려 준 그 정리다 — 여백이 모자라다던 마지막 정리보다야
사정이 낫지만, 증명 미루는 버릇은 한결같다. 확장 유클리드 호제법이 더
빠르지만, 예순 걸음 남짓의 거듭제곱도 충분히 싸고 코드가 곧다.
@<유한체@>=
func (f *Fp) mul(a, b uint64) uint64 { return mulmod(a, b, f.p) }
func (f *Fp) pow(a, e uint64) uint64 { return powmod(a, e, f.p) }
func (f *Fp) inv(a uint64) uint64    { return powmod(a, f.p-2, f.p) }

@ 시험은 메르센 소수 $2^{61}-1$ 위에서 |math/big|과 대조한다. 답안지 채점을
채점 대상보다 믿음직한 자에게 맡기는 셈이다. 역원은 $a\cdot a^{-1}=1$만
확인하면 된다.
@(elliptic_test.go@>=
const p61 = 2305843009213693951 // 메르센 소수 $2^{61}-1$

var rng = mrand.New(mrand.NewPCG(42, 2026))

func TestFp(t *testing.T) {
	f := NewFp(p61)
	P := new(big.Int).SetUint64(p61)
	for range 200 {
		a, b := rng.Uint64()%p61, rng.Uint64()%p61
		A := new(big.Int).SetUint64(a)
		B := new(big.Int).SetUint64(b)
		check := func(op string, want *big.Int, got uint64) {
			if want.Mod(want, P).Uint64() != got {
				t.Fatalf("%s(%d, %d) = %d이 되었다", op, a, b, got)
			}
		}
		check("add", new(big.Int).Add(A, B), f.add(a, b))
		check("sub", new(big.Int).Sub(A, B), f.sub(a, b))
		check("mul", new(big.Int).Mul(A, B), f.mul(a, b))
		if a != 0 && f.mul(a, f.inv(a)) != 1 {
			t.Fatalf("inv(%d)가 역원이 아니다", a)
		}
	}
}
