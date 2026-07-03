# 타원 곡선 답사기

타원 곡선과 그 암호를 스스로 익히려고 쓴 [GWEB](https://github.com/sjnam/gweb)
문학적 프로그램이다. `.w` 원본 하나에서 Go 패키지 `elliptic`과 한글
문서(PDF)가 함께 나온다. 자매편으로는 볼록 껍질 트릭 모음
[cht](https://github.com/sjnam/cht), 동시 큐
[go-lcrq](https://github.com/sjnam/go-lcrq)가 있다.

목적지는 두 곳이다. 하나는 곡선 위의 점을 세는 **Schoof 알고리즘**, 다른
하나는 그 곡선을 무기로 쓰는 **타원 곡선 암호**(ECDSA)다. 가는 길에 유한체,
다항식과 NTT 빠른 곱셈, 나눗셈 다항식, 이산 로그 공격 삼종을 지난다.

## 층 구성

| 장(`.w`) | 내용 |
| --- | --- |
| `fp.w` | `uint64` 위의 유한체 $F_p$ 산술 ($p<2^{63}$) |
| `poly.w` | $F_p[x]$ 다항식과 세 소수 NTT + CRT 빠른 곱셈 |
| `curve.w` | `big.Int` 위의 군 법칙 (MetaPost 그림과 함께) |
| `divpoly.w` | 등분점을 붙드는 나눗셈 다항식 $\psi_n$ |
| `schoof.w` | 프로베니우스 특성 방정식으로 점 세기 |
| `dlp.w` | Shanks, Pollard $\rho$, Pohlig–Hellman |
| `ecdsa.w` | secp256k1 위의 서명·검증 |

낮은 층(유한체·다항식)은 속도가 생명이라 힙을 꺼리는 `uint64`로, 높은
층(곡선·암호)은 256비트 수를 다뤄야 하니 `big.Int`로 짰다. Schoof가 두
세계에 다리를 놓아 큰 곡선의 문제를 작은 체의 빠른 산술로 푼다 — 옛
`big.Int` 전용 구현이 60비트 곡선에 며칠을 쓰던 것을 초 단위로 줄인 대목이다.

## 빌드

```bash
make          # tangle + 문서(PDF) 조판
make tangle   # elliptic.w -> elliptic.go, elliptic_test.go
make test     # go test ./...
make doc      # elliptic.pdf (한글이라 luatex, 그림은 MetaPost)
make clean    # 생성물 삭제 (.w 원본은 남김)
```

문서 조판에는 `gweave`/`gtangle`, `luatex`, `mpost`가 필요하고, GWEB 매크로
(`kotexgweb.tex`, `gwebmac.tex`)는 설치된 texmf 트리에서 자동으로 찾는다.
