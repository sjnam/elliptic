@s Int int
@s Reader int
@s WaitGroup int

\input kotexgweb.tex
\def\title{타원 곡선 답사기}

@* 들어가며.
이 문서는 타원 곡선을 스스로 익히려고 쓴 답사기다. 목적지는 두 곳이다.
하나는 곡선 위의 점을 세는 {\it Schoof 알고리즘} — 순수 수학이 어떻게 다항
시간 마술을 부리는지 보여 주는 봉우리이고, 다른 하나는 그 곡선을 무기로
쓰는 {\it 타원 곡선 암호} — 이산 로그의 어려움 위에 세운 서명 ECDSA다.
가는 길에 유한체, 다항식과 그 빠른 곱셈(NTT), 나눗셈 다항식, 이산 로그
공격 삼종을 차례로 지난다.

답사에는 사연이 있다. 전에 같은 것을 |math/big|만으로 짰더니 Schoof가
너무 느려, 60비트 소수 하나 세는 데 밤을 새웠다. 그래서 이번 판은 점 세기
엔진을 통째로 |uint64| 산술 위에 다시 세우고 다항식 곱셈에 수론적
변환(NTT)을 얹었다 — 학습용 장난감치고는 제법 이빨이 있다. 하지만 이 글의
진짜 목적은 속도가 아니라 이해다. 그래서 함수 하나하나에 ``무엇을''과 함께
``왜''를 적으려 했고, 곡선의 역사와 뒷이야기도 곳곳에 끼워 넣었다.

@ 프로그램은 한 \GO/ 패키지 |elliptic|으로 나온다. \.{gtangle}하면 주
출력 \.{elliptic.go}와 시험 파일 \.{elliptic\_test.go}가 함께 떨어진다.
층은 아래에서 위로 쌓인다.

$$\vbox{\halign{\indent#\hfil&\quad#\hfil\cr
{\bf 유한체}&    |uint64| 위의 $F_p$ 산술\cr
{\bf 다항식}&    $F_p[x]$와 NTT 빠른 곱셈\cr
{\bf 타원 곡선}& |big.Int| 위의 군 법칙 (그림과 함께)\cr
{\bf 나눗셈 다항식}& 등분점을 붙드는 $\psi_n$\cr
{\bf Schoof}&    프로베니우스로 점 세기\cr
{\bf 이산 로그}& Shanks, Pollard $\rho$, Pohlig--Hellman\cr
{\bf ECDSA}&     secp256k1 위의 서명\cr}}$$

낮은 두 층(유한체·다항식)은 속도가 생명이라 힙을 꺼리는 |uint64|로 짜고,
높은 층(곡선·암호)은 256비트 수를 다뤄야 하니 |big.Int|로 짠다. Schoof는
이 두 세계에 다리를 놓아, 큰 곡선의 문제를 작은 체의 빠른 산술로 푼다.
@c
package elliptic

import (
	"errors"
	"io"
	"math/big"
	"math/bits"
	"sort"
	"sync"

	"crypto/rand"
)

@<유한체@>@;
@<다항식@>@;
@<타원 곡선@>@;
@<나눗셈 다항식@>@;
@<Schoof 알고리즘@>@;
@<이산 로그 문제@>@;
@<ECDSA@>@;

@ 시험 파일은 패키지 안에서 돌며 내부 함수까지 들여다본다. 각 장이 제
시험을 \.{elliptic\_test.go}에 조금씩 보태므로, 여기서는 임포트만 모아
둔다. 무작위성은 재현 가능하도록 씨앗을 고정한 |math/rand/v2|의 PCG를 쓰되,
암호 연산(키 생성·서명)만은 진짜 |crypto/rand|를 쓴다.
@(elliptic_test.go@>=
package elliptic

import (
	"crypto/rand"
	"crypto/sha256"
	"math/big"
	mrand "math/rand/v2"
	"testing"
)

@i fp.w
@i poly.w
@i curve.w
@i divpoly.w
@i schoof.w
@i dlp.w
@i ecdsa.w

@* 찾아보기.
