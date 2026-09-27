# 타원 곡선 문학적 프로그램 빌드용 Makefile.
#
#   make            # tangle + 문서 조판
#   make tangle     # elliptic.w -> elliptic.go, elliptic_test.go
#   make doc        # elliptic.pdf 조판 (한글이라 luatex)
#   make test       # go test ./...
#   make clean      # 생성물 삭제 (.w 원본은 남김)
#
# 매크로(gwebmac.tex, kotexgweb.tex)는 설치된 texmf 트리에서 자동으로 찾는다.

GTANGLE ?= gtangle
GWEAVE  ?= gweave

WSRC := elliptic.w fp.w poly.w curve.w divpoly.w schoof.w dlp.w ecdsa.w

.PHONY: all tangle doc test clean
.DEFAULT_GOAL := all

all: tangle doc

tangle: elliptic.go

elliptic.go: $(WSRC)
	$(GTANGLE) elliptic.w

doc: elliptic.pdf

# 그림은 따로 만들 것이 없다. .w 안에 적힌 MetaPost 코드를 luamplib 이
# 조판 중에 직접 그리므로 mpost 를 돌릴 일이 없다.
elliptic.pdf: $(WSRC)
	$(GWEAVE) elliptic.w && luatex elliptic.tex </dev/null

test: tangle
	go test ./...

clean:
	rm -f elliptic.go elliptic_test.go
	rm -f elliptic.tex elliptic.log elliptic.toc elliptic.scn elliptic.idx
