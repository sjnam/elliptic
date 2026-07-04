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

# MetaPost 그림: ecfig.mp -> ecfig-1.pdf, ecfig-2.pdf
ecfig-1.pdf ecfig-2.pdf: ecfig.mp
	mptopdf ecfig.mp

elliptic.pdf: $(WSRC) ecfig-1.pdf ecfig-2.pdf
	$(GWEAVE) elliptic.w && luatex elliptic.tex </dev/null

test: tangle
	go test ./...

clean:
	rm -f elliptic.go elliptic_test.go
	rm -f elliptic.tex elliptic.log elliptic.toc elliptic.scn elliptic.idx
	rm -f ecfig.1 ecfig-1.pdf ecfig.2 ecfig-2.pdf ecfig.log ecfig.mpx mptextmp.mp mpxerr.tex
