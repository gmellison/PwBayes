bin_PROGRAMS=	pwb

pwbb_SOURCES=	\
	bayes.c \
	bayes.h \
	best.c \
	best.h \
	command.c \
	command.h \
	likelihood.c \
	likelihood.h \
	mbbeagle.c \
	mbbeagle.h \
	mcmc.c \
	mcmc.h \
	model.c \
	model.h \
	pairwise.c \
	pairwise.h \
	proposal.c \
	proposal.h \
	sumpt.c \
	sumpt.h \
	utils.c \
	utils.h 

pwb_CFLAGS=	$(SIMD_FLAGS) $(CPUEXT_FLAGS)

pwb_LDADD=	$(BEAGLE_LIBS)
pwb_CFLAGS+=	$(BEAGLE_CFLAGS)
